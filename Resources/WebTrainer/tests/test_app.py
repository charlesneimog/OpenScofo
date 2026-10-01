"""Regression checks for untrusted ZIP uploads and reviewed dataset state."""

from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import app


class DatasetTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.base = Path(self.work.name)

    def archive(self, entries):
        path = self.base / "uploaded-hash.zip"
        with zipfile.ZipFile(path, "w") as archive:
            for name in entries:
                archive.writestr(
                    name, b"recording" if not str(name).endswith("/") else b""
                )
        return path

    def extract(self, entries):
        return app.extract_dataset(self.archive(entries), self.base / "dataset")

    def test_direct_classes_and_all_extensions(self):
        root, counts = self.extract(
            ["normal/a.wav", "normal/b.aiff", "key-click/c.aif"]
        )
        self.assertEqual(root, self.base / "dataset")
        self.assertEqual(counts, [("key-click", 1), ("normal", 2)])

    def test_optional_wrapper_and_ignored_files(self):
        root, counts = self.extract(
            [
                "Flute/normal/a.WAV",
                "Flute/key-click/b.AIFF",
                "Flute/key-click/c.AIF",
                "README.md",
                "Flute/README.txt",
                "Flute/.DS_Store",
                "Flute/.hidden/audio.wav",
                "__MACOSX/Flute/._normal.wav",
            ]
        )
        self.assertEqual(root.name, "Flute")
        self.assertEqual(counts, [("key-click", 2), ("normal", 1)])
        self.assertTrue((root / "normal/a.wav").exists())

    def test_unsafe_paths_even_for_ignored_files(self):
        for path in [
            "../../file",
            "/absolute.wav",
            "normal/../../escape.wav",
            "normal\\a.wav",
            "C:/normal/a.wav",
        ]:
            with self.subTest(path=path), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(app.DatasetError, "unsafe"):
                    app.extract_dataset(
                        self.archive([path]), Path(directory) / "dataset"
                    )
        self.assertFalse((self.base / "escape.wav").exists())

    def test_symlink(self):
        entry = zipfile.ZipInfo("normal/link.wav")
        entry.create_system = 3
        entry.external_attr = (stat.S_IFLNK | 0o777) << 16
        with self.assertRaisesRegex(app.DatasetError, "unsafe"):
            self.extract([entry])

    def test_invalid_zip(self):
        path = self.base / "invalid.zip"
        path.write_text("not a ZIP")
        with self.assertRaisesRegex(app.DatasetError, "valid, unencrypted ZIP"):
            app.extract_dataset(path, self.base / "dataset")

    def test_one_class(self):
        with self.assertRaisesRegex(app.DatasetError, "two technique folders"):
            self.extract(["normal/a.wav"])

    def test_empty_class(self):
        with self.assertRaisesRegex(app.DatasetError, "empty folders"):
            self.extract(["normal/a.wav", "click/"])

    def test_unsupported_only_class(self):
        with self.assertRaisesRegex(app.DatasetError, "empty folders"):
            self.extract(["normal/a.wav", "click/README.txt"])

    def test_no_audio(self):
        with self.assertRaisesRegex(app.DatasetError, "No supported audio"):
            self.extract(["normal/README.md", "click/example.mp3"])

    def test_nested_audio(self):
        with self.assertRaisesRegex(app.DatasetError, "nested folders"):
            self.extract(["normal/session/a.wav", "click/b.wav"])

    def test_root_audio(self):
        with self.assertRaisesRegex(app.DatasetError, "audio outside"):
            self.extract(["a.wav", "normal/b.wav", "click/c.wav"])

    def test_case_normalization_collision(self):
        with self.assertRaisesRegex(app.DatasetError, "duplicate"):
            self.extract(["normal/a.WAV", "normal/a.wav", "click/b.wav"])

    def test_limits(self):
        with (
            patch.object(app, "MAX_BYTES", 1),
            self.assertRaisesRegex(app.DatasetError, "too large"),
        ):
            self.extract(["normal/a.wav", "click/b.wav"])

    def test_preview_cleans_files_and_enables_train(self):
        path = self.archive(["normal/a.wav", "click/b.wav"])
        with patch.object(app, "WORK_DIR", self.base / "work"):
            rows, summary, reviewed, button, download, status = app.preview_dataset(
                str(path)
            )
            self.assertEqual(rows, [("click", 1), ("normal", 1)])
            self.assertIn("2 audio files", summary)
            self.assertEqual(reviewed["digest"], app.archive_digest(path))
            self.assertTrue(button.interactive)
            self.assertEqual(list((self.base / "work").iterdir()), [])
        self.assertIsNone(download)

    def test_invalid_preview_disables_train(self):
        result = app.preview_dataset(str(self.archive(["normal/a.wav"])))
        self.assertIsNone(result[2])
        self.assertFalse(result[3].interactive)

    def test_training_requires_review(self):
        result = list(app.train_model(None, None, app.DESCRIPTORS, 10, 0.05, 2, 5))
        self.assertIn("inspect", result[-1][0])

    def test_changed_upload_requires_new_review(self):
        path = self.archive(["normal/a.wav", "click/b.wav"])
        with patch.object(app, "WORK_DIR", self.base / "work"):
            reviewed = app.preview_dataset(str(path))[2]
            self.archive(["different/a.wav", "click/b.wav"])
            result = list(
                app.train_model(str(path), reviewed, app.DESCRIPTORS, 10, 0.05, 2, 5)
            )
            self.assertIn("upload changed", result[-1][0])
            self.assertEqual(list((self.base / "work").iterdir()), [])


if __name__ == "__main__":
    unittest.main()
