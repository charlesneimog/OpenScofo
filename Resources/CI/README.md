# OpenScofo CI commands

`ci.py` detects Linux, macOS, or Windows with Python's `platform` module and
runs the matching implementation. The GitHub Actions workflow only selects an
operation; it does not contain platform-specific shell scripts.

Examples:

```sh
python Resources/CI/ci.py --install-dependencies
python Resources/CI/ci.py --install-qt-framework
python Resources/CI/ci.py --build --precision 32
python Resources/CI/ci.py --files-test
python Resources/CI/ci.py --stage-puredata-package
```

Run `python Resources/CI/ci.py --help` for the complete list of operations.
