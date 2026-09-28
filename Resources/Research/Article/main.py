import OpenScofo
import librosa
from pathlib import Path
import math


class ExtendedTechniqueGate:
    def __init__(self):
        print("Init...")
        self.sr = 48000
        self.fft_size = 2048
        self.hop_size = 512
        self.scofo = OpenScofo.OpenScofo(self.sr, self.fft_size, self.hop_size)
        self.orchidea_sol = "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea"
        self.orchidea_sol_files = self.list_SOL_files()

        self.db_silence_gate = -60
        self.first_x_ms_after_silence = 200

        self.scofo.activate_all_descriptors()
        self.check_gate(
            {
                "Tongue Ram": "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/tongue_ram/Fl-tng_ram-A3-mf-N-N.wav",
                "Pizzicato": "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/pizzicato/Fl-pizz-D4-f-N-N.wav",
                "Jet Whistle": "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/jet_whistle/Fl-jet_wh-N-N-N-N.wav",
            }
        )

    def list_SOL_files(self):
        return list(Path(self.orchidea_sol).rglob("*.wav"))

    def check_gate(self, audio_paths):
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(14, 10),
            sharex=False,
            sharey=True,
        )

        for ax, (name, audio_path) in zip(axes, audio_paths.items()):
            audio, _ = librosa.load(audio_path, sr=self.sr, mono=True)

            times = []
            silence_gate = []
            tech_gate = []
            pitch_gate = []

            for i in range(0, len(audio), self.hop_size):
                block = audio[i : i + self.hop_size]

                if len(block) < self.hop_size:
                    block = librosa.util.fix_length(
                        block,
                        size=self.hop_size,
                    )

                if not self.scofo.process_block(block):
                    raise RuntimeError("Failed to process block")

                desc = self.scofo.get_description()

                silence_prob = desc.silence
                ext_prob = desc.ext

                sound_prob = 1.0 - silence_prob
                tech_weight = ext_prob
                pitch_weight = 1.0 - ext_prob

                times.append(i / self.sr)

                silence_gate.append(silence_prob)
                tech_gate.append(tech_weight * sound_prob)
                pitch_gate.append(pitch_weight * sound_prob)

            ax.plot(times, silence_gate, label="Silence")
            ax.plot(times, tech_gate, label="Extended Technique")
            ax.plot(times, pitch_gate, label="Pitch / Note")

            ax.set_title(name)
            ax.set_ylabel("Gate")
            ax.set_ylim(0.0, 1.0)
            ax.grid(alpha=0.3)
            ax.legend(loc="upper right")

        axes[-1].set_xlabel("Time (s)")

        fig.suptitle("OpenScofo Gate Analysis", fontsize=16)

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    ExtendedTechniqueGate()
