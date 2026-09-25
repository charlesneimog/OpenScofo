import OpenScofo
import librosa
from pathlib import Path


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
            # "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/tongue_ram/Fl-tng_ram-A3-mf-N-N.wav"
            # "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/pizzicato/Fl-pizz-D4-f-N-N.wav"
            "/home/neimog/Nextcloud/Music-Resources/Samples/Orchidea/Winds/Flute/jet_whistle/Fl-jet_wh-N-N-N-N.wav"
        )

    def list_SOL_files(self):
        return list(Path(self.orchidea_sol).rglob("*.wav"))

    def check_gate(self, audio_path):
        import matplotlib.pyplot as plt

        audio, _ = librosa.load(audio_path, sr=self.sr, mono=True)

        times = []
        silence_gate = []
        sound_gate = []
        tech_gate = []
        pitch_gate = []

        for i in range(0, len(audio), self.hop_size):
            block = audio[i : i + self.hop_size]

            if len(block) < self.hop_size:
                block = librosa.util.fix_length(block, size=self.hop_size)

            if not self.scofo.process_block(block):
                raise RuntimeError("Failed to process block")

            desc = self.scofo.get_description()

            silence_prob = desc.silence
            ext_prob = desc.ext

            sound_prob = max(0.0, 1.0 - silence_prob)
            tech_weight = ext_prob
            pitch_weight = 1.0 - ext_prob

            times.append(i / self.sr)

            silence_gate.append(silence_prob)
            sound_gate.append(sound_prob)
            tech_gate.append(tech_weight * sound_prob)
            pitch_gate.append(desc.pitch_confidence * pitch_weight * sound_prob)

        plt.figure(figsize=(14, 5))

        plt.plot(times, silence_gate, label="Silence")
        plt.plot(times, tech_gate, label="Extended Technique")
        plt.plot(times, pitch_gate, label="Pitch / Note")

        plt.xlabel("Time (s)")
        plt.ylabel("Gate")
        plt.ylim(0.0, 1.0)
        plt.legend()
        plt.grid(alpha=0.3)

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    ExtendedTechniqueGate()
