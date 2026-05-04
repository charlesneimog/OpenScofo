import OpenScofo
import librosa

scofo = OpenScofo.OpenScofo(48000, 2048, 512)


def process_score(score_file, audio_file):
    ok = scofo.load_score(score_file)
    if not ok:
        raise Exception("Failed to load score")

    x, _ = librosa.load(audio_file, sr=48000)
    current_pos = 0
    for i in range(0, len(x), 512):
        if i + 512 > len(x):
            break

        block = x[i : i + 512]
        scofo.process_block(block)
        pos = scofo.get_current_score_position()
        if pos != current_pos:
            current_pos = pos
            print("New event " + str(current_pos))


process_score(
    "../Resources/Extended-Techniques-Random-Tests/score-1.txt",
    "../Resources/Extended-Techniques-Random-Tests/score-1.wav",
)

process_score(
    "../Resources/Extended-Techniques-Random-Tests/score-2.txt",
    "../Resources/Extended-Techniques-Random-Tests/score-2.wav",
)
