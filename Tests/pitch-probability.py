import OpenScofo
import librosa
import numpy as np
import json
from pathlib import Path

scofo = OpenScofo.OpenScofo(48000, 2048, 512)


def collect_pitch_data(audio_path, expected_freq, score_positions):
    """
    Coleta dados de pitch para calibração.
    expected_freq: frequência esperada (ex: 440 para A4)
    score_positions: lista de posições na partitura onde a nota ocorre
    """
    x, sr = librosa.load(audio_path, sr=48000)
    frame_size = 2048

    data = []
    for pos in score_positions:
        # Pega o frame na posição esperada
        start = int(pos * sr)
        if start + frame_size > len(x):
            continue

        frame = x[start : start + frame_size]
        scofo.process_block(frame)

        # Obtém probabilidade para a frequência correta
        prob_correct = scofo.get_pitch_prob(expected_freq)

        # Obtém probabilidade para uma frequência errada (ex: 5 semitons acima)
        wrong_freq = expected_freq * (2 ** (5 / 12))
        prob_wrong = scofo.get_pitch_prob(wrong_freq)

        data.append(
            {
                "position": pos,
                "prob_correct": prob_correct,
                "prob_wrong": prob_wrong,
                "ratio": prob_correct / (prob_wrong + 1e-6),
            }
        )

    return data


# Exemplo de uso
data = collect_pitch_data(
    audio_path="/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/ordinario/Fl-ord-A4-mf-N-N.wav",
    expected_freq=440.0,
    score_positions=[0.5, 1.0, 1.5, 2.0],  # Ajuste conforme sua partitura
)

print(json.dumps(data, indent=2))
