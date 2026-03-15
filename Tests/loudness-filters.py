import numpy as np


def kweighting_coeffs_for_fs(b0_48k, b1_48k, b2_48k, a1_48k, a2_48k, fs):
    """
    Calcula os coeficientes digitais do biquad (b0,b1,b2,a1,a2)
    para uma taxa de amostragem fs diferente de 48000 Hz.
    Baseado na fórmula do Klangfreund / BS.1770.
    """
    # Derivações a partir dos coeficientes de 48 kHz
    KoverQ = (2.0 - 2.0 * a2_48k) / (a2_48k - a1_48k + 1.0)
    K = np.sqrt((a1_48k + a2_48k + 1.0) / (a2_48k - a1_48k + 1.0))
    Q = K / KoverQ
    arctanK = np.arctan(K)
    VB = (b0_48k - b2_48k) / (1.0 - a2_48k)
    VH = (b0_48k - b1_48k + b2_48k) / (a2_48k - a1_48k + 1.0)
    VL = (b0_48k + b1_48k + b2_48k) / (a1_48k + a2_48k + 1.0)

    # Recalcula K para a nova taxa de amostragem
    K_new = np.tan(arctanK * 48000.0 / fs)
    commonFactor = 1.0 / (1.0 + K_new / Q + K_new * K_new)

    # Coeficientes digitais ajustados
    b0 = (VH + VB * K_new / Q + VL * K_new * K_new) * commonFactor
    b1 = 2.0 * (VL * K_new * K_new - VH) * commonFactor
    b2 = (VH - VB * K_new / Q + VL * K_new * K_new) * commonFactor
    a1 = 2.0 * (K_new * K_new - 1.0) * commonFactor
    a2 = (1.0 - K_new / Q + K_new * K_new) * commonFactor

    return b0, b1, b2, a1, a2


# Exemplo: stage 1 do K-weighting
b0, b1, b2, a1, a2 = (
    1.53512485958697,
    -2.69169618940638,
    1.19839281085285,
    -1.69065929318241,
    0.73248077421585,
)
fs = 44100
coeffs = kweighting_coeffs_for_fs(b0, b1, b2, a1, a2, fs)
print("Coeficientes para fs =", fs)
print("b0,b1,b2,a1,a2 =", coeffs)
