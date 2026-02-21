import os
import time
import numpy as np
import matplotlib.pyplot as plt
import random
import OScofo

# --- Configuration ---
sr = 48000
fft_size = 4096
hop = 1024
eps = 1e-12  # Safety epsilon to avoid log(0) or division by zero

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# --- Initialize OScofo ---
scofo = OScofo.OScofo(sr, fft_size, hop)
ok = scofo.parse_score("../Resources/tests/canticos.txt")

if not ok:
    raise Exception("Failed to parse score.")

# --- Reference Template (The "Observation") ---
# We treat Event 5 at a specific time as our 'ground truth' for this test
ref_pos = 5
t_in_event = 128
time_template_ref = scofo.get_time_coherence_template(ref_pos, t_in_event)

# Prepare P (Observation Distribution)
P = np.asarray(time_template_ref, dtype=np.float64)
P = np.clip(P, eps, None)
P /= P.sum()  # Normalize to PMF

# --- Search for Best Match (Hypothesis Testing) ---
best_kl = float("inf")  # Mathematically, we want to MINIMIZE divergence
best_temp = None
best_r = 0
best_t = 0

print(f"Searching for best match to Event {ref_pos}...")

for i in range(20):  # Increased iterations for a better search
    r = random.randint(0, 10)
    t = random.randint(0, 200)
    if r == 5:
        print(r, t)

    # Get Hypothesis Template (Q)
    time_template_hypo = scofo.get_time_coherence_template(r, t)
    Q = np.asarray(time_template_hypo, dtype=np.float64)
    Q = np.clip(Q, eps, None)
    Q /= Q.sum()  # Normalize to PMF

    # Calculate Forward Kullback-Leibler Divergence: D_KL(P || Q)
    # This measures the "information gain" or "extra bits" required
    # if we use Q to represent P.
    kl_PQ = np.sum(P * np.log(P / Q))

    # Update if this hypothesis is a better fit (lower divergence)
    if kl_PQ < best_kl:
        best_kl = kl_PQ
        best_temp = time_template_hypo
        best_r = r
        best_t = t

# --- Results ---
# Likelihood estimation: exp(-D_KL) is often used as a similarity weight
likelihood = np.exp(-best_kl)

print("-" * 30)
print(f"Target: Event {ref_pos}, Time {t_in_event}")
print(f"Best Match: Event {best_r}, Time {best_t}")
print(f"KL Divergence: {best_kl:.6f} (Lower is better)")
print(f"Similarity Score: {likelihood:.4f} (1.0 is identical)")
print("-" * 30)

if best_temp is None:
    raise Exception("No valid template found.")

# --- Visualization ---
plt.figure(figsize=(10, 6))
plt.plot(time_template_ref, linewidth=2, label=f"Reference (Ev:{ref_pos})", alpha=0.8)
plt.plot(
    best_temp, linewidth=2, linestyle="--", label=f"Best Match (Ev:{best_r})", alpha=0.8
)

plt.title("Time Coherence Template Matching via KL Divergence")
plt.xlabel("Analysis Frames (Time)")
plt.ylabel("Probability Density")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.show()
