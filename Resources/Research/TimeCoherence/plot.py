import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/home/neimog/Documents/Git/OScofo/Resources/Research/TimeCoherence/bwv-1013.csv")

plt.figure(figsize=(10, 4))
plt.plot(df["time"], df["probability"])
plt.xlabel("Time (s)")
plt.ylabel("Probability")
plt.title("Time Coherence Template")
plt.ylim(0, 1.05)
plt.tight_layout()
plt.show()

