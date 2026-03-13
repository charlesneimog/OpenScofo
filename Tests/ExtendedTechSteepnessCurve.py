import numpy as np
import matplotlib.pyplot as plt

# Generate points between 0 and 1
x = np.linspace(0, 1, 500)

# Calculate the curves
linear = x
exponential = x**2
sigmoid_12 = 1 / (1 + np.exp(-15 * (x - 0.5)))
sigmoid_50 = 1 / (1 + np.exp(-50 * (x - 0.5)))
tanh_curve = 0.5 * (1 + np.tanh(15 * (x - 0.5)))  # tanh version
step = np.where(x > 0.5, 1.0, 0.0)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x, linear, label="Original (Linear)", linestyle="--", color="gray")
plt.plot(x, exponential, label="Exponential (x^2)", color="blue")
plt.plot(x, sigmoid_12, label="Sigmoid (steepness=15)", color="green")
plt.plot(x, sigmoid_50, label="Sigmoid (steepness=50)", color="orange")
plt.plot(x, tanh_curve, label="tanh (steepness=15)", color="purple")
plt.plot(x, step, label="Hard Step (> 0.5)", linestyle=":", color="red", linewidth=2)

# Styling
plt.axvline(0.5, color="black", alpha=0.3, linestyle="-.")  # Center line
plt.title("Probability Shaping Curves")
plt.xlabel("Input Probability")
plt.ylabel("Modified Probability")
plt.legend()
plt.grid(True, alpha=0.5)

# Save or show
# plt.savefig("probability_curves.png")
plt.show()
