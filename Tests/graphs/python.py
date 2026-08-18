import math
import numpy as np
import matplotlib.pyplot as plt


# Expected event duration in frames
mean = 100.0

# Negative Binomial dispersion values
r_values = [4, 8, 16, 32]

# Maximum duration to display
u_max = 220
u = np.arange(0, u_max + 1)


def negative_binomial_pmf(mu, r, x_values):
    """
    Negative Binomial parameterization:

        mean = r * (1 - p) / p

    therefore:

        p = r / (r + mean)

    and:

        variance = mean + mean^2 / r
    """

    p = r / (r + mu)

    pmf = np.zeros_like(x_values, dtype=float)

    for i, x in enumerate(x_values):
        log_pmf = (
            math.lgamma(x + r)
            - math.lgamma(r)
            - math.lgamma(x + 1)
            + r * math.log(p)
            + x * math.log1p(-p)
        )

        pmf[i] = math.exp(log_pmf)

    return pmf


def poisson_pmf(mu, x_values):
    pmf = np.zeros_like(x_values, dtype=float)

    for i, x in enumerate(x_values):
        log_pmf = (
            -mu
            + x * math.log(mu)
            - math.lgamma(x + 1)
        )

        pmf[i] = math.exp(log_pmf)

    return pmf


def survivor_from_pmf(pmf):
    """
    D(u) = P(U >= u)
    """
    return np.flip(np.cumsum(np.flip(pmf)))


# ---------------------------------------------------------
# Occupancy / duration PMF
# ---------------------------------------------------------

plt.figure(figsize=(10, 6))

for r in r_values:
    pmf = negative_binomial_pmf(mean, r, u)

    plt.plot(
        u,
        pmf,
        label=f"Negative Binomial, r={r}"
    )

poisson = poisson_pmf(mean, u)

plt.plot(
    u,
    poisson,
    "--",
    label="Poisson"
)

plt.axvline(
    mean,
    linestyle=":",
    label=f"Expected duration = {mean:.0f}"
)

plt.xlabel("Elapsed frames $u$")
plt.ylabel("Probability")
plt.title("Occupancy PMF $d_j(u)$")
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ---------------------------------------------------------
# Survivor function
# ---------------------------------------------------------

plt.figure(figsize=(10, 6))

for r in r_values:
    pmf = negative_binomial_pmf(mean, r, u)
    survivor = survivor_from_pmf(pmf)

    plt.plot(
        u,
        survivor,
        label=f"Negative Binomial, r={r}"
    )

poisson_survivor = survivor_from_pmf(poisson)

plt.plot(
    u,
    poisson_survivor,
    "--",
    label="Poisson"
)

plt.axvline(
    mean,
    linestyle=":",
    label=f"Expected duration = {mean:.0f}"
)

plt.xlabel("Elapsed frames $u$")
plt.ylabel("Probability")
plt.title("Survivor function $D_j(u) = P(U \\geq u)$")
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
