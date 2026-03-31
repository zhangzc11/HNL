import numpy as np
import matplotlib.pyplot as plt

# =========================
# Parameters
# =========================
mW = 80.4  # GeV
mN_list = [5, 10, 15, 20, 30, 50]  # GeV

# =========================
# Analytic integral:
# I(smax) = \int_0^{smax} ds [(mN^2-s)(mN^2+2s)]/(s-mW^2)^2
# =========================
def integral_analytic(smax, mN, mW):
    mN2 = mN**2
    mW2 = mW**2

    term1 = -2.0 * smax
    term2 = (4.0 * mW2 - mN2) * np.log(mW2 / (mW2 - smax))
    term3 = - (mN2 + 2.0 * mW2) * (mW2 - mN2) * (1.0 / (mW2 - smax) - 1.0 / mW2)

    return term1 + term2 + term3

# =========================
# Unnormalized dGamma/dM
# Since dGamma/d(M^2) ∝ I(smax),
# then dGamma/dM = 2M * dGamma/d(M^2) ∝ 2M * I(smax)
# =========================
def dGammadM(M, mN, mW):
    mW2 = mW**2
    mN2 = mN**2
    M2 = M**2

    # kinematic endpoint from smax >= 0
    Mmax2 = mW2 - mN2
    if Mmax2 <= 0:
        return np.zeros_like(M)

    val = np.zeros_like(M, dtype=float)

    mask = (M2 >= 0) & (M2 <= Mmax2)
    if not np.any(mask):
        return val

    smax = mN2 * (1.0 - M2[mask] / (mW2 - mN2))
    I = integral_analytic(smax, mN, mW)

    val[mask] = 2.0 * M[mask] * I
    return val

# =========================
# Plot
# =========================
plt.figure(figsize=(8, 6))

for mN in mN_list:
    Mmax = np.sqrt(mW**2 - mN**2)
    M = np.linspace(1e-6, Mmax, 1200)

    y = dGammadM(M, mN, mW)

    # numerical normalization
    norm = np.trapz(y, M)
    if norm > 0:
        y /= norm

    plt.plot(M, y, lw=2, label=fr"$m_N={mN}\,\mathrm{{GeV}}$")

plt.xlabel(r"$M_{\mu\mu}\ \mathrm{[GeV]}$", fontsize=13)
plt.ylabel(r"Normalized distribution", fontsize=13)
plt.title(r"Normalized $M_{\mu\mu}$ distributions", fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()