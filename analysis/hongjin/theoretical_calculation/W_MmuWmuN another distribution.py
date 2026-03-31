import numpy as np
import matplotlib.pyplot as plt

# =========================
# Parameters
# =========================
mW = 80.4  # GeV
mN_list = [5, 10, 15, 20, 30, 50]  # GeV

# =========================
# Distribution p(m)
# p(m) = 4m(C-m^2)/C^2,  C = mW^2 - mN^2
# valid for 0 <= m <= sqrt(C)
# =========================
def p_m(m, mN, mW):
    C = mW**2 - mN**2
    y = np.zeros_like(m, dtype=float)
    if C <= 0:
        return y

    mask = (m >= 0) & (m <= np.sqrt(C))
    y[mask] = 4.0 * m[mask] * (C - m[mask]**2) / C**2
    return y

# =========================
# Plot
# =========================
plt.figure(figsize=(8, 6))

for mN in mN_list:
    C = mW**2 - mN**2
    if C <= 0:
        continue

    mmax = np.sqrt(C)
    m = np.linspace(0, mmax, 1200)
    y = p_m(m, mN, mW)

    plt.plot(m, y, lw=2, label=fr"$m_N={mN}\,\mathrm{{GeV}}$")

plt.xlabel(r"$M_{\mu\mu}\ \mathrm{[GeV]}$", fontsize=13)
plt.ylabel(r"$p(M_{\mu\mu})$", fontsize=13)
plt.title(r"Distribution $p(m)=\frac{4m(C-m^2)}{C^2}$", fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()