import numpy as np
import matplotlib.pyplot as plt

# =========================
# 基本参数
# =========================
mW = 80.4  # GeV
mN_list = [5, 10, 15, 20, 30, 50]

# =========================
# 原函数 F(s)
# =========================
def F_of_s(s, mN, mW):
    mN2 = mN**2
    mW2 = mW**2
    term1 = -2.0 * s
    term2 = (4.0 * mW2 - mN2) * np.log(mW2 / (mW2 - s))
    term3 = -(mN2 + 2.0 * mW2) * (mW2 - mN2) * (1.0 / (mW2 - s) - 1.0 / mW2)
    return term1 + term2 + term3

def s_plus(p, mN, mW):
    return mN**2 - 2.0 * mN**2 * p / mW

def s_minus(p, mN, mW):
    return np.maximum(0.0, mN**2 - 2.0 * mW * p)

# =========================
# dΓ/dp
# 只保留形状
# =========================
def dGamma_dp(p, mN, mW):
    p = np.asarray(p)
    sp = s_plus(p, mN, mW)
    sm = s_minus(p, mN, mW)

    val = np.zeros_like(p, dtype=float)
    mask = (p >= 0.0) & (p <= mW / 2.0) & (sp >= 0.0)
    val[mask] = F_of_s(sp[mask], mN, mW) - F_of_s(sm[mask], mN, mW)
    return np.maximum(val, 0.0)

# =========================
# dΓ/dpT
# =========================
def dGamma_dpT(pT_array, mN, mW, n_int=1200):
    pT_array = np.asarray(pT_array)
    result = np.zeros_like(pT_array, dtype=float)

    p_max = mW / 2.0

    for i, pT in enumerate(pT_array):
        if pT <= 0.0 or pT >= p_max:
            result[i] = 0.0
            continue

        umax = np.sqrt(p_max**2 - pT**2)
        u = np.linspace(0.0, umax, n_int)

        p = np.sqrt(pT**2 + u**2)
        dp_du = np.zeros_like(u)
        dp_du[1:] = u[1:] / p[1:]
        dp_du[0] = 0.0

        kernel = np.zeros_like(u)
        mask = p > pT
        kernel[mask] = pT / (p[mask]**2 * np.sqrt(1.0 - pT**2 / p[mask]**2))

        integrand = dGamma_dp(p, mN, mW) * kernel * dp_du
        result[i] = np.trapz(integrand, u)

    return result

# =========================
# 以 lg pT 为横轴作图
# 横轴 x = log10(pT / GeV)
# 若要画成概率密度 dΓ/d(log10 pT)，需乘上 Jacobian: pT * ln(10)
# =========================
plt.figure(figsize=(8, 5.5))

pT_min = 1e-2          # 避免 log(0)
pT_max = mW / 2.0 - 1e-2
pT_grid = np.logspace(np.log10(pT_min), np.log10(pT_max), 320)

for mN in mN_list:
    y_pT = dGamma_dpT(pT_grid, mN, mW, n_int=1400)

    # 转成对 log10(pT) 的分布:
    # dΓ/d(log10 pT) = (dΓ/dpT) * dpT/d(log10 pT) = (dΓ/dpT) * pT * ln(10)
    y_log = y_pT * pT_grid * np.log(10.0)

    # 对 log10(pT) 归一化
    x_log = np.log10(pT_grid)+9
    area = np.trapz(y_log, x_log)
    if area > 0:
        y_log = y_log / area

    plt.plot(x_log, y_log, label=fr"$m_N={mN}\,\mathrm{{GeV}}$")

plt.xlabel(r"$\log_{10}\!\left(p_T(\mu_N)/\mathrm{GeV}\right)$")
plt.ylabel("Normalized distribution")
plt.title(r"Normalized distribution of $\mu_N$ versus $\log_{10} p_T$")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()