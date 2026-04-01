import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) # 忽略解析解在边界的除零警告

# ==========================================
# 严格的理论数值积分 (Semi-Analytical Integration)
# 完全不使用随机数，求解极其平滑的理论曲线
# ==========================================

m_W = 80.4
mass_points =[5, 10, 15, 20, 30, 50]
colors =['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']

# --- 1. W_PT 理论概率密度函数 ---
def theory_W_pt_pdf(pt, peak_val=8.0):
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    tail = 0.05 / (pt + 5)**2
    dist = sudakov + tail
    return np.where(pt > 0, dist, 0)

# --- 2. 积分网格构建 (Grid Generation) ---
# 定义高密度网格以替代随机蒙特卡洛，等效于完美的三重积分
# N_PTW * N_cos * N_phi = 200 * 200 * 200 = 800万个确定性积分点
grid_PT_W = np.linspace(0.01, 80.0, 200)
grid_cos  = np.linspace(-1.0, 1.0, 200)
grid_phi  = np.linspace(0, 2*np.pi, 200)

# 创建 3D 网格
PT_W_3d, Cos_3d, Phi_3d = np.meshgrid(grid_PT_W, grid_cos, grid_phi, indexing='ij')

# 计算每个格点的积分权重 (f(P_TW) * dPT * dcos * dphi)
dPT = grid_PT_W[1] - grid_PT_W[0]
dCos = grid_cos[1] - grid_cos[0]
dPhi = grid_phi[1] - grid_phi[0]

# 权重 = P(P_TW) * P(cos) * P(phi) * 积分微元
# 因为 cos 均匀分布在[-1,1] 密度为1/2; phi 均匀分布在[0,2pi] 密度为1/(2pi)
Weights_3d = theory_W_pt_pdf(PT_W_3d) * (1.0/2.0) * (1.0/(2*np.pi)) * (dPT * dCos * dPhi)
# 将权重归一化（理论上积分应为1）
Weights_3d /= np.sum(Weights_3d)

# 将网格展平以便计算
PT_W_flat = PT_W_3d.ravel()
Cos_flat = Cos_3d.ravel()
Sin_flat = np.sqrt(1 - Cos_flat**2)
Phi_flat = Phi_3d.ravel()
W_flat = Weights_3d.ravel()

# 横向洛伦兹因子
M_T_flat = np.sqrt(m_W**2 + PT_W_flat**2)
gamma_T = M_T_flat / m_W
beta_T = PT_W_flat / M_T_flat

# ================= 开始绘图 =================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

bins_pt = np.linspace(0, 60, 200) # 非常精细的 200 个 bin，展示平滑曲线
bin_centers = 0.5 * (bins_pt[:-1] + bins_pt[1:])

for m_N, color in zip(mass_points, colors):
    print(f"正在进行解析与数值积分计算: m_N = {m_N} GeV ...")
    
    # 物理常数
    p_star = (m_W**2 - m_N**2) / (2 * m_W)
    E_N_star = (m_W**2 + m_N**2) / (2 * m_W)
    E_mu_star = p_star
    
    # ----------------------------------------------------
    # 情景 A: 纯解析解 (Idealized, P_T,W = 0)
    # 利用刚刚推导的公式： f(Pt) = Pt / (p* \sqrt{p*^2 - Pt^2})
    # ----------------------------------------------------
    pt_ideal = np.linspace(0, p_star - 0.001, 1000) # 避开奇点
    pdf_ideal = pt_ideal / (p_star * np.sqrt(p_star**2 - pt_ideal**2))
    
    ax1.plot(pt_ideal, pdf_ideal, color=color, linestyle=':', lw=2.0, alpha=0.8)
    ax2.plot(pt_ideal, pdf_ideal, color=color, linestyle=':', lw=2.0, alpha=0.8) # 对于N也是同样的公式形状，只是P_T上限也是p_star

    # ----------------------------------------------------
    # 情景 B: 数值积分 (Realistic, P_T,W > 0)
    # 使用展平的网格计算实验室系下的横动量，并按权重累加
    # ----------------------------------------------------
    # 1. 瞬时缪子 (mu1)
    px_mu_lab = gamma_T * (p_star * Sin_flat * np.cos(Phi_flat) + beta_T * E_mu_star)
    py_mu_lab = p_star * Sin_flat * np.sin(Phi_flat)
    pt_mu_lab = np.sqrt(px_mu_lab**2 + py_mu_lab**2)
    
    hist_mu, _ = np.histogram(pt_mu_lab, bins=bins_pt, weights=W_flat)
    hist_mu_pdf = hist_mu / (bins_pt[1] - bins_pt[0]) # 转为概率密度
    ax1.plot(bin_centers, hist_mu_pdf, color=color, linestyle='-', lw=2.5, label=f'$m_N={m_N}$ GeV')

    # 2. 重中微子 (N) - 注意静止系中 N 与 mu1 动量反向 (cos 和 phi 有负号)
    px_N_lab = gamma_T * (-p_star * Sin_flat * np.cos(Phi_flat) + beta_T * E_N_star)
    py_N_lab = -p_star * Sin_flat * np.sin(Phi_flat)
    pt_N_lab = np.sqrt(px_N_lab**2 + py_N_lab**2)
    
    hist_N, _ = np.histogram(pt_N_lab, bins=bins_pt, weights=W_flat)
    hist_N_pdf = hist_N / (bins_pt[1] - bins_pt[0]) # 转为概率密度
    ax2.plot(bin_centers, hist_N_pdf, color=color, linestyle='-', lw=2.5, label=f'$m_N={m_N}$ GeV')

# ================= 图表修饰 =================
ax1.axvline(20, color='black', linestyle='--', lw=2.5, label='LHCb Cut (20 GeV)')

# 图 1 设置
ax1.set_title(r"【解析与数值积分】瞬时缪子横动量 $P_{T,\mu_1}$", fontsize=15, pad=10)
ax1.set_xlabel(r"$P_{T,\mu_1}$ [GeV]", fontsize=14)
ax1.set_ylabel("Probability Density", fontsize=14)
ax1.set_xlim(0, 60)
ax1.set_ylim(0, 0.28)

# 图 2 设置
ax2.set_title(r"【解析与数值积分】重中微子横动量 $P_{T,N}$", fontsize=15, pad=10)
ax2.set_xlabel(r"$P_{T,N}$ [GeV]", fontsize=14)
ax2.set_ylabel("Probability Density", fontsize=14)
ax2.set_xlim(0, 60)
ax2.set_ylim(0, 0.28)

# 自定义完美图例
custom_lines =[Line2D([0], [0], color=c, lw=3) for c in colors]
custom_lines.extend([
    Line2D([0], [0], color='gray', lw=2, linestyle='-'),
    Line2D([0], [0], color='gray', lw=2, linestyle=':')
])
legend_labels =[f'$m_N = {m}$ GeV' for m in mass_points] +[r'数值积分 (含 $W_{P_T}$)', r'解析解 ($P_{T,W}=0$)']

ax1.legend(custom_lines, legend_labels, fontsize=12, loc='upper left')
ax2.legend(custom_lines, legend_labels, fontsize=12, loc='upper right')

for ax in [ax1, ax2]:
    ax.grid(alpha=0.3, linestyle='-.')

plt.suptitle(r"基于 $W_{P_T}$ 理论的子粒子横动量理论曲线 (网格数值积分法，无统计涨落)", fontsize=18, y=1.02)
plt.tight_layout()
plt.savefig("NuPT2_Analytical_Smooth.png", dpi=300, bbox_inches='tight')
plt.show()
print("计算完成，图像已保存。")