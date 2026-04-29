import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import uproot
import glob
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================
# 1. 配置路径与参数
# ==========================================
# 请确保此路径指向你截图中的目录
data_dir = "/home/ben/HNL/analysis/zeyu/HNLdata" 
mass_points = [5, 10, 15, 20, 30, 50]
colors = ['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
m_W = 80.4

# ==========================================
# 2. 理论数值积分
# ==========================================
def theory_W_pt_pdf(pt, peak_val=8.0):
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    tail = 0.05 / (pt + 5)**2
    return np.where(pt > 0, sudakov + tail, 0)

# 构建积分网格 
grid_size = 200
g_pt = np.linspace(0.01, 80.0, grid_size)
g_cos = np.linspace(-1.0, 1.0, grid_size)
g_phi = np.linspace(0, 2*np.pi, grid_size)
PT_W_3d, Cos_3d, Phi_3d = np.meshgrid(g_pt, g_cos, g_phi, indexing='ij')

Weights_3d = theory_W_pt_pdf(PT_W_3d) * (0.5) * (1.0/(2*np.pi))
Weights_3d /= np.sum(Weights_3d)

PT_W_f, Cos_f, Phi_f, W_f = PT_W_3d.ravel(), Cos_3d.ravel(), Phi_3d.ravel(), Weights_3d.ravel()
Sin_f = np.sqrt(1 - Cos_f**2)
gamma_T = np.sqrt(m_W**2 + PT_W_f**2) / m_W
beta_T = PT_W_f / (m_W * gamma_T)

# ==========================================
# 3. ROOT 数据读取 
# ==========================================
def get_mc_data(m_N, branch_name):
    # 匹配截图中的 429120xx_... 格式
    file_pattern = f"{data_dir}/*_{m_N}GeVCtau0ps_*.root"
    files = glob.glob(file_pattern)
    if not files:
        return np.array([])
    
    data = []
    try:
        with uproot.open(files[0]) as f:
            for tree_name in ["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]:
                if tree_name in f and branch_name in f[tree_name].keys():
                    data.extend(f[tree_name][branch_name].array(library="np"))
        return np.array(data) / 1000.0  # MeV -> GeV
    except:
        return np.array([])

# ==========================================
# 4. 绘图主程序
# ==========================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
bins_pt = np.linspace(0, 60, 100)
bin_centers = 0.5 * (bins_pt[:-1] + bins_pt[1:])

for m, c in zip(mass_points, colors):
    print(f"[*] Processing m_N = {m} GeV...")
    
    # --- A. 理论部分 (虚线) ---
    p_star = (m_W**2 - m**2) / (2 * m_W)
    E_N_star = (m_W**2 + m**2) / (2 * m_W)
    
    # Muon 理论
    pt_mu_th = np.sqrt((gamma_T*(p_star*Sin_f*np.cos(Phi_f) + beta_T*p_star))**2 + (p_star*Sin_f*np.sin(Phi_f))**2)
    h_mu_th, _ = np.histogram(pt_mu_th, bins=bins_pt, weights=W_f, density=True)
    ax1.plot(bin_centers, h_mu_th, color=c, linestyle='--', lw=1.8, alpha=0.7)

    # Neutrino 理论
    pt_N_th = np.sqrt((gamma_T*(-p_star*Sin_f*np.cos(Phi_f) + beta_T*E_N_star))**2 + (-p_star*Sin_f*np.sin(Phi_f))**2)
    h_N_th, _ = np.histogram(pt_N_th, bins=bins_pt, weights=W_f, density=True)
    ax2.plot(bin_centers, h_N_th, color=c, linestyle='--', lw=1.8, alpha=0.7)

    # --- B. ROOT 数据 (实线) ---
    # 使用你确认的分支名 MuR_PT 和 NuR_PT
    mc_mu = get_mc_data(m, "Mu_PT")
    if len(mc_mu) > 0:
        ax1.hist(mc_mu, bins=bins_pt, density=True, histtype='step', color=c, lw=2.2, label=f'$m_N={m}$ GeV')
        
    mc_N = get_mc_data(m, "NuR_PT")
    if len(mc_N) > 0:
        ax2.hist(mc_N, bins=bins_pt, density=True, histtype='step', color=c, lw=2.2)

# ==========================================
# 5. 图表修饰
# ==========================================
# 统一图例
custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='-'),
                Line2D([0], [0], color='black', lw=2, linestyle='--')]
legend_labels = ['MC Full Sim (Solid)', 'Numerical Theory (Dashed)']

for ax, title in zip([ax1, ax2], [r"Prompt Muon $P_{T,\mu_1}$", r"Heavy Neutrino $P_{T,N}$"]):
    ax.set_title(title, fontsize=15)
    ax.set_xlabel("$P_T$ [GeV]", fontsize=12)
    ax.set_ylabel("Normalized Density", fontsize=12)
    ax.set_xlim(0, 60)
    ax.grid(alpha=0.2, linestyle=':')
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + custom_lines, l + legend_labels, loc='upper right')

plt.tight_layout()
plt.savefig("Validation_Mu_NuR_PT.png", dpi=300)

print("[+] 保存成功：Validation_Mu_NuR_PT.png")
