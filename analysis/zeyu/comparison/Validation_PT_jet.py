import numpy as np
import matplotlib.pyplot as plt
import uproot
import glob

# ========================================================
# 1. 理论网格积分与数据读取
# ========================================================
M_W = 80.4

def get_pt_theoretical_distribution_jet(m_N, bins=100, pt_max=60):
    """利用四维张量网格进行平滑的理论边缘积分投影 (针对 Jet)"""
    N_grid = 50
    x = np.linspace(0.01, 0.99, N_grid)
    c_th_N = np.linspace(-0.99, 0.99, N_grid)
    c_th_S = np.linspace(-0.99, 0.99, N_grid)
    phi = np.linspace(0.01, 2*np.pi-0.01, N_grid)
    
    X, CT_N, CT_S, PHI = np.meshgrid(x, c_th_N, c_th_S, phi, indexing='ij')
    
    weight = 2 * X**2 * (3 - 2*X)  # Michel Spectrum
    weight = weight / np.sum(weight)
    
    gamma_N = (M_W**2 + m_N**2) / (2 * M_W * m_N)
    beta_N = (M_W**2 - m_N**2) / (M_W**2 + m_N**2)
    
    ST_N = np.sqrt(1 - CT_N**2)
    ST_S = np.sqrt(1 - CT_S**2)
    
    # 核心改变：对于 Jet 的能量 E*
    E_star = m_N * (1 - X / 2.0)
        
    term1 = gamma_N * E_star * (CT_S + beta_N) * ST_N
    term2 = E_star * ST_S * np.cos(PHI) * CT_N
    p_lab_x = term1 + term2
    p_lab_y = E_star * ST_S * np.sin(PHI)
    PT = np.sqrt(p_lab_x**2 + p_lab_y**2)
    
    hist, bin_edges = np.histogram(PT.flatten(), bins=bins, range=(0, pt_max), weights=weight.flatten())
    density = hist / (bin_edges[1] - bin_edges[0])
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    return bin_centers, density

def get_mc_data(file_pattern, branch_name):
    """读取 LHCb MC ROOT 数据"""
    files = glob.glob(file_pattern)
    if not files: return np.array([])
    data =[]
    try:
        with uproot.open(files[0]) as f:
            for tree_name in ["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]:
                if tree_name in f and branch_name in f[tree_name].keys():
                    data.extend(f[tree_name][branch_name].array(library="np"))
        return np.array(data) / 1000.0  # MeV to GeV
    except Exception: return np.array([])

# ========================================================
# 2. 绘图主程序
# ========================================================
def main():
    masses =[5, 10, 15, 20, 30, 50]
    colors =['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
    data_dir = "../HNLdata"
    
    plt.figure(figsize=(9, 7))
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
    
    print("[*] 正在计算理论四重积分生成 Jet 横动量验证图，请稍候(约5秒)...")
    
    for m, c in zip(masses, colors):
        file_pattern = f"{data_dir}/*_{m}GeVCtau0ps_*.root"
        
        # 1. 读取 MC
        mc_pt_jet = get_mc_data(file_pattern, "Jet1_PT")
        if len(mc_pt_jet) > 0:
            plt.hist(mc_pt_jet, bins=60, range=(0, 60), density=True, 
                     histtype='step', color=c, lw=2.0)
            
        # 2. 纯理论计算
        b_jet, d_jet = get_pt_theoretical_distribution_jet(m, pt_max=60)
        plt.plot(b_jet, d_jet, color=c, lw=2.5, linestyle='--', label=f'$m_N={m}$ GeV')

    plt.title(r"Validation: Jet $P_T$ ($P_{T, \text{jet}}$)", fontsize=16)
    plt.xlabel(r"$P_{T, \text{jet}}$ [GeV]", fontsize=14)
    plt.ylabel("Normalized Density", fontsize=14)
    plt.xlim(0, 60)
    plt.grid(alpha=0.3)
    
    from matplotlib.lines import Line2D
    handles, labels = plt.gca().get_legend_handles_labels()
    custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='-'),
                    Line2D([0], [0], color='black', lw=2, linestyle='--')]
    plt.legend(handles + custom_lines, labels +['MC Full Sim (Solid)', 'Pure Theory (Dashed)'], 
               loc='upper right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig("Validation_PT_jet.png", dpi=300)
    print("[+] 保存成功：Validation_PT_jet.png")

if __name__ == "__main__":
    main()