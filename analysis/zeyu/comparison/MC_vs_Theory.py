import numpy as np
import matplotlib.pyplot as plt
import uproot
import glob
import os

# ========================================================
# 1. 纯理论解析公式定义 (Quantum V-A 理论)
# ========================================================
M_W = 80.4

def pdf_M_mumu_quantum(M, m_N):
    """理论双缪子不变质量解析公式 (V-A理论下的Michel谱积分)"""
    C = M_W**2 - m_N**2
    max_M = np.sqrt(C)
    y = M**2 / C
    pdf = (4 * M / C) * (5./6. - 1.5 * y**2 + (2./3.) * y**3)
    return np.where(M < max_M, pdf, 0)

def get_pt_theoretical_distribution(m_N, particle='mu', bins=100, pt_max=60):
    """利用四维张量网格进行平滑的理论边缘积分投影 (绝对光滑的理论横动量曲线)"""
    N_grid = 50
    x = np.linspace(0.01, 0.99, N_grid)
    c_th_N = np.linspace(-0.99, 0.99, N_grid)
    c_th_S = np.linspace(-0.99, 0.99, N_grid)
    phi = np.linspace(0.01, 2*np.pi-0.01, N_grid)
    
    X, CT_N, CT_S, PHI = np.meshgrid(x, c_th_N, c_th_S, phi, indexing='ij')
    
    # 量子 V-A 理论下的 Michel Spectrum 权重
    weight = 2 * X**2 * (3 - 2*X)  
    weight = weight / np.sum(weight)
    
    gamma_N = (M_W**2 + m_N**2) / (2 * M_W * m_N)
    beta_N = (M_W**2 - m_N**2) / (M_W**2 + m_N**2)
    
    ST_N = np.sqrt(1 - CT_N**2)
    ST_S = np.sqrt(1 - CT_S**2)
    
    if particle == 'mu':
        E_star = X * m_N / 2.0
    else: # jet
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

# ========================================================
# 2. 读取本地 LHCb MC ROOT 数据的函数 (修改了相对路径)
# ========================================================
def get_mc_data(file_pattern, branch_name):
    """使用 uproot 从对应的 root 文件中读取分支数据 (自动除以 1000 转为 GeV)"""
    files = glob.glob(file_pattern)
    if not files:
        print(f"[-] 警告: 未找到匹配的文件: {file_pattern}")
        return np.array([])
    
    target_file = files[0]
    data =[]
    try:
        with uproot.open(target_file) as f:
            # 兼容读取 OS1J 和 OS2J 目录下的树
            for tree_name in ["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]:
                if tree_name in f:
                    tree = f[tree_name]
                    if branch_name in tree.keys():
                        arr = tree[branch_name].array(library="np")
                        data.extend(arr)
        
        # LHCb 的动量单位通常是 MeV，这里除以 1000 转换为 GeV
        return np.array(data) / 1000.0  
    except Exception as e:
        print(f"[-] 读取 {target_file} 报错: {e}")
        return np.array([])

# ========================================================
# 3. 主程序：绘制同框验证图
# ========================================================
def main():
    masses =[5, 10, 15, 20, 30, 50]
    colors =['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
    
    # 核心修改：利用相对路径指向上一级的 HNLdata 目录
    data_dir = "../HNLdata"
    
    # 将质量对应到具体文件名的通配符
    file_map = {
        5: f"{data_dir}/*_5GeVCtau0ps_*.root",
        10: f"{data_dir}/*_10GeVCtau0ps_*.root",
        15: f"{data_dir}/*_15GeVCtau0ps_*.root",
        20: f"{data_dir}/*_20GeVCtau0ps_*.root",
        30: f"{data_dir}/*_30GeVCtau0ps_*.root",
        50: f"{data_dir}/*_50GeVCtau0ps_*.root"
    }

    # 创建 1x3 的并排大图
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plt.rcParams['font.sans-serif'] =['DejaVu Sans', 'Arial']

    print("[*] 开始读取MC数据与理论积分，这大约需要 20-30 秒，请稍候...")
    
    for m, c in zip(masses, colors):
        print(f"  --> 正在处理 m_N = {m} GeV ...")
        
        # --- 读取 MC 数据 ---
        file_pattern = file_map[m]
        mc_mumu = get_mc_data(file_pattern, "W_MmuWmuN")  # 双缪子不变质量
        mc_pt_mu = get_mc_data(file_pattern, "MuNuR_PT")   # 次级缪子横动量
        mc_pt_jet = get_mc_data(file_pattern, "Jet1_PT")   # 喷注横动量
        
        # --- 图 1: M_mumu ---
        if len(mc_mumu) > 0:
            axes[0].hist(mc_mumu, bins=80, range=(0, 85), density=True, 
                         histtype='step', color=c, lw=2.0, linestyle='-')
        # 理论 M_mumu (虚线)
        x_m = np.linspace(0, 85, 500)
        axes[0].plot(x_m, pdf_M_mumu_quantum(x_m, m), color=c, lw=2.5, linestyle='--')

        # --- 图 2: P_{T, \mu_N} ---
        if len(mc_pt_mu) > 0:
            axes[1].hist(mc_pt_mu, bins=60, range=(0, 60), density=True, 
                         histtype='step', color=c, lw=2.0, linestyle='-')
        # 理论 P_{T, \mu_N} (虚线)
        b_mu, d_mu = get_pt_theoretical_distribution(m, particle='mu', pt_max=60)
        axes[1].plot(b_mu, d_mu, color=c, lw=2.5, linestyle='--')

        # --- 图 3: P_{T, jet} ---
        if len(mc_pt_jet) > 0:
            axes[2].hist(mc_pt_jet, bins=60, range=(0, 60), density=True, 
                         histtype='step', color=c, lw=2.0, linestyle='-')
        # 理论 P_{T, jet} (虚线)
        b_jet, d_jet = get_pt_theoretical_distribution(m, particle='jet', pt_max=60)
        axes[2].plot(b_jet, d_jet, color=c, lw=2.5, linestyle='--')

    # --- 统一装饰图表 ---
    titles =[r"Di-muon Mass $M_{\mu\mu}$", r"Secondary Muon $P_{T, \mu_N}$", r"Jet $P_{T, \text{jet}}$"]
    
    for i, ax in enumerate(axes):
        ax.set_title(titles[i], fontsize=15)
        ax.set_xlabel("[GeV]", fontsize=13)
        ax.set_ylabel("Normalized Density", fontsize=13)
        ax.grid(alpha=0.3)
        
        # 添加手动图例以说明实线和虚线的含义
        from matplotlib.lines import Line2D
        custom_lines = [
            Line2D([0], [0], color='black', lw=2, linestyle='-'),
            Line2D([0],[0], color='black', lw=2, linestyle='--')
        ]
        ax.legend(custom_lines, ['MC Full Sim (Solid)', 'Pure Theory (Dashed)'], loc='upper right', fontsize=11)

    axes[0].set_xlim(0, 85)
    axes[1].set_xlim(0, 60)
    axes[2].set_xlim(0, 60)

    plt.suptitle("Validation: LHCb Full Simulation vs. Quantum V-A Theoretical Model", fontsize=18, y=1.02)
    plt.tight_layout()
    
    output_filename = "MC_vs_Theory_Validation.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n[+] 大功告成！MC 与纯理论同框叠加的验证图已保存为 {output_filename}")

if __name__ == "__main__":
    main()