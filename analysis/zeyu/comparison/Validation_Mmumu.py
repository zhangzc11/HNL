import numpy as np
import matplotlib.pyplot as plt
import uproot
import glob

# ========================================================
# 1. 理论公式与数据读取
# ========================================================
M_W = 80.4

def pdf_M_mumu_quantum(M, m_N):
    """理论双缪子不变质量解析公式 (Quantum V-A Theory)"""
    C = M_W**2 - m_N**2
    max_M = np.sqrt(C)
    y = M**2 / C
    pdf = (4 * M / C) * (5./6. - 1.5 * y**2 + (2./3.) * y**3)
    return np.where(M < max_M, pdf, 0)

def get_mc_data(file_pattern, branch_name):
    """读取 LHCb MC ROOT 数据"""
    files = glob.glob(file_pattern)
    if not files:
        return np.array([])
    target_file = files[0]
    data =[]
    try:
        with uproot.open(target_file) as f:
            for tree_name in["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]:
                if tree_name in f:
                    tree = f[tree_name]
                    if branch_name in tree.keys():
                        data.extend(tree[branch_name].array(library="np"))
        return np.array(data) / 1000.0  # MeV to GeV
    except Exception as e:
        print(f"读取报错: {e}")
        return np.array([])

# ========================================================
# 2. 绘图主程序
# ========================================================
def main():
    masses =[5, 10, 15, 20, 30, 50]
    colors =['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
    data_dir = "../HNLdata"
    
    plt.figure(figsize=(9, 7))
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
    
    print("[*] 正在生成双缪子不变质量 M_mumu 验证图...")
    
    for m, c in zip(masses, colors):
        file_pattern = f"{data_dir}/*_{m}GeVCtau0ps_*.root"
        
        # 1. 读取 MC 数据并画实线直方图
        mc_mumu = get_mc_data(file_pattern, "W_MmuWmuN")
        if len(mc_mumu) > 0:
            plt.hist(mc_mumu, bins=80, range=(0, 85), density=True, 
                     histtype='step', color=c, lw=2.0)
            
        # 2. 计算纯理论并画虚线
        x_m = np.linspace(0, 85, 500)
        plt.plot(x_m, pdf_M_mumu_quantum(x_m, m), color=c, lw=2.5, linestyle='--', label=f'$m_N={m}$ GeV')

    # 图表修饰
    plt.title(r"Validation: Di-muon Invariant Mass ($M_{\mu\mu}$)", fontsize=16)
    plt.xlabel(r"$M_{\mu\mu}$ [GeV]", fontsize=14)
    plt.ylabel("Normalized Density", fontsize=14)
    plt.xlim(0, 85)
    plt.grid(alpha=0.3)
    
    # 自定义图例 (区分实线和虚线)
    from matplotlib.lines import Line2D
    handles, labels = plt.gca().get_legend_handles_labels()
    custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='-'),
                    Line2D([0], [0], color='black', lw=2, linestyle='--')]
    
    # 合并质量点图例和线型说明
    plt.legend(handles + custom_lines, labels + ['MC Full Sim (Solid)', 'Pure Theory (Dashed)'], 
               loc='upper right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig("Validation_Mmumu.png", dpi=300)
    print("[+] 保存成功：Validation_Mmumu.png")

if __name__ == "__main__":
    main()