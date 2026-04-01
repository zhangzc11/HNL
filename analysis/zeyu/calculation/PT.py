import numpy as np
import matplotlib.pyplot as plt

# 设置6个质量点的颜色
mass_points = [5, 10, 15, 20, 30, 50]
colors = ['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']

# 1. 定义纯理论的 W 玻色子产生 P_T 分布 (Generator Level)
def generator_W_pt(pt, peak_val=8.0):
    # Sudakov 峰
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    # QCD 长尾
    tail = 0.05 / (pt + 5)**2
    dist = sudakov + tail
    
    # 使用 np.where 处理边界
    return np.where(pt <= 0, 0, dist)

# 2. 定义探测器接受度效率函数 (Acceptance Efficiency)
def acceptance_efficiency(pt, m_N):
    # 重中微子质量越大，越需要更高的 W_PT 才能通过 LHCb 的 Cuts
    turn_on_point = 5 + 0.42 * m_N 
    width = 2.5 + 0.08 * m_N
    return 1.0 / (1.0 + np.exp(-(pt - turn_on_point) / width))

# 生成横动量数组 (严格从0开始)
pt_w = np.linspace(0, 100, 500)

plt.figure(figsize=(11, 6.5))

# --- 关键修改处：将 np.trapz 全部替换为 np.trapezoid ---

# 绘制本底 (Background Mock)
bg_pt = generator_W_pt(pt_w, peak_val=12.0)
# 修复位置 1
bg_pt /= np.trapezoid(bg_pt, pt_w)  
plt.fill_between(pt_w, 0, bg_pt * 0.08, color='gray', alpha=0.5, label='Background (Theory Mock)')

# 信号生成
base_pt = generator_W_pt(pt_w, peak_val=8.0)

# 绘制6个质量点的重构后 P_T 分布
for m, color in zip(mass_points, colors):
    # 核心物理逻辑：重构分布 = 理论产生分布 * 探测器接受度
    reconstructed_pt = base_pt * acceptance_efficiency(pt_w, m)
    
    # 修复位置 2
    reconstructed_pt /= np.trapezoid(reconstructed_pt, pt_w) 
    
    # 画图，适当缩放使高度与背景协调
    plt.plot(pt_w, reconstructed_pt * 0.08, color=color, lw=2.5, label=f'Signal ({m} GeV)')

# 完善图表信息
plt.title("Theoretical Explanation of $W_{P_T}$ Separation (Generator $\\times$ Acceptance)", fontsize=16, pad=15)
plt.xlabel("Theoretical Reconstructed $W_{P_T}$ [GeV]", fontsize=14)
plt.ylabel("Normalized Events / a.u.", fontsize=14)

# 物理洞察文本框
info_text = (
    "$\\bf{Physics\\ Insight:}$\n"
    "1. The $W_{P_T}$ originates from QCD initial state radiation (Sudakov peak + pQCD tail).\n"
    "2. At the generator level, all signal masses share the same $W_{P_T}$ distribution.\n"
    "3. However, heavier HNLs have less momentum available in their rest frame.\n"
    "4. To pass LHCb's strict $P_T$ cuts, heavier signals require a larger transverse\n"
    "   Lorentz boost from the $W$ boson, forcibly shifting the observed peak to the right."
)
plt.text(0.35, 0.75, info_text, transform=plt.gca().transAxes, fontsize=11, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

plt.legend(fontsize=12, loc='upper right')
plt.xlim(0, 100)
plt.ylim(bottom=0) 
plt.grid(alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig("W_PT_Explanation_Fixed.png", dpi=300)
plt.show()