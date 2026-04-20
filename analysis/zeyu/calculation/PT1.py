import numpy as np
import matplotlib.pyplot as plt

# 设置 6 个质量点和颜色
mass_points =[5, 10, 15, 20, 30, 50]
colors =['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']

# 1. 定义纯理论的 W 玻色子产生 P_T 概率密度函数 (PDF)
def theory_W_pt_pdf(pt, peak_val=8.0):
    # Sudakov 峰 (Rayleigh 分布形态契合软胶子重和特征)
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    # 高 P_T 的微扰 QCD 长尾 (1/P_T^2 衰减)
    tail = 0.05 / (pt + 5)**2
    dist = sudakov + tail
    dist = np.where(pt > 0, dist, 0)
    return dist

# 2. 核心改进：使用舍选法 (Accept-Reject) 真正“生成” MC 事件
def generate_W_PT_events(N, peak_val=8.0):
    print(f"正在生成 {N} 个 W_PT 物理事件 (Peak={peak_val})...")
    pt_events = np.zeros(N)
    count = 0
    # 估算 PDF 的最大值 (在 peak_val 附近)
    max_pdf = theory_W_pt_pdf(np.array([peak_val]), peak_val)[0] * 1.2 
    
    while count < N:
        x = np.random.uniform(0.01, 100.0, N)
        y = np.random.uniform(0, max_pdf, N)
        valid = x[y < theory_W_pt_pdf(x, peak_val)]
        n_v = len(valid)
        if n_v > 0:
            fill = min(n_v, N - count)
            pt_events[count:count+fill] = valid[:fill]
            count += fill
    return pt_events

# 3. 探测器接受度概率 (Acceptance Probability)
def acceptance_prob(pt, m_N):
    turn_on_point = 5 + 0.42 * m_N 
    width = 2.5 + 0.08 * m_N
    return 1.0 / (1.0 + np.exp(-(pt - turn_on_point) / width))

# ================= 开始绘图 =================
plt.figure(figsize=(12, 7))
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

N_total = 500000  # 模拟五十万个真实事件
bins = np.linspace(0, 100, 100)

# 生成本底 Mock 数据 (Peak=12)
bg_events = generate_W_PT_events(N_total, peak_val=12.0)
plt.hist(bg_events, bins=bins, density=True, color='gray', alpha=0.4, label='Background (Theory Mock, Peak=12)')

# 生成基础信号数据 (Peak=8)
signal_base_events = generate_W_PT_events(N_total, peak_val=8.0)

# 遍历 6 个质量点，施加物理 Cut
for m, color in zip(mass_points, colors):
    # 计算每个事件通过探测器的概率
    probs = acceptance_prob(signal_base_events, m)
    # 抛骰子决定该事件是否被探测器记录
    random_dice = np.random.uniform(0, 1, N_total)
    passed_events = signal_base_events[random_dice < probs]
    
    # 画出通过探测器的事件分布
    plt.hist(passed_events, bins=bins, density=True, histtype='step', color=color, lw=2.5, label=f'Signal ($m_N={m}$ GeV)')

# 图表修饰
plt.title(r"基于 MC 抽样的 $W_{P_T}$ 分离演化 (Generator Level $\to$ Detector Accept)", fontsize=16, pad=15)
plt.xlabel(r"Reconstructed $W_{P_T}$ [GeV]", fontsize=14)
plt.ylabel("Probability Density", fontsize=14)

info_text = (
    "这里我们使用了 $\mathbf{Accept-Reject}$ 抽样法，\n"
    "创造了 50 万个底层的物理碰撞事件，\n"
)
plt.text(0.35, 0.70, info_text, transform=plt.gca().transAxes, fontsize=12, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))

plt.legend(fontsize=12, loc='upper right')
plt.xlim(0, 100)
plt.grid(alpha=0.3, linestyle='--')
plt.tight_layout()
save_path = "W_PT_MC_Simulation.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"图像已成功保存")
plt.show()