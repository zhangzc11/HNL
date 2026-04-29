import numpy as np
import matplotlib.pyplot as plt
import uproot
import glob

# =============================
# 基本设置
# =============================
mass_points = [5, 10, 15, 20, 30, 50]
colors = ['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
data_dir = "../HNLdata"

# =============================
# PT.py: Generator-level
# =============================
def generator_W_pt(pt, peak_val=8.0):
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    tail = 0.05 / (pt + 5)**2
    return np.where(pt <= 0, 0, sudakov + tail)

# =============================
# PT.py: Acceptance
# =============================
def acceptance_efficiency(pt, m_N):
    turn_on_point = 5 + 0.42 * m_N 
    width = 2.5 + 0.08 * m_N
    return 1.0 / (1.0 + np.exp(-(pt - turn_on_point) / width))

# =============================
# 理论分布（PT方法 + 归一化 + 截断）
# =============================
def get_theory_pt(m_N, pt_max=80):

    pt = np.linspace(0, pt_max, 1000)

    dist = generator_W_pt(pt) * acceptance_efficiency(pt, m_N)

    # ✔ 归一化（关键）
    dist /= np.trapezoid(dist, pt)

    return pt, dist

# =============================
# 读取 MC 数据
# =============================
def get_mc_data(file_pattern, branch_name):
    files = glob.glob(file_pattern)
    if not files:
        return np.array([])

    data = []
    try:
        with uproot.open(files[0]) as f:
            for tree_name in ["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]:
                if tree_name in f and branch_name in f[tree_name].keys():
                    data.extend(f[tree_name][branch_name].array(library="np"))
        return np.array(data) / 1000.0
    except Exception:
        return np.array([])

# =============================
# 主程序
# =============================
print("[*] 计算 W_PT 理论积分 (带 20GeV 截断推导)...")
plt.figure(figsize=(9, 7))
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']

for m, c in zip(mass_points, colors):

    file_pattern = f"{data_dir}/*_{m}GeVCtau0ps_*.root"

    # =============================
    # MC 数据（实线）
    # =============================
    mc_data = get_mc_data(file_pattern, "W_PT")

    if len(mc_data) > 0:
        plt.hist(
            mc_data,
            bins=60,
            range=(0, 80),   # ✔ 保留 W_PT 截断
            density=True,
            histtype='step',
            color=c,
            lw=2.0
        )

    # =============================
    # 理论（虚线）
    # =============================
    pt, dist = get_theory_pt(m)

    plt.plot(
        pt,
        dist,
        linestyle='--',   # ✔ 虚线
        color=c,
        lw=2.5,
        label=f'$m_N={m}$ GeV'
    )

# =============================
# 画图设置
# =============================
plt.title(r"Validation: $W$ Boson Transverse Momentum ($W_{P_T}$)", fontsize=16)
plt.xlabel(r"$W_{P_T}$ [GeV]", fontsize=14)
plt.ylabel("Normalized Density", fontsize=14)

plt.xlim(0, 80)  # ✔ 截断
plt.grid(alpha=0.3)

from matplotlib.lines import Line2D
handles, labels = plt.gca().get_legend_handles_labels()

custom_lines = [
    Line2D([0], [0], color='black', lw=2, linestyle='-'),
    Line2D([0], [0], color='black', lw=2, linestyle='--')
]

plt.legend(
    handles + custom_lines,
    labels + ['MC Full Sim', 'Theory (Generator × Acceptance)'],
    loc='upper right'
)

plt.tight_layout()
plt.savefig("Validation_W_PT.png", dpi=300)

print("[+] 保存成功：Validation_W_PT.png")