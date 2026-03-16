import uproot
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. 配置参数
# ==========================================
# 您提供的 branch 列表（示例，请根据您的 filtered_branches 替换）
filtered_branches = [
    "W_MINIPCHI2",                            # 核心几何：W 主顶点冲击参数卡方
    "MuNuR_MINIP",                            # 核心几何：重建系统最小冲击参数
    "NuR_MAX_SDOCACHI2",                      # 顶点质量：径迹间最近距离卡方
    "NuR_ALV",                                # 顶点辅助变量
    "Mu_PT",                                  # 动力学：缪子横向动量
    "NuR_MIN_PT",                             # 动力学门槛：衰变子粒子最小横向动量
    "W_PT",                                   # 动力学：系统总横向动量
    "W_ETA",                                  # 角度：系统赝快度
    "Mu_PROBNN_MU",                           # PID：缪子神经网络概率
    "Mu_PROBNN_GHOST",                        # 径迹质量：幽灵径迹概率
    "Mu_PROBNN_PI",                           # 背景排查：pion 神经网络概率
    "Mu_HCALEOP",                             # 鉴别辅助：HCAL 能量参数
    "Mu_POSITION_STATEAT_LastMeasurement_Z",  # 空间位置：最后测量点 Z 坐标
    "nTTracks",                               # 事件复杂度：总径迹数
    "Mu_BREMTRACKBASEDENERGY",                # 能量修正：轫致辐射能量
    "MuNuR_E",                                # 能量量级：系统总能量
    "NuR_log_MIN_P",                          # 辅助动力学：最小动量对数
    "W_THETA"                                 # 角度：极角
]

# 文件路径
data_path = "/home/steven/HNL/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
signal_path = "/home/steven/HNL/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"

# 假设 ROOT 文件中的 Tree 名称为 "DecayTree"，请根据实际情况修改
tree_name = "myTupleSS1J/DecayTree"

def load_data(path, branches):
    """从 ROOT 文件加载数据并转换为 Pandas DataFrame"""
    with uproot.open(path) as file:
        tree = file[tree_name]
        # 只读取存在的 branch，避免报错
        available_branches = [b for b in branches if b in tree.keys()]
        df = tree.arrays(available_branches, library="pd")
    return df

# ==========================================
# 2. 读取数据
# ==========================================
print("正在读取信号与数据文件...")
df_sig = load_data(signal_path, filtered_branches)
df_data = load_data(data_path, filtered_branches)

# 处理潜在的空值
df_sig.dropna(inplace=True)
df_data.dropna(inplace=True)

print(f"信号事件数: {len(df_sig)}, 数据事件数: {len(df_data)}")

# ==========================================
# 3. 计算 KS 统计量 (Kolmogorov-Smirnov Test)
# ==========================================
print("--- KS 统计量分析 (区分度排序) ---")
ks_results = []

for branch in df_sig.columns:
    if branch in df_data.columns:
        # 计算 KS 统计量和 P 值
        # ks_stat 越大，表示信号与背景分布差异越明显
        ks_stat, p_val = ks_2samp(df_sig[branch], df_data[branch])
        ks_results.append({"Branch": branch, "KS_Stat": ks_stat, "P_Value": p_val})

# 转换为 DataFrame 并按 KS 统计量降序排列
ks_df = pd.DataFrame(ks_results).sort_values(by="KS_Stat", ascending=False)
print(ks_df.to_string(index=False))

# ==========================================
# 4. 计算相关矩阵 (Correlation Matrix)
# ==========================================
print("正在计算相关矩阵...")

def plot_corr(df, title):
    plt.figure(figsize=(10, 8))
    corr = df.corr()
    mask = np.triu(np.ones_like(corr, dtype=bool)) # 只显示下三角
    sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap='coolwarm', center=0)
    plt.title(title)
    plt.tight_layout()

# 计算信号的相关性（了解物理变量间的内在联系）
plot_corr(df_sig, "Correlation Matrix: Signal (MC)")
plt.savefig("corr_signal.png")

# 计算数据的相关性（了解背景变量间的相关性，防止过拟合冗余变量）
plot_corr(df_data, "Correlation Matrix: Data (Background)")
plt.savefig("corr_data.png")

print("相关矩阵图已保存为 corr_signal.png 和 corr_data.png")

# ==========================================
# 5. 综合建议输出
# ==========================================
top_features = ks_df[ks_df['KS_Stat'] > 0.1]['Branch'].tolist()
print(f"[建议]: 基于 KS 统计量，以下变量对信号与背景有较好区分度: \n{top_features}")
print("注意：若某两个变量相关性 > 0.9，建议在机器学习中只保留 KS 值更高的那一个。")
