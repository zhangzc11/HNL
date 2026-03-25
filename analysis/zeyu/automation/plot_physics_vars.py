import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ==========================================
# 1. 设定你要观察的变量
# ==========================================
# 你可以每次换一个变量来画图观察，比如 'Jet1_PT', 'MuNuR_PT', 'W_PT'
var_name = 'MuNuR_PT'  
tree_name = "myTupleOS1J/DecayTree"

# ==========================================
# 2. 定义文件路径 (Data + 几个有代表性的MC)
# ==========================================
data_path = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"

# 挑选三个不同质量的 MC 进行对比 (5GeV, 15GeV, 50GeV)
mc_paths = {
    "5GeV_10ps": "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912011_5GeVCtau10ps_options_20260228_111712_skim.root",
    "15GeV_10ps": "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912015_15GeVCtau10ps_options_20260228_111712_skim.root",
    "50GeV_10ps": "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912021_50GeVCtau10ps_options_20260228_111712_skim.root"
}

# ==========================================
# 3. 读取数据
# ==========================================
print(f"[*] 正在读取变量: {var_name}")
with uproot.open(data_path) as f:
    df_data = f[tree_name].arrays([var_name], library="pd")

mc_dfs = {}
for label, path in mc_paths.items():
    with uproot.open(path) as f:
        mc_dfs[label] = f[tree_name].arrays([var_name], library="pd")

# ==========================================
# 4. 绘图 (画在同一张图上)
# ==========================================
plt.figure(figsize=(10, 7))

# 确定统一的画图范围 (排除极个别异常值，取 1% 到 99% 的分位数比较美观)
min_val = np.percentile(df_data[var_name], 1)
max_val = np.percentile(mc_dfs["50GeV_10ps"][var_name], 99)
bins = np.linspace(min_val, max_val, 60)

# 画 Data 本底 (用灰色填充)
plt.hist(df_data[var_name], bins=bins, density=True, alpha=0.3, color='black', label='Background (Data)')

# 画 3 种 MC 信号 (用不同颜色的线条)
colors = ['blue', 'green', 'red']
for (label, df), color in zip(mc_dfs.items(), colors):
    plt.hist(df[var_name], bins=bins, density=True, histtype='step', linewidth=2.5, color=color, label=f'Signal MC ({label})')

plt.title(f'Physics Comparison: {var_name} across different HNL Masses', fontsize=14)
plt.xlabel(var_name, fontsize=12)
plt.ylabel('Normalized Events / a.u.', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)

out_name = f"Physics_Compare_{var_name}.png"
plt.savefig(out_name, dpi=300)
print(f"[+] 物理对比图已保存为: {out_name}")
