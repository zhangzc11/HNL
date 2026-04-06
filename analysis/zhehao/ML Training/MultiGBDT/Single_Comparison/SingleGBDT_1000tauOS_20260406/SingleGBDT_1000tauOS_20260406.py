import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split

# ==============================
# 1. 路径与常量设置
# ==============================
DATA_DIR = "/home/steven/HNL" 
DATA_FILE = os.path.join(DATA_DIR, "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root")
file_names = [
    "42912031_5GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912033_10GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912035_15GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912037_20GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912039_30GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912041_50GeVCtau1000ps_options_20260228_111712_skim.root"
]
MC_FILES = sorted([os.path.join(DATA_DIR, f) for f in file_names])

# 这里定义 1J 和 2J 两个 Tree
TREE_1J = "myTupleOS1J/DecayTree" 
TREE_2J = "myTupleOS2J/DecayTree" 

OUTPUT_DIR = "./SingleGBDT_1000tauOS_20260406"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# W 玻色子 PDG 质量 (MeV)
W_MASS_PDG = 80360.0 

# ==============================
# 2. 变量定义
# ==============================
VARIABLES = [
    "W_log_P","W_ETA","W_THETA","W_MINIP", 
    "Mu_log_PT", "Mu_MINIP","Mu_log_P","Mu_ETA","Mu_THETA",
    "MuNuR_E", "MuNuR_log_PT", "NuR_log_PT","MuNuR_log_P","NuR_log_P","MuNuR_MINIP",
]

PRESEL_VARS = ["Mu_PROBNN_MU", "MuNuR_PROBNN_MU", "Mu_PROBNN_GHOST", "MuNuR_PROBNN_GHOST", "W_M"]


EVENT_VARS = ["EVENTNUMBER"] 

LOAD_VARS = list(set(VARIABLES + PRESEL_VARS + EVENT_VARS))

# ==============================
# 3 & 4. 读取与合并 Tree 的封装函数
# ==============================
def load_and_combine_trees(file_path, is_mc=True, max_entries=200000):
    """同时读取 1J 和 2J，合并并打上标签"""
    df_list = []
    
    with uproot.open(file_path) as f:
        # 读取 1J
        if TREE_1J in f:
            df_1j = f[TREE_1J].arrays(LOAD_VARS, entry_stop=max_entries, library="pd")
            df_1j["topology"] = "1J" # 标记来源
            df_list.append(df_1j)
            
        # 读取 2J
        if TREE_2J in f:
            df_2j = f[TREE_2J].arrays(LOAD_VARS, entry_stop=max_entries, library="pd")
            df_2j["topology"] = "2J"
            df_list.append(df_2j)
            
    if not df_list:
        return pd.DataFrame() # 空 df
        
    df_combined = pd.concat(df_list, ignore_index=True)
    df_combined["label"] = 1 if is_mc else 0
    return df_combined

# ==============================
# 5 & 6. 批量读取数据 (1J & 2J)
# ==============================
print("开始读取数据 (包含 1J 和 2J)...")
mc_list = []
for file_path in MC_FILES:
    try:
        df_mc = load_and_combine_trees(file_path, is_mc=True)
        if not df_mc.empty:
            mc_list.append(df_mc)
    except Exception as e:
        print(f"跳过文件 {file_path}: {e}")

mc = pd.concat(mc_list, ignore_index=True)

data = load_and_combine_trees(DATA_FILE, is_mc=False)

# ==============================
# 7. Best Candidate Selection & Pre-selection
# ==============================
dataset = pd.concat([mc, data], ignore_index=True)
dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()
print(f"清理前(包含1J+2J全部候选体)数据量: {len(dataset)}")

# 【核心逻辑 A：跨目录 Best Candidate Selection】
# 1. 计算每个候选体的质量与 W_PDG 的偏差
dataset["W_M_diff"] = abs(dataset["W_M"] - W_MASS_PDG)

# 2. 排序：优先按 runNumber, eventNumber 排在一起，同一个事件内按 W_M_diff 从小到大排
dataset = dataset.sort_values(by=["EVENTNUMBER", "W_M_diff"])

# 3. 去重：对同一个事件，只保留第一行（即 W_M_diff 最小的那个，实现了二选一 / 择优录取）
dataset = dataset.drop_duplicates(subset=["EVENTNUMBER"], keep="first")

print(f"执行 Best Candidate Selection (每个事件只留1个最接近W质量的候选体) 后，剩余数据量: {len(dataset)}")

# 【逻辑 B：普通的通过性 Pre-selection】
mask = (
    (dataset["Mu_PROBNN_MU"] > 0.5) & 
    (dataset["MuNuR_PROBNN_MU"] > 0.5) & 
    (dataset["Mu_PROBNN_GHOST"] < 0.1) & 
    (dataset["MuNuR_PROBNN_GHOST"] < 0.1) & 
    (dataset["W_M"] > 50000.0) & 
    (dataset["W_M"] < 100000.0) 
)
dataset = dataset[mask].copy()

print(f"最后执行 Pre-selection 质量窗口等截断后，最终数据量 = {len(dataset)}")
print(f"  --> 信号(Signal)数 = {(dataset['label'] == 1).sum()}")
print(f"  --> 本底(Background)数 = {(dataset['label'] == 0).sum()}")

# 准备训练数据
X = dataset[VARIABLES]
y = dataset["label"]

# ==============================
# 8. 划分集
# ==============================
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

# ==============================
# 9 & 10. 训练与评估
# ==============================
print("正在训练单 GBDT 模型...")
model = GradientBoostingClassifier(
    n_estimators=100, 
    learning_rate=0.1,
    max_depth=3,
    random_state=42
)
model.fit(X_train, y_train)

y_score_test = model.predict_proba(X_test)[:, 1]
y_score_train = model.predict_proba(X_train)[:, 1]

auc_test = roc_auc_score(y_test, y_score_test)
auc_train = roc_auc_score(y_train, y_score_train)

print(f"训练集 AUC = {auc_train:.6f}")
print(f"测试集 AUC = {auc_test:.6f}")

# ==============================
# 11. 变量重要性
# ==============================

importances = model.feature_importances_
ranking = sorted(
    zip(VARIABLES, importances),
    key=lambda x: x[1],
    reverse=True
)

print("\n区分度较高的变量（按 importance 排序）:")
for var, score in ranking:
    print(f"{var:15s}  {score:.6f}")

# 保存文字结果
with open(os.path.join(OUTPUT_DIR, "feature_ranking.txt"), "w", encoding="utf-8") as f:
    f.write(f"AUC = {auc_test:.6f}\n\n")
    f.write("Feature ranking:\n")
    for var, score in ranking:
        f.write(f"{var:15s}  {score:.6f}\n")

# ==============================
# 12. 画 feature importance 图
# ==============================

plt.figure(figsize=(8, 5))
plt.barh(VARIABLES, importances)
plt.xlabel("Importance")
plt.title("Feature Importance")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "feature_importance.png"), dpi=150)
plt.close()

# ==============================
# 13. 画 ROC 曲线
# ==============================

fpr, tpr, _ = roc_curve(y_test, y_score_test)

plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f"AUC = {auc_test:.4f}")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "roc_curve.png"), dpi=150)
plt.close()


# ==============================
# 13.5 绘制 BDT 分数分布图 (Signal vs Background)
# ==============================

plt.figure(figsize=(8, 6))

# 设置 bin 的范围和数量
bins = np.linspace(0, 1, 51)

# 提取测试集中 信号 和 本底 的 BDT 分数
score_sig = y_score_test[y_test == 1]
score_bkg = y_score_test[y_test == 0]

# 绘制直方图
# alpha 是透明度，density=False 表示显示原始数量
plt.hist(score_sig, bins=bins, alpha=0.5, label='Signal (MC)', color='blue', histtype='stepfilled')
plt.hist(score_bkg, bins=bins, alpha=0.5, label='Background (Data)', color='red', histtype='stepfilled')

# 高能物理中，由于本底和信号数量级可能相差很大，通常建议开启对数纵轴
plt.yscale('log') 

plt.xlabel('BDT Score')
plt.ylabel('Events (Log Scale)')
plt.title('BDT Response Distribution')
plt.legend(loc='upper center')
plt.grid(axis='y', alpha=0.3)

# 保存图片
plt.savefig(os.path.join(OUTPUT_DIR, "bdt_distribution.png"), dpi=150)
print(f"  - bdt_distribution.png (信号/本底分数分布图)")
plt.close()

# 额外建议：如果想看形状对比（归一化到单位面积），可以再画一张
plt.figure(figsize=(8, 6))
plt.hist(score_sig, bins=bins, alpha=0.3, label='Signal (Normalized)', color='blue', density=True, lw=2, histtype='step')
plt.hist(score_bkg, bins=bins, alpha=0.3, label='Background (Normalized)', color='red', density=True, lw=2, histtype='step')
plt.xlabel('BDT Score')
plt.ylabel('Probability Density')
plt.title('Normalized BDT Shape Comparison')
plt.legend(loc='upper center')
plt.savefig(os.path.join(OUTPUT_DIR, "bdt_shape_comparison.png"), dpi=150)
plt.close()

# ==============================
# 14. 输出结束信息
# ==============================

print("\n运行完成。输出文件保存在：", OUTPUT_DIR)
print("包括：")
print("  - feature_importance.png")
print("  - roc_curve.png")
print("  - feature_ranking.txt")