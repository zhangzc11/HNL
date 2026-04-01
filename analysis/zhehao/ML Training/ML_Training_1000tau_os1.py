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
# 1. 路径设置（WSL / Linux）
# ==============================

DATA_DIR = "/home/steven/HNL" 

DATA_FILE = os.path.join(
    DATA_DIR,
    "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
)
file_names = [
    "42912031_5GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912033_10GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912035_15GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912037_20GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912039_30GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912041_50GeVCtau1000ps_options_20260228_111712_skim.root"
]

# 使用列表推导式生成完整路径
MC_FILES = sorted([os.path.join(DATA_DIR, f) for f in file_names])

# ROOT 文件里的树名
TREE_NAME = "myTupleOS1J/DecayTree" # 输入你想分析的branch类别

# 输出目录
OUTPUT_DIR = "./ML_Training_1000tau_os1"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# 2. 选择要使用的变量
#    这些变量名需要和 ROOT 文件一致
# ==============================

VARIABLES = [
    "Mu_ELECTRONENERGY", 
    "MuNuR_ELECTRONENERGY",

    "Mu_ELECTRONSHOWERDLL", 
    "MuNuR_ELECTRONSHOWERDLL",

    #"Mu_ELECTRONSHOWEREOP", 
    "MuNuR_ELECTRONSHOWEREOP",



    #"Mu_HCALEOP", "MuNuR_HCALEOP", 



    #"W_MmuWmuN",

    "MuNuR_log_PT",
    "NuR_log_PT",
    "Mu_log_PT", 
    "W_log_PT",


    "MuNuR_log_P",                            
    "NuR_log_P",                           
    "Mu_log_P",                               
    "W_log_P",

    #"W_M",
    #"W_E",

    #"MuNuR_M",
    "MuNuR_E", 

    "W_ETA",
    "W_THETA",
    "Mu_ETA",
    "Mu_THETA",
    
    #"Mu_PROBNN_PI", "Mu_PROBNN_MU", "Mu_PROBNN_GHOST", "Mu_PROBNN_E","Mu_PROBNN_P","Mu_PROBNN_K",

    #"MuNuR_PROBNN_PI","MuNuR_PROBNN_K", "MuNuR_PROBNN_GHOST","MuNuR_PROBNN_E","MuNuR_PROBNN_P","MuNuR_PROBNN_MU",

    #"nTTracks", 
    #"nVeloTracks", 

    "MuNuR_MINIP",          
    "Mu_MINIP",                                 
    "W_MINIP"
]

# ==============================
# 3. 读取 ROOT 的函数
# ==============================

def load_root(file_path, branches, tree_name="DecayTree", max_entries=None):
    """
    从 ROOT 文件中读取指定 branches，返回 pandas.DataFrame
    max_entries=None 表示读全部；否则只读前 max_entries 个事件
    """
    with uproot.open(file_path) as f:
        tree = f[tree_name]
        df = tree.arrays(
            branches,
            entry_stop=max_entries,
            library="pd"
        )
    return df

# ==============================
# 4. 检查文件
# ==============================

if not MC_FILES:
    raise FileNotFoundError(f"没有在 {DATA_DIR} 下找到 MC 文件（429*.root）")

if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(f"没有找到真实数据文件：{DATA_FILE}")

print("找到的 MC 文件数：", len(MC_FILES))
print("真实数据文件：", DATA_FILE)

# ==============================
# 5. 读取 MC（signal）
# ==============================

print("\n开始读取 MC 文件...")

mc_list = []
for file_path in MC_FILES:
    try:
        df_mc = load_root(
            file_path,
            VARIABLES,
            tree_name=TREE_NAME,
            max_entries=200000  # 可改成 None 读取全部
        )
        df_mc["label"] = 1
        df_mc["source_file"] = os.path.basename(file_path)
        mc_list.append(df_mc)
        print(f"读取成功: {os.path.basename(file_path)}  事件数 = {len(df_mc)}")
    except Exception as e:
        print(f"跳过文件: {file_path}")
        print(f"原因: {e}")

if not mc_list:
    raise RuntimeError("所有 MC 文件都读取失败了，请检查 TREE_NAME 和变量名。")

mc = pd.concat(mc_list, ignore_index=True)

# ==============================
# 6. 读取真实数据（background）
# ==============================

print("\n开始读取真实数据文件...")

data = load_root(
    DATA_FILE,
    VARIABLES,
    tree_name=TREE_NAME,
    max_entries=200000  # 可改成 None 读取全部
)
data["label"] = 0
data["source_file"] = os.path.basename(DATA_FILE)

print(f"真实数据事件数 = {len(data)}")

# ==============================
# 7. 合并数据集
# ==============================

dataset = pd.concat([mc, data], ignore_index=True)

print("--- 检查缺失值 ---")
nan_counts = dataset.isna().sum()
print(nan_counts[nan_counts > 0]) # 只打印有缺失值的变量

print("--- 检查无穷值 ---")
inf_counts = (dataset == np.inf).sum() + (dataset == -np.inf).sum()
print(inf_counts[inf_counts > 0])

# 看看删掉之前有多少行
print(f"清理前数据量: {len(dataset)}")
# 去掉缺失值和无穷大

dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()

#pre-selsection
#dataset = dataset[
 #   (dataset["Mu_PT"] > 25000) &          # 砍掉极低动量的软背景
#    (dataset["Mu_PROBNN_MU"] > 0.5) &    # 要求它在物理上至少算个缪子，砍掉纯强子背景
 #   (dataset["Mu_PROBNN_GHOST"] < 0.1) & # 砍掉明显的假径迹
 #   (dataset["W_PT"] < 100000)&             # 系统动量截断
 #   (dataset["NuR_MIN_PT"] > 5000)          # 横向动量截断
#]


print("\n合并后总事件数 =", len(dataset))
print("signal 事件数 =", (dataset["label"] == 1).sum())
print("background 事件数 =", (dataset["label"] == 0).sum())

X = dataset[VARIABLES]
y = dataset["label"]

# ==============================
# 8. 划分训练集 / 测试集
# ==============================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.3,
    random_state=42,
    stratify=y
)

print("\n训练集大小 =", len(X_train))
print("测试集大小 =", len(X_test))

# ==============================
# 9. 训练模型
# ==============================

print("\n开始训练模型...")

model = GradientBoostingClassifier(
    random_state=42
)

model.fit(X_train, y_train)

# ==============================
# 10. 模型评估
# ==============================
print("开始评估模型...")

# 1. 评估测试集 (Test Set) - 这是衡量模型泛化能力的核心指标
y_score_test = model.predict_proba(X_test)[:, 1]
auc_test = roc_auc_score(y_test, y_score_test)

# 2. 评估训练集 (Training Set) - 用于与测试集结果对比，检查是否过拟合
y_score_train = model.predict_proba(X_train)[:, 1]
auc_train = roc_auc_score(y_train, y_score_train)

# 3. 在终端打印两个 AUC 值
print("--- 模型评估结果 ---")
print(f"训练集 AUC (Training AUC) = {auc_train:.6f}")
print(f"测试集 AUC (Test AUC)     = {auc_test:.6f}")

# 注意：后续的 ROC 曲线和特征重要性等，仍然使用测试集的结果，因为这才是对模型泛化能力的真实评估。
# 为了让后续代码能继续运行，我们将测试集的结果赋给原来的变量名。
y_score = y_score_test
auc = auc_test

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
    f.write(f"AUC = {auc:.6f}\n\n")
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

fpr, tpr, _ = roc_curve(y_test, y_score)

plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f"AUC = {auc:.4f}")
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
score_sig = y_score[y_test == 1]
score_bkg = y_score[y_test == 0]

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