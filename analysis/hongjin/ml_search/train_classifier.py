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

DATA_DIR = "/home/data" # 数据路径

DATA_FILE = os.path.join(
    DATA_DIR,
    "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
)

MC_FILES = sorted(glob.glob(os.path.join(DATA_DIR, "429*.root")))

# ROOT 文件里的树名
TREE_NAME = "myTupleSS1J/DecayTree" # 输入你想分析的branch类别

# 输出目录
OUTPUT_DIR = "./ml_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# 2. 选择要使用的变量
#    这些变量名需要和 ROOT 文件一致
# ==============================

VARIABLES = [
    "W_PT",
    "Mu_PT",
    "MuNuR_PT",
    "nPVs",
    "Mu_ETA",
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

# 去掉缺失值和无穷大
dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna()

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

y_score = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_test, y_score)

print(f"\nAUC = {auc:.6f}")

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
# 14. 输出结束信息
# ==============================

print("\n运行完成。输出文件保存在：", OUTPUT_DIR)
print("包括：")
print("  - feature_importance.png")
print("  - roc_curve.png")
print("  - feature_ranking.txt")