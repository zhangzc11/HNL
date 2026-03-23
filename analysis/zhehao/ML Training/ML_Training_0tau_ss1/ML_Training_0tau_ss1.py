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
    "42912010_5GeVCtau0ps_options_20260228_111712_skim.root",
    "42912012_10GeVCtau0ps_options_20260228_111712_skim.root",
    "42912014_15GeVCtau0ps_options_20260228_111712_skim.root",
    "42912016_20GeVCtau0ps_options_20260228_111712_skim.root",
    "42912018_30GeVCtau0ps_options_20260228_111712_skim.root",
    "42912020_50GeVCtau0ps_options_20260228_230639_skim.root"
]

# 使用列表推导式生成完整路径
MC_FILES = sorted([os.path.join(DATA_DIR, f) for f in file_names])

# ROOT 文件里的树名
TREE_NAME = "myTupleSS1J/DecayTree" # 输入你想分析的branch类别

# 输出目录
OUTPUT_DIR = "./ML_Training_0tau_ss1"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================
# 2. 选择要使用的变量
#    这些变量名需要和 ROOT 文件一致
# ==============================

VARIABLES = [
    #"Mu_BREMTRACKBASEDENERGY",                # 缪子重建基于径迹的辐射能量（有重复）
    #"Mu_ELECTRONCOL",                        # 缪子重建电子簇射的列坐标(相差不是很大？)
    "Mu_ELECTRONENERGY",                      # 缪子重建电子在探测器中沉积的能量 （在0附近堆积严重）
    #"Mu_ELECTRONSHOWERDLL",                   # 缪子重建电子簇射的对数似然差（不是很清楚具体含义，但是图像上确实相差较大）
    #"Mu_ELECTRONSHOWEREOP",                    # 缪子重建电子簇射的EOP参数，在0附近堆积严重（EOP更加像是一个用于判断是否是μ子的量而非判断HNL的谬子有何性质特殊性）
    "Mu_HCALEOP",                             # 缪子在强子量能器（HCAL）的EOP参数
    #"MuNuR_POSITION_STATEAT_LastMeasurement_Z", #munur最后测量点位置Z坐标，缪子径迹最后测量点的束流方向位置（在图像上data表现为孤立的点簇？）
    #"Mu_POSITION_STATEAT_LastMeasurement_Z",  #缪子最后测量点位置Z坐标，缪子径迹最后测量点的束流方向位置（在图像上data表现为孤立的点簇？）
    #"W_MmuWmuN",                               #W玻色子缪子缪子W缪子N质量，用于W玻色子研究的特定质量组合变量，不是很清楚？
    "MuNuR_log_PT",                             #munur的对数横向动量，有明显偏移
    #"NuR_log_MIN_PT",                           #NUR对数横向动量（有重复）
    #"MuNuR_TRACKPT",                            #munur的对数横向动量，有明显偏移，与上文似有重复
    "W_M",
    "MuNuR_QOVERP",                               #在该（缪子中微子）事件中，重建出的缪子径迹的电荷与动量之比（物理含义不是很清晰）
    #"Mu_PROBNN_MU",                              #缪子神经网络缪子概率，神经网络算法判断缪子为缪子的概率值
    #"MuNuR_log_P",                              #munur的对数动量，有明显偏移(有重合)
    "MuNuR_E",                                 #munur能量
    #"MuNuR_PZ",                                 #munur的pz，有分离，但是不确定是否有重复（有强相关的量）
    "Mu_PROBNN_PI",                            #有一个明显的峰的偏差？但是不清楚为何
    #"MuNuR_PX",                                 #munur的px，有一个峰，但是其是否具有实际物理分辨意义存疑
    #"MuNuR_PY",                                 #munur的py，有一个峰，但是其是否具有实际物理分辨意义存疑
    "Mu_log_PT",                               #μ子横向动量的对数，有显著的很好的峰的偏差
    "W_MIN_PT",
    "MuNuR_HCALEOP",                            #缪子中微子重建强子量能器簇射的端点位置，有堆积但分离相对明显
    #"MuNuR_PROBNN_PI",
    #"MuNuR_M",                                  #munur的质量，有一个峰，但是其是否具有实际物理分辨意义存疑
    #"MuNuR_PROBNN_K",                           #重建的概率，但是分离相对较好
    #"Mu_PROBNN_GHOST",                          #缪子神经网络幽灵径迹概率，神经网络算法判断缪子径迹为虚假径迹的概率值（看上去更加像是一个重建过程中用到的量）
    #"MuNuR_PROBNN_GHOST",
    "W_ETA",
    #"W_THETA",                                    #与eta强相关
    #"nTTracks",                                  #T站径迹数量，记录事件中在T站径迹探测器中重建出的径迹总数（不是很明白其含义）
    #"MuNuR_ELECTRONROW",                        #缪子中微子重建电子在探测器中的行位置（差别并不显著。且有堆积）
    #"MuNuR_ELECTRONCOL",                        #缪子中微子重建电子在探测器中的列位置（差别并不显著，且有堆积）
    #"Mu_MINIP",                                 #缪子中微子重建最小冲击参数，缪子中微子重建径迹到主顶点的最小垂直距离，堆积严重
    #"nVeloTracks",                              #顶点定位器径迹数量，记录事件中在顶点定位器中重建出的径迹总数（不是很明白其物理含义）
    "Mu_QOVERP",                                #在该缪子事件中，缪子径迹的电荷与动量之比
    "MuNuR_MINIP",                              #缪子中微子重建最小冲击参数，缪子中微子重建径迹到主顶点的最小垂直距离 （0点附件堆积严重）           #重建的概率，但是分离相对较好
    #"MuNuR_PROBNN_PI"                      #重建的概率，但是分离相对较好
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