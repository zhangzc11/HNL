import glob
import os
import re
import warnings
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ============================================================
# 1. User settings
# ============================================================
DATA_DIR = "/home/data"
DATA_FILE = os.path.join(
    DATA_DIR,
    "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root",
)
MC_FILES = sorted(glob.glob(os.path.join(DATA_DIR, "429*.root")))

# 改成你想分析的树
TREE_NAME = "myTupleOS1J/DecayTree"

# 变量列表：TOP28
TOP30_VARIABLES = [
    "W_M",
    "W_log_MIN_PT",
    "W_MIN_PT",
    "Mu_TRACKPT",
    "Mu_PT",
    "Mu_log_PT",
    "W_log_SUM_PT",
    "W_SUM_PT",
    "W_MmuWmuN",
    "NuR_ALV",
    "W_log_MAX_PT",
    "W_MAX_PT",
    "Mu_HCALEOP",
    "Mu_ELECTRONSHOWEREOP",
    "Mu_PROBNN_GHOST",
    "MuNuR_log_PT",
    "MuNuR_TRACKPT",
    "MuNuR_PT",
    "NuR_MIN_PT",
    "NuR_log_MIN_PT",
    "NuR_log_PT",
    "NuR_PT",
    "Mu_PROBNN_MU",
    "nFTClusters",
    "Mu_BREMTRACKBASEDENERGY",
    "NuR_M",
    "MuNuR_PROBNN_MU",
    "MuNuR_PROBNN_GHOST",
]

# 只做 tau=0ps 的机器学习
TARGET_TAU_PS = 0

# 读取上限。先调小快速试跑，确认无误后可改成 None
MAX_ENTRIES_PER_FILE = None

# 为了避免程序“像卡住一样”，这里用更轻的 GBDT 参数
GBDT_PARAMS = dict(
    n_estimators=60,
    learning_rate=0.08,
    max_depth=2,
    min_samples_leaf=100,
    subsample=0.7,
    random_state=42,
)

# 为了减少训练时间，可在每个 fold 中把 background 采样到不超过 signal 的若干倍
BACKGROUND_TO_SIGNAL_RATIO_CAP = 1.5

OUTPUT_DIR = "./ml_tau0_grouped_fast_output"
ROC_DIR = os.path.join(OUTPUT_DIR, "roc_curves")
SCORE_DIR = os.path.join(OUTPUT_DIR, "score_distributions")
FOLD_DIR = os.path.join(OUTPUT_DIR, "fold_tables")
os.makedirs(ROC_DIR, exist_ok=True)
os.makedirs(SCORE_DIR, exist_ok=True)
os.makedirs(FOLD_DIR, exist_ok=True)

# ============================================================
# 2. Helpers
# ============================================================
FILE_PATTERN = re.compile(r"(?P<id>429\d+?)_(?P<mass>\d+)GeVCtau(?P<tau>\d+)ps_.*\.root$")

def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)

def parse_mc_filename(path: str):
    name = os.path.basename(path)
    m = FILE_PATTERN.match(name)
    if not m:
        return None
    return {
        "file": path,
        "name": name,
        "mass": int(m.group("mass")),
        "tau": int(m.group("tau")),
    }

def load_root(file_path, branches, tree_name, max_entries=None):
    with uproot.open(file_path) as f:
        tree = f[tree_name]
        return tree.arrays(branches, entry_stop=max_entries, library="pd")

def clean_df(df: pd.DataFrame, vars_to_use):
    out = df.replace([np.inf, -np.inf], np.nan).dropna(subset=vars_to_use).copy()
    return out

def normalized_weights(n):
    if n <= 0:
        return None
    return np.ones(n, dtype=float) / float(n)

# ============================================================
# 3. Discover files
# ============================================================
if not MC_FILES:
    raise FileNotFoundError(f"没有在 {DATA_DIR} 下找到 MC 文件（429*.root）")
if not os.path.exists(DATA_FILE):
    raise FileNotFoundError(f"没有找到本底数据文件：{DATA_FILE}")

mc_meta = []
for path in MC_FILES:
    info = parse_mc_filename(path)
    if info is not None:
        mc_meta.append(info)

mc_tau0 = sorted([x for x in mc_meta if x["tau"] == TARGET_TAU_PS], key=lambda z: z["mass"])
if not mc_tau0:
    raise RuntimeError(f"没有找到 tau={TARGET_TAU_PS}ps 的 MC 文件。")

# ============================================================
# 4. Read tau=0 signal and background once
# ============================================================
print("[1/2] 构建 tau=0ps 机器学习数据集 ...")

signal_parts = []
for info in mc_tau0:
    print(f"  读取 tau=0 signal: {info['name']}")
    df = load_root(info["file"], TOP30_VARIABLES, TREE_NAME, MAX_ENTRIES_PER_FILE)
    df["mass"] = info["mass"]
    df["label"] = 1
    df["source_file"] = info["name"]
    signal_parts.append(df)

signal_df = pd.concat(signal_parts, ignore_index=True)
signal_df = clean_df(signal_df, TOP30_VARIABLES)

background_df = load_root(DATA_FILE, TOP30_VARIABLES, TREE_NAME, MAX_ENTRIES_PER_FILE)
background_df["label"] = 0
background_df["source_file"] = os.path.basename(DATA_FILE)
background_df = clean_df(background_df, TOP30_VARIABLES)

masses = sorted(signal_df["mass"].unique().tolist())
print("tau=0ps signal 质量点 =", masses)
print("tau=0ps signal 总事件数 =", len(signal_df))
print("background 事件数     =", len(background_df))
print("ROOT 树 =", TREE_NAME)

# 为 background 预先切成 len(masses) 份，后面和 leave-one-mass-out 对应
rng = np.random.RandomState(42)
bkg_perm = rng.permutation(len(background_df))
bkg_splits = np.array_split(bkg_perm, len(masses))

# ============================================================
# 5. Grouped leave-one-mass-out ML
# ============================================================
print("\n[2/2] 训练 29 个 GBDT 模型并输出 ROC/AUC ...")
experiments = [("all_30_variables", TOP30_VARIABLES.copy(), None)]
for var in TOP30_VARIABLES:
    experiments.append((f"drop_{var}", [v for v in TOP30_VARIABLES if v != var], var))

results = []

for exp_i, (tag, used_vars, removed_var) in enumerate(experiments, start=1):
    print(f"\n[{exp_i:02d}/31] {tag}")
    all_y_true = []
    all_y_score = []
    fold_rows = []

    for fold_i, heldout_mass in enumerate(masses, start=1):
        sig_test = signal_df.loc[signal_df["mass"] == heldout_mass, used_vars + ["mass"]].copy()
        sig_train = signal_df.loc[signal_df["mass"] != heldout_mass, used_vars + ["mass"]].copy()

        bkg_test = background_df.iloc[bkg_splits[fold_i - 1]][used_vars].copy()
        bkg_train = background_df.drop(background_df.index[bkg_splits[fold_i - 1]])[used_vars].copy()

        # 限制 background 训练样本量，避免太慢
        max_bkg_train = int(BACKGROUND_TO_SIGNAL_RATIO_CAP * len(sig_train))
        if len(bkg_train) > max_bkg_train and max_bkg_train > 0:
            choose_idx = rng.choice(len(bkg_train), size=max_bkg_train, replace=False)
            bkg_train = bkg_train.iloc[choose_idx].copy()

        max_bkg_test = int(BACKGROUND_TO_SIGNAL_RATIO_CAP * len(sig_test))
        if len(bkg_test) > max_bkg_test and max_bkg_test > 0:
            choose_idx = rng.choice(len(bkg_test), size=max_bkg_test, replace=False)
            bkg_test = bkg_test.iloc[choose_idx].copy()

        X_train = pd.concat([sig_train[used_vars], bkg_train[used_vars]], ignore_index=True)
        y_train = np.concatenate([
            np.ones(len(sig_train), dtype=int),
            np.zeros(len(bkg_train), dtype=int),
        ])

        X_test = pd.concat([sig_test[used_vars], bkg_test[used_vars]], ignore_index=True)
        y_test = np.concatenate([
            np.ones(len(sig_test), dtype=int),
            np.zeros(len(bkg_test), dtype=int),
        ])

        model = GradientBoostingClassifier(**GBDT_PARAMS)
        model.fit(X_train, y_train)
        y_score = model.predict_proba(X_test)[:, 1]
        fold_auc = roc_auc_score(y_test, y_score)

        print(
            f"    fold {fold_i}/{len(masses)} | test mass = {heldout_mass:>2} GeV | "
            f"Ntrain = {len(X_train):>6} | Ntest = {len(X_test):>5} | AUC = {fold_auc:.8f}"
        )

        all_y_true.append(y_test)
        all_y_score.append(y_score)
        fold_rows.append({
            "fold": fold_i,
            "heldout_mass_GeV": heldout_mass,
            "n_train_signal": int(len(sig_train)),
            "n_train_background": int(len(bkg_train)),
            "n_test_signal": int(len(sig_test)),
            "n_test_background": int(len(bkg_test)),
            "fold_auc": float(fold_auc),
        })

    y_true_all = np.concatenate(all_y_true)
    y_score_all = np.concatenate(all_y_score)

    auc = roc_auc_score(y_true_all, y_score_all)
    fpr, tpr, _ = roc_curve(y_true_all, y_score_all)

    print(f"    ==> overall grouped AUC = {auc:.16f}")

    # ROC
    title = "Grouped ROC (all 30 variables)" if removed_var is None else f"Grouped ROC (drop {removed_var})"
    plt.figure(figsize=(6.5, 6.5))
    plt.plot(fpr, tpr, linewidth=2, label=f"AUC = {auc:.16f}")
    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    roc_path = os.path.join(ROC_DIR, f"{exp_i:02d}_{safe_filename(tag)}.png")
    plt.savefig(roc_path, dpi=160)
    plt.close()

    # score distributions
    sig_scores = y_score_all[y_true_all == 1]
    bkg_scores = y_score_all[y_true_all == 0]
    plt.figure(figsize=(7, 5.5))
    plt.hist(
        sig_scores,
        bins=50,
        range=(0, 1),
        weights=normalized_weights(len(sig_scores)),
        histtype="step",
        linewidth=1.8,
        label="Signal",
    )
    plt.hist(
        bkg_scores,
        bins=50,
        range=(0, 1),
        weights=normalized_weights(len(bkg_scores)),
        histtype="stepfilled",
        alpha=0.35,
        edgecolor="black",
        linewidth=1.0,
        label="Background",
    )
    plt.xlabel("GBDT score")
    plt.ylabel("Normalized Events")
    plt.title(title.replace("ROC", "Score distribution"))
    plt.legend(frameon=False)
    plt.tight_layout()
    score_path = os.path.join(SCORE_DIR, f"{exp_i:02d}_{safe_filename(tag)}_score.png")
    plt.savefig(score_path, dpi=160)
    plt.close()

    fold_df = pd.DataFrame(fold_rows)
    fold_csv = os.path.join(FOLD_DIR, f"{exp_i:02d}_{safe_filename(tag)}_fold_auc.csv")
    fold_df.to_csv(fold_csv, index=False, encoding="utf-8-sig")

    results.append({
        "tag": tag,
        "removed_variable": "NONE" if removed_var is None else removed_var,
        "n_variables_used": len(used_vars),
        "auc": float(auc),
        "roc_curve_file": os.path.basename(roc_path),
        "score_file": os.path.basename(score_path),
        "fold_auc_file": os.path.basename(fold_csv),
    })

# ============================================================
# 6. Save summary
# ============================================================
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(by="auc", ascending=False).reset_index(drop=True)
results_df.to_csv(os.path.join(OUTPUT_DIR, "auc_summary_sorted.csv"), index=False, encoding="utf-8-sig")

baseline_auc = float(results_df.loc[results_df["removed_variable"] == "NONE", "auc"].iloc[0])
leave_one_out_df = results_df[results_df["removed_variable"] != "NONE"].copy()
leave_one_out_df["auc_drop_vs_all30"] = baseline_auc - leave_one_out_df["auc"]
leave_one_out_df = leave_one_out_df.sort_values(by="auc_drop_vs_all30", ascending=False)
leave_one_out_df.to_csv(
    os.path.join(OUTPUT_DIR, "leave_one_out_importance_by_auc_drop.csv"),
    index=False,
    encoding="utf-8-sig",
)

with open(os.path.join(OUTPUT_DIR, "README_results.txt"), "w", encoding="utf-8") as f:
    f.write("输出说明\n")
    f.write("- roc_curves/: 29 张 grouped ROC 图\n")
    f.write("- score_distributions/: 29 张 score 分布图\n")
    f.write("- fold_tables/: 每个实验的 fold-by-fold AUC\n")
    f.write("- auc_summary_sorted.csv: 29 个实验按 AUC 排序\n")
    f.write("- leave_one_out_importance_by_auc_drop.csv: 去掉某变量后 AUC 的下降量\n")
    f.write(f"- baseline all-28 grouped AUC = {baseline_auc:.12f}\n")

print("\n完成。")
print("输出文件夹:", OUTPUT_DIR)
print("- ROC 曲线图:", ROC_DIR)
print("- score 分布图:", SCORE_DIR)
print("- fold AUC 表:", FOLD_DIR)
print("- AUC 汇总: auc_summary_sorted.csv")
print("- 去掉单变量后的 AUC 下降量: leave_one_out_importance_by_auc_drop.csv")
