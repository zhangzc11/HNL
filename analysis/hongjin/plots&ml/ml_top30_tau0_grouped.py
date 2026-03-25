
import glob
import os
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import KFold

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

# Change here if needed:
TREE_NAME = "myTupleOS1J/DecayTree"   # e.g. "myTupleSS1J/DecayTree"
MAX_ENTRIES_PER_FILE = None           # e.g. 50000 for a fast test
RANDOM_STATE = 42
N_ESTIMATORS = 200
LEARNING_RATE = 0.05
MAX_DEPTH = 3

OUTPUT_DIR = "./top30_ml_output"
ROC_DIR = os.path.join(OUTPUT_DIR, "roc_curves")
SCORE_DIR = os.path.join(OUTPUT_DIR, "score_distributions")
os.makedirs(ROC_DIR, exist_ok=True)
os.makedirs(SCORE_DIR, exist_ok=True)

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

FILE_PATTERN = re.compile(
    r"(?P<id>429\d+?)_(?P<mass>\d+)GeVCtau(?P<tau>\d+)ps_.*\.root$"
)


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
        "sample_id": m.group("id"),
    }


def load_root(file_path, branches, tree_name, max_entries=None):
    with uproot.open(file_path) as f:
        tree = f[tree_name]
        return tree.arrays(branches, entry_stop=max_entries, library="pd")


def sanitize_dataframe(df: pd.DataFrame, columns):
    work = df.copy()
    for c in columns:
        work[c] = pd.to_numeric(work[c], errors="coerce")
    work = work.replace([np.inf, -np.inf], np.nan)
    work = work.dropna(subset=columns)
    return work


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def build_tau0_datasets():
    mc_meta = []
    for path in MC_FILES:
        info = parse_mc_filename(path)
        if info is not None:
            mc_meta.append(info)

    if not mc_meta:
        raise RuntimeError("没有找到可解析的 MC 文件名，请检查路径和命名格式。")

    tau0_infos = sorted(
        [x for x in mc_meta if x["tau"] == 0],
        key=lambda x: (x["mass"], x["name"]),
    )
    if not tau0_infos:
        raise RuntimeError("没有 tau=0ps 的 MC 文件，无法进行 tau=0ps 分类训练。")

    signal_by_mass = {}
    for info in tau0_infos:
        print(f"  读取 tau=0 signal: {info['name']}")
        df = load_root(info["file"], TOP30_VARIABLES, TREE_NAME, MAX_ENTRIES_PER_FILE)
        df = sanitize_dataframe(df, TOP30_VARIABLES)
        df["mass"] = info["mass"]
        df["tau"] = info["tau"]
        df["source_file"] = info["name"]
        df["label"] = 1
        signal_by_mass[info["mass"]] = df

    background_df = load_root(DATA_FILE, TOP30_VARIABLES, TREE_NAME, MAX_ENTRIES_PER_FILE)
    background_df = sanitize_dataframe(background_df, TOP30_VARIABLES)
    background_df["label"] = 0
    background_df["source_file"] = os.path.basename(DATA_FILE)

    return signal_by_mass, background_df


def cross_validated_auc(signal_by_mass, background_df, used_vars, tag):
    """Use leave-one-mass-out for signal + KFold on background.
    This is stricter than random event-level splitting and helps avoid
    unrealistically perfect AUC from same-source leakage.
    """
    masses = sorted(signal_by_mass.keys())
    n_folds = len(masses)
    if n_folds < 2:
        raise RuntimeError("tau=0ps 质量点少于 2 个，无法做按质量分组验证。")

    bg_index = np.arange(len(background_df))
    bg_kfold = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_STATE)
    bg_splits = list(bg_kfold.split(bg_index))

    y_true_all = []
    y_score_all = []
    fold_rows = []

    for fold_id, mass in enumerate(masses):
        sig_test = signal_by_mass[mass]
        sig_train_parts = [signal_by_mass[m] for m in masses if m != mass]
        sig_train = pd.concat(sig_train_parts, ignore_index=True)

        bg_train_idx, bg_test_idx = bg_splits[fold_id]
        bg_train = background_df.iloc[bg_train_idx].copy()
        bg_test = background_df.iloc[bg_test_idx].copy()

        train_df = pd.concat([sig_train, bg_train], ignore_index=True)
        test_df = pd.concat([sig_test, bg_test], ignore_index=True)

        X_train = train_df[used_vars]
        y_train = train_df["label"]
        X_test = test_df[used_vars]
        y_test = test_df["label"]

        model = GradientBoostingClassifier(
            random_state=RANDOM_STATE,
            n_estimators=N_ESTIMATORS,
            learning_rate=LEARNING_RATE,
            max_depth=MAX_DEPTH,
        )
        model.fit(X_train, y_train)
        y_score = model.predict_proba(X_test)[:, 1]

        fold_auc = roc_auc_score(y_test, y_score)
        y_true_all.append(np.asarray(y_test))
        y_score_all.append(np.asarray(y_score))
        fold_rows.append(
            {
                "fold": fold_id + 1,
                "heldout_signal_mass_GeV": mass,
                "n_signal_test": int((y_test == 1).sum()),
                "n_background_test": int((y_test == 0).sum()),
                "auc": float(fold_auc),
            }
        )

    y_true_all = np.concatenate(y_true_all)
    y_score_all = np.concatenate(y_score_all)

    auc = roc_auc_score(y_true_all, y_score_all)
    fpr, tpr, _ = roc_curve(y_true_all, y_score_all)

    # ROC
    plt.figure(figsize=(6.5, 6.5))
    plt.plot(fpr, tpr, linewidth=2, label=f"AUC = {auc:.12f}")
    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve ({tag})\nSignal: leave-one-mass-out, Background: KFold")
    plt.legend()
    plt.tight_layout()
    roc_png = os.path.join(ROC_DIR, f"{safe_filename(tag)}.png")
    plt.savefig(roc_png, dpi=160)
    plt.close()

    # score distribution
    plt.figure(figsize=(7.2, 5.8))
    sig_scores = y_score_all[y_true_all == 1]
    bg_scores = y_score_all[y_true_all == 0]
    plt.hist(sig_scores, bins=50, density=True, histtype="step", linewidth=2, label="Signal test score")
    plt.hist(bg_scores, bins=50, density=True, histtype="stepfilled", alpha=0.35, label="Background test score")
    plt.xlabel("GBDT score")
    plt.ylabel("Normalized events")
    plt.title(f"Score distribution ({tag})")
    plt.legend()
    plt.tight_layout()
    score_png = os.path.join(SCORE_DIR, f"{safe_filename(tag)}.png")
    plt.savefig(score_png, dpi=160)
    plt.close()

    fold_df = pd.DataFrame(fold_rows)
    fold_csv = os.path.join(OUTPUT_DIR, f"{safe_filename(tag)}_fold_auc.csv")
    fold_df.to_csv(fold_csv, index=False, encoding="utf-8-sig")

    return {
        "auc": float(auc),
        "roc_curve_file": os.path.basename(roc_png),
        "score_plot_file": os.path.basename(score_png),
        "fold_auc_file": os.path.basename(fold_csv),
        "n_test_signal_total": int((y_true_all == 1).sum()),
        "n_test_background_total": int((y_true_all == 0).sum()),
        "auc_mean_over_folds": float(fold_df["auc"].mean()),
        "auc_std_over_folds": float(fold_df["auc"].std(ddof=0)),
    }


def main():
    if not os.path.exists(DATA_FILE):
        raise FileNotFoundError(f"没有找到本底数据文件: {DATA_FILE}")

    print("ROOT 树名:", TREE_NAME)
    print("输出目录:", OUTPUT_DIR)
    print("\n[1/2] 构建 tau=0ps 机器学习数据集 ...")
    signal_by_mass, background_df = build_tau0_datasets()

    print("  tau=0ps signal 质量点 =", sorted(signal_by_mass.keys()))
    print("  tau=0ps signal 总事件数 =", sum(len(df) for df in signal_by_mass.values()))
    print("  background 事件数      =", len(background_df))

    print("\n[2/2] 训练 29 个 GBDT 模型并输出 ROC/AUC ...")
    experiments = [("all_30_variables", TOP30_VARIABLES.copy(), None)]
    for var in TOP30_VARIABLES:
        experiments.append((f"drop_{var}", [v for v in TOP30_VARIABLES if v != var], var))

    results = []
    for idx, (tag, used_vars, removed_var) in enumerate(experiments, start=1):
        metrics = cross_validated_auc(signal_by_mass, background_df, used_vars, tag)

        results.append(
            {
                "experiment_index": idx,
                "tag": tag,
                "removed_variable": "NONE" if removed_var is None else removed_var,
                "n_variables_used": len(used_vars),
                "auc": metrics["auc"],
                "auc_mean_over_folds": metrics["auc_mean_over_folds"],
                "auc_std_over_folds": metrics["auc_std_over_folds"],
                "n_test_signal_total": metrics["n_test_signal_total"],
                "n_test_background_total": metrics["n_test_background_total"],
                "roc_curve_file": metrics["roc_curve_file"],
                "score_plot_file": metrics["score_plot_file"],
                "fold_auc_file": metrics["fold_auc_file"],
            }
        )
        print(
            f"  [{idx:02d}/29] {tag:45s} "
            f"AUC = {metrics['auc']:.12f} "
            f"(fold mean = {metrics['auc_mean_over_folds']:.12f}, std = {metrics['auc_std_over_folds']:.12f})"
        )

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="auc", ascending=False).reset_index(drop=True)
    results_df.to_csv(
        os.path.join(OUTPUT_DIR, "auc_summary_sorted.csv"),
        index=False,
        encoding="utf-8-sig",
    )

    with open(os.path.join(OUTPUT_DIR, "auc_summary_sorted.txt"), "w", encoding="utf-8") as f:
        f.write("tau=0ps GBDT AUC 排序结果\n")
        f.write("评估方式：signal 按质量点留一法，background 使用 KFold\n")
        f.write("=" * 100 + "\n")
        for _, row in results_df.iterrows():
            f.write(
                f"AUC={row['auc']:.12f} | fold_mean={row['auc_mean_over_folds']:.12f} "
                f"| fold_std={row['auc_std_over_folds']:.12f} | used={int(row['n_variables_used']):2d} | "
                f"removed={row['removed_variable']} | roc={row['roc_curve_file']} | "
                f"score={row['score_plot_file']} | fold_csv={row['fold_auc_file']}\n"
            )

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
        f.write("- roc_curves/: 共29张ROC曲线图\n")
        f.write("- score_distributions/: 共29张测试集打分分布图\n")
        f.write("- 每个实验各有一个 *_fold_auc.csv，记录6个质量点留一法的AUC\n")
        f.write("- auc_summary_sorted.csv/txt: 29个实验按AUC排序\n")
        f.write("- leave_one_out_importance_by_auc_drop.csv: 去掉某变量后AUC相对all_30_variables的下降量\n")
        f.write("- 这版不再使用随机事件级 train_test_split，而使用更严格的按质量点分组评估\n")
        f.write(f"- baseline all-28 AUC = {baseline_auc:.12f}\n")

    print("\n完成。")
    print("输出文件夹:", OUTPUT_DIR)
    print("- ROC 曲线图:", ROC_DIR)
    print("- 打分分布图:", SCORE_DIR)
    print("- 29个实验AUC排序: auc_summary_sorted.csv / .txt")
    print("- 去掉单变量后的AUC下降量: leave_one_out_importance_by_auc_drop.csv")


if __name__ == "__main__":
    main()
