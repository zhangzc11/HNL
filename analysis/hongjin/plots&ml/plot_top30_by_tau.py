
import glob
import os
import re
import warnings
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

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
N_BINS = 60

OUTPUT_DIR = "./top30_plot_output"
PLOT_DIR = os.path.join(OUTPUT_DIR, "variable_plots")
os.makedirs(PLOT_DIR, exist_ok=True)

TOP30_VARIABLES = [
    "W_M",
    "Mu_POSITION_STATEAT_LastMeasurement_Z",
    "W_log_MIN_PT",
    "W_MIN_PT",
    "MuNuR_POSITION_STATEAT_LastMeasurement_Z",
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

MASS_ORDER = [5, 10, 15, 20, 30, 50]
MASS_COLORS = {
    5: "red",
    10: "blue",
    15: "green",
    20: "magenta",
    30: "c",
    50: "orange",
}

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


def sanitize_series(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    return s


def union_range(arrays):
    valid = []
    for arr in arrays:
        if arr is None or len(arr) == 0:
            continue
        x = np.asarray(arr)
        x = x[np.isfinite(x)]
        if x.size:
            valid.append(x)
    if not valid:
        return None
    lo = min(np.min(x) for x in valid)
    hi = max(np.max(x) for x in valid)
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    if lo == hi:
        pad = 0.5 if lo == 0 else abs(lo) * 0.05
        lo -= pad
        hi += pad
    return (lo, hi)


def normalized_weights(n):
    if n <= 0:
        return None
    return np.ones(n, dtype=float) / float(n)


def safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main():
    mc_meta = []
    for path in MC_FILES:
        info = parse_mc_filename(path)
        if info is not None:
            mc_meta.append(info)

    if not mc_meta:
        raise RuntimeError("没有找到可解析的 MC 文件名，请检查路径和命名格式。")

    if not os.path.exists(DATA_FILE):
        raise FileNotFoundError(f"没有找到本底数据文件: {DATA_FILE}")

    mc_by_tau = defaultdict(list)
    for info in mc_meta:
        mc_by_tau[info["tau"]].append(info)
    for tau in mc_by_tau:
        mc_by_tau[tau] = sorted(mc_by_tau[tau], key=lambda x: (x["mass"], x["name"]))

    tau_values = sorted(mc_by_tau.keys())

    print("找到的 tau 类别:", tau_values)
    print("ROOT 树名:", TREE_NAME)
    print("输出目录:", OUTPUT_DIR)
    print("\n开始画变量分布图 ...")

    for variable in TOP30_VARIABLES:
        print(f"  变量: {variable}")

        try:
            bg_df = load_root(DATA_FILE, [variable], TREE_NAME, MAX_ENTRIES_PER_FILE)
            bg_values = sanitize_series(bg_df[variable]).to_numpy()
        except Exception as e:
            print(f"    跳过变量 {variable}：读取本底失败 -> {e}")
            continue

        for tau in tau_values:
            mass_arrays = {}
            for info in mc_by_tau[tau]:
                try:
                    df = load_root(info["file"], [variable], TREE_NAME, MAX_ENTRIES_PER_FILE)
                    vals = sanitize_series(df[variable]).to_numpy()
                    mass_arrays[info["mass"]] = vals
                except Exception as e:
                    print(f"    读取失败 {info['name']} / {variable}: {e}")

            xrange = union_range(list(mass_arrays.values()) + [bg_values])
            if xrange is None:
                print(f"    tau={tau}ps, 变量 {variable} 无有效数据，跳过绘图")
                continue

            plt.figure(figsize=(10, 7))

            for mass in MASS_ORDER:
                if mass not in mass_arrays:
                    continue
                vals = mass_arrays[mass]
                if vals.size == 0:
                    continue
                plt.hist(
                    vals,
                    bins=N_BINS,
                    range=xrange,
                    weights=normalized_weights(len(vals)),
                    histtype="step",
                    linewidth=1.8,
                    label=f"{mass}GeV",
                    color=MASS_COLORS.get(mass, None),
                )

            if bg_values.size > 0:
                plt.hist(
                    bg_values,
                    bins=N_BINS,
                    range=xrange,
                    weights=normalized_weights(len(bg_values)),
                    histtype="stepfilled",
                    alpha=0.35,
                    edgecolor="black",
                    linewidth=1.0,
                    label="Data",
                )

            plt.text(
                0.06,
                0.93,
                f"Fixed Tau = {tau}ps",
                transform=plt.gca().transAxes,
                fontsize=16,
                fontweight="bold",
                va="top",
            )
            plt.text(
                0.06,
                0.88,
                f"Dir: {TREE_NAME.split('/')[0]}",
                transform=plt.gca().transAxes,
                fontsize=15,
                fontweight="bold",
                va="top",
            )
            plt.xlabel(variable)
            plt.ylabel("Normalized Events")
            plt.legend(frameon=False, fontsize=10)
            plt.tight_layout()
            out_path = os.path.join(PLOT_DIR, f"tau{tau}ps_{safe_filename(variable)}.png")
            plt.savefig(out_path, dpi=160)
            plt.close()

    with open(os.path.join(OUTPUT_DIR, "README_plots.txt"), "w", encoding="utf-8") as f:
        f.write("输出说明\n")
        f.write("- variable_plots/: 每个变量、每个tau一张分布图\n")
        f.write(f"- ROOT tree: {TREE_NAME}\n")
        f.write(f"- variables: {len(TOP30_VARIABLES)}\n")

    print("\n完成。")
    print("输出文件夹:", OUTPUT_DIR)
    print("- 变量分布图:", PLOT_DIR)


if __name__ == "__main__":
    main()
