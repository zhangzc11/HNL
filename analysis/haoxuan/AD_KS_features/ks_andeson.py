import uproot
import os
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from utils import *

mc_file = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"
data_file = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"

# directory = "myTupleOS2J"
# tree_name_data  = directory  + "/DecayTree"
# tree_name_mc = directory + "MC/DecayTree"

args = sys.argv

directory = "myTupleOS2J"
if len(args) > 1:
    directory = args[1]
if len(args) > 2:
    mc_file = args[2]
if len(args) > 3:
    data_file = args[3]

tree_name_data = f"{directory}/DecayTree"
tree_name_mc = f"{directory}/DecayTree"
mc = load_file(mc_file, tree_name_mc)
data = load_file(data_file, tree_name_data)

all_keys = mc.keys()

i = 1
ks = []
branck = []
andeson = []
for key in all_keys:
    if key not in mc or key not in data:
        # print(f"Skipping {key}: Not found in both datasets.")
        continue

    mc_a = mc[key]
    data_a = data[key]
    mc_a = np.array(mc_a).ravel() 
    data_a = np.array(data_a).ravel()

    try:
        mc_a = np.array(mc_a, dtype=float)
        data_a = np.array(data_a, dtype=float)
    except ValueError as e:
        # print(f"Skipping {key}: Error converting to float - {e}")
        continue

    mc_a = mc_a[~np.isnan(mc_a)]
    data_a = data_a[~np.isnan(data_a)]
    mc_arr = np.asarray(mc_a).ravel()
    data_arr = np.asarray(data_a).ravel()

    if mc_arr.ndim != 1 or data_arr.ndim != 1:
        # print(f"Skipping {key}: Data is not 1D vector (MC shape: {mc_arr.shape}, Data shape: {data_arr.shape}).")
        continue
    if mc_arr.size == 0 or data_arr.size == 0:
        # print(f"Skipping {key}: Empty array encountered.")
        continue
    if not (np.issubdtype(mc_arr.dtype, np.number) and np.issubdtype(data_arr.dtype, np.number)):
        # print(f"Skipping {key}: Non-numeric data encountered.")
        continue

    # 使用展平后的数组进行 KS 检验
    ks_stat, p_value = stats.ks_2samp(mc_arr, data_arr)
    para = 0
    try:
        andeson_result = stats.anderson_ksamp([mc_arr, data_arr])
        para = andeson_result[0]
    except ValueError as e:
        continue
    
    print(f"|{i:03d} | KS = {ks_stat:.4f} | Andeson = {para: .4f} |Branch: {key}")
    ks.append(ks_stat)
    branck.append(key)
    andeson.append(para)
    i += 1

ks_m = max(ks)
and_max = max(andeson)
weight = []

for i in range(len(ks)):
    weight.append(0.5*ks[i]/ks_m + 0.5*andeson[i]/and_max)

csv_file = "./csv/branch_ks.csv"
pd.DataFrame({"Branch": branck, "KS": ks, "Andeson": andeson, "Weight": weight}).to_csv(csv_file, index=False)
print(f"KS values saved to {csv_file}")
