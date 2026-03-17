import pandas as pd
import uproot
import matplotlib.pyplot as plt
from utils import *
import sys

ks = pd.read_csv('./csv/branch_ks.csv')
ks_sorted = ks.sort_values(by='KS', ascending=False)
ks_sorted.to_csv('./csv/ks_sorted.csv', index=False)
andeson_sorted = ks.sort_values(by='Andeson', ascending=False)
andeson_sorted.to_csv('./csv/andeson_sorted.csv', index=False)
weight_sorted = ks.sort_values(by='Weight', ascending=False)
weight_sorted.to_csv('./csv/weight_sorted.csv', index=False)

mc_file = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"
data_file = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"

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

def load_file(file_path, tree_name):
    file = uproot.open(file_path)
    tree = file[tree_name]
    physics_keys = ks_sorted['Branch'].values
    # 提取物理意义的数据
    physics_data = tree.arrays(physics_keys, library="np")
    return physics_data

mc = load_file(mc_file, tree_name_mc)
data = load_file(data_file, tree_name_data)

# 前20个绘图
branches = weight_sorted['Branch'].values[:20]
for i in range(len(branches)):
    branch = branches[i]
    mc_a = mc[branch]
    data_a = data[branch]
    plt.hist(mc_a, bins=50, alpha=0.5, label='MC', density=True)
    plt.hist(data_a, bins=50, alpha=0.5, label='Data', density=True)
    plt.title( f'{directory}_{branch}')
    plt.legend()
    plt.xlabel(branch)
    plt.text(0.05, 0.95, 
             f'KS = {weight_sorted[weight_sorted["Branch"] == branch]["KS"].values[0]:.4f}, Andersen = {weight_sorted[weight_sorted["Branch"] == branch]["Andeson"].values[0]:.0f}',
               transform=plt.gca().transAxes)
    plt.savefig(f'./png/{str(i+1)}_{directory}_{branch}.png')
    plt.close()
