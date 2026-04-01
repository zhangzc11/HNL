import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

data_path = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_base_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"
masses = [5, 10, 15, 20, 30, 50]
lifetimes = [0, 10, 100, 1000]

def find_mc_file(mass, lifetime):
    pattern = f"*_{mass}GeVCtau{lifetime}ps_options_*_skim.root"
    matches = glob.glob(os.path.join(mc_base_dir, pattern))
    return matches[0] if matches else None

# 【新增：智能推断物理单位】
def get_unit(var_name):
    if any(var_name.endswith(s) for s in["_PT", "_P", "_E", "_M", "_PX", "_PY", "_PZ", "MmuWmuN"]):
        return "[MeV]"
    elif any(var_name.endswith(s) for s in ["_VX", "_VY", "_VZ"]):
        return "[mm]"
    elif "LTIME" in var_name:
        return "[ps]"
    elif "ETA" in var_name or "PHI" in var_name or "THETA" in var_name:
        return "[rad]"
    else:
        return "[a.u.]"

print("="*60)
# 【新增：交互式选择物理区域】
print("请选择要分析的物理区域 (Channel Topology):")
print("  [1] OS1J (异号单喷流 - 适合轻质量)")
print("  [2] OS2J (异号双喷流 - 适合高质量)")
print("  [3] SS1J (同号单喷流 - 极低本底)")
print("  [4] SS2J (同号双喷流 - 极低本底)")
region_map = {"1": "OS1J", "2": "OS2J", "3": "SS1J", "4": "SS2J"}
reg_choice = input("👉 请选择 1/2/3/4: ").strip()
region_name = region_map.get(reg_choice, "OS1J")
tree_name = f"myTuple{region_name}/DecayTree"

var_name = input("\n👉 请输入你想绘制的变量名 (例如 W_MmuWmuN): ").strip()
unit = get_unit(var_name)

mode = input("\n👉 请选择模式 [1]扫描质量 [2]扫描寿命: ").strip()

output_dir = f"Plots_{region_name}_{var_name}"
if not os.path.exists(output_dir): os.makedirs(output_dir)

print(f"\n[*] 正在读取 {region_name} 区域的数据...")
with uproot.open(data_path) as f:
    df_data = f[tree_name].arrays([var_name], library="pd")

def draw_and_save(df_bkg, dict_sig, title_prefix, out_prefix):
    if not dict_sig: return
    all_vals = [df[var_name].values for df in dict_sig.values() if len(df) > 0]
    if len(df_bkg) > 0: all_vals.append(df_bkg[var_name].values)
    
    global_min = min([np.percentile(v, 1) for v in all_vals if len(v) > 0])
    global_max = max([np.percentile(v, 99) for v in all_vals if len(v) > 0])
    bins = np.linspace(global_min, global_max, 50)
    colors = plt.cm.jet(np.linspace(0.1, 0.9, len(dict_sig)))

    for use_log in [False, True]:
        plt.figure(figsize=(10, 7))
        plt.hist(df_bkg[var_name], bins=bins, density=True, alpha=0.3, color='black', label='Background (Data)')
        for (label, df), color in zip(dict_sig.items(), colors):
            plt.hist(df[var_name], bins=bins, density=True, histtype='step', lw=2.5, color=color, label=f'Signal ({label})')

        scale_type = "Log" if use_log else "Linear"
        if use_log: plt.yscale('log')

        # 【新增：标题带上区域，X轴带上单位】
        plt.title(f"[{region_name}] {title_prefix} - {var_name} ({scale_type})", fontsize=15, fontweight='bold')
        plt.xlabel(f"{var_name} {unit}", fontsize=13)
        plt.ylabel('Normalized Events', fontsize=13)
        
        # 【新增：在图内添加 LHCb 文本框】
        plt.text(0.05, 0.95, f"LHCb Internal\nRegion: {region_name}", transform=plt.gca().transAxes, 
                 fontsize=14, fontweight='bold', va='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        plt.legend(fontsize=11, loc='upper right')
        plt.grid(True, alpha=0.3, ls="--")
        plt.savefig(f"{output_dir}/{out_prefix}_{scale_type}.png", dpi=300)
        plt.close()
    print(f"[+] 成功生成图像: {out_prefix}")

# (下方保持原来的模式 1 和模式 2 收集逻辑不变...)
if mode == '1':
    for tau in lifetimes:
        sig_dict = {}
        for m in masses:
            fpath = find_mc_file(m, tau)
            if fpath:
                try:
                    with uproot.open(fpath) as f_mc:
                        sig_dict[f"{m} GeV"] = f_mc[tree_name].arrays([var_name], library="pd")
                except: pass
        draw_and_save(df_data, sig_dict, f"Fixed Lifetime {tau}ps", f"MassScan_{var_name}_Tau{tau}ps")
