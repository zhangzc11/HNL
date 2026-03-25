import uproot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# ==========================================
# 0. 基础路径和参数设置
# ==========================================
data_path = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_base_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"
tree_name = "myTupleOS1J/DecayTree"

masses =[5, 10, 15, 20, 30, 50]
lifetimes =[0, 10, 100, 1000]


# ==========================================
# 1. 辅助函数：寻找对应的 MC 文件
# ==========================================
def find_mc_file(mass, lifetime):
    # 使用通配符寻找匹配的文件，防止时间戳不同导致找不到
    pattern = f"*_{mass}GeVCtau{lifetime}ps_options_*_skim.root"
    matches = glob.glob(os.path.join(mc_base_dir, pattern))
    return matches[0] if matches else None

# ==========================================
# 2. 交互式选择变量
# ==========================================
print("="*60)
print("      欢迎使用 HNL 物理特征随参数演化绘图工具")
print("="*60)

csv_file = "global_feature_ranking.csv"
if os.path.exists(csv_file):
    df_rank = pd.read_csv(csv_file)
    print("\n[*] 根据之前的 KS 检验，为您推荐区分度排名前 20 的变量：")
    # 打印成多列，比较美观
    top_vars = df_rank['feature'].head(20).to_list()
    for i in range(0, len(top_vars), 2):
        col1 = top_vars[i]
        col2 = top_vars[i+1] if i+1 < len(top_vars) else ""
        print(f"  - {col1:<30} {col2}")
else:
    print("\n[-] 未找到 global_feature_ranking.csv 文件，请根据记忆手动输入。")

print("-" * 60)
var_name = input("👉 请输入你想绘制的变量名 (例如 Jet1_PT): ").strip()

# ==========================================
# 3. 交互式选择绘图模式
# ==========================================
print("\n请选择对比模式：")
print("  [1] 固定飞行时间，观察信号随【质量 (Mass)】的变化")
print("  [2] 固定质量，观察信号随【飞行时间 (Lifetime)】的变化")
mode = input("👉 请输入 1 或 2: ").strip()

# --- 自动化文件夹命名 ---
mode_suffix = "ScanMass" if mode == "1" else "ScanLifetime"
output_dir = f"Plots_{var_name}_{mode_suffix}"
if not os.path.exists(output_dir): os.makedirs(output_dir)

# ==========================================
# 4. 读取本底数据 (Data)
# ==========================================
print(f"\n[*] 正在读取本底数据 (Data) 的 {var_name} 变量...")
try:
    with uproot.open(data_path) as f:
        df_data = f[tree_name].arrays([var_name], library="pd")
except Exception as e:
    print(f"[-] 读取本底数据失败，请确认变量名 {var_name} 是否拼写正确！")
    exit()

# ==========================================
# 5. 核心绘图函数 (包含正常坐标和对数坐标)
# ==========================================
def draw_and_save(df_bkg, dict_sig, title_prefix, out_prefix):
    if not dict_sig:
        print(f"[-] 警告: 没有找到信号数据，跳过绘制 {title_prefix}")
        return

    # 确定统一的 X 轴画图范围 (为了美观，排除了极端的 1% 的长尾)
    all_vals = [df[var_name].values for df in dict_sig.values() if len(df) > 0]
    if len(df_bkg) > 0:
        all_vals.append(df_bkg[var_name].values)
    
    global_min = min([np.percentile(v, 1) for v in all_vals if len(v) > 0])
    global_max = max([np.percentile(v, 99) for v in all_vals if len(v) > 0])
    bins = np.linspace(global_min, global_max, 50)

    # 用渐变色来画不同的信号线，更美观
    colors = plt.cm.jet(np.linspace(0.1, 0.9, len(dict_sig)))

    for use_log in [False, True]:
        plt.figure(figsize=(10, 7))
        
        # 画本底 (灰色填充)
        plt.hist(df_bkg[var_name], bins=bins, density=True, alpha=0.3, color='black', label='Background (Data)')

        # 画信号 (渐变色线条)
        for (label, df), color in zip(dict_sig.items(), colors):
            plt.hist(df[var_name], bins=bins, density=True, histtype='step', lw=2.5, color=color, label=f'Signal ({label})')

        scale_type = "Log" if use_log else "Linear"
        if use_log:
            plt.yscale('log')

        plt.title(f"{title_prefix} - {var_name} ({scale_type} Scale)", fontsize=15)
        plt.xlabel(var_name, fontsize=13)
        plt.ylabel('Normalized Events / a.u.', fontsize=13)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3, which="both", ls="--")

        # 保存
        out_filename = f"{output_dir}/{out_prefix}_{scale_type}.png"
        plt.savefig(out_filename, dpi=300)
        plt.close()
    
    print(f"[+] 成功生成图像: {output_dir}/{out_prefix}_Linear/Log.png")

# ==========================================
# 6. 根据模式收集数据并生成图像
# ==========================================
if mode == '1':
    # 模式 1：固定寿命，遍历质量
    for tau in lifetimes:
        print(f"\n[*] 正在处理: 固定寿命 = {tau}ps ...")
        sig_dict = {}
        for m in masses:
            fpath = find_mc_file(m, tau)
            if fpath:
                try:
                    with uproot.open(fpath) as f_mc:
                        df = f_mc[tree_name].arrays([var_name], library="pd")
                        sig_dict[f"{m} GeV"] = df
                except: pass
        
        draw_and_save(df_data, sig_dict, 
                      title_prefix=f"Fixed Lifetime {tau}ps", 
                      out_prefix=f"Mass_Scan_{var_name}_Tau{tau}ps")

elif mode == '2':
    # 模式 2：固定质量，遍历寿命
    for m in masses:
        print(f"\n[*] 正在处理: 固定质量 = {m}GeV ...")
        sig_dict = {}
        for tau in lifetimes:
            fpath = find_mc_file(m, tau)
            if fpath:
                try:
                    with uproot.open(fpath) as f_mc:
                        df = f_mc[tree_name].arrays([var_name], library="pd")
                        sig_dict[f"{tau} ps"] = df
                except: pass
        
        draw_and_save(df_data, sig_dict, 
                      title_prefix=f"Fixed Mass {m}GeV", 
                      out_prefix=f"Lifetime_Scan_{var_name}_Mass{m}GeV")

else:
    print("[-] 输入错误，程序退出。")

print("\n[✔] 所有绘图任务完成！请前往 Output_Physics_Plots 文件夹查看。")
