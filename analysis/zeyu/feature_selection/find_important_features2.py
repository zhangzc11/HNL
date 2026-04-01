import uproot
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import os
import glob

data_path = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"

def get_real_key(file_path):
    with uproot.open(file_path) as f:
        for k in f.keys():
            if "DecayTree" in k: return k.split(';')[0]
    return None

def main():
    d_key = get_real_key(data_path)
    with uproot.open(data_path) as f_data:
        data_keys = f_data[d_key].keys()
        
        # 【修改点 1：精准过滤逻辑】
        # 坚决剔除的无用后缀（完全匹配，防止误杀 PID 等变量）
        bad_suffixes = ["_KEY", "_ID", "EVENTNUMBER", "RUNNUMBER", "GPSTIME"]
        # 必须排除 TRUE 相关的蒙卡真值（Data中没有）
        
        physics_branches =[]
        for b in data_keys:
            if "TRUE" in b.upper() or "MC_" in b.upper(): continue
            if any(b.endswith(bad) for bad in bad_suffixes): continue
            if "Hlt" in b or "Spruce" in b: continue # 暂时排除繁杂的触发变量
            physics_branches.append(b)
            
    print(f"[*] 经过精准过滤，最终提取 {len(physics_branches)} 个核心物理变量。")

    print("[*] 正在加载 Data 数据...")
    with uproot.open(data_path) as f_data:
        data_dict = f_data[d_key].arrays(physics_branches, library="np")

    mc_files = glob.glob(os.path.join(mc_dir, "*_skim.root"))
    all_mc_results =[]

    for i, mc_file in enumerate(mc_files):
        fname = os.path.basename(mc_file)
        m_key = get_real_key(mc_file)
        if not m_key: continue

        try:
            with uproot.open(mc_file) as f_mc:
                t_mc = f_mc[m_key]
                common_vars =[v for v in physics_branches if v in t_mc.keys()]
                mc_dict = t_mc.arrays(common_vars, library="np")
                
                print(f"[{i+1}/{len(mc_files)}] 分析中: {fname}")
                for var in common_vars:
                    # 【修改点 2：忽略常数变量，比如 Mu_M】
                    if np.std(mc_dict[var]) < 1e-5: continue 
                    score, _ = ks_2samp(data_dict[var], mc_dict[var])
                    all_mc_results.append({'mc_file': fname, 'feature': var, 'ks_d': score})
        except:
            continue

    df_all = pd.DataFrame(all_mc_results)
    df_avg = df_all.groupby('feature')['ks_d'].agg(['mean', 'std']).reset_index()
    df_avg = df_avg.sort_values(by='mean', ascending=False)

    print("\n" + "="*60)
    print("   全样本综合区分度 Top 30 (已剔除常数质量)")
    print("="*60)
    print(df_avg.head(30).to_string(index=False))
    df_avg.to_csv("global_feature_ranking_v2.csv", index=False)

if __name__ == "__main__":
    main()

