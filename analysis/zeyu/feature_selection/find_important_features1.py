import uproot
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import os
import glob

data_path = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"

try:
    from branch_config import all_branches
except ImportError:
    print("[-] 错误: 找不到 branch_config.py")
    exit()

def get_real_key(file_path):
    with uproot.open(file_path) as f:
        for k in f.keys():
            if "DecayTree" in k:
                return k.split(';')[0]
    return None

def main():
    # 1. 探测 Data 的 Tree 和 变量
    d_key = get_real_key(data_path)
    with uproot.open(data_path) as f_data:
        t_data = f_data[d_key]
        data_keys = t_data.keys()
        
        # 2. 过滤变量：
        #    a) 必须在 Data 里存在
        #    b) 排除 Key, ID, Type (索引类)
        #    c) 排除 TRUE (真值类，因为 Data 没有)
        physics_branches = [
            b for b in all_branches 
            if b in data_keys 
            and not any(x in b.lower() for x in ["_key", "_id", "_type", "true"])
        ]
        
    print(f"[*] 经过过滤，最终选取 {len(physics_branches)} 个 Data/MC 共有的物理变量。")

    # 3. 预加载 Data
    print("[*] 正在加载 Data 数据...")
    with uproot.open(data_path) as f_data:
        t_data = f_data[d_key]
        data_dict = t_data.arrays(physics_branches, library="np")

    # 4. 遍历 MC
    mc_files = glob.glob(os.path.join(mc_dir, "*_skim.root"))
    all_mc_results = []

    for i, mc_file in enumerate(mc_files):
        fname = os.path.basename(mc_file)
        m_key = get_real_key(mc_file)
        if not m_key: continue

        try:
            with uproot.open(mc_file) as f_mc:
                t_mc = f_mc[m_key]
                # 确保只读取当前 MC 也拥有的变量
                common_vars = [v for v in physics_branches if v in t_mc.keys()]
                mc_dict = t_mc.arrays(common_vars, library="np")
                
                print(f"[{i+1}/{len(mc_files)}] 分析中: {fname}")
                for var in common_vars:
                    score, _ = ks_2samp(data_dict[var], mc_dict[var])
                    all_mc_results.append({'mc_file': fname, 'feature': var, 'ks_d': score})
        except:
            continue

    # 5. 汇总
    df_all = pd.DataFrame(all_mc_results)
    df_avg = df_all.groupby('feature')['ks_d'].agg(['mean', 'std']).reset_index()
    df_avg = df_avg.sort_values(by='mean', ascending=False)

    print("\n" + "="*50)
    print("   全样本（不同质量/寿命）综合区分度 Top 15")
    print("="*50)
    print(df_avg.head(15).to_string(index=False))
    print("="*50)
    df_avg.to_csv("global_feature_ranking.csv", index=False)

if __name__ == "__main__":
    main()
