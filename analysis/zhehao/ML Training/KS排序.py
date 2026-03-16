import uproot
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import warnings
from filtered_branches import filtered_branches 

# 忽略计算中的小警告（如某些分布全为0时的情况）
warnings.filterwarnings("ignore")

# ==========================================
# 1. 配置参数
# ==========================================
# 文件路径（请确保路径正确）
data_path = "/home/steven/HNL/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
signal_path = "/home/steven/HNL/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"

# ROOT 树的名称
tree_name =  "myTupleSS1J/DecayTree"

# 设置你想要保留的最高 KS 值的变量数量
TOP_N = 50 

# ==========================================
# 2. 自动化 KS 筛选逻辑
# ==========================================
def get_high_ks_branches(sig_path, dat_path, t_name):
    ks_list = []
    
    # 1. 获取所有可用的 Branch 列表
    all_branches = filtered_branches
    
    print(f"检测到总 Branch 数量: {len(all_branches)}")
    print("正在开始逐一扫描变量区分度，请稍候...")

    # 2. 逐个 Branch 进行分析
    with uproot.open(sig_path) as f_sig, uproot.open(dat_path) as f_dat:
        t_sig = f_sig[t_name]
        t_dat = f_dat[t_name]
        
        for i, branch in enumerate(all_branches):
            try:
                # 读取单个 Branch 数据
                sig_arr = t_sig.arrays([branch], library="np")[branch]
                dat_arr = t_dat.arrays([branch], library="np")[branch]
                
                # 数据预处理：剔除 NaN 和 Inf
                sig_arr = sig_arr[np.isfinite(sig_arr)]
                dat_arr = dat_arr[np.isfinite(dat_arr)]
                
                # 检查是否为有效数值且不为空
                if len(sig_arr) == 0 or len(dat_arr) == 0:
                    continue
                if not np.issubdtype(sig_arr.dtype, np.number):
                    continue
                
                # 计算 KS 统计量
                # ks_stat 越大，区分度越好
                ks_stat, _ = ks_2samp(sig_arr, dat_arr)
                
                ks_list.append((branch, ks_stat))
                
                # 打印进度（每50个变量提示一次）
                if (i + 1) % 50 == 0:
                    print(f"已处理 {i + 1} / {len(all_branches)} 个变量...")
                    
            except Exception as e:
                # 跳过无法处理的特殊变量（如复杂的嵌套对象）
                continue

    # 3. 排序
    # 根据 KS 统计量从大到小排序
    ks_list.sort(key=lambda x: x[1], reverse=True)
    
    return ks_list

# ==========================================
# 3. 执行分析与输出结果
# ==========================================
if __name__ == "__main__":
    final_results = get_high_ks_branches(signal_path, data_path, tree_name)
    
    print("" + "="*50)
    print(f" 筛选完成！以下是 KS 值最高的前 {TOP_N} 个 Branch 列表 ")
    print("="*50)
    
    # 输出格式化后的列表，方便复制到训练程序中
    top_branches_list = []
    for name, score in final_results[:TOP_N]:
        top_branches_list.append(name)
        print(f"KS: {score:.4f} | Branch: {name}")
    
    print("--- 建议用于机器学习训练的列表 (Python 格式) ---")
    print(top_branches_list)
    print("="*50)
    print("注意：在使用上述变量前，请手动剔除诸如 'EventNumber'、'runNumber' 等物理无关变量。")
