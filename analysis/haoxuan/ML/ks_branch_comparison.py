import uproot
import awkward as ak
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import os
import sys

candidate_directories = ["myTupleOS1J", "myTupleOS2J", "myTupleSS1J", "myTupleSS2J"]
# ----------------- 配置 -----------------
# 文件路径（MC 与 data），以及输入/输出文件名
mc_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"
data_dir = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/"
mc_file = "42912010_5GeVCtau0ps_options_20260228_111712_skim.root"
data_file = "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_file_path = os.path.join(mc_dir, mc_file)
local_data_file = os.path.join(data_dir, data_file)

# 本地 data 文件路径优先
data_file_path = local_data_file
if not os.path.exists(data_file_path):
    print(f"错误：本地数据文件 {data_file_path} 未找到。")
    sys.exit(1)

if not os.path.exists(mc_file_path):
    print(f"错误：MC 文件 {mc_file_path} 未找到。")
    sys.exit(1)

top = 50
branch_list_file = "./csv/filtered_branches.csv"
output_csv = f"./csv/top{top}_ks_branches.csv"

# 候选目录列表（将在两个文件中寻找共同存在的目录）
candidate_directories = ["myTupleOS1J", "myTupleOS2J", "myTupleSS1J", "myTupleSS2J"]

def main():
    print("开始进行 KS 检验分析...")

    # 1) 读取分支列表（来自 filtered_branches.csv）
    if not os.path.exists(branch_list_file):
        print(f"错误：分支列表文件 {branch_list_file} 未找到。")
        sys.exit(1)

    df_branches = pd.read_csv(branch_list_file)
    if "Branch" not in df_branches.columns:
        print("错误：CSV 必须包含 'Branch' 列。")
        sys.exit(1)

    target_branches = df_branches["Branch"].tolist()
    print(f"从 {branch_list_file} 加载了 {len(target_branches)} 个分支")

    # 2) 在 data 与 MC 文件中选择一个共同存在的目录
    selected_dir = None

    try:
        with uproot.open(data_file_path) as f_data, uproot.open(mc_file_path) as f_mc:
            # 在候选目录中查找第一个两文件均包含的目录
            for dirname in candidate_directories:
                tree_path = f"{dirname}/DecayTree"
                if tree_path in f_data and tree_path in f_mc:
                    selected_dir = dirname
                    break

            if selected_dir is None:
                print("错误：在候选目录中未找到两个文件均包含的共同目录。")
                sys.exit(1)

            # 选定的树路径
            tree_name = f"{selected_dir}/DecayTree"
            tree_data = f_data[tree_name]
            tree_mc = f_mc[tree_name]

            # 验证分支在两个 tree 中是否存在（去掉可能的 ";x" 后缀）
            available_data_branches = set(k.split(";")[0] for k in tree_data.keys())
            available_mc_branches = set(k.split(";")[0] for k in tree_mc.keys())

            # 取交集并排序，作为有效比较分支
            target_set = set(target_branches)
            valid_branches = sorted(list(target_set.intersection(available_data_branches).intersection(available_mc_branches)))
            print(f"用于比较的有效分支数量: {len(valid_branches)}（总 {len(target_branches)}）")

            # 3) 对每个有效分支执行 KS 检验
            print("正在计算 KS 统计量...")
            ks_results = []

            total = len(valid_branches)
            len_ks = 0
            for i, branch in enumerate(valid_branches):
                # 每 10 个分支更新一次进度
                if i % 10 == 0:
                    print(f"Processing [{i+1}/{total}] {branch}...", end="\r")

                try:
                    # 读取分支数组（使用 awkward 以兼容 jagged/矢量分支）
                    b_data = tree_data[branch].array(library="ak")
                    b_mc = tree_mc[branch].array(library="ak")

                    # 展平 jagged 数组
                    data_flat = ak.flatten(b_data, axis=None)
                    mc_flat = ak.flatten(b_mc, axis=None)

                    # 转为 numpy 并剔除 NaN
                    data_np = ak.to_numpy(data_flat)
                    mc_np = ak.to_numpy(mc_flat)
                    data_np = data_np[~np.isnan(data_np)]
                    mc_np = mc_np[~np.isnan(mc_np)]

                    if len(data_np) == 0 or len(mc_np) == 0:
                        continue

                    # 计算 KS 统计量（返回 statistic, pvalue）
                    res = ks_2samp(data_np, mc_np)
                    ks_stat = res.statistic
                    ks_results.append({"Branch": branch, "KS_Value": ks_stat})
                    len_ks += 1
                except Exception:
                    # 若读取或计算失败则跳过该分支
                    continue

    except Exception as e:
        print(f"严重错误: {e}")
        sys.exit(1)

    print(f"\nKS 检验完成。成功计算的分支数: {len_ks}")

    # 4) 排序并输出 top 30
    if ks_results:
        df_res = pd.DataFrame(ks_results)
        df_res = df_res.sort_values(by="KS_Value", ascending=False)
        df_top = df_res.head(top)
        df_top.to_csv(output_csv, index=False)
        print(f"已将 top {top} KS 分支保存为 {output_csv}")
        print(df_top)
    else:
        print("无结果可保存。")

if __name__ == "__main__":
    main()
