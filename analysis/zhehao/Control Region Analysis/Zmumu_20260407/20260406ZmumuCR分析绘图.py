import uproot
import pandas as pd
import matplotlib.pyplot as plt
import hist
import numpy as np
from scipy.stats import ks_2samp
import os
from tqdm.auto import tqdm # For progress bars
import re # 导入 re 模块

# --- 配置部分 ---
# MC 模拟数据文件路径
MC_ROOT_FILE_PATH = "42112001_Zmumu_20260329_162457.root"
MC_TREE_PATH = "myTupleZmumu/DecayTree"

# Data 实际数据文件路径
DATA_ROOT_FILE_PATH = "data_MagUp_Sprucing24r1_WZControlRegion_all_disk_20260329_162828.root"
DATA_TREE_PATH = "myTupleZmumu/DecayTree"

OUTPUT_DIR = "20260406_Control_Region_Zmumu" # 输出文件夹


# 待绘制的branch列表，每个元素是一个字典，包含：
# 'name': branch名称 (字符串)
# 'bins': 直方图的bin设置 (整数表示bin数量，或一个数组表示bin边缘)
# 'range': 绘制的数据范围 (元组 (min, max)，可选，如果提供则优先使用)
# 'label': 绘图时x轴的标签 (字符串，可选，默认为name)
BRANCHES_TO_PLOT = [
    {'name': 'Z_M', 'bins': 50, 'range': [80000, 96000],'label': r'$Z$ $m$ (MeV)'}, # 聚焦Z峰
    {'name': 'Mup_PROBNN_D', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mup\_PROBNN\_D$'},
    {'name': 'Mup_PROBNN_E', 'bins': 50, 'range': [0, 0.01], 'label': r'$Mup\_PROBNN\_E$'},
    {'name': 'Mup_PROBNN_GHOST', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mup\_PROBNN\_GHOST$'},
    {'name': 'Mup_PROBNN_K', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mup\_PROBNN\_K$'},
    {'name': 'Mup_PROBNN_MU', 'bins': 50, 'range': [0.8, 1], 'label': r'$Mup\_PROBNN\_MU$'},
    {'name': 'Mup_PROBNN_P', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mup\_PROBNN\_P$'},
    {'name': 'Mup_PROBNN_PI', 'bins': 50, 'range': [0,0.1], 'label': r'$Mup\_PROBNN\_PI$'},
    {'name': 'Mum_PROBNN_D', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mum\_PROBNN\_D$'},
    {'name': 'Mum_PROBNN_E', 'bins': 50, 'range': [0, 0.01], 'label': r'$Mum\_PROBNN\_E$'},
    {'name': 'Mum_PROBNN_GHOST', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mum\_PROBNN\_GHOST$'},
    {'name': 'Mum_PROBNN_K', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mum\_PROBNN\_K$'},
    {'name': 'Mum_PROBNN_MU', 'bins': 50, 'range': [0.8, 1], 'label': r'$Mum\_PROBNN\_MU$'},
    {'name': 'Mum_PROBNN_P', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mum\_PROBNN\_P$'},
    {'name': 'Mum_PROBNN_PI', 'bins': 50, 'range': [0, 0.1], 'label': r'$Mum\_PROBNN\_PI$'}
]

# Pre-selection 条件 (字符串形式，可用于pandas query)
# 这是一个示例，请根据你的物理需求进行修改
# 注意：query 方法不支持多行字符串直接包含换行符，所以这里移除了换行符，用空格连接
# 同时确保运算符两侧有空格，避免解析错误
PRESELECTION_CONDITION = (
    "( (Mup_Hlt1SingleHighPtMuonDecision_TOS == 1) | (Mum_Hlt1SingleHighPtMuonDecision_TOS == 1) ) & "
    "( (Mup_Hlt2QEE_SingleHighPtMuonFullDecision_TOS == 1) | (Mum_Hlt2QEE_SingleHighPtMuonFullDecision_TOS == 1) ) & "
    "(Mup_PT > 20000) & (Mum_PT > 20000) & "
    "(Mup_ETA > 2.2) & (Mup_ETA < 4.4) & (Mum_ETA > 2.2) & (Mum_ETA < 4.4) & "
    "(Mup_P < 2000000) & (Mum_P < 2000000) & "
    "(Mup_TRCHI2DOF < 2.5) & (Mum_TRCHI2DOF < 2.5) & "
    "(Mup_MINIPCHI2 < 16) & (Mum_MINIPCHI2 < 16) & "   # 原先文档中没有，进一步确保来自原初顶点
    #这里处理了 Mup_OWNPVNDOF 可能为 0 的情况，避免除以零错误
    "( (Mup_OWNPVNDOF > 0 & Mup_OWNPVIPCHI2 / Mup_OWNPVNDOF < 100) & Mup_OWNPVNDOF != 0 ) & " 
    "( (Mum_OWNPVNDOF > 0 & Mum_OWNPVIPCHI2 / Mum_OWNPVNDOF < 100) & Mum_OWNPVNDOF != 0 ) & "    
    "(Z_M > 80000) & (Z_M < 96000)"
)


# --- 主函数 ---
def analyze_zmumu_data():
    """
    读取ROOT文件中的Zmumu数据和MC数据，进行预选，
    绘制MC-data对比图和data/MC比率图，并计算KS值。
    """
    print(f"开始分析 MC 文件: {MC_ROOT_FILE_PATH}")
    print(f"开始分析 Data 文件: {DATA_ROOT_FILE_PATH}")

    # 确保输出目录存在
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"图表和KS值将保存到目录: {OUTPUT_DIR}")

    # 收集需要读取的所有 branch 名称
    branches_to_read = []
    # 从 BRANCHES_TO_PLOT 中添加
    for b_info in BRANCHES_TO_PLOT:
        if b_info['name'] not in branches_to_read: # 避免重复添加
            branches_to_read.append(b_info['name'])

    # 从 PRESELECTION_CONDITION 中提取变量并添加
    # 查找所有像 'Z_M' 这样的单词，排除数字和运算符
    preselection_vars_from_condition = set(re.findall(r'\b[a-zA-Z_]\w*\b', PRESELECTION_CONDITION))
    # 排除 Python 关键字和数字
    keywords_to_exclude = {'True', 'False', 'and', 'or', 'not'}
    preselection_vars_from_condition = {v for v in preselection_vars_from_condition if v not in keywords_to_exclude and not v.isdigit()}
    
    for var in preselection_vars_from_condition:
        if var not in branches_to_read: # 避免重复添加
            branches_to_read.append(var)
            
    # 获取所有可用的 branch 名称，以避免读取不存在的 branch
    try:
        with uproot.open(DATA_ROOT_FILE_PATH) as data_file:
            data_tree_keys = data_file[DATA_TREE_PATH].keys()
        with uproot.open(MC_ROOT_FILE_PATH) as mc_file:
            mc_tree_keys = mc_file[MC_TREE_PATH].keys()
        
        # 为了确保两个文件中都存在，取交集
        all_available_branches = set(data_tree_keys) & set(mc_tree_keys)
    except Exception as e:
        print(f"无法打开 ROOT 文件以获取 branch 列表: {e}")
        return

    # 过滤掉不存在的 branch，只读取共同存在的
    final_branches_to_read = []
    for branch_name in branches_to_read:
        if branch_name in all_available_branches:
            final_branches_to_read.append(branch_name)
        else:
            print(f"警告: Branch '{branch_name}' 在一个或两个 ROOT 文件中不存在，将不会被读取。")

    print(f"将要读取的唯一 Branch 列表 ({len(final_branches_to_read)} 个): {final_branches_to_read}")

    # 1. 读取数据
    print("1. 读取数据...")
    try:
        # 读取数据 (实际数据)
        print(f"  从 '{DATA_ROOT_FILE_PATH}' 读取 {DATA_TREE_PATH}，分支: {len(final_branches_to_read)} 个...")
        with uproot.open(DATA_ROOT_FILE_PATH) as data_file:
            data_tree = data_file[DATA_TREE_PATH]
            data_df = data_tree.arrays(final_branches_to_read, library="pd")

        # 读取MC (模拟数据)
        print(f"  从 '{MC_ROOT_FILE_PATH}' 读取 {MC_TREE_PATH}，分支: {len(final_branches_to_read)} 个...")
        with uproot.open(MC_ROOT_FILE_PATH) as mc_file:
            mc_tree = mc_file[MC_TREE_PATH]
            mc_df = mc_tree.arrays(final_branches_to_read, library="pd")

    except Exception as e:
        print(f"错误：无法读取 ROOT 文件或指定的 Tree: {e}")
        return

    print(f"  数据样本原始事件数: {len(data_df)}")
    print(f"  MC样本原始事件数: {len(mc_df)}")

    # 2. 应用Pre-selection
    print("2. 应用Pre-selection...")
    try:
        data_df_preselected = data_df.query(PRESELECTION_CONDITION)
        mc_df_preselected = mc_df.query(PRESELECTION_CONDITION)
    except Exception as e:
        print(f"应用 PRESELECTION_CONDITION 失败。请检查条件中的变量名是否与 DataFrame 列名匹配，以及条件语法是否正确。错误: {e}")
        # 打印出 DataFrame 列名，帮助调试
        print(f"Data DataFrame columns: {data_df.columns.tolist()}")
        print(f"MC DataFrame columns: {mc_df.columns.tolist()}")
        return

    print(f"  数据样本预选后事件数: {len(data_df_preselected)}")
    print(f"  MC样本预选后事件数: {len(mc_df_preselected)}")

    if len(data_df_preselected) == 0 or len(mc_df_preselected) == 0:
        print("警告：预选后数据或MC样本为空，无法继续绘图和分析。")
        return

    # 确定MC归一化因子
    # 简单的归一化方式：让MC的总事件数与数据一致
    mc_norm_factor = 1.0
    #简单归一化：len(data_df_preselected) / len(mc_df_preselected)
    print(f"  MC 归一化因子 (使MC总事件数与数据一致): {mc_norm_factor:.4f}")

    # 3. 绘制图表并计算KS值
    print("3. 绘制图表并计算KS值...")
    # 遍历 BRANCHES_TO_PLOT，而不是 branches_to_read，因为你只打算绘制这些
    for branch_info in tqdm(BRANCHES_TO_PLOT, desc="处理Branch"):
        branch_name = branch_info['name']
        
        # 仅绘制存在于预选DataFrame中的branch
        if branch_name not in data_df_preselected.columns or branch_name not in mc_df_preselected.columns:
            print(f"  跳过 Branch '{branch_name}' (预选DataFrame中不存在或已过滤)。")
            continue

        bins = branch_info['bins']
        plot_range = branch_info.get('range', None)
        x_label = branch_info.get('label', branch_name)

        print(f"  处理 Branch: {branch_name}")

        data_values = data_df_preselected[branch_name].values
        mc_values = mc_df_preselected[branch_name].values

        actual_min = 0.0
        actual_max = 0.0
        if plot_range: # 如果用户指定了范围，就用指定的
            actual_min, actual_max = plot_range
            # 同时裁剪 values
            data_values = data_values[(data_values >= actual_min) & (data_values <= actual_max)]
            mc_values = mc_values[(mc_values >= actual_min) & (mc_values <= actual_max)]
        else: # 否则，自动从数据和MC中确定
            # 确保数组非空，避免 np.min/max 错误
            all_values = np.concatenate([data_values, mc_values])
            if all_values.size == 0:
                print(f"    警告：Branch '{branch_name}' 在预选后没有数据，无法确定自动范围，跳过。")
                continue
            
            actual_min, actual_max = np.min(all_values), np.max(all_values)
            
            # 对自动确定的范围做一些填充，防止min==max导致绘图失败
            if actual_min == actual_max:
                if actual_min == 0:
                    actual_max = 1.0 # 随便给个范围
                else:
                    actual_min *= 0.9
                    actual_max *= 1.1
            # 确保 min < max, 否则 hist 库会报错
            if actual_min >= actual_max:
                actual_max = actual_min + 1.0 # 强制给一个非零宽度

        if len(data_values) == 0 or len(mc_values) == 0:
            print(f"    警告：Branch '{branch_name}' 在指定范围预选后没有有效数据，跳过。")
            continue

        # 使用 hist 库创建直方图
        if isinstance(bins, int):
            h_data = hist.new.Reg(bins, actual_min, actual_max, name=branch_name, label=x_label).Double().fill(data_values)
            h_mc = hist.new.Reg(bins, actual_min, actual_max, name=branch_name, label=x_label).Double().fill(mc_values)
        elif isinstance(bins, (list, np.ndarray)): # 如果是 bin 边缘数组，则范围由数组决定
            h_data = hist.new.Var(bins, name=branch_name, label=x_label).Double().fill(data_values)
            h_mc = hist.new.Var(bins, name=branch_name, label=x_label).Double().fill(mc_values)
        else:
            raise ValueError(f"'{branch_name}' 的 bin 设置 '{bins}' 类型不支持。请使用整数或 numpy 数组。")
        
        # 归一化MC直方图
        h_mc_norm = h_mc * mc_norm_factor

        # 绘制 MC-data 对比图
        fig, (ax_main, ax_ratio) = plt.subplots(
            nrows=2,
            ncols=1,
            gridspec_kw={"height_ratios": (3, 1)},
            sharex=True,
            figsize=(10, 8)
        )
        fig.subplots_adjust(hspace=0.05)

        # 主图
        h_mc_norm.plot(ax=ax_main, histtype='step', label='MC (Normalized)', color='red', density=False)
        h_data.plot(ax=ax_main, histtype='errorbar', label='Data', color='black', marker='.', linestyle='None', density=False)

        ax_main.set_ylabel("Events")
        ax_main.set_title(f"Data vs MC Comparison: {branch_name}")
        ax_main.legend()
        ax_main.grid(axis='y', linestyle='--', alpha=0.7)

        # 比率图 (Data / MC)
        mc_counts = h_mc_norm.values()
        data_counts = h_data.values()

        ratio_vals = np.zeros_like(data_counts, dtype=float)
        valid_indices = mc_counts > 0
        ratio_vals[valid_indices] = data_counts[valid_indices] / mc_counts[valid_indices]

        data_errors = np.sqrt(data_counts)
        ratio_yerr = np.zeros_like(data_errors, dtype=float)
        error_valid_indices = (mc_counts > 0) & (data_counts > 0)
        ratio_yerr[error_valid_indices] = data_errors[error_valid_indices] / mc_counts[error_valid_indices]

        bin_centers = h_data.axes[0].centers

        ax_ratio.errorbar(bin_centers[valid_indices], ratio_vals[valid_indices],
                          yerr=ratio_yerr[valid_indices],
                          fmt='k.',
                          label='Data/MC')
        ax_ratio.axhline(1, color='gray', linestyle='--', alpha=0.7)
        ax_ratio.set_ylabel("Data/MC")
        ax_ratio.set_xlabel(x_label)
        ax_ratio.set_ylim(0.5, 1.5) # 通常比率图的Y轴范围
        ax_ratio.grid(axis='y', linestyle='--', alpha=0.7)

        plot_filename = os.path.join(OUTPUT_DIR, f"Z_{branch_name}.png")
        plt.savefig(plot_filename)
        plt.close(fig)
        print(f"    图表保存到: {plot_filename}")

        # 计算 KS 值
        ks_statistic, ks_pvalue = ks_2samp(data_values, mc_values)

        # 将 KS 值记录到 txt 文件
        ks_filename = os.path.join(OUTPUT_DIR, f"Z_{branch_name}.txt")
        with open(ks_filename, 'w') as f:
            f.write(f"Branch: {branch_name}")
            f.write(f"Pre-selection condition: {PRESELECTION_CONDITION}")
            f.write(f"MC Normalization Factor: {mc_norm_factor:.4f}")
            f.write(f"Kolmogorov-Smirnov Statistic: {ks_statistic:.4f}")
            f.write(f"Kolmogorov-Smirnov p-value: {ks_pvalue:.4e}")
            if plot_range:
                f.write(f"Analysis Range: [{plot_range[0]}, {plot_range[1]}]")

        print(f"    KS 统计量: {ks_statistic:.4f}, p-value: {ks_pvalue:.4e}")
        print(f"    KS 值保存到: {ks_filename}")

    print("分析完成！")

# --- 运行分析 ---
if __name__ == "__main__":
    analyze_zmumu_data()

