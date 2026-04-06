import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler # 引入StandardScaler用于特征标准化

# ==============================
# 1. 路径设置与全局变量
# ==============================
DATA_DIR = "/home/steven/HNL"  # 数据文件根目录
# 真实数据文件 (Background)
DATA_FILE = os.path.join(
    DATA_DIR, "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
)

# 输出目录
OUTPUT_DIR = "./ML_Training_HNL_CascadeGBDT"
os.makedirs(OUTPUT_DIR, exist_ok=True) # 确保输出目录存在

# ==============================
# 2. 读取 ROOT 的函数
# ==============================
def load_root(file_path, branches, tree_name="DecayTree", max_entries=None):
    """
    从 ROOT 文件中读取指定 branches，返回 pandas.DataFrame
    max_entries=None 表示读全部；否则只读前 max_entries 个事件
    """
    try:
        with uproot.open(file_path) as f:
            # 检查文件是否包含 tree_name，并处理子目录情况
            # uproot v5 行为可能有所不同，这里采取兼容性写法
            if '/' in tree_name:
                top_dir, sub_tree = tree_name.split('/', 1)
                if top_dir not in f:
                    print(f"Warning: Top directory '{top_dir}' not found in {file_path}. Skipping.")
                    return pd.DataFrame()
                if sub_tree not in f[top_dir]:
                    print(f"Warning: Subtree '{sub_tree}' not found in {file_path}/{top_dir}. Skipping.")
                    return pd.DataFrame()
                tree = f[top_dir][sub_tree]
            elif tree_name not in f:
                print(f"Warning: Tree '{tree_name}' not found in {file_path}. Skipping.")
                return pd.DataFrame()
            else:
                tree = f[tree_name]

            # 检查是否有所有需要的 branches
            missing_branches = [b for b in branches if b not in tree.keys()]
            if missing_branches:
                print(f"Warning: Missing branches {missing_branches} in tree '{tree_name}' in {file_path}. Skipping.")
                return pd.DataFrame()
            
            df = tree.arrays(branches, entry_stop=max_entries, library="pd")
            return df
    except Exception as e:
        print(f"Error reading {file_path} with tree '{tree_name}': {e}")
        return pd.DataFrame() # 发生错误时返回空DataFrame

# ==============================
# 3. 交互式选择 MC 文件
# ==============================
# 请在这里填充你的 MC 文件列表
# 这是一个示例列表，请替换为你的真实 MC 文件名
MC_FILE_OPTIONS = [
    "42912010_5GeVCtau0ps_options_20260228_111712_skim.root",
    "42912011_5GeVCtau10ps_options_20260228_111712_skim.root",
    "42912012_10GeVCtau0ps_options_20260228_111712_skim.root",
    "42912013_10GeVCtau10ps_options_20260228_111712_skim.root",
    "42912014_15GeVCtau0ps_options_20260228_111712_skim.root",
    "42912015_15GeVCtau10ps_options_20260228_111712_skim.root",
    "42912016_20GeVCtau0ps_options_20260228_111712_skim.root",
    "42912017_20GeVCtau10ps_options_20260228_111712_skim.root",
    "42912018_30GeVCtau0ps_options_20260228_111712_skim.root",
    "42912019_30GeVCtau10ps_options_20260228_230635_skim.root",
    "42912020_50GeVCtau0ps_options_20260228_230639_skim.root",
    "42912021_50GeVCtau10ps_options_20260228_111712_skim.root",
    "42912030_5GeVCtau100ps_options_20260228_111712_skim.root",
    "42912031_5GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912032_10GeVCtau100ps_options_20260228_111712_skim.root",
    "42912033_10GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912034_15GeVCtau100ps_options_20260228_111712_skim.root",
    "42912035_15GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912036_20GeVCtau100ps_options_20260228_111712_skim.root",
    "42912037_20GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912038_30GeVCtau100ps_options_20260228_111712_skim.root",
    "42912039_30GeVCtau1000ps_options_20260228_111712_skim.root",
    "42912040_50GeVCtau100ps_options_20260228_111712_skim.root",
    "42912041_50GeVCtau1000ps_options_20260228_111712_skim.root",
]

def select_mc_files_interactive(data_dir, mc_file_options):
    """
    交互式选择 MC 文件，用户可以从列表中选择多个。
    """
    if not mc_file_options:
        print("\nMC 文件列表为空。请手动输入文件名（一行一个，空行结束）：")
        selected_filenames = []
        while True:
            line = input().strip()
            if not line:
                break
            selected_filenames.append(line)
    else:
        print("\n请选择要使用的 MC 文件（输入序号，多个用逗号分隔，'all' 选择全部）：")
        for i, fname in enumerate(mc_file_options):
            print(f"   {i+1}. {fname}")
        
        while True:
            choice = input("请输入你的选择: ").strip().lower()
            if choice == 'all':
                selected_filenames = mc_file_options
                break
            
            try:
                indices = [int(idx.strip()) - 1 for idx in choice.split(',')]
                selected_filenames = [mc_file_options[i] for i in indices if 0 <= i < len(mc_file_options)]
                if selected_filenames:
                    break
                else:
                    print("无效的输入。请重新输入有效的序号。")
            except ValueError:
                print("无效的输入格式。请重新输入序号或 'all'。")

    if not selected_filenames:
        raise ValueError("未选择任何 MC 文件。")

    mc_files_full_paths = sorted([os.path.join(data_dir, f) for f in selected_filenames])
    
    # 检查文件是否存在
    existing_mc_files = []
    for f in mc_files_full_paths:
        if os.path.exists(f):
            existing_mc_files.append(f)
        else:
            print(f"警告: 文件 {f} 不存在于 {data_dir} 目录中，已跳过。")
            
    if not existing_mc_files:
        raise FileNotFoundError(f"没有找到任何存在的 MC 文件。请检查输入和 {data_dir} 目录。")
            
    print(f"已选择 {len(existing_mc_files)} 个 MC 文件：")
    for f in existing_mc_files:
        print(f" - {os.path.basename(f)}")
    return existing_mc_files

# ==============================
# 4. 交互式选择 TREE_NAME
# ==============================
def select_tree_name_interactive():
    """
    交互式选择 ROOT 文件中的树名，支持单选或合并选择。
    """
    available_trees = [
        "myTupleOS1J/DecayTree",              
        "myTupleOS2J/DecayTree",              
        "myTupleSS1J/DecayTree",              
        "myTupleSS2J/DecayTree"
    ]
    
    print("请选择 ROOT 文件中的 DecayTree 读取模式：")
    print("1. 单选模式：从以下列表中选择一个 DecayTree。")
    for i, tree in enumerate(available_trees):
        print(f"   {i+1}. {tree}")
    print("2. 合并模式：选择 SS 或 OS，系统将合并两个相应的 DecayTree 并选择更接近真实 W_M 的事件。")
    
    while True:
        choice = input("请输入你的选择 (1-4 代表单选，'SS' 或 'OS' 代表合并)：").strip().upper()
        
        if choice in ['1', '2', '3', '4']:
            idx = int(choice) - 1
            selected_tree_names = [available_trees[idx]]
            print(f"已选择单选模式：{selected_tree_names[0]}")
            return selected_tree_names, False # False 表示非合并模式
        elif choice == 'SS':
            selected_tree_names = ["myTupleSS1J/DecayTree", "myTupleSS2J/DecayTree"]
            print(f"已选择合并模式 (SS)：{selected_tree_names}")
            return selected_tree_names, True # True 表示合并模式
        elif choice == 'OS':
            selected_tree_names = ["myTupleOS1J/DecayTree", "myTupleOS2J/DecayTree"]
            print(f"已选择合并模式 (OS)：{selected_tree_names}")
            return selected_tree_names, True # True 表示合并模式
        else:
            print("无效的输入。请重新输入。")

# ==============================
# 5. 选择要使用的变量 (VARIABLES)
# ==============================
# 你的原始 VARIABLES 列表，我将保留并取消注释部分，以提供更全面的特征。
# 注意：确保这些变量在你的ROOT文件中都存在！
VARIABLES = [
    "Mu_ELECTRONENERGY", "MuNuR_ELECTRONENERGY",
    "Mu_ELECTRONSHOWERDLL", "MuNuR_ELECTRONSHOWERDLL",
    "Mu_ELECTRONSHOWEREOP", "MuNuR_ELECTRONSHOWEREOP",
    "Mu_HCALEOP", "MuNuR_HCALEOP",
    "W_MmuWmuN",
    "MuNuR_log_PT", "NuR_log_PT", "Mu_log_PT", "W_log_PT",
    "MuNuR_log_P", "NuR_log_P", "Mu_log_P", "W_log_P",
    "W_M",
    "W_E", "MuNuR_M", "MuNuR_E",
    "W_ETA", "W_THETA", "Mu_ETA", "Mu_THETA",
    "Mu_PROBNN_PI", "Mu_PROBNN_MU", "Mu_PROBNN_GHOST", "Mu_PROBNN_E","Mu_PROBNN_P","Mu_PROBNN_K",
    "MuNuR_PROBNN_PI","MuNuR_PROBNN_K", "MuNuR_PROBNN_GHOST","MuNuR_PROBNN_E","MuNuR_PROBNN_P","MuNuR_PROBNN_MU",
    "nTTracks", "nVeloTracks",
    "MuNuR_MINIP", "Mu_MINIP", "W_MINIP"
]

# 在合并模式下需要额外的变量
ADDITIONAL_MERGE_VARS = ["EVENTNUMBER", "W_M"] 
TRUE_W_M_FOR_SELECTION = 80360  # 真实的 W 质量，单位 MeV/c^2

# ==============================
# 6. Pre-selection 函数
# 这是你未来需要填充切割条件的地方
# ==============================
def apply_pre_selection(dataset):
    """
    对数据集应用预筛选（pre-selection）切割。
    请在这里根据物理需求添加具体的切割条件。
    """
    original_len = len(dataset)
    original_sig_len = (dataset['label'] == 1).sum()
    original_bkg_len = (dataset['label'] == 0).sum()

    print(f"--- 原始数据量 (Pre-selection前) ---")
    print(f"信号事件数: {original_sig_len}")
    print(f"本底事件数: {original_bkg_len}")
    
    # --- 示例切割条件 (请根据你的物理需求修改或添加) ---
    dataset = dataset[
        (dataset["Mu_PROBNN_MU"] > 0.5) & (dataset["MuNuR_PROBNN_MU"] > 0.5) & # 保证是muon
        (dataset["Mu_PROBNN_GHOST"] < 0.1) & (dataset["Mu_PROBNN_GHOST"] < 0.1) & # cut明显的假径迹
        (dataset["W_M"] > 50000) & (dataset["W_M"] < 100000) # W质量窗口 
     ]
    # --------------------------------------------------

    # 重要的：数据清洗，去除缺失值和无穷大
    # 这应该在任何切割条件之前进行，以确保数据质量
    dataset = dataset.replace([np.inf, -np.inf], np.nan).dropna(subset=VARIABLES) # 只对使用的特征进行清洗
    
    cut_sig_len = (dataset['label'] == 1).sum()
    cut_bkg_len = (dataset['label'] == 0).sum()

    print(f"--- Pre-selection (包含清洗) 后数据量 ---")
    print(f"信号事件数: {cut_sig_len}")
    print(f"本底事件数: {cut_bkg_len}")
    
    pre_sel_sig_eff = cut_sig_len / original_sig_len if original_sig_len > 0 else 0
    pre_sel_bkg_eff = cut_bkg_len / original_bkg_len if original_bkg_len > 0 else 0
    
    return dataset, pre_sel_sig_eff, pre_sel_bkg_eff

# ==============================
# 7. 级联 GBDT 的配置
# ==============================
# 定义级联 GBDT 的数量和每个模型的特征子集
# 这里的 FEATURES_FOR_GBDT_STAGES 是一个列表的列表，每个子列表包含一个 GBDT 阶段使用的特征
# 你可以根据物理考量来划分特征，例如：
# Stage 1: W
# Stage 2:第一个mu
# Stage 3: 第二个mu以及nu
FEATURES_FOR_GBDT_STAGES = [
    ["W_log_P","W_ETA","W_THETA","W_MINIP"], # W branches
    ["Mu_log_PT", "Mu_MINIP","Mu_log_P","Mu_ETA","Mu_THETA"], # First Muon branches
    ["MuNuR_E", "MuNuR_log_PT", "NuR_log_PT","MuNuR_log_P","NuR_log_P","MuNuR_MINIP"], # Munu&nu branches
]

# 确保所有 FEATURES_FOR_GBDT_STAGES 中的特征都包含在 VARIABLES 中
for stage_features in FEATURES_FOR_GBDT_STAGES:
    for feature in stage_features:
        if feature not in VARIABLES:
            raise ValueError(f"Feature '{feature}' in GBDT stage configuration is not in the main VARIABLES list.")

NUM_GBDTS = len(FEATURES_FOR_GBDT_STAGES)
GBDT_MODELS = [] # 存储训练好的 GBDT 模型
SCALERS = [] # 存储每个 GBDT 阶段的 StandardScaler

# ==============================
# Main Script Logic
# ==============================

if __name__ == "__main__":
    # --- 1. 获取 MC 文件列表 ---
    try:
        MC_FILES = select_mc_files_interactive(DATA_DIR, MC_FILE_OPTIONS)
    except (ValueError, FileNotFoundError) as e:
        print(f"错误: {e}")
        exit()

    # --- 2. 选择 TREE_NAME ---
    SELECTED_TREE_NAMES, IS_MERGE_MODE = select_tree_name_interactive()
    
    # 如果是合并模式，需要添加额外的变量到 VARIABLES 列表中，确保能读取
    current_variables = list(VARIABLES) # 复制一份，避免修改原始 VARIABLES 列表
    if IS_MERGE_MODE:
        for var in ADDITIONAL_MERGE_VARS:
            if var not in current_variables:
                current_variables.append(var)
    
    # --- 3. 检查 Data 文件 ---
    if not os.path.exists(DATA_FILE):
        raise FileNotFoundError(f"没有找到真实数据文件：{DATA_FILE}")
    print("真实数据文件：", DATA_FILE)

    # --- 4. 读取 MC (signal) ---
    print("开始读取 MC 文件...")
    mc_list = []

    for file_path in MC_FILES:
        if IS_MERGE_MODE:
            file_event_candidates = {} # {EVENTNUMBER: []}
            for tree_name_candidate in SELECTED_TREE_NAMES:
                df_mc_tree = load_root(file_path, current_variables, tree_name=tree_name_candidate, max_entries=200000)
                if not df_mc_tree.empty:
                    df_mc_tree["label"] = 1
                    df_mc_tree["source_file"] = os.path.basename(file_path)
                    df_mc_tree["original_tree"] = tree_name_candidate # 记录来自哪个树
                    
                    for _, row in df_mc_tree.iterrows():
                        event_num = row["EVENTNUMBER"]
                        if event_num not in file_event_candidates:
                            file_event_candidates[event_num] = []
                        file_event_candidates[event_num].append(row)
            
            # 为当前文件选择最佳候选事件
            selected_mc_rows = []
            for event_num, candidates in file_event_candidates.items():
                if not candidates:
                    continue
                # 选择 W_M 最接近 TRUE_W_M_FOR_SELECTION 的事件
                best_candidate = min(candidates, key=lambda x: abs(x["W_M"] - TRUE_W_M_FOR_SELECTION))
                selected_mc_rows.append(best_candidate.drop(["EVENTNUMBER", "original_tree"])) # 移除辅助变量
            
            if selected_mc_rows:
                mc_list.append(pd.DataFrame(selected_mc_rows))
                print(f"读取并合并MC文件: {os.path.basename(file_path)} 事件数 = {len(selected_mc_rows)}")
            else:
                print(f"跳过文件: {os.path.basename(file_path)} (读取失败或合并后为空)")

        else: # 非合并模式
            df_mc = load_root(file_path, current_variables, tree_name=SELECTED_TREE_NAMES[0], max_entries=200000)
            if not df_mc.empty:
                df_mc["label"] = 1
                df_mc["source_file"] = os.path.basename(file_path)
                mc_list.append(df_mc)
                print(f"读取MC文件: {os.path.basename(file_path)} 事件数 = {len(df_mc)}")
            else:
                print(f"跳过文件: {os.path.basename(file_path)} (读取失败或为空)")

    if not mc_list:
        raise RuntimeError("所有 MC 文件都读取失败了，请检查 TREE_NAME 和变量名。")
    
    mc = pd.concat(mc_list, ignore_index=True)
    print(f"总计读取 MC 事件数 = {len(mc)}")


    # --- 5. 读取真实数据 (background) ---
    print("开始读取真实数据文件...")
    data_df_list = []
    # 对于真实数据，如果选择了合并模式，同样尝试读取多个树
    for tree_name_candidate in SELECTED_TREE_NAMES:
        df_data_tree = load_root(DATA_FILE, current_variables, tree_name=tree_name_candidate, max_entries=200000)
        if not df_data_tree.empty:
            data_df_list.append(df_data_tree)
            print(f"读取Data文件 (树: {tree_name_candidate}) 事件数 = {len(df_data_tree)}")
        else:
            print(f"跳过Data文件 (树: {tree_name_candidate}) (读取失败或为空)")

    data = pd.concat([df for df in data_df_list if not df.empty], ignore_index=True)
    if data.empty:
         raise RuntimeError("真实数据文件读取失败或为空。")

    data["label"] = 0
    data["source_file"] = os.path.basename(DATA_FILE)
    
    # 如果Data也是合并模式，且有EVENTNUMBER，则也需要进行合并筛选（通常Data没有TRUE_W_M，这里简化）
    # 如果你的Data合并有特殊逻辑，需要在这里添加
    if IS_MERGE_MODE and "EVENTNUMBER" in data.columns and len(SELECTED_TREE_NAMES) > 1:
        print("Data 也在合并模式下，尝试根据 EVENTNUMBER 和 W_M 进行去重/选择...")
        # 针对Data的合并逻辑可能需要更具体，这里假设如果有重复EVENTNUMBER，只取W_M最接近真实的事件
        # 实际中Data的去重逻辑可能更复杂，例如基于 Run/Lumi/Event number。
        # 这里为了演示，我们假设对于重复的 EVENTNUMBER，也选择 W_M 最接近真实值的。
        # 请根据你的实际数据和物理需求调整此部分！
        
        # 假设 Data 的合并逻辑和 MC 类似 (选择W_M最接近 TRUE_W_M_FOR_SELECTION 的事件)
        # 注意：这在物理上可能不完全合理，因为Data通常不应该假设有“真实”的W_M
        # 更常见的做法是对Data只读取一个默认树，或采用其他启发式合并策略。
        
        data_processed_rows = []
        # Group by EVENTNUMBER and select the best candidate
        for event_num, group in data.groupby("EVENTNUMBER"):
            if len(group) > 1:
                best_candidate = group.iloc[
                    (group["W_M"] - TRUE_W_M_FOR_SELECTION).abs().argsort()[:1]
                ]
                data_processed_rows.append(best_candidate.drop("EVENTNUMBER", axis=1)) # 移除辅助变量
            else:
                data_processed_rows.append(group.drop("EVENTNUMBER", axis=1))
        
        data = pd.concat(data_processed_rows, ignore_index=True)
        print(f"真实数据合并及去重后事件数 = {len(data)}")

    print(f"真实数据事件数 = {len(data)}")


    # --- 6. 合并数据集 ---
    dataset = pd.concat([mc, data], ignore_index=True)
    print(f"合并后总事件数 = {len(dataset)}")
    print(f"Signal 事件数 = {(dataset['label'] == 1).sum()}")
    print(f"Background 事件数 = {(dataset['label'] == 0).sum()}")

    # --- 7. 应用 Pre-selection ---
    print("--- 应用 Pre-selection ---")
    
    dataset, pre_sel_sig_eff, pre_sel_bkg_eff = apply_pre_selection(dataset)
    
    print(f"Pre-selection 信号通过率: {pre_sel_sig_eff:.4f}")
    print(f"Pre-selection 本底通过率: {pre_sel_bkg_eff:.4f}")

    # 检查经过预选后是否还有数据
    if dataset.empty or (dataset['label'] == 1).sum() == 0 or (dataset['label'] == 0).sum() == 0:
        print("预选后信号或本底数据为空，无法进行 GBDT 训练。程序退出。")
        exit()

    X = dataset[VARIABLES] # 确保 X 使用的是统一的 VARIABLES 列表
    y = dataset["label"]

    # --- 8. 划分训练集 / 测试集 ---
    # 使用 stratify=y 确保训练集和测试集中的类别比例与原始数据集一致
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    print(f"训练集大小 = {len(X_train)}")
    print(f"测试集大小 = {len(X_test)}")

    # ==================================================
    # 9. 训练级联 GBDT 模型
    # ==================================================
    print("开始训练级联 GBDT 模型...")
    
    # 存储每个 GBDT 的 ROC 曲线数据和 AUC
    all_gbdt_roc_data = {}
    
    # 用于记录筛选效率
    selection_efficiency_log = {
        "Pre-selection": {
            "signal_eff": pre_sel_sig_eff,
            "background_eff": pre_sel_bkg_eff
        }
    }
    
    # 记录每个GBDT的初始事件数 (测试集，用于计算效率)
    # 这里的 initial_test_sig_count 和 initial_test_bkg_count 应该是 Pre-selection 后的总量
    initial_test_sig_count_total = (y_test == 1).sum()
    initial_test_bkg_count_total = (y_test == 0).sum()
    
    # 用于计算累积通过率的测试集 (从 Pre-selection 后的测试集开始)
    cumulative_X_test = X_test.copy() 
    cumulative_y_test = y_test.copy()

    for i in range(NUM_GBDTS):
        stage_name = f"GBDT_{i+1}"
        features_for_this_stage = FEATURES_FOR_GBDT_STAGES[i]
        
        # 确保当前阶段的特征都在 X_train 和 X_test 中
        missing_feats_train = [f for f in features_for_this_stage if f not in X_train.columns]
        missing_feats_test = [f for f in features_for_this_stage if f not in X_test.columns]
        if missing_feats_train or missing_feats_test:
            print(f"警告: {stage_name} 阶段缺少特征。训练集缺少: {missing_feats_train}, 测试集缺少: {missing_feats_test}. 跳过此阶段。")
            # 可以选择跳过此阶段或用默认值填充
            continue

        print(f"--- 训练 {stage_name} (使用特征: {features_for_this_stage}) ---")

        # 对当前GBDT阶段的特征进行标准化
        scaler = StandardScaler()
        # 训练Scaler应该在完整的训练集上
        X_train_scaled = scaler.fit_transform(X_train[features_for_this_stage])
        X_test_scaled = scaler.transform(X_test[features_for_this_stage]) # 对原始测试集进行缩放
        SCALERS.append(scaler)

        model = GradientBoostingClassifier(random_state=42, n_estimators=100, learning_rate=0.1, max_depth=3) # 可以调整参数
        model.fit(X_train_scaled, y_train) # 在完整的训练集上训练
        GBDT_MODELS.append(model)

        # 评估当前 GBDT (在原始测试集上)
        y_score_test = model.predict_proba(X_test_scaled)[:, 1]
        auc_test = roc_auc_score(y_test, y_score_test)
        
        y_score_train = model.predict_proba(scaler.transform(X_train[features_for_this_stage]))[:, 1]
        auc_train = roc_auc_score(y_train, y_score_train)

        print(f"  {stage_name} 训练集 AUC = {auc_train:.6f}")
        print(f"  {stage_name} 测试集 AUC = {auc_test:.6f}")
        
        # 存储当前GBDT的ROC数据和AUC（用原始测试集）
        fpr, tpr, _ = roc_curve(y_test, y_score_test)
        all_gbdt_roc_data[stage_name] = {"fpr": fpr, "tpr": tpr, "auc": auc_test}

        # ==================================================
        # 计算级联筛选效率
        # 这里我们使用一个简单的阈值（例如0.5）来模拟级联筛选
        # 并在cumulative_X_test上进行，以反映逐步筛选的效果
        # ==================================================
        
        # 对当前的 cumulative_X_test 应用当前GBDT
        if not cumulative_X_test.empty:
            cumulative_X_test_scaled_for_this_stage = SCALERS[i].transform(cumulative_X_test[features_for_this_stage])
            cumulative_y_score = GBDT_MODELS[i].predict_proba(cumulative_X_test_scaled_for_this_stage)[:, 1]
            
            # 选择一个阈值进行切割
            bdt_threshold = 0.5 # 示例阈值，实际应用中会优化
            
            # 筛选通过的事件
            passed_indices_in_cumulative = (cumulative_y_score >= bdt_threshold)
            
            # 更新通过事件
            cumulative_X_test = cumulative_X_test[passed_indices_in_cumulative]
            cumulative_y_test = cumulative_y_test[passed_indices_in_cumulative]
        else:
            # 如果cumulative_X_test已经为空，则后续效率均为0
            cumulative_y_test = pd.Series([]) # 设置为空
            
        # 计算通过率 (相对于最初的 Pre-selection 后的测试集总量)
        current_sig_count = (cumulative_y_test == 1).sum()
        current_bkg_count = (cumulative_y_test == 0).sum()
        
        sig_eff_stage = current_sig_count / initial_test_sig_count_total if initial_test_sig_count_total > 0 else 0
        bkg_eff_stage = current_bkg_count / initial_test_bkg_count_total if initial_test_bkg_count_total > 0 else 0
        
        selection_efficiency_log[stage_name] = {
            "signal_eff": sig_eff_stage,
            "background_eff": bkg_eff_stage
        }
        print(f"  {stage_name} (累积) 信号通过率: {sig_eff_stage:.4f}")
        print(f"  {stage_name} (累积) 本底通过率: {bkg_eff_stage:.4f}")
        
    # ==================================================
    # 10. 输出结果
    # ==================================================
    print("--- 输出结果 ---")

    # --- 10.1 绘制 ROC 曲线 (每个 GBDT 独立) ---
    plt.figure(figsize=(8, 8))
    for stage_name, data in all_gbdt_roc_data.items():
        plt.plot(data["fpr"], data["tpr"], label=f"{stage_name} AUC = {data['auc']:.4f}")
    
    # 绘制对角线作为随机分类器的参考
    plt.plot([0, 1], [0, 1], linestyle="--", color='gray', label="Random Classifier")
    
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curves for Cascade GBDT Stages (evaluated on original test set)")
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    roc_curve_filename = os.path.join(OUTPUT_DIR, "cascade_roc_curves.png")
    plt.savefig(roc_curve_filename, dpi=150)
    plt.close()
    print(f" - ROC 曲线图已保存至: {roc_curve_filename}")

    # --- 10.2 绘制 BDT 分数分布图 (Signal vs Background) 和特征重要性图 ---
    for i in range(NUM_GBDTS):
        stage_name = f"GBDT_{i+1}"
        model = GBDT_MODELS[i]
        scaler = SCALERS[i]
        features_for_this_stage = FEATURES_FOR_GBDT_STAGES[i]

        # 确保特征存在
        if not all(f in X_test.columns for f in features_for_this_stage):
            print(f"跳过 {stage_name} 的绘图和重要性分析，因为缺少必要的特征。")
            continue

        # 使用原始测试集进行预测
        X_test_scaled_for_plot = scaler.transform(X_test[features_for_this_stage])
        y_score_plot = model.predict_proba(X_test_scaled_for_plot)[:, 1]

        score_sig = y_score_plot[y_test == 1]
        score_bkg = y_score_plot[y_test == 0]

        bins = np.linspace(0, 1, 51) # BDT 分数通常在 0 到 1 之间

        # BDT Response Distribution (带对数纵轴)
        plt.figure(figsize=(8, 6))
        plt.hist(score_sig, bins=bins, alpha=0.5, label='Signal (MC)', color='blue', histtype='stepfilled')
        plt.hist(score_bkg, bins=bins, alpha=0.5, label='Background (Data)', color='red', histtype='stepfilled')
        plt.yscale('log')
        plt.xlabel(f'{stage_name} Score')
        plt.ylabel('Events (Log Scale)')
        plt.title(f'{stage_name} Response Distribution')
        plt.legend(loc='upper center')
        plt.grid(axis='y', alpha=0.3)
        bdt_dist_filename = os.path.join(OUTPUT_DIR, f"{stage_name}_bdt_distribution.png")
        plt.savefig(bdt_dist_filename, dpi=150)
        plt.close()
        print(f" - {stage_name} BDT 分布图已保存至: {bdt_dist_filename}")

        # Normalized BDT Shape Comparison
        plt.figure(figsize=(8, 6))
        # histtype='step' 更好用于归一化图
        plt.hist(score_sig, bins=bins, alpha=0.7, label='Signal (Normalized)', color='blue', density=True, histtype='step', lw=2)
        plt.hist(score_bkg, bins=bins, alpha=0.7, label='Background (Normalized)', color='red', density=True, histtype='step', lw=2)
        plt.xlabel(f'{stage_name} Score')
        plt.ylabel('Probability Density')
        plt.title(f'Normalized {stage_name} Shape Comparison')
        plt.legend(loc='upper center')
        bdt_shape_filename = os.path.join(OUTPUT_DIR, f"{stage_name}_bdt_shape_comparison.png")
        plt.savefig(bdt_shape_filename, dpi=150)
        plt.close()
        print(f" - {stage_name} BDT 形状对比图已保存至: {bdt_shape_filename}")
        
        # --- Feature Importance 图 ---
        importances = model.feature_importances_
        feature_names = features_for_this_stage
        
        # 按照重要性排序
        feature_importance_ranking = sorted(
            zip(feature_names, importances),
            key=lambda x: x[1],
            reverse=True
        )
        
        # 打印到终端
        print(f"{stage_name} 区分度较高的变量（按 importance 排序）:")
        for var, score in feature_importance_ranking:
            print(f"  {var:25s} {score:.6f}")

        # 绘制图
        plt.figure(figsize=(10, max(5, len(feature_names) * 0.4))) # 动态调整图高
        sorted_features = [f[0] for f in feature_importance_ranking]
        sorted_importances = [f[1] for f in feature_importance_ranking]
        
        plt.barh(sorted_features, sorted_importances)
        plt.xlabel("Importance")
        plt.title(f"{stage_name} Feature Importance")
        plt.gca().invert_yaxis() # 重要性最高的在最上面
        plt.tight_layout()
        feature_importance_filename = os.path.join(OUTPUT_DIR, f"{stage_name}_feature_importance.png")
        plt.savefig(feature_importance_filename, dpi=150)
        plt.close()
        print(f" - {stage_name} 特征重要性图已保存至: {feature_importance_filename}")
        
        # 保存文字结果 (特征重要性)
        feature_ranking_txt_filename = os.path.join(OUTPUT_DIR, f"{stage_name}_feature_ranking.txt")
        with open(feature_ranking_txt_filename, "w", encoding="utf-8") as f:
            f.write(f"--- {stage_name} Feature Ranking ---")
            f.write(f"AUC = {all_gbdt_roc_data[stage_name]['auc']:.6f}")
            f.write("Feature ranking:")
            for var, score in feature_importance_ranking:
                f.write(f"{var:25s} {score:.6f}")
        print(f" - {stage_name} 特征重要性排名已保存至: {feature_ranking_txt_filename}")

     # ==================================================
    # 11. 级联 GBDT 联合筛选结果
    # ==================================================
    print("--- 计算并绘制级联 GBDT 联合筛选结果 ---")

    # 重新初始化测试集，用于计算联合分数
    X_test_combined = X_test.copy()
    y_test_combined = y_test.copy()

    # 存储每个阶段的累积分数
    all_cumulative_scores = np.ones(len(X_test_combined)) # 初始化为1，代表未筛选

    for i in range(NUM_GBDTS):
        stage_name = f"GBDT_{i+1}"
        model = GBDT_MODELS[i]
        scaler = SCALERS[i]
        features_for_this_stage = FEATURES_FOR_GBDT_STAGES[i]

        # 确保特征存在
        if not all(f in X_test_combined.columns for f in features_for_this_stage):
            print(f"跳过 {stage_name} 的联合分数计算，因为缺少必要的特征。")
            continue

        # 对整个测试集进行预测
        X_scaled_for_stage = scaler.transform(X_test_combined[features_for_this_stage])
        stage_scores = model.predict_proba(X_scaled_for_stage)[:, 1]

        # 联合分数：这里可以有多种策略。
        # 常见策略是取所有 GBDT 分数的乘积，或对数加权和，或直接取最后一个 GBDT 的分数。
        # 暂时我们采用乘积，因为分数越高代表信号可能性越高，乘积会增强这种效应。
        #all_cumulative_scores *= stage_scores
        # 另一种常见方法是只用最后一个 GBDT 作为最终判别器，或者对其进行加权。
        all_cumulative_scores = stage_scores # 如果只关心最后一个GBDT的输出
        
    final_bdt_score = all_cumulative_scores

    # --- 11.1 联合 GBDT 的 ROC 曲线 ---
    fpr_final, tpr_final, thresholds_final = roc_curve(y_test_combined, final_bdt_score)
    auc_final = roc_auc_score(y_test_combined, final_bdt_score)

    plt.figure(figsize=(8, 8))
    plt.plot(fpr_final, tpr_final, label=f"Combined GBDT AUC = {auc_final:.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--", color='gray', label="Random Classifier")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Combined GBDT ROC Curve")
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    roc_curve_combined_filename = os.path.join(OUTPUT_DIR, "combined_gbdt_roc_curve.png")
    plt.savefig(roc_curve_combined_filename, dpi=150)
    plt.close()
    print(f" - 联合 GBDT ROC 曲线图已保存至: {roc_curve_combined_filename}")

    # --- 11.2 联合 BDT 分数分布 ---
    score_sig_final = final_bdt_score[y_test_combined == 1]
    score_bkg_final = final_bdt_score[y_test_combined == 0]

    bins_final = np.linspace(0, 1, 51) # 联合分数可能不再是0-1，但我们可以归一化或者直接使用。
                                       # 如果分数是乘积，它的范围会更小，可能需要调整 bins。
                                       # 为了绘图方便，暂时仍使用 0-1 的区间，但需注意解读。
    
    # 调整 bin 范围以适应乘积分数
    min_score = np.min(final_bdt_score) if len(final_bdt_score) > 0 else 0
    max_score = np.max(final_bdt_score) if len(final_bdt_score) > 0 else 1
    bins_final = np.linspace(min_score, max_score if max_score > min_score else 1, 51)
    if len(bins_final) <= 1: bins_final = np.linspace(0,1,51) # Fallback if all scores are same

    plt.figure(figsize=(8, 6))
    plt.hist(score_sig_final, bins=bins_final, alpha=0.5, label='Signal (MC)', color='blue', histtype='stepfilled')
    plt.hist(score_bkg_final, bins=bins_final, alpha=0.5, label='Background (Data)', color='red', histtype='stepfilled')
    plt.yscale('log')
    plt.xlabel("Combined GBDT Score (Product)")
    plt.ylabel('Events (Log Scale)')
    plt.title("Combined GBDT Score Distribution")
    plt.legend(loc='upper center')
    plt.grid(axis='y', alpha=0.3)
    bdt_dist_combined_filename = os.path.join(OUTPUT_DIR, "combined_gbdt_distribution.png")
    plt.savefig(bdt_dist_combined_filename, dpi=150)
    plt.close()
    print(f" - 联合 GBDT 分布图已保存至: {bdt_dist_combined_filename}")

    # --- 11.3 信号效率 vs 本底抑制率 (Signal Efficiency vs Background Rejection) ---
    # 这本质上是ROC曲线的另一种展现形式，y轴是TPR，x轴是1-FPR（本底抑制率）
    # 但通常物理分析中更关注在某个高信号效率下的本底抑制率
    
    # 计算不同 BDT 阈值下的信号效率和本底效率
    sig_efficiencies = []
    bkg_rejections = []
    
    # 对阈值进行排序，从高到低
    thresholds_final_sorted = np.sort(np.unique(final_bdt_score))[::-1]
    
    for threshold in thresholds_final_sorted:
        selected_sig = final_bdt_score[(y_test_combined == 1) & (final_bdt_score >= threshold)]
        selected_bkg = final_bdt_score[(y_test_combined == 0) & (final_bdt_score >= threshold)]
        
        sig_eff = len(selected_sig) / initial_test_sig_count_total if initial_test_sig_count_total > 0 else 0
        bkg_eff = len(selected_bkg) / initial_test_bkg_count_total if initial_test_bkg_count_total > 0 else 0
        
        sig_efficiencies.append(sig_eff)
        bkg_rejections.append(1 - bkg_eff) # 本底抑制率 = 1 - 本底通过率

    plt.figure(figsize=(8, 8))
    plt.plot(sig_efficiencies, bkg_rejections, label="Combined GBDT")
    plt.xlabel("Signal Efficiency (TPR)")
    plt.ylabel("Background Rejection (1 - FPR)")
    plt.title("Combined GBDT: Signal Efficiency vs Background Rejection")
    plt.legend(loc='lower left')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    eff_rej_combined_filename = os.path.join(OUTPUT_DIR, "combined_gbdt_eff_rejection.png")
    plt.savefig(eff_rej_combined_filename, dpi=150)
    plt.close()
    print(f" - 联合 GBDT 效率 vs 抑制率图已保存至: {eff_rej_combined_filename}")

    # --- 11.4 级联筛选效率条形图 ---
    stage_labels = list(selection_efficiency_log.keys())
    sig_effs = [selection_efficiency_log[s]["signal_eff"] for s in stage_labels]
    bkg_effs = [selection_efficiency_log[s]["background_eff"] for s in stage_labels]

    x = np.arange(len(stage_labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width/2, sig_effs, width, label='Signal Efficiency', color='skyblue')
    rects2 = ax.bar(x + width/2, bkg_effs, width, label='Background Efficiency', color='salmon')

    ax.set_ylabel('Cumulative Efficiency')
    ax.set_title('Cumulative Selection Efficiency per Stage')
    ax.set_xticks(x)
    ax.set_xticklabels(stage_labels, rotation=45, ha="right")
    ax.legend()
    ax.set_ylim(0, 1.05) # 效率通常在0-1之间

    fig.tight_layout()
    cascade_eff_bar_filename = os.path.join(OUTPUT_DIR, "cascade_cumulative_efficiency_bar_chart.png")
    plt.savefig(cascade_eff_bar_filename, dpi=150)
    plt.close()
    print(f" - 级联累积筛选效率条形图已保存至: {cascade_eff_bar_filename}")


    # --- 10.3 筛选效率报告 (更新为包含联合结果) ---
    efficiency_report_filename = os.path.join(OUTPUT_DIR, "cascade_selection_efficiency.txt")
    with open(efficiency_report_filename, "w", encoding="utf-8") as f:
        f.write("--- Cascade GBDT Selection Efficiency Report ---")
        f.write(f"Pre-selection 后测试集初始信号事件数: {initial_test_sig_count_total}")
        f.write(f"Pre-selection 后测试集初始本底事件数: {initial_test_bkg_count_total}")
        
        for step, eff_data in selection_efficiency_log.items():
            f.write(f"{step}:")
            f.write(f"  信号通过率 (相对于 Pre-selection 后测试集初始值): {eff_data['signal_eff']:.6f}")
            f.write(f"  本底通过率 (相对于 Pre-selection 后测试集初始值): {eff_data['background_eff']:.6f}")

        f.write(f"--- Combined GBDT Score Evaluation (Final Stage) ---")
        f.write(f"Combined GBDT AUC on test set: {auc_final:.6f}")
    print(f" - 级联筛选效率报告已保存至: {efficiency_report_filename}")

    print(f"所有输出文件保存在：{OUTPUT_DIR}")
    print("运行完成。")