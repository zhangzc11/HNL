import uproot
import pandas as pd
import numpy as np
import awkward as ak
import os
import re
import sys
import errno
import ROOT
import array
from concurrent.futures import ProcessPoolExecutor

# ROOT Batch Mode
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0) # 不显示统计框

# ----------------- 配置 -----------------
MC_DIR = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"
DATA_DIR = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/"
LOCAL_DATA_FILE = "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
LOCAL_DATA_FILE = os.path.join(DATA_DIR, LOCAL_DATA_FILE)
CSV_FILE = "./csv/top30_ks_branches.csv"
CANDIDATE_DIRS = ["myTupleOS1J", "myTupleOS2J", "myTupleSS1J", "myTupleSS2J"]
FILE_PATTERN = r"_(\d+GeV)Ctau(\d+ps)_.*skim\.root"

# ----------------- 辅助函数 -----------------

def parse_filename(filename):
    """
    解析文件名获取质量(Mass)和寿命(Lifetime/Tau)信息。
    文件名示例: 42912010_5GeVCtau0ps_options_20260228_111712_skim.root
    返回: (mass_str, tau_str) 例如 ("5GeV", "0ps")
    """
    match = re.search(FILE_PATTERN, filename)
    if match:
        return match.group(1), match.group(2)
    return None, None

def get_tree(file_path):
    """
    打开ROOT文件并获取第一个有效的TTree。
    """
    try:
        with uproot.open(file_path) as f:
            for dirname in CANDIDATE_DIRS:
                path = f"{dirname}/DecayTree"
                if path in f:
                    return f[path]
            # Fallback
            if "DecayTree" in f:
                return f["DecayTree"]
    except Exception as e:
        print(f"无法打开文件 {file_path}: {e}")
    return None

def load_data(file_paths_and_labels, branch_name):
    """
    加载指定文件列表中某分支的数据。
    file_paths_and_labels: list of (path, label) tuples.
    返回 (字典 {label: numpy_array}, found_directory_name)
    """
    data_dict = {}
    found_dir_name = None
    
    # 强制要求输入为列表
    # 如果只传了一个路径字符串，包装一下
    if isinstance(file_paths_and_labels, str):
        file_paths_and_labels = [(file_paths_and_labels, "Data")]

    if not isinstance(file_paths_and_labels, list):
         print(f"Error: expected list, got {type(file_paths_and_labels)}")
         return {}, None

    for item in file_paths_and_labels:
        if isinstance(item, tuple) or isinstance(item, list):
            path, label = item
        else:
            print(f"Error: Item {item} is not a tuple")
            continue

        try:
             with uproot.open(path) as f:
                tree = None
                for dirname in CANDIDATE_DIRS:
                    if f"{dirname}/DecayTree" in f:
                        tree = f[f"{dirname}/DecayTree"]
                        if found_dir_name is None:
                            found_dir_name = dirname
                        break
                if tree is None and "DecayTree" in f:
                    tree = f["DecayTree"]
                    if found_dir_name is None:
                        found_dir_name = "DecayTree"
                
                if tree is None:
                    continue
                
                try:
                    # 使用 array(library='ak') 读取
                    # 注意 branch_name 可能需要精确匹配
                    if branch_name not in tree:
                        continue
                        
                    arr = tree[branch_name].array(library="ak")
                    # 展平 Jagged Array
                    flat_arr = ak.to_numpy(ak.flatten(arr, axis=None))
                    # 将 NaN 替换为 0
                    flat_arr = np.nan_to_num(flat_arr, nan=0)
                    data_dict[label] = flat_arr
                except Exception as e:
                     print(f"Warning: Failed to read branch {branch_name} in {os.path.basename(path)}: {e}")
                     continue
                     
        except Exception as e:
            print(f"读取数据失败 {os.path.basename(path)}: {e}")

    return data_dict, found_dir_name

def fill_th1(name, title, data_arr, bins, color, style=None):
    """
    创建一个 TH1D 并用 numpy 数组填充。
    """
    if data_arr is None or len(data_arr) == 0:
        return None
        
    # data_arr = data_arr[np.isfinite(data_arr)]
    # if len(data_arr) == 0: return None
        
    # 转换为标准的 float64 C-contiguous 数组，以便 ROOT 使用
    data_arr = np.ascontiguousarray(data_arr, dtype=np.float64)
    weights = np.ones(len(data_arr), dtype=np.float64)
    
    # 创建直方图
    # bins: array of edges
    h = ROOT.TH1D(name, title, len(bins)-1, array.array('d', bins))
    h.Sumw2()
    
    # 快速填充 (FillN)
    # 注意: PyROOT 中可以直接传递 numpy 数组给 FillN，但在某些旧版本中可能需要 array.array
    # 这里尝试直接传递
    h.FillN(len(data_arr), data_arr, weights)
    
    # 归一化 (Integtral -> 1)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
        
    # 样式设置
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetStats(0)
    
    if style == 'filled':
        # Data 样式 (灰色填充)
        h.SetFillColor(ROOT.kGray)
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)
        h.SetFillStyle(1001) # Solid
    
    return h

def plot_branch_task(args):
    """
    Multiprocessing worker function.
    args: (branch_name, data_info, mc_file_list, output_dir, file_prefix, title_info)
    """
    branch_name, data_path, mc_files_info_tuples, output_dir, file_prefix, title_info = args
    
    # print(f"Processing {branch_name} ...")
    
    # 1. Load Data
    # Reuse load_data logic but simplified for single call
    data_dict, data_dir = load_data([(data_path, "Data")], branch_name)
    data_arr = data_dict.get("Data", None)
    
    # 2. Load MCs
    # mc_files_info_tuples: list of (path, label)
    mc_dict, mc_dir = load_data(mc_files_info_tuples, branch_name)
    
    # Prioritize data_dir, if not found use mc_dir
    active_dir = data_dir if data_dir else mc_dir
    
    if data_arr is None and not mc_dict:
        return f"Skipped {branch_name} (No data)"

    # 3. Determine Range and Bins
    all_arrays = []
    if data_arr is not None and len(data_arr) > 0:
        data_arr = data_arr[np.isfinite(data_arr)]
        if len(data_arr) > 0: all_arrays.append(data_arr)
        
    valid_mc = {}
    for k, v in mc_dict.items():
        v = v[np.isfinite(v)]
        if len(v) > 0:
            valid_mc[k] = v
            all_arrays.append(v)
            
    if not all_arrays:
        return f"Skipped {branch_name} (Empty arrays)"
        
    combined = np.concatenate(all_arrays)
    min_val = np.min(combined)
    max_val = np.max(combined)
    
    # 50 bins
    bins = np.linspace(min_val, max_val, 51) # 51 edges for 50 bins
    
    # 4. Plot with ROOT
    canvas = ROOT.TCanvas(f"c_{branch_name}", f"c_{branch_name}", 800, 600)
    
    max_y = 0
    
    # Create Data Hist
    h_data = None
    if data_arr is not None and len(data_arr) > 0:
        h_data = fill_th1(f"h_data_{branch_name}", "", data_arr, bins, ROOT.kBlack, style='filled')
        if h_data:
            if h_data.GetMaximum() > max_y: max_y = h_data.GetMaximum()
            
    # Create MC Hists
    h_mcs = []
    # ROOT Colors: kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, ...
    root_colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1, ROOT.kOrange+1, ROOT.kViolet, ROOT.kAzure+1, ROOT.kPink, ROOT.kSpring]
    
    legend = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    
    # Sort MCs by label for consistent legend
    # Use the order from mc_files_info_tuples which is sorted numerically in main
    sorted_mc_items = []
    seen_labels = set()
    for _, label in mc_files_info_tuples:
        if label in valid_mc and label not in seen_labels:
            sorted_mc_items.append((label, valid_mc[label]))
            seen_labels.add(label)
    
    for i, (label, arr) in enumerate(sorted_mc_items):
        color = root_colors[i % len(root_colors)]
        h = fill_th1(f"h_mc_{i}_{branch_name}", "", arr, bins, color)
        if h:
            h_mcs.append(h)
            legend.AddEntry(h, label, "l")
            if h.GetMaximum() > max_y: max_y = h.GetMaximum()

    if h_data:
        legend.AddEntry(h_data, "Data", "f")

    # Draw
    # Draw Data first if filled? No, draw axes first.
    # Create a frame histogram for axes
    frame = ROOT.TH1D(f"frame_{branch_name}", f"{branch_name}", len(bins)-1, array.array('d', bins))
    frame.SetStats(0)
    frame.GetXaxis().SetTitle(branch_name)
    frame.GetYaxis().SetTitle("Normalized Events")
    frame.SetMaximum(max_y * 1.3) # Add MORE headroom for text
    frame.Draw("AXIS")
    
    # Draw Data (filled)
    if h_data:
        h_data.Draw("HIST SAME")
        
    # Draw MCs (lines)
    for h in h_mcs:
        h.Draw("HIST SAME")
        
    legend.Draw()

    # Add Info Text
    # Must draw AFTER other histograms to be visible on top
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextAlign(13) # Align Top-Left
    
    info_text = f"{title_info}"
    if active_dir:
        dir_text = f"Dir: {active_dir}"
        latex.DrawLatex(0.15, 0.88, info_text)
        latex.DrawLatex(0.15, 0.83, dir_text)
    else:
        latex.DrawLatex(0.15, 0.88, info_text)
    
    # Save
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError:
            pass
            
    filename = f"{file_prefix}_{branch_name}.png"
    filename = re.sub(r'[\\/*?:"<>|]', "", filename)
    save_path = os.path.join(output_dir, filename)
    
    canvas.SaveAs(save_path)
    # prevent memory leak in ROOT
    if h_data: h_data.Delete()
    for h in h_mcs: h.Delete()
    frame.Delete()
    canvas.Close()
    
    return f"Saved {save_path}"

def main():
    print("=== 开始作图脚本 (ROOT并行版) ===")
    
    # 1. 检查必要文件
    if not os.path.exists(LOCAL_DATA_FILE):
        print(f"错误: 找不到数据文件 {LOCAL_DATA_FILE}")
        return
    if not os.path.exists(CSV_FILE):
        print(f"错误: 找不到 CSV 文件 {CSV_FILE}")
        return
        
    # 2. 读取 CSV 获取前 30 个分支
    try:
        df_branches = pd.read_csv(CSV_FILE)
        top_branches = df_branches['Branch'].head(30).tolist()
        print(f"已加载 {len(top_branches)} 个分支进行作图。")
    except Exception as e:
        print(f"读取 CSV 失败: {e}")
        return

    # 3. 扫描 MC 目录
    mc_files = [] 
    if not os.path.exists(MC_DIR):
        print(f"错误: MC 目录 {MC_DIR} 不存在")
        return
        
    print(f"正在扫描 MC 目录: {MC_DIR}")
    for f in os.listdir(MC_DIR):
        full_path = os.path.join(MC_DIR, f)
        if not os.path.isfile(full_path):
            continue
        m, t = parse_filename(f)
        if m and t:
            mc_files.append({'path': full_path, 'mass': m, 'tau': t, 'filename': f})
    
    print(f"找到 {len(mc_files)} 个符合条件的 MC 文件。")
    
    # 4. 模式选择 (argv 或 input)
    mode = None
    target_val = None
    
    if len(sys.argv) >= 3:
        mode = sys.argv[1]
        target_val = sys.argv[2]
        print(f"使用命令行参数: Mode={mode}, Value={target_val}")
    else:
        print("\n请选择作图模式 (输入 1 或 2):")
        print("1. 固定质量 (Mass)，变化寿命 (Lifetime)")
        print("2. 固定寿命 (Lifetime)，变化质量 (Mass)")
        try:
            mode = input("请输入: ").strip()
        except EOFError:
             return 

    # Handle input arguments logic here to get available values 
    if mode == '1': # Fixed Mass
         available_masses = sorted(list(set(f['mass'] for f in mc_files)))
         if not target_val:
             print(f"\n检测到的可选质量值 (Mass): {available_masses}")
             target_val = input("请输入目标质量值 (例如 5GeV): ").strip()
         
         target_files = [f for f in mc_files if f['mass'] == target_val]
         # Sort by tau
         try:
            target_files.sort(key=lambda x: int(re.search(r'\d+', x['tau']).group()))
         except:
            target_files.sort(key=lambda x: x['tau'])
            
         subdir_name = f"M{target_val}"
         file_prefix = f"M{target_val}"
         title_context = f"Fixed Mass = {target_val}"
         
    elif mode == '2': # Fixed Tau
         available_taus = sorted(list(set(f['tau'] for f in mc_files)))
         if not target_val:
             print(f"\n检测到的可选寿命值 (Lifetime/Tau): {available_taus}")
             target_val = input("请输入目标寿命值 (例如 0ps): ").strip()

         target_files = [f for f in mc_files if f['tau'] == target_val]
         # Sort by Mass
         try:
            target_files.sort(key=lambda x: int(re.search(r'\d+', x['mass']).group()))
         except:
            target_files.sort(key=lambda x: x['mass'])
            
         subdir_name = f"tau{target_val}"
         file_prefix = f"tau{target_val}"
         title_context = f"Fixed Tau = {target_val}"
    
    else:
        print("无效模式。")
        return

    if not target_files:
        print(f"未找到匹配文件。Val={target_val}")
        return

    # 5. Output Dir
    output_dir = os.path.join("./png", subdir_name)
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError:
            pass
    print(f"输出目录: {output_dir}")

    # 6. Prepare Parallel Tasks
    print(f"开始生成图片 (共 {len(top_branches)} 张)... 使用多进程并行。")
    
    # Construct task arguments
    tasks = []
    
    for branch in top_branches:
        # Prepare MC paths and labels
        mc_paths_labels = []
        for f in target_files:
            if mode == '1':
                label = f['tau']
            else:
                label = f['mass']
            mc_paths_labels.append((f['path'], label))
            
        args = (branch, LOCAL_DATA_FILE, mc_paths_labels, output_dir, file_prefix, title_context)
        tasks.append(args)
        
    # Execute in parallel
    # Determine number of workers
    n_workers = min(len(tasks), os.cpu_count() or 4)
    print(f"启动 {n_workers} 个工作进程...")
    
    # Use ProcessPoolExecutor to map tasks
    # Note: ProcessPoolExecutor works best with functions defined at module level
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(plot_branch_task, tasks))
        
    for res in results:
        pass

    print(f"\n全部完成！图片保存在: {os.path.abspath(output_dir)}")


if __name__ == "__main__":
    main()
