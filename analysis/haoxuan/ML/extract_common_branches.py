import uproot
import os
import re
import pandas as pd
import sys

# ----------------- 配置 -----------------
# 文件路径和模式配置
mc_dir = "/pool1/lhcb/zhicaiz/mc/2024/HNL_Wmumuqq/ntuple/ntuple_20260228/"
data_dir = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/"
local_data_file = "data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"

file_pattern = r".*skim\.root.*"
candidate_directories = ["myTupleOS1J", "myTupleOS2J", "myTupleSS1J", "myTupleSS2J"]

# ----------------- 辅助函数 -----------------

def find_files_with_regex(root_dir, pattern):
    """
    在给定根目录下查找匹配正则模式的文件。

    参数:
        root_dir (str): 根目录路径
        pattern (str): 正则表达式模式，用于匹配文件名

    返回:
        list: 匹配文件的完整路径列表
        """
    matched_files = []
    regex = re.compile(pattern)
    if not os.path.exists(root_dir):
        print(f"Warning: Directory {root_dir} does not exist or is not accessible.")
        return []
        
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if regex.search(filename):
                full_path = os.path.join(dirpath, filename)
                matched_files.append(full_path)
    return matched_files

def is_numeric_type_string(typename):
    typename = typename.lower()
    
    # Non-numeric types
    if 'string' in typename or 'char*' in typename or 'char *' in typename:
        return False
    numeric_keywords = ['int', 'float', 'double', 'short', 'long', 'bool', 'unsigned', 'char', 'byte']
    if any(k in typename for k in numeric_keywords):
        return True
    return False 

def is_physical_meaningful(name):
    """
    判断分支名是否具有分析上的物理意义，排除常见的元数据/索引/ID 字段。

    规则（简要）：
        - 排除以 ID/Key/index/Number/TCK/Time 等表示元数据或索引的字段（避免事件/轨迹 ID 等）
        - 保留 PID 相关（粒子鉴别）字段
        - 排除显式列入排除列表的字段（如 RunNumber、EventNumber 等）

    参数:
        name (str): 分支名

    返回:
        bool: True 表示保留该分支（可能有物理意义），False 表示剔除
    """
    if name in "ID":
        if not name in "PID": 
            return False
    if name in "Id":
        if not name in "Pid":
            return False
    bad_suffixes = ("Key", "KEY", "indx", "INDX", "index", 
                    "NUMBER", "Number", "number",
                    "Key", "ID", "indx", "Decision", "ALLPV", "NUMBER","TYPE","TCK","GPSTIME")
    if any(keyword in name for keyword in bad_suffixes):
        return False


    excluded_names = {
        "odin", "ODIN_TCK", 
        "RunNumber", "EventNumber", "FillNumber", 
        "GpsTime", "GPSTIME",
        "BUNCHCROSSING_ID", "BUNCHCROSSING_TYPE"
    }
    if name in excluded_names:
        return False
        
    if "_ID_" in name or "_Id_" in name:
        return False
        
    return True


# ----------------- 主函数 -----------------

def main():
    print("开始提取分支...")
    
    # 1. 收集文件列表
    files_to_process = []
    
    # 检查本地数据文件优先级
    # 优先顺序：本地文件 -> data_dir 中的文件
    if os.path.exists(local_data_file):
        files_to_process.append(os.path.abspath(local_data_file))
        print(f"已添加本地数据文件: {local_data_file}")
    elif os.path.exists(os.path.join(data_dir, local_data_file)):
        files_to_process.append(os.path.join(data_dir, local_data_file))
        print(f"已从 data_dir 添加数据文件: {os.path.join(data_dir, local_data_file)}")
    else:
        print(f"警告：数据文件 {local_data_file} 在本地或 {data_dir} 中未找到")

    # MC Files
    print(f"开始读取MC文件夹: {mc_dir}")
    mc_files = find_files_with_regex(mc_dir, file_pattern)
    print(f"找到 {len(mc_files)} MC 文件.")
    files_to_process.extend(mc_files)
    
    files_to_process = sorted(list(set(files_to_process)))
    
    if not files_to_process:
        print("错误：未找到任何待处理文件。")
        return

    # 开始处理文件以提取公共分支
    common_branches = None
    processed_count = 0
    
    print("\n开始处理文件以提取公共分支...")
    
    for i, file_path in enumerate(files_to_process):
        try:
            with uproot.open(file_path) as f:
                # 查找有效的树（候选目录之一）
                tree = None
                valid_tree_name = None
                
                for dirname in candidate_directories:
                    # 构建 TTree 路径
                    # 典型结构：目录/DecayTree
                    path = f"{dirname}/DecayTree" 
                    if path in f:
                        tree = f[path]
                        valid_tree_name = path
                        break
                
                if tree is None:
                    # 若上面未找到，尝试在顶层寻找以防结构不同
                    if "DecayTree" in f:
                        tree = f["DecayTree"]
                        valid_tree_name = "DecayTree"
                    else:
                        print(f"警告：在 {os.path.basename(file_path)} 中未找到有效的树，已跳过。")
                        continue
                
                # 获取所有分支及其类型
                all_branches = tree.typenames() # Dict[name, type_string]
                
                current_file_branches = set()
                
                for b_name, b_type in all_branches.items():
                    # 应用物理意义过滤
                    if not is_physical_meaningful(b_name):
                        continue
                    
                    # 应用数值类型过滤
                    if not is_numeric_type_string(b_type):
                        continue
                        
                    current_file_branches.add(b_name)
                
                # Update intersection
                if common_branches is None:
                    common_branches = current_file_branches
                else:
                    common_branches = common_branches.intersection(current_file_branches)
                
                processed_count += 1
                print(f"[{i+1}/{len(files_to_process)}] {os.path.basename(file_path)}: {len(current_file_branches)} 个有效分支。当前公共集合大小: {len(common_branches) if common_branches else 0}")
                
                if common_branches is not None and len(common_branches) == 0:
                    print("警告：公共分支集合已为空！")
                    # 如果交集为空，可以提前停止
                    break
                    
        except Exception as e:
            print(f"Error processing {os.path.basename(file_path)}: {e}")
            continue

    # Result
    print("-" * 50)
    if common_branches and len(common_branches) > 0:
        print(f"Final Result: {len(common_branches)} common branches found.")
        sorted_branches = sorted(list(common_branches))
        
        output_csv = "./csv/filtered_branches.csv"
        df = pd.DataFrame(sorted_branches, columns=["Branch"])
        df.to_csv(output_csv, index=False)
        print(f"Successfully saved to {output_csv}")
    else:
        print("No common branches found.")

if __name__ == "__main__":
    main()
