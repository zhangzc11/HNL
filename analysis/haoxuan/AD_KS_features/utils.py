import uproot

def load_file(file_path, tree_name):
    file = uproot.open(file_path)
    tree = file[tree_name]
    # 定义要过滤的关键词列表
    exclude_keywords = ["Key", "ID", "indx", "Decision", "ALLPV", "NUMBER","TYPE","TCK","GPSTIME"]
    physics_keys = [key for key in tree.keys() 
                   if not any(keyword in key for keyword in exclude_keywords)]
    # 提取物理意义的数据
    physics_data = tree.arrays(physics_keys, library="np")
    return physics_data