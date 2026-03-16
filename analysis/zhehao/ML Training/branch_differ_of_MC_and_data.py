import uproot

# 替换成你的真实路径
data_file_path = "/home/steven/HNL/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_file_path = "/home//steven/HNL/42912010_5GeVCtau0ps_options_20260228_111712_skim.root"
tree_name = "myTupleSS1J/DecayTree"

with uproot.open(data_file_path) as f:
    data_keys = set(f[tree_name].keys())

with uproot.open(mc_file_path) as f:
    mc_keys = set(f[tree_name].keys())

# 找出只在 MC 中存在的变量
mc_only_keys = mc_keys - data_keys

# 找出公共的branch
common_keys = []
for i  in mc_keys:
    if i in data_keys:
        common_keys.append(i)
#按照一定的要求筛选branch（排除某些关键词）
def filter_branches(keys):
    exclude_keywords = [
        'ID', 'KEY', 'DECISION', 'TRUE', 'MC', 'BKGCAT', 
        'TIME', 'NUMBER', 'EVENT', 'INDX', 'TCK', 'TIS', 'TOS', 
        'GPS', 'BUNCH', 'TYPE', 'OWNPVIP', 'OBJECT', 'CLUSTER', 'MATCH'
    ]
    return [b for b in keys if not any(k in b.upper() for k in exclude_keywords)]

filtered_keys = filter_branches(common_keys)
print("只在 MC 文件中存在的变量示例:")
print("common_keys length:")
print(len(common_keys))
print("filtered_keys length:")
print(filtered_keys)
print(len(filtered_keys))