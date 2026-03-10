#!/usr/bin/env python3
"""
ROOT TTree分支提取脚本
功能：遍历ROOT文件中所有TTree，提取分支名称和类型信息，输出到CSV文件
"""

import ROOT
import csv
import sys
import os

def extract_tree_branches(root_file_path, output_csv_path):
    """
    从ROOT文件中提取所有TTree的分支信息
    
    参数:
        root_file_path: ROOT文件路径
        output_csv_path: 输出CSV文件路径
    """
    
    # 打开ROOT文件
    if not os.path.exists(root_file_path):
        print(f"错误：文件 {root_file_path} 不存在！")
        return False
    
    root_file = ROOT.TFile.Open(root_file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"错误：无法打开ROOT文件 {root_file_path}")
        return False
    
    print(f"成功打开ROOT文件: {root_file_path}")
    print("=" * 60)
    
    # 存储所有分支信息
    all_branches = []
    
    # 获取文件中的所有keys
    keys = root_file.GetListOfKeys()
    
    # 遍历所有顶级对象
    for key in keys:
        obj_name = key.GetName()
        obj_class = key.GetClassName()
        
        print(f"\n处理对象: {obj_name} (类型: {obj_class})")
        
        # 如果是TDirectoryFile，进入目录查找TTree
        if "TDirectoryFile" in obj_class:
            dir_obj = root_file.Get(obj_name)
            if dir_obj:
                # 获取目录中的所有keys
                dir_keys = dir_obj.GetListOfKeys()
                for dir_key in dir_keys:
                    inner_name = dir_key.GetName()
                    inner_class = dir_key.GetClassName()
                    
                    if "TTree" in inner_class:
                        tree = dir_obj.Get(inner_name)
                        if tree:
                            print(f"  找到TTree: {inner_name}")
                            branches_info = get_tree_branches(tree, obj_name, inner_name)
                            all_branches.extend(branches_info)
                            print(f"    分支数量: {len(branches_info)}")
                    break
                break
        
        # 如果直接是TTree
        elif "TTree" in obj_class:
            tree = root_file.Get(obj_name)
            if tree:
                print(f"  找到TTree: {obj_name}")
                branches_info = get_tree_branches(tree, "", obj_name)
                all_branches.extend(branches_info)
                print(f"    分支数量: {len(branches_info)}")
                break
    
    # 关闭ROOT文件
    root_file.Close()
    
    # 写入CSV文件
    if all_branches:
        write_to_csv(all_branches, output_csv_path)
        print("\n" + "=" * 60)
        print(f"成功！共提取 {len(all_branches)} 个分支信息")
        print(f"输出文件: {output_csv_path}")
        return True
    else:
        print("警告：未找到任何TTree分支信息")
        return False


def get_tree_branches(tree, directory_name, tree_name):
    """
    获取TTree中所有分支的信息
    
    参数:
        tree: TTree对象
        directory_name: 所在目录名
        tree_name: TTree名称
    
    返回:
        分支信息列表
    """
    branches_info = []
    
    # 获取所有分支
    branches = tree.GetListOfBranches()
    
    for branch in branches:
        branch_name = branch.GetName()
        branch_title = branch.GetTitle()
        branch_class = branch.GetClassName()
        
        # 获取叶子(leaf)信息
        leaves = branch.GetListOfLeaves()
        leaf_types = []
        leaf_names = []
        
        for leaf in leaves:
            leaf_names.append(leaf.GetName())
            leaf_types.append(leaf.GetTypeName() if leaf.GetTypeName() else leaf.GetType())
        
        # 如果分支有子分支（如vector等）
        sub_branches = branch.GetListOfBranches()
        if sub_branches and sub_branches.GetEntries() > 0:
            # 递归处理子分支
            for sub_branch in sub_branches:
                sub_info = {
                    'Directory': directory_name,
                    'Tree': tree_name,
                    'Branch': f"{branch_name}.{sub_branch.GetName()}",
                    'BranchTitle': sub_branch.GetTitle(),
                    'Type': get_leaf_type(sub_branch),
                    'LeafName': sub_branch.GetName()
                }
                branches_info.append(sub_info)
        else:
            branch_info = {
                'Directory': directory_name,
                'Tree': tree_name,
                'Branch': branch_name,
                'BranchTitle': branch_title,
                'Type': ', '.join(leaf_types) if leaf_types else branch_class,
                'LeafName': ', '.join(leaf_names) if leaf_names else branch_name
            }
            branches_info.append(branch_info)
    
    return branches_info


def get_leaf_type(branch):
    """
    获取分支的叶子类型
    """
    leaves = branch.GetListOfLeaves()
    if leaves and leaves.GetEntries() > 0:
        leaf = leaves.First()
        type_name = leaf.GetTypeName()
        if type_name:
            return type_name
        return leaf.GetType()
    return "Unknown"


def write_to_csv(branches_info, output_path):
    """
    将分支信息写入CSV文件
    
    参数:
        branches_info: 分支信息列表
        output_path: 输出文件路径
    """
    fieldnames = ['Directory', 'Tree', 'Branch', 'BranchTitle', 'Type', 'LeafName']
    
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(branches_info)


def print_tree_info(root_file_path):
    """
    打印ROOT文件中所有TTree的详细信息（类似Print()输出）
    """
    root_file = ROOT.TFile.Open(root_file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"错误：无法打开ROOT文件 {root_file_path}")
        return
    
    print("\n" + "=" * 80)
    print("ROOT文件TTree详细信息 (类似Print()输出)")
    print("=" * 80)
    
    keys = root_file.GetListOfKeys()
    
    for key in keys:
        obj_name = key.GetName()
        obj_class = key.GetClassName()
        
        if "TDirectoryFile" in obj_class:
            dir_obj = root_file.Get(obj_name)
            if dir_obj:
                dir_keys = dir_obj.GetListOfKeys()
                for dir_key in dir_keys:
                    if "TTree" in dir_key.GetClassName():
                        tree = dir_obj.Get(dir_key.GetName())
                        if tree:
                            print(f"\n{'='*60}")
                            print(f"目录: {obj_name} | TTree: {dir_key.GetName()}")
                            print(f"{'='*60}")
                            tree.Print()
                            print(f"\n条目数: {tree.GetEntries()}")
    
    root_file.Close()


if __name__ == "__main__":
    # 默认文件路径
    default_root_file = "/pool1/lhcb/zhicaiz/data/Collision24/ntuple_HNL_Wmumuqq/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229.root"
    default_output_csv = "./data/tree_branches_info.csv"
    
    # 解析命令行参数
    if len(sys.argv) >= 2:
        root_file = sys.argv[1]
    else:
        root_file = default_root_file
    
    if len(sys.argv) >= 3:
        output_csv = sys.argv[2]
    else:
        output_csv = default_output_csv
    
    print(f"输入ROOT文件: {root_file}")
    print(f"输出CSV文件: {output_csv}")
    
    # 提取分支信息并保存到CSV
    success = extract_tree_branches(root_file, output_csv)
    
    # 可选：打印详细的TTree信息
    if success:
        print("\n是否打印详细的TTree信息？(y/n): ", end="")
        try:
            choice = input().strip().lower()
            if choice == 'y':
                print_tree_info(root_file)
        except:
            pass
