import ROOT
import re
import os
from ROOT import TCanvas, TH1F, TLegend, TLatex, gStyle

try:
    from branch_config import all_branches
    print(f"[+] 成功从外部配置文件读取了 {len(all_branches)} 个 Branch 名称。")
except ImportError:
    print("[-] 错误: 找不到 branch_config.py 文件，请确保它与本脚本在同一目录下。")
    exit()
# ==========================================
# 0. 基础设置与身份验证
# ==========================================
print("==================================================")
print("[*] 欢迎使用批量处理自动机 - Code by Steven with the help of Gemini 3.1 Pro")
print("==================================================")

# 后台静默模式，防止循环作图时弹出大量窗口卡死系统
ROOT.gROOT.SetBatch(True)

gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# ==========================================
# 1. 数据源与遍历列表定义
# ==========================================
data_file_path = "/home/steven/HNL/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root" # 这是skim过的data文件
mc_base_dir = "/home/steven/HNL/" # 请确保此路径正确，如果需要本地使用请修改路径

# 原始MC文件列表
raw_mc_files = """
42912010_5GeVCtau0ps_options_20260228_111712_skim.root
42912011_5GeVCtau10ps_options_20260228_111712_skim.root
42912012_10GeVCtau0ps_options_20260228_111712_skim.root
42912013_10GeVCtau10ps_options_20260228_111712_skim.root
42912014_15GeVCtau0ps_options_20260228_111712_skim.root
42912015_15GeVCtau10ps_options_20260228_111712_skim.root
42912016_20GeVCtau0ps_options_20260228_111712_skim.root
42912017_20GeVCtau10ps_options_20260228_111712_skim.root
42912018_30GeVCtau0ps_options_20260228_111712_skim.root
42912019_30GeVCtau10ps_options_20260228_230635_skim.root
42912020_50GeVCtau0ps_options_20260228_230639_skim.root
42912021_50GeVCtau10ps_options_20260228_111712_skim.root
42912030_5GeVCtau100ps_options_20260228_111712_skim.root
42912031_5GeVCtau1000ps_options_20260228_111712_skim.root
42912032_10GeVCtau100ps_options_20260228_111712_skim.root
42912033_10GeVCtau1000ps_options_20260228_111712_skim.root
42912034_15GeVCtau100ps_options_20260228_111712_skim.root
42912035_15GeVCtau1000ps_options_20260228_111712_skim.root
42912036_20GeVCtau100ps_options_20260228_111712_skim.root
42912037_20GeVCtau1000ps_options_20260228_111712_skim.root
42912038_30GeVCtau100ps_options_20260228_111712_skim.root
42912039_30GeVCtau1000ps_options_20260228_111712_skim.root
42912040_50GeVCtau100ps_options_20260228_111712_skim.root
42912041_50GeVCtau1000ps_options_20260228_111712_skim.root
"""
mc_file_list = [f.strip().replace(" ", "_") for f in raw_mc_files.splitlines() if f.strip()]
#这一行实际上

directories = [
    "myTupleOS1J/DecayTree",              
    "myTupleOS2J/DecayTree",              
    "myTupleSS1J/DecayTree",              
    "myTupleSS2J/DecayTree"
]

# ==========================================
# 2. 交互式界面：正则过滤 Branch
# ==========================================
print("[交互选项] 当前检测到总共 1163 个 Branch。")
print("请输入正则表达式以选择想要遍历的 Branch。")
print("示例 1: .*_PT$   (只画所有以 _PT 结尾的变量)")
print("示例 2: ^W_      (只画所有以 W_ 开头的变量)")
print("示例 3: .*       (匹配全部变量)")
regex_input = input("请输入正则表达式 (按回车默认使用 '.*_PT$'): ")

if not regex_input.strip():
    regex_input = ".*_PT$"

try:
    pattern = re.compile(regex_input)
    target_branches = [b for b in all_branches if pattern.match(b)]
    print(f"[✔] 正则匹配成功！共筛选出 {len(target_branches)} 个符合条件的 Branch 进行遍历。")
except Exception as e:
    print(f"[-] 正则表达式错误: {e}")
    exit()

if len(target_branches) == 0:
    print("[-] 没有找到匹配的 Branch，程序退出。")
    exit()

use_logy = input("[交互选项] 纵轴(Y-axis)是否开启对数坐标？(y/n, 默认 n): ").lower().strip() == 'y'
use_logx = input("[交互选项] 横轴(X-axis)是否开启对数坐标？(y/n, 默认 n): ").lower().strip() == 'y'

if use_logy: print("[↑] Y 轴已设为对数模式。")
if use_logx: print("[↔] X 轴已设为对数模式。")

# 创建保存图片的文件夹
output_dir = "Output_Plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================
# 3. 核心获取与绘图函数
# ==========================================
def draw_and_save(mc_path, data_path, tree_path, branch_name):
    # 解析输出文件名所需的要素
    full_mc_filename = mc_path.split("/")[-1].replace(".root", "")
    short_mc_name = full_mc_filename.split("_")[1]
    dir_name_safe = tree_path.split("/")[0] # 例如取出 myTupleOS1J
    dir_name_safe =dir_name_safe.replace("myTuple", "")

    output_filename = f"{output_dir}/{short_mc_name}_{dir_name_safe}_{branch_name}.png"

    # 打开文件
    f_mc = ROOT.TFile.Open(mc_path)
    f_data = ROOT.TFile.Open(data_path)
    
    if not f_mc or f_mc.IsZombie() or not f_data or f_data.IsZombie():
        print(f"[-] 错误: 无法打开文件 (MC 或 Data)。跳过。")
        return

    tree_mc = f_mc.Get(tree_path)
    tree_data = f_data.Get(tree_path)

    if not tree_mc or not tree_data:
        print(f"[-] 警告: {tree_path} 不存在。跳过。")
        f_mc.Close()
        f_data.Close()
        return

    # 动态获取坐标轴极值范围 (Min/Max)
    mc_min = tree_mc.GetMinimum(branch_name)
    mc_max = tree_mc.GetMaximum(branch_name)
    data_min = tree_data.GetMinimum(branch_name)
    data_max = tree_data.GetMaximum(branch_name)

    global_min = min(mc_min, data_min)
    global_max = max(mc_max, data_max)

    # 剔除极端无意义的极值 (如果某些Tree全空会导致GetMinimum返回极大极小值)
    if global_max < global_min or global_max > 1e10:
        print(f"[!] 警告: {branch_name} 极值异常，跳过绘图。")
        return

    if use_logx:
        # 如果开启 X 对数，强制最小值必须大于 0，通常设为一个极小的正数（如 1e-3）
        # 否则 ROOT 的 LogX 无法正常显示
        global_min = max(global_min, 1e-3)
        if global_max <= 0:
            print(f"[!] 警告: {branch_name} 最大值 <= 0，无法绘制 LogX。跳过。")
            return

    # 创建直方图
    h_name_mc = f"h_mc_{branch_name}_{short_mc_name}_{dir_name_safe}"
    h_name_data = f"h_data_{branch_name}_{short_mc_name}_{dir_name_safe}"
    
    h_mc = TH1F(h_name_mc, f";{branch_name};Normalized Yields", 100, global_min, global_max)
    h_data = TH1F(h_name_data, f";{branch_name};Normalized Yields", 100, global_min, global_max)

    # 核心修复：>> 后面必须匹配上面定义的 h_name
    tree_mc.Draw(f"{branch_name} >> {h_name_mc}", "", "goff")
    tree_data.Draw(f"{branch_name} >> {h_name_data}", "", "goff")


    # 归一化
    if h_mc.Integral() > 0: h_mc.Scale(1.0 / h_mc.Integral())
    if h_data.Integral() > 0: h_data.Scale(1.0 / h_data.Integral())

    h_mc.SetLineColor(ROOT.kBlue)
    h_mc.SetLineWidth(2)
    h_data.SetLineColor(ROOT.kRed)
    h_data.SetMarkerColor(ROOT.kRed)
    h_data.SetMarkerStyle(20) # Data 用实心点

    # 绘图
    c1 = TCanvas("c1", "c1", 800, 600)

    if use_logy:
        c1.SetLogy(1)  # 开启 Y 轴对数
    if use_logx:
        c1.SetLogx(1)  # 开启 X 轴对数

    # 动态 Y 轴高度调整
    max_y = max(h_mc.GetMaximum(), h_data.GetMaximum())
    
    if use_logy:
        # 对数模式：起点设为 1e-4（适合归一化后的产额），终点留出 10 倍空间给图例
        h_mc.GetYaxis().SetRangeUser(1e-4, max_y * 10)
    else:
        # 线性模式：保持原样
        h_mc.GetYaxis().SetRangeUser(0, max_y * 1.3)
    
    h_mc.Draw("HIST")
    h_data.Draw("PE SAME")

    # 图例
    leg = TLegend(0.55, 0.7, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_mc, f"MC ({short_mc_name})", "l")
    leg.AddEntry(h_data, "Data Background", "lp")
    leg.Draw()

    # LHCb 文本说明
    tex = TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.DrawLatex(0.15, 0.82, "#bf{LHCb} internal")
    tex.DrawLatex(0.15, 0.75, f"Dir: {dir_name_safe}")
    tex.DrawLatex(0.15, 0.68, f"Var: {branch_name}")

    c1.SaveAs(output_filename)
    print(f"[+] 保存成功: {output_filename}")

    # 清理内存，防止循环卡顿
    c1.Close()
    h_mc.Delete()
    h_data.Delete()
    f_mc.Close()
    f_data.Close()

# ==========================================
# 4. 启动三层嵌套自动机遍历
# ==========================================
print("[*] 自动机启动，开始遍历任务...")

total_tasks = len(mc_file_list) * len(directories) * len(target_branches)
current_task = 0

for mc_file in mc_file_list:
    full_mc_path = os.path.join(mc_base_dir, mc_file)
    for directory in directories:
        for branch in target_branches:
            current_task += 1
            print(f"[{current_task}/{total_tasks}] 处理 -> MC: {mc_file[:20]}... | Dir: {directory.split('/')[0]} | Branch: {branch}")
            draw_and_save(full_mc_path, data_file_path, directory, branch)

print("[✔] 所有遍历绘制任务执行完毕！图片已保存在 Output_Plots 文件夹下。")
