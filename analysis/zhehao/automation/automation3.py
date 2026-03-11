import ROOT
import re
import os
from ROOT import TCanvas, TH1F, TLegend, TLatex, gStyle

# 设置身份识别
# 我是 Gemini 3.1 Pro，一个 AI 智能助手。

try:
    from branch_config import all_branches
    print(f"[+] 成功从外部配置文件读取了 {len(all_branches)} 个 Branch 名称。")
except ImportError:
    print("[-] 错误: 找不到 branch_config.py 文件，请确保它与本脚本在同一目录下。")
    exit()

# ==========================================
# 0. 基础设置
# ==========================================
print("==================================================")
print("[*] 欢迎使用批量处理自动机 - Code by Steven with Gemini 3.1 Pro")
print("==================================================")

ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

# ==========================================
# 1. 数据源定义
# ==========================================
data_file_path = "/home/steven/HNL/data_MagUp_Sprucing24r1_all_disk_options_20260228_011229_skim.root"
mc_base_dir = "/home/steven/HNL/"

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
mc_file_list = [f.strip() for f in raw_mc_files.splitlines() if f.strip()]

directories = [
    "myTupleOS1J/DecayTree",              
    "myTupleOS2J/DecayTree",              
    "myTupleSS1J/DecayTree",              
    "myTupleSS2J/DecayTree"
]

# 定义绘图样式（支持最多6个MC + 1个Data）
# 颜色：红, 蓝, 绿, 紫, 橙, 青; Data用黑色
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+7, ROOT.kCyan+2]
# Data 样式固定为黑色点
data_color = ROOT.kBlack

# ==========================================
# 2. 交互式选择逻辑
# ==========================================

# 1. 选择 Directory
print("[步骤 1] 请选择要分析的 Directory (输入编号):")
for i, d in enumerate(directories):
    print(f"  {i+1}. {d}")
dir_choice = int(input("选择 (1-4): ")) - 1
selected_directory = directories[dir_choice]

# 2. 选择 MC 文件
print("[步骤 2] 请选择要参与绘图的 MC 文件编号 (最多6个, 用下划线分隔, 如 1_5_12):")
for i, f in enumerate(mc_file_list):
    print(f"  {i+1:2d}. {f}")
mc_choices_raw = input("选择编号: ")
selected_mc_indices = [int(x)-1 for x in mc_choices_raw.split('_') if x.strip()]
if len(selected_mc_indices) > 6:
    print("[!] 警告: 选择了超过6个MC，将只取前6个。")
    selected_mc_indices = selected_mc_indices[:6]
selected_mc_files = [mc_file_list[i] for i in selected_mc_indices]

# 3. 正则过滤 Branch
print(f"[步骤 3] 当前全集共 {len(all_branches)} 个 Branch。请输入正则过滤条件。")
regex_input = input("正则 (默认 '.*_PT$'): ")
if not regex_input.strip(): regex_input = ".*_PT$"
pattern = re.compile(regex_input)
target_branches = [b for b in all_branches if pattern.match(b)]
print(f"[✔] 筛选出 {len(target_branches)} 个符合条件的 Branch。")

use_logy = input("Y轴对数？(y/n): ").lower().strip() == 'y'
use_logx = input("X轴对数？(y/n): ").lower().strip() == 'y'

output_dir = "Combined_Plots"
if not os.path.exists(output_dir): os.makedirs(output_dir)

# ==========================================
# 3. 核心联合绘图函数
# ==========================================
def draw_combined(branch_name, mc_files, data_path, tree_path):
    # 容器：存储所有直方图和文件引用以防内存提前释放
    hist_list = []
    file_handles = []
    
    # 1. 处理 Data
    f_data = ROOT.TFile.Open(data_path)
    if not f_data or f_data.IsZombie(): return
    t_data = f_data.Get(tree_path)
    if not t_data: return
    
    d_min = t_data.GetMinimum(branch_name)
    d_max = t_data.GetMaximum(branch_name)
    
    # 2. 处理所有选中的 MC
    g_min, g_max = d_min, d_max
    mc_info = [] # 存储 (tree, short_name)

    for i, mc_f in enumerate(mc_files):
        full_path = os.path.join(mc_base_dir, mc_f)
        f_mc = ROOT.TFile.Open(full_path)
        if not f_mc or f_mc.IsZombie(): continue
        t_mc = f_mc.Get(tree_path)
        if not t_mc: continue
        
        g_min = min(g_min, t_mc.GetMinimum(branch_name))
        g_max = max(g_max, t_mc.GetMaximum(branch_name))
        
        short_name = mc_f.split("_")[1] # 提取如 5GeVCtau0ps
        mc_info.append((t_mc, short_name))
        file_handles.append(f_mc)

    # 范围异常检查
    if g_max <= g_min or g_max > 1e12: return
    if use_logx: g_min = max(g_min, 1e-3)

    # 3. 创建直方图并填充
    # Data Histogram
    h_data = TH1F(f"h_data_{branch_name}", f";{branch_name};Normalized Yields", 100, g_min, g_max)
    t_data.Draw(f"{branch_name} >> h_data_{branch_name}", "", "goff")
    if h_data.Integral() > 0: h_data.Scale(1.0 / h_data.Integral())
    h_data.SetLineColor(data_color)
    h_data.SetMarkerColor(data_color)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(0.8)

    # MC Histograms
    for i, (tree, name) in enumerate(mc_info):
        h_name = f"h_mc_{i}_{branch_name}"
        h_temp = TH1F(h_name, "", 100, g_min, g_max)
        tree.Draw(f"{branch_name} >> {h_name}", "", "goff")
        if h_temp.Integral() > 0: h_temp.Scale(1.0 / h_temp.Integral())
        
        h_temp.SetLineColor(colors[i % len(colors)])
        h_temp.SetLineWidth(2)
        hist_list.append((h_temp, name))

    # 4. 画布渲染
    c = TCanvas("c", "c", 900, 700)
    if use_logy: c.SetLogy(1)
    if use_logx: c.SetLogx(1)

    # 计算 Y 轴最大值
    max_y = h_data.GetMaximum()
    for h, _ in hist_list:
        max_y = max(max_y, h.GetMaximum())
    
    # 设置坐标轴
    if use_logy:
        h_data.GetYaxis().SetRangeUser(1e-4, max_y * 15)
    else:
        h_data.GetYaxis().SetRangeUser(0, max_y * 1.4)

    h_data.Draw("PE") # 先画 Data
    for h, _ in hist_list:
        h.Draw("HIST SAME")

    # 5. 图例与标签
    leg = TLegend(0.50, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(h_data, "Data Background", "lp")
    for h, name in hist_list:
        leg.AddEntry(h, f"MC: {name}", "l")
    leg.Draw()

    tex = TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.DrawLatex(0.15, 0.84, "#bf{LHCb} internal")
    dir_tag = tree_path.split("/")[0].replace("myTuple", "")
    tex.DrawLatex(0.15, 0.79, f"Region: {dir_tag}")
    tex.DrawLatex(0.15, 0.74, f"Var: {branch_name}")

    # 保存
    save_path = f"{output_dir}/Combined_{dir_tag}_{branch_name}.png"
    c.SaveAs(save_path)
    
    # 清理
    c.Close()
    f_data.Close()
    for f in file_handles: f.Close()

# ==========================================
# 4. 运行
# ==========================================
print(f"[*] 开始任务：处理 {len(target_branches)} 个变量...")
for i, branch in enumerate(target_branches):
    print(f"[{i+1}/{len(target_branches)}] 正在绘制: {branch}")
    draw_combined(branch, selected_mc_files, data_file_path, selected_directory)

print(f"[✔] 全部完成！结果已存入 {output_dir} 文件夹。")
