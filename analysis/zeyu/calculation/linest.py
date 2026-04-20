import numpy as np
import matplotlib
matplotlib.use('Agg') # 服务器端必须使用 Agg 模式绘图，否则会报错
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. 定义 Sigmoid 函数
def sigmoid(x, x0, k, max_eff):
    return max_eff / (1.0 + np.exp(-(x - x0) / k))

# 2. 数据读取函数（在这里接入你服务器上的真实数据）
def get_data_from_server(m_N):
    """
    TODO: 替换为你的真实读取逻辑
    你需要返回：
    - pt_bins: P_T 的分箱中心点 (ndarray)
    - efficiency: 每个分箱对应的真实效率值 (ndarray)
    """
    # 临时占位逻辑（为了让代码能跑通，演示用）
    # 实际操作时，你应该从你的 ROOT 或 txt 文件读取这些点
    x_test = np.linspace(5, 75, 15)
    y_test = sigmoid(x_test, 5.0 + 0.42 * m_N, 2.5 + 0.08 * m_N, 0.05)
    return x_test, y_test

# 3. 主分析逻辑
mass_points = [5, 10, 15, 20, 30, 50]
extracted_x0 = []
extracted_k = []

plt.figure(figsize=(8, 6))

for m in mass_points:
    x_data, y_data = get_data_from_server(m)
    
    # 执行拟合
    try:
        # p0 是初始猜测值 [x0, k, max_eff]
        popt, _ = curve_fit(sigmoid, x_data, y_data, p0=[20, 5, 0.01])
        extracted_x0.append(popt[0])
        extracted_k.append(popt[1])
        
        # 绘图校验
        plt.plot(x_data, y_data, 'o', alpha=0.5)
        x_fine = np.linspace(0, 80, 100)
        plt.plot(x_fine, sigmoid(x_fine, *popt), label=f'm={m} GeV')
    except:
        print(f"Mass {m} fit failed!")

plt.legend()
plt.savefig('acceptance_check.png') # 结果保存为图片

# 4. 核心：线性回归得到你向导师汇报的公式
slope_x, inter_x = np.polyfit(mass_points, extracted_x0, 1)
slope_k, inter_k = np.polyfit(mass_points, extracted_k, 1)

print("\n" + "="*30)
print("拟合完成！你的线性公式参数如下：")
print(f"开启点位置: x0 = {slope_x:.4f} * m_N + {inter_x:.4f}")
print(f"开启宽度:   k  = {slope_k:.4f} * m_N + {inter_k:.4f}")
print("="*30)
