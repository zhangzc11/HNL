import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ==========================================
# 1. 定义你图片中的 W_PT 理论 PDF 与生成器
# ==========================================
def theory_W_pt_pdf(pt, peak_val=8.0):
    # Sudakov 峰 + QCD 长尾
    sudakov = (pt / peak_val**2) * np.exp(-(pt**2) / (2 * peak_val**2))
    tail = 0.05 / (pt + 5)**2
    dist = sudakov + tail
    return np.where(pt > 0, dist, 0)

def generate_W_PT(N, peak_val=8.0):
    """基于舍选法生成随机 W_PT"""
    pt = np.zeros(N)
    count = 0
    max_pdf = theory_W_pt_pdf(np.array([peak_val]), peak_val)[0] * 1.2
    while count < N:
        x = np.random.uniform(0.01, 100.0, N)
        y = np.random.uniform(0, max_pdf, N)
        valid = x[y < theory_W_pt_pdf(x, peak_val)]
        n_v = len(valid)
        if n_v > 0:
            fill = min(n_v, N - count)
            pt[count:count+fill] = valid[:fill]
            count += fill
    return pt

# ==========================================
# 2. 相对论运动学工具
# ==========================================
def boost(p, beta):
    b2 = np.sum(beta**2, axis=1)
    b2 = np.clip(b2, 0, 0.99999)
    gamma = 1.0 / np.sqrt(1 - b2)
    bp = np.sum(p[:, 1:] * beta, axis=1)
    b2_safe = np.where(b2 == 0, 1e-10, b2)
    gamma2 = (gamma - 1.0) / b2_safe
    
    p_prime = np.zeros_like(p)
    p_prime[:, 0] = gamma * (p[:, 0] + bp)
    p_prime[:, 1] = p[:, 1] + (gamma2 * bp + gamma * p[:, 0]) * beta[:, 0]
    p_prime[:, 2] = p[:, 2] + (gamma2 * bp + gamma * p[:, 0]) * beta[:, 1]
    p_prime[:, 3] = p[:, 3] + (gamma2 * bp + gamma * p[:, 0]) * beta[:, 2]
    return p_prime

def get_pt(p):
    return np.sqrt(p[:, 1]**2 + p[:, 2]**2)

# ==========================================
# 3. 模拟主流程
# ==========================================
m_W = 80.4
mass_points = [5, 10, 15, 20, 30, 50]
colors = ['blue', 'dodgerblue', 'turquoise', 'greenyellow', 'orange', 'red']
N_events = 500000 

# --- 场景 A & B 的 W 运动学基础 ---
Y_W_A = np.random.uniform(-6, 6, N_events) 
beta_W_A = np.column_stack((np.zeros(N_events), np.zeros(N_events), np.tanh(Y_W_A)))

PT_W_B = generate_W_PT(N_events, peak_val=8.0) 
Phi_W = np.random.uniform(0, 2*np.pi, N_events)
Y_W_B = np.random.uniform(-6, 6, N_events) # 【核心修改】：同步扩展场景 B 的快度范围
M_T = np.sqrt(m_W**2 + PT_W_B**2)
E_W_B = M_T * np.cosh(Y_W_B)
beta_W_B = np.column_stack((PT_W_B * np.cos(Phi_W), 
                            PT_W_B * np.sin(Phi_W), 
                            M_T * np.sinh(Y_W_B))) / E_W_B[:, np.newaxis]

# 准备绘图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

for m_N, color in zip(mass_points, colors):
    # --- 步骤 1: W 静止系动量生成 ---
    E_mu1_W = (m_W**2 - m_N**2) / (2 * m_W)
    p_star = E_mu1_W
    E_N_W = (m_W**2 + m_N**2) / (2 * m_W)

    cos_th1 = np.random.uniform(-1, 1, N_events)
    phi1 = np.random.uniform(0, 2*np.pi, N_events)
    sin_th1 = np.sqrt(1 - cos_th1**2)

    p_mu1_W = np.column_stack((np.full(N_events, E_mu1_W), 
                               p_star*sin_th1*np.cos(phi1), 
                               p_star*sin_th1*np.sin(phi1), 
                               p_star*cos_th1))
    p_N_W = np.column_stack((np.full(N_events, E_N_W), -p_mu1_W[:, 1:]))

    # --- 步骤 2: 洛伦兹变换 ---
    p_mu1_lab_A = boost(p_mu1_W, beta_W_A)
    p_N_lab_A = boost(p_N_W, beta_W_A)
    p_mu1_lab_B = boost(p_mu1_W, beta_W_B)
    p_N_lab_B = boost(p_N_W, beta_W_B)

    # --- 步骤 3: 提取并绘图 ---
    pt_mu_A, pt_mu_B = get_pt(p_mu1_lab_A), get_pt(p_mu1_lab_B)
    pt_N_A, pt_N_B = get_pt(p_N_lab_A), get_pt(p_N_lab_B)

    ax1.hist(pt_mu_A, bins=100, range=(0, 60), density=True, histtype='step', color=color, linestyle=':', lw=1.2, alpha=0.6)
    ax1.hist(pt_mu_B, bins=100, range=(0, 60), density=True, histtype='step', color=color, linestyle='-', lw=1.8)
    
    ax2.hist(pt_N_A, bins=100, range=(0, 60), density=True, histtype='step', color=color, linestyle=':', lw=1.2, alpha=0.6)
    ax2.hist(pt_N_B, bins=100, range=(0, 60), density=True, histtype='step', color=color, linestyle='-', lw=1.8)

# 装饰
ax1.axvline(20, color='black', linestyle='--', lw=2, label='Trigger Cut')
ax1.set_title("Prompt Muon $P_T$")
ax2.set_title("HNL $P_T$ ")

custom_lines = [Line2D([0], [0], color='gray', linestyle='-'), Line2D([0], [0], color='gray', linestyle=':')]
ax2.legend(custom_lines, ['Realistic (With $P_{T,W}$)', 'Idealized ($P_{T,W}=0$)'], loc='upper right')

plt.tight_layout()

save_path = "NuPT_Simulation.png"
plt.savefig(save_path, dpi=300)
print(f"图像已成功保存")
plt.show()