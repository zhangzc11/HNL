## 一、简介

```train_classifier.py``` 程序是一个基于机器学习的 HNL 信号与本底区分程序（使用ChatGPT 5.4 辅助完成）。本程序用于对 **HNL物理分析**中的 **Monte Carlo 信号样本** 与 **LHCb 实验真实数据** 进行比较，并利用机器学习方法自动寻找能够有效区分 **信号（signal）** 与 **本底（background）** 的变量。由

程序主要完成以下任务：

1. 读取 ROOT 文件中的物理变量（branch）
2. 将 MC 样本作为 **信号**
3. 将真实数据作为 **本底**
4. 构建训练数据集
5. 训练机器学习分类器（Gradient Boosting）
6. 自动评估模型区分能力
7. 输出最具有区分能力的物理变量

---

## 二、程序功能

程序主要包含以下功能：

### 1 读取 ROOT 数据

利用 `uproot` 读取 ROOT 文件中的 `TTree`，提取指定的变量，例如原程序默认分析的变量为（**可根据需要更改**）：

* `W_PT`
* `Mu_PT`
* `MuNuR_PT`
* `Mu_ETA`
* `nPVs`

ROOT 文件中的树结构类似：

```
myTupleSS1J/DecayTree
myTupleSS2J/DecayTree
myTupleOS1J/DecayTree
myTupleOS2J/DecayTree
```

程序默认读取：

```
myTupleSS1J/DecayTree
```

---

### 2 构建机器学习数据集

程序将数据分为两类：

| 数据来源       | 标签             |
| ---------- | -------------- |
| MC ROOT 文件 | signal = 1     |
| LHCb 实验数据  | background = 0 |

最终构建的数据结构类似：

| W_PT | Mu_PT | MuNuR_PT | Mu_ETA | nPVs | label |
| ---- | ----- | -------- | ------ | ---- | ----- |
| ...  | ...   | ...      | ...    | ...  | 1     |
| ...  | ...   | ...      | ...    | ...  | 0     |

---

### 3 训练机器学习模型

程序使用 **Gradient Boosting Decision Tree (GBDT)** 进行训练。该算法在高能物理数据分析中常用于 **信号-本底分类**。

---

### 4 评估分类性能

程序计算模型的 **ROC AUC 值**：

```
AUC = 0.5  → 完全无法区分
AUC ≈ 0.7 → 有一定区分能力
AUC ≈ 0.8 → 区分能力较好
AUC > 0.9 → 区分能力很强
```

---

### 5 自动寻找最重要变量

机器学习模型会输出每个变量的重要性：

```
Feature Importance
```

例如本程序会输出：

```
Mu_PT      0.41
W_PT       0.25
MuNuR_PT   0.20
Mu_ETA     0.09
nPVs       0.05
```

这些变量通常是区分信号与本底的关键物理量。

---

# 三、项目结构

示例目录结构：

```
project/
│
├── train_classifier.py
├── README.md
│
├── data/
│   ├── 42912010_5GeVCtau0ps.root
│   ├── 42912011_5GeVCtau10ps.root
│   ├── ...
│   └── data_MagUp_..._skim.root
│
└── ml_output/
    ├── feature_importance.png
    ├── roc_curve.png
    └── feature_ranking.txt
```

---

# 四、运行环境

推荐在 **WSL Linux 环境**中运行。

需要安装以下 Python 库：

```
uproot
pandas
numpy
scikit-learn
matplotlib
```

Linux终端的安装命令为：

```
python3 -m pip install uproot pandas numpy scikit-learn matplotlib
```

# 五、数据路径设置

程序默认数据路径为：

```
/home/data
```

代码中：

```python
DATA_DIR = "/home/data"
```

其中：

| 文件类型             | 说明          |
| ---------------- | ----------- |
| 429*.root        | MC 信号样本     |
| data_MagUp*.root | LHCb 实验真实数据 |

要将其修改为实际运行代码时的数据路径与数据文件名称。

---

# 六、运行程序

在 Linux 终端运行：

```bash
python3 train_classifier.py
```

程序运行流程：
1. 搜索所有 MC ROOT 文件
2. 读取 ROOT 变量
3. 构建训练数据
4. 训练机器学习模型
5. 输出结果

---

# 七、程序输出

程序会在 `ml_output` 目录中生成以下文件：

### 1 变量重要性图

```
feature_importance.png
```

展示每个变量对信号-本底区分的贡献。

---

### 2 ROC 曲线

```
roc_curve.png
```

展示模型分类性能。

---

### 3 变量排名

```
feature_ranking.txt
```

列出变量区分能力排序。

若已存在名称为 `ml_output` 目录，则会覆盖前一次的结果。因此注意及时修改上一次的结果目录名称。

---

# 八、可修改参数

### 修改读取的 ROOT 树

```python
TREE_NAME = "myTupleSS1J/DecayTree"
```

可改为：

```
myTupleSS2J/DecayTree
myTupleOS1J/DecayTree
myTupleOS2J/DecayTree
```

---

### 修改机器学习变量

```python
VARIABLES = [
    "W_PT",
    "Mu_PT",
    "MuNuR_PT",
    "Mu_ETA",
    "nPVs"
]
```

可以根据 ROOT 文件中的 branch 进行扩展。

---
