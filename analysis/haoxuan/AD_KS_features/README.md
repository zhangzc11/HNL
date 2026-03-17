# 基于KS检验和AD检验的特征选择工具

本项目通过比较蒙特卡洛模拟（MC）与真实对撞数据（Data）在ROOT树（TTree）每个分支上的分布差异，计算Kolmogorov-Smirnov统计量（KS）和Anderson-Darling准则（AD），并将两者组合为加权分数，最终对分支进行排序并可视化差异最大的前若干个分支。

工作流由`Makefile`自动化，包含两个主要Python脚本：`ks_andeson.py`和`weight_sort.py`，以及一个辅助模块`utils.py`。

## 文件结构

```
.
├── Makefile               # 自动构建脚本
├── ks_andeson.py          # 对所有分支执行KS和Anderson-Darling检验
├── weight_sort.py         # 对检验结果排序并绘制前20个分支的直方图
├── utils.py               # 工具函数：加载ROOT文件并过滤无关分支
├── csv/                   # 存放CSV输出文件的目录（自动创建）
└── png/                   # 存放直方图文件的目录（自动创建）
```

---

## 使用方法

### Makefile方式

`Makefile`中定义了三个可覆盖的变量：

| 变量名      | 含义                                  | 示例值                                                       |
| ----------- | ------------------------------------- | ------------------------------------------------------------ |
| `directory` | ROOT文件中TTree所在目录名（如`myTupleSS1J`） | `myTupleOS2J`                                                |
| `mc_file`   | MC模拟的ROOT文件完整路径               | `/pool1/.../mc.root`                                         |
| `data_file` | 真实数据的ROOT文件完整路径             | `/pool1/.../data.root`                                       |

**重要：** `Makefile`中的默认路径仅为示例，使用时必须根据实际情况覆盖这些变量。

#### 目标说明
- **`sort`**：执行完整流程。  
  1. 运行`ks_andeson.py`，生成`csv/branch_ks.csv`。  
  2. 运行`weight_sort.py`，基于该CSV生成排序后的副本，并绘制按组合权重排序的前20个分支的直方图。  
- **`clean`**：删除`csv/`和`png/`目录下所有生成的文件。

#### 执行示例
```bash
make sort directory="myTupleOS2J" mc_file="/path/to/your/mc.root" data_file="/path/to/your/data.root"
```

### 手动执行

如果不使用Makefile，可以依次运行两个脚本：

1. **生成分支统计结果**  
   ```bash
   python3 ks_andeson.py <directory> <mc_file> <data_file>
   ```
   输出：`csv/branch_ks.csv`

2. **排序并绘图**  
   ```bash
   python3 weight_sort.py <directory> <mc_file> <data_file>
   ```
   输出：排序后的CSV文件位于`csv/`，直方图位于`png/`

---

## 参数说明

两个Python脚本均接受三个命令行参数（顺序一致）：

1. **`directory`**：ROOT文件内TTree的路径前缀，例如`"myTupleOS2J"`，最终构造的树名为`f"{directory}/DecayTree"`。  
2. **`mc_file`**：MC模拟ROOT文件的完整路径。  
3. **`data_file`**：数据ROOT文件的完整路径。

> 注意：脚本内部也定义了默认值（`directory = "myTupleOS2J"`，以及两个硬编码的ROOT路径），但建议总是通过命令行参数指定，以确保使用正确的文件。

---

## 输出结果

### CSV文件（位于`csv/`目录）

- **`branch_ks.csv`**：原始结果，包含四列：`Branch`（分支名）、`KS`（KS统计量）、`Andeson`（Anderson-Darling统计量）、`Weight`（组合权重）。
- **`ks_sorted.csv`**：按KS统计量降序排列。
- **`andeson_sorted.csv`**：按Anderson-Darling统计量降序排列。
- **`weight_sorted.csv`**：按组合权重降序排列（用于绘图的排序依据）。

### 绘图文件（位于`png/`目录）

对于`weight_sorted.csv`中排名前20的分支，分别生成直方图对比图。  
文件名格式：`{排名}_{directory}_{分支名}.png`  
每张图包含：
- MC（蓝色半透明）与Data（橙色半透明）的归一化直方图。
- 标题为`{directory}_{分支名}`。
- 图内文本显示该分支的KS统计量和Anderson-Darling统计量。

---

## 自定义配置

### 1. 过滤不需要的分支
### 2. 修改组合权重公式

`ks_andeson.py`中计算权重的公式为：
```python
weight = 0.5 * ks[i]/ks_m + 0.5 * andeson[i]/and_max
```
您可以调整系数或替换为其他组合方式，例如只使用KS、只使用AD，或采用不同的归一化方法。

### 3. 调整绘图数量

`weight_sort.py`中通过切片`[:20]`控制绘制的分支数量。若要绘制更多或更少，修改该数字即可。

---

## 注意事项

1. **数据类型**  
   - 非数值型分支会被跳过。
   - 包含NaN值的条目会被剔除。
   - 若分支数据为空，也会跳过。

2. **Anderson-Darling检验异常**  
   当两个样本完全相同或样本量过小时，`scipy.stats.anderson_ksamp`可能抛出`ValueError`，此时该分支会被静默跳过（不输出到CSV）。
---
