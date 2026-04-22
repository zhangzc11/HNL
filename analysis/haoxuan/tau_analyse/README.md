# tau_feature 筛选说明（M5GeV, avg tau）

## 数据来源
- 输入相关系数矩阵: coree/M5GeV/M_5GeV_tau_avgps_corr.csv
- 矩阵规模: 50 x 50（50 个 branch）

## 目标
- 选出“可用于独立训练”的代表变量集合。
- 要求所选变量之间线性相关性弱。

## 筛选方法
1. 以绝对相关系数阈值 |r| >= 0.7 建图：
   - 每个 branch 是一个节点。
   - 若两变量 |r| >= 0.7，则连接一条边。
2. 对图做连通分量分组：
   - 每个分组内变量彼此存在较强相关链路。
3. 每组保留 1 个代表变量：
   - 组内选择“最大相关最小”的变量（先最小化组内 max|r|，再最小化组内 mean|r|）。
4. 输出代表变量表到 tau_feature.csv。

## 类型定义
- type 0: 与其他变量的最大绝对相关系数 max|r| > 0.9
- type 1: 0.8 < max|r| <= 0.9
- type 2: 0.7 < max|r| <= 0.8
- type 3: max|r| <= 0.7（补充类型，表示整体更独立）

## 结果概览
- 原始变量数: 50
- 分组数: 26
- 保留变量数: 26（每组 1 个）
- 保留变量两两之间最大绝对相关系数: 0.596760
- 最相关的一对保留变量: NuR_END_VZ 与 NuR_BPV_LTIME

说明：最大值 0.596760 < 0.7，满足“保留变量之间相关性不强”的目标。

## 强相关分组（group_size > 1）
- G01 (11): MuNuR_ISDOWNSTREAM, MuNuR_ISLONG, MuNuR_NHITS, MuNuR_NVPHITS, MuNuR_NVPLAYERS, MuNuR_POSITION_STATEAT_FirstMeasurement_Z, MuNuR_TRACKHASVELO, MuNuR_TRACKHISTORY, MuNuR_TRCHI2, MuNuR_TRCHI2DOF, MuNuR_TRNDOF
  - 代表变量: MuNuR_TRCHI2DOF
- G02 (3): NuR_END_VX, NuR_OWNPV_FDVECX, NuR_OWNPV_VDX
  - 代表变量: NuR_END_VX
- G03 (3): NuR_END_VZ, NuR_OWNPV_FDVECZ, NuR_OWNPV_VDZ
  - 代表变量: NuR_END_VZ
- G04 (3): NuR_OWNPV_DIRA, NuR_OWNPV_ETA, NuR_OWNPV_FDIRZ
  - 代表变量: NuR_OWNPV_ETA
- G05 (2): Jet1_ConSumTRACKHASVELO, Jet1_ConSumTRACKISLONG
  - 代表变量: Jet1_ConSumTRACKHASVELO
- G06 (2): MuNuR_MINIP, MuNuR_OWNPVIP
  - 代表变量: MuNuR_MINIP
- G07 (2): MuNuR_OWNPVIPCHI2, NuR_MIN_OWNPVIPCHI2
  - 代表变量: MuNuR_OWNPVIPCHI2
- G08 (2): NuR_BPV_IPCHI2, NuR_OWNPVIPCHI2
  - 代表变量: NuR_BPV_IPCHI2
- G09 (2): NuR_BPV_LTIME, NuR_OWNPV_LTIME
  - 代表变量: NuR_BPV_LTIME
- G10 (2): NuR_MINIP, NuR_OWNPVIP
  - 代表变量: NuR_MINIP
- G11 (2): NuR_OWNPV_FD, NuR_OWNPV_VDRHO
  - 代表变量: NuR_OWNPV_FD
- G12 (2): NuR_OWNPV_FDVECY, NuR_OWNPV_VDY
  - 代表变量: NuR_OWNPV_FDVECY

## 从强相关分组提取的学习特征（12个）
- MuNuR_TRCHI2DOF
- NuR_END_VX
- NuR_END_VZ
- NuR_OWNPV_ETA
- Jet1_ConSumTRACKHASVELO
- MuNuR_MINIP
- MuNuR_OWNPVIPCHI2
- NuR_BPV_IPCHI2
- NuR_BPV_LTIME
- NuR_MINIP
- NuR_OWNPV_FD
- NuR_OWNPV_FDVECY


MuNuR_TRCHI2DOF, NuR_END_VX, NuR_END_VZ, NuR_OWNPV_ETA, Jet1_ConSumTRACKHASVELO, MuNuR_MINIP, MuNuR_OWNPVIPCHI2, NuR_BPV_IPCHI2, NuR_BPV_LTIME, NuR_MINIP, NuR_OWNPV_FD, NuR_OWNPV_FDVECY
