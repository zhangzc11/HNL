# tau=0ps GBDT AUC 结果整理

## 结果概览

- 本次读取到的 `auc_summary_sorted.csv` 共包含 **29** 条实验结果。
- 其中基准实验 `all_30_variables` 的 AUC 为 **0.979383**。
- 当前这份结果里，最高 AUC 的实验是 **drop_NuR_M**，AUC = **0.980373**。
- 最低 AUC 的实验是 **drop_MuNuR_PROBNN_MU**，AUC = **0.975503**。
- 最高与最低之间的 AUC 差值为 **0.004870**。

> 说明：你原先希望得到 31 个实验，但这份 CSV 中实际只有 **29** 条记录，因此当前文档只基于这 **29** 条结果整理。

## 基准实验

- 实验名：`all_30_variables`
- 使用变量数：**28**
- AUC：**0.979383**
- ROC 文件：`01_all_30_variables.png`
- Score 分布文件：`01_all_30_variables_score.png`

## AUC 排名前 5

|tag|removed_variable|n_variables_used|auc|
|---|---|---|---|
|drop_NuR_M|NuR_M|27|0.980373|
|drop_NuR_log_MIN_PT|NuR_log_MIN_PT|27|0.979390|
|all_30_variables|NONE|28|0.979383|
|drop_Mu_TRACKPT|Mu_TRACKPT|27|0.979362|
|drop_W_MIN_PT|W_MIN_PT|27|0.979358|

## AUC 排名后 5

|tag|removed_variable|n_variables_used|auc|
|---|---|---|---|
|drop_Mu_PROBNN_MU|Mu_PROBNN_MU|27|0.977670|
|drop_Mu_HCALEOP|Mu_HCALEOP|27|0.977615|
|drop_Mu_PROBNN_GHOST|Mu_PROBNN_GHOST|27|0.977515|
|drop_W_M|W_M|27|0.977108|
|drop_MuNuR_PROBNN_MU|MuNuR_PROBNN_MU|27|0.975503|

## 去掉某个变量后，AUC 下降最多的前 10 个变量

这部分可以理解为“对模型最重要的变量候选”。因为一旦去掉它，AUC 相比 `all_30_variables` 下降得更多。

|removed_variable|auc|auc_drop_vs_all30|
|---|---|---|
|MuNuR_PROBNN_MU|0.975503|0.003880|
|W_M|0.977108|0.002275|
|Mu_PROBNN_GHOST|0.977515|0.001868|
|Mu_HCALEOP|0.977615|0.001768|
|Mu_PROBNN_MU|0.977670|0.001713|
|nFTClusters|0.978340|0.001043|
|MuNuR_PROBNN_GHOST|0.978387|0.000996|
|W_MmuWmuN|0.978436|0.000947|
|Mu_ELECTRONSHOWEREOP|0.978773|0.000610|
|Mu_log_PT|0.979107|0.000276|

## 去掉某个变量后，AUC 反而上升的实验

如果去掉某个变量后 AUC 上升，通常说明这个变量对当前模型没有帮助，甚至可能带来冗余或噪声。

|removed_variable|auc|delta_auc_vs_all30|
|---|---|---|
|NuR_M|0.980373|0.000990|
|NuR_log_MIN_PT|0.979390|0.000007|

## 简要结论

1. `all_30_variables` 的基准 AUC 已经达到 **0.979383**，说明当前这组变量整体区分能力很强。
2. 从 “AUC 下降最多” 的结果看，比较关键的变量包括：`MuNuR_PROBNN_MU`、`W_M`、`Mu_PROBNN_GHOST`、`Mu_HCALEOP`、`Mu_PROBNN_MU`。
3. 从 “去掉变量后 AUC 上升” 的结果看，去掉 `NuR_M` 后 AUC 升到 **0.980373**，去掉 `NuR_log_MIN_PT` 后也略有上升，说明这两个变量在当前配置下可能存在一定冗余。
4. 当前最高 AUC 和基准 AUC 的差只有 **0.000990**，所以虽然存在轻微改进，但总体变化不算特别大。
5. 当前最低 AUC 仍有 **0.975503**，说明即使去掉部分关键变量，模型整体仍保持了较强区分能力。

## 原始结果文件

- 源文件：`auc_summary_sorted.csv`
