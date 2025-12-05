# sequence_similarity_score2 修订日志（相对 sequence_similarity_score2-old）

概览
- 文件：utils/sequence_similarity_score2.py 相对 utils/sequence_similarity_score2-old.py 的变更摘要
- 目标：纠正未配对位点的评分语义、减少未配对位点对总分的“拉高”效应，并提供更灵活的加权与更高的性能。

修复
- 未配对位点处理语义修复：
  - 旧版在无预测结构矩阵时，会基于“序列互补性”推断伪配对，导致参考结构标记为未配对的 (i,j) 位置被误罚（相同序列得分 < 1）。
  - 新版在“仅序列”评分场景下，参考结构的未配对位点一律计为正确（每个未配对位置 +1），不再因互补性被误罚，确保相同序列基础配对分可达 1.0。

改动
- 归一化策略调整：
  - base_pair_score 默认聚焦“成对位点”，按 4 × paired_positions 归一化。
  - 未配对位点分数独立统计（normalized_unpaired_score），默认不并入主分，避免未配对位点数量过多导致得分失真。
- 加权选项：
  - 新增 unpaired_weight（默认 0.0），允许按需引入未配对组件到 base_pair_score 中。
  - CLI 增加 --unpaired-weight 参数，便于命令行控制该权重。
- 输出指标与日志：
  - 新增 normalized_paired_score、normalized_unpaired_score、raw_score_paired、raw_score_unpaired、max_possible_score_paired/UNPAIRED 等统计项，提升可解释性。
- 性能与鲁棒性：
  - 新增全向量化评分路径（_score_loop_vectorized），无需 Numba 也有较高性能。
  - 引入上三角索引缓存（按 n 与 min_separation 缓存）与序列编码缓存，降低训练循环重复开销。
  - Numba JIT 失败时自动回退至向量化路径；统一大小写处理避免无效碱基索引。

行为变化与影响
- base_pair_score 更强调结构正确性（成对位点），在未配对位点占比高的数据上，分数可能低于旧版的“未配对拉高”结果。
- 在无预测结构矩阵的场景下：
  - 未配对位点默认记为正确（unpaired_only = 1.0），主分仍只看成对位点；如参考结构无成对位点，应结合序列相似度（如编辑距离/恢复率）评估。

使用建议
- 若需兼顾未配对位点贡献，可设置 --unpaired-weight（如 0.1）以小幅并入未配对组件。
- 根据任务选择 min_separation（默认 4），避免近邻伪配对对结构评分的干扰。

兼容性
- 保持原有评分接口与结果字典键值的主体结构；新增键不影响既有逻辑。
- CLI 向后兼容，未指定 --unpaired-weight 时行为为“仅成对位点”评分并附带未配对统计输出。