# Experiment Audit Report

本次审查为只读审查。没有恢复旧图，没有新增实验，也没有修改论文、代码、CSV 或 PNG。

## Overall Verdict

第 6 章实验主体可以支撑论文当前的核心结论：六个核心实验均有对应脚本、CSV 和正文图表或表格支撑；正文表格中的代表数值与 CSV 可对应；第 6 章已明确多个“不证明什么”的边界。

需要处理的是少数表述边界问题：不要把无预处理 GMRES 观察写成完整性能研究，不要把 exp05 的 near-resonance 失败直接写成已经验证 shifted-Laplacian 等预处理器的必要性，也不要让“高效”一类措辞看起来来自未做的计时实验。

## A-Class Issues

无 A 类问题。未发现核心结论完全缺少脚本、CSV、PNG 或正文表格支撑的情况。

## B-Class Issues

### B-01: Exp04 GMRES 图需要标注 capped or failed rows

- 位置：`thesis/chapters/6_experiments.tex:217-218`；`code/python/experiments/exp04_modified_vs_true.py:331-349`
- 证据：`exp04_modified_vs_true.csv` 中存在 `gmres_success=False` 且迭代数达到上限的行；当前图 `exp04_gmres_iters_vs_sigma.png` 没有显式区分未收敛点。
- 风险：读者可能把达到上限的点理解成成功收敛的迭代次数。
- 建议：不新增实验。只在图注/正文中说明 capped or not converged markers，或后续若改图则用不同 marker 区分。

### B-02: “通常需要 shifted-Laplacian 等预处理技术”应降级

- 位置：`thesis/chapters/6_experiments.tex:287-288`；`thesis/chapters/7_conclusion.tex:78`
- 证据：exp05 支撑的是 near-resonance 下无预处理 GMRES 在给定上限内未收敛；它没有测试 shifted-Laplacian、MINRES 或任何预处理器。
- 风险：把“未预处理基线失败”过度升级为“本文实验验证了某类预处理技术必要”。
- 建议：改为“提示实际大规模或近共振 true Helmholtz 求解通常需要考虑预处理；本文未测试预处理器”。

## C-Class Issues

### C-01: 第 7 章“高效”措辞不应暗示第 6 章做了计时结论

- 位置：`thesis/chapters/7_conclusion.tex:58-60`
- 证据：第 6 章明确不保留无统计重复计时和旧 timing/complexity 图作为核心证据。
- 风险：轻微。算法章节可以有复杂度分析，但实验章没有独立计时证据。
- 建议：写成“具备 FFT 快速直接求解结构；第 6 章验证离散一致性与精度”。

### C-02: Neumann 方法的“优雅、无缝、高效性”可更实验化

- 位置：`thesis/chapters/7_conclusion.tex:83-85`
- 证据：exp03 支撑的是测试的五点 Neumann/mixed 情形下二阶收敛和实现一致性。
- 风险：轻微。形容词不影响数学结论，但可以更稳。
- 建议：改为“在测试情形中支持统一求解并呈现预期二阶收敛”。

### C-03: Exp05 CSV baseline rows 有未来误引风险

- 位置：`code/python/experiments/results/exp05_resonance.csv`
- 证据：正文 near-resonance 表格只使用 `delta` 非空的扫描行，这些行支撑结论；baseline rows 不参与正文结论。
- 风险：如果未来引用 baseline rows，需要确认 `gmres_success` 字段是否真实反映求解器结果。
- 建议：正文继续只引用 near-resonance rows；若修代码，修 CSV metadata 即可，不新增实验。

## Six-Experiment Inventory Status

- exp00: `exp00_fft_vs_sparse.py` + `exp00_fft_vs_sparse.csv` + 正文表格。支撑离散 FFT-vs-sparse 一致性。
- exp01: `exp01_convergence.py` + `exp01_convergence.csv` + `exp01_convergence.png` + 正文表格。支撑二阶与四阶收敛阶。
- exp02: `exp02_nonhom_bc.py` + `exp02_nonhom_bc.csv` + `exp02_nonhom_bc.png` + 正文表格。支撑非齐次 Dirichlet 验证。
- exp03: `exp03_neumann_mixed.py` + `exp03_neumann_mixed.csv` + `exp03_neumann_mixed.png` + 正文表格。支撑五点 Neumann/mixed 验证。
- exp04: `exp04_modified_vs_true.py` + `exp04_modified_vs_true.csv` + 三张 PNG + 正文表格。支撑 Gaussian RHS 下无预处理 GMRES 观察和谱指标。
- exp05: `exp05_true_helmholtz_resonance.py` + `exp05_resonance.csv` + `exp05_resonance.png` + 正文表格。支撑 true Helmholtz near-resonance 扫描。

结论：第 6 章当前确认为六个核心实验结构。

## Legacy Structure And Figure Status

- 未发现第 6 章正文引用“实验七”到“实验十六”作为实验结构。
- 未发现第 6 章正文引用旧 heatmap、3D surface、restart plot、FACR parameter plot 或 old timing plot 作为证据。
- `thesis/figures/` 中仍存在若干历史 PNG，例如 heatmap、surface、restart、FACR parameter、complexity/timing-like figures；第 6 章明确说明这些不参与最终结论。按用户限制，不建议恢复或删除它们。

## CSV And Table Checks

- exp00: 正文 `n=65` 的误差上界与 CSV 对应行一致，例如 Dirichlet Poisson 最大约 `1.91e-14`，Dirichlet true Helmholtz 最大约 `2.39e-14`。
- exp01: 正文多项式轨道 `sigma=10` 表格与 CSV 一致，五点约二阶，FFT9 约四阶。
- exp02: 正文 `n=65` 表格与 CSV 一致，且脚本使用非齐次 Dirichlet shifted manufactured solution。
- exp03: 正文 FA 代表行与 CSV 中 `n=65` FA 行一致；CSV 还包含 CR/FACR 行。
- exp04: 正文 `n=33` 表格与 CSV 一致，包括 `sigma=-50` 的 near-resonance capped/failed 行。
- exp05: 正文 representative rows 与 CSV 中 `delta=-0.1,-0.01,-0.001,-0.0001` 的 near-resonance 行一致。

## Final Assessment

第 6 章实验整体证据链完整。终稿前优先做文字降级，不需要新增实验，也不需要恢复旧图。若继续打磨，最值得处理的是 exp04 未收敛点的图注说明，以及第 6/7 章关于预处理和效率的保守表述。
