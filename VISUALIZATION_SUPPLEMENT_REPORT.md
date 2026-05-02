# Visualization Supplement Report

生成时间：2026-04-29

## 1. Summary

本报告汇总当前第 6 章的可视化补充图。这些图用于解释实验一中的非齐次 Dirichlet 空间结构和多模态 Dirichlet sanity check、实验二的 Neumann/mixed 代表性解场与误差分布，以及实验五的 true Helmholtz near-resonance 小分母机制。

这些内容未新增核心实验、未恢复旧图，也未改变第 6 章 5 个核心实验结构。

本次优化内容：

- 温度场图中 analytical/exact field 与 FFT9 numerical solution 使用相同 `vmin/vmax`。
- 温度场误差图标题加入最大误差 `max=...`。
- 实验五正文主图保留 `exp05_near_resonance_summary.png`，同时展示解放大、GMRES capped/not-converged 标记和解场对比。
- near-resonance heatmap 改为 resonance strength `-\log10|d_{p,q}|`。
- near-resonance heatmap 标题改为 selected modal pair `(2,3)/(3,2), delta=1e-3`。
- near-resonance heatmap 使用统一文本框列出 `min modes: (2,3), (3,2)`，并加入 `p,q<=10` zoomed inset。
- 多模态 manufactured solution 补充图已新增并验证，脚本为 `exp06_complex_manufactured_visualization.py`，但不计入核心实验编号。
- 实验四 Neumann/mixed 验证的 `exp03_neumann_mixed_summary.png` 已优化为 1x3 横版三子图，并新增 `exp03_neumann_mixed_fields.png` 展示代表性解场与误差分布；收敛阶与边界残差结论仍由 summary 图和 CSV 支撑。

## 2. Generated figures

已生成并同步复制以下 PNG：

| Figure | Code-side path | Thesis path | Purpose |
|---|---|---|---|
| Nonhomogeneous Dirichlet temperature field | `code/python/experiments/figures/exp02_temperature_field_comparison.png` | `thesis/figures/exp02_temperature_field_comparison.png` | exact 与 FFT9 数值解使用相同色标，展示解析温度场、误差分布最大值与 `y=0.5` 截线 |
| Neumann/Mixed fields | `code/python/experiments/figures/exp03_neumann_mixed_fields.png` | `thesis/figures/exp03_neumann_mixed_fields.png` | 展示 Pure Neumann 与 Mixed `(N,D)` 的 representative exact/numerical fields 和 log10 误差分布；不替代 1x3 summary 图和 CSV 收敛证据 |
| True Helmholtz near-resonance summary | `code/python/experiments/figures/exp05_near_resonance_summary.png` | `thesis/figures/exp05_near_resonance_summary.png` | 正文实验五主图，展示简并模态对 `(2,3)/(3,2)` 下的解放大曲线、GMRES capped/not-converged 标记，以及同一 detuning branch 下的 raw field 幅值对比 |
| True Helmholtz denominator heatmap | `code/python/experiments/figures/exp05_denominator_heatmap.png` | `thesis/figures/exp05_denominator_heatmap.png` | 备份可视化，展示 selected modal pair `(2,3)/(3,2)` 下的 resonance strength `-log10|d_{p,q}|`；不作为正文主图 |
| Complex manufactured fields | `code/python/experiments/figures/exp06_complex_manufactured_fields.png` | `thesis/figures/exp06_complex_manufactured_fields.png` | 展示多模态 manufactured solution 的 exact field、FA/FFT9 数值解、两幅相同对数色标的误差图和 `y=0.5` 截线 |
| Complex manufactured convergence | `code/python/experiments/figures/exp06_complex_manufactured_convergence.png` | `thesis/figures/exp06_complex_manufactured_convergence.png` | 展示同一多模态补充问题中 FA 二阶与 FFT9 四阶收敛 |

生成命令：

```powershell
cd code/python
python -m experiments.generate_visualization_supplements
python -m experiments.exp03_neumann_mixed
python -m experiments.exp05_true_helmholtz_resonance
```

## 3. Thesis changes

已修改 `thesis/chapters/6_experiments.tex`：

- 在实验一内部加入非齐次 Dirichlet 子情形和 `fig:exp02_temperature_field_comparison`。
- 实验二正文主图使用 `fig:exp03_neumann_mixed` 1x3 summary 图，并加入 `fig:exp03_neumann_mixed_fields` 作为代表性解场可视化补充。
- 实验五正文主图使用 `fig:exp05_near_resonance_summary`；`fig:exp05_denominator_heatmap` 保留为备份图像但不再作为正文主图。
- 在实验一内部加入“多模态 Dirichlet 可视化检验”段落，以及 `fig:exp06_complex_fields` 和 `fig:exp06_complex_convergence`。
- 图注均说明其用途是空间结构或频域机制可视化，不作为新增收敛阶、性能比较或预处理器证据。
- 实验五 summary figure 的正文说明和图注已同步为“解放大 + GMRES capped/not-converged + raw field 幅值对比”，并明确绝对幅值由 panel (a) 与子图标题中的 `max|u|` 给出。

未新增 `实验七`，也未恢复旧 heatmap、restart、FACR 参数图或旧 timing 图。

同步更新：

- `README.md`：记录核心实验产物与可视化补充产物，并同步最终 PDF 为 62 pages。
- `FINAL_POLISH_REPORT.md`：说明补充图不改变 5 个核心实验结构，并同步最终 PDF 为 62 pages。

## 4. Verification

| Check | Command | Result |
|---|---|---|
| Figure generation | `cd code/python; python -m experiments.generate_visualization_supplements` | PASS |
| EXP03 V2 script | `cd code/python; python -m experiments.exp03_neumann_mixed` | PASS: regenerated CSV, summary figure, and fields figure |
| Complex manufactured script | `cd code/python; python -m experiments.exp06_complex_manufactured_visualization` | PASS: FA order 1.991/2.004, FFT9 order 3.988/3.989 |
| pytest | `cd code/python; python -m pytest -q` | PASS: 103 passed, 2 skipped, 2 warnings |
| XeLaTeX/BibTeX | `cd thesis; xelatex; bibtex; xelatex; xelatex` | PASS: `main.pdf` generated, 62 pages |
| LaTeX log scan | `Select-String thesis\main.log ...` | PASS: no fatal error, no undefined citation/reference, no multiply-defined label |

Non-blocking warnings:

- pytest 仍有既有 Neumann Poisson compatibility warnings，属于预期测试行为。
- XeLaTeX 仍有 SimSun 字体替代 warning 和 hyperref bookmark token warnings，不影响 PDF 生成、交叉引用或正文渲染。

## 5. Five-experiment structure

第 6 章整理为以下 5 个核心实验；exp00 作为 implementation validation 保留，不计入核心实验编号：

1. Dirichlet 离散格式收敛性验证；
2. Neumann 与 mixed 边界处理；
3. 精度--成本对比与复杂度实证；
4. Modified/True Helmholtz 的谱结构与 GMRES 行为；
5. True Helmholtz 近共振模态放大。

所有补充图均归入实验一、实验二或实验五的可视化/机制说明，不构成新增实验。

## 6. Final status

READY TO SUBMIT
