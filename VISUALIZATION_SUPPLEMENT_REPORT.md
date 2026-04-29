# Visualization Supplement Report

生成时间：2026-04-29

## 1. Summary

本报告汇总当前第 6 章的可视化补充图。这些图用于解释实验三的非齐次 Dirichlet 空间结构、实验六的 true Helmholtz near-resonance 小分母机制，以及多模态 manufactured solution 上 FA/FFT9 的可视化验证。

这些内容未新增核心实验、未恢复旧图，也未改变第 6 章 6 个核心实验结构。

本次优化内容：

- 温度场图中 analytical/exact field 与 FFT9 numerical solution 使用相同 `vmin/vmax`。
- 温度场误差图标题加入最大误差 `max=...`。
- near-resonance heatmap 改为 resonance strength `-\log10|d_{p,q}|`。
- near-resonance heatmap 标题改为 selected modal pair `(2,3)/(3,2), delta=1e-3`。
- near-resonance heatmap 使用统一文本框列出 `min modes: (2,3), (3,2)`，并加入 `p,q<=10` zoomed inset。
- 多模态 manufactured solution 补充图已新增并验证，脚本为 `exp06_complex_manufactured_visualization.py`，但不计入核心实验编号。

## 2. Generated figures

已生成并同步复制以下 PNG：

| Figure | Code-side path | Thesis path | Purpose |
|---|---|---|---|
| Nonhomogeneous Dirichlet temperature field | `code/python/experiments/figures/exp02_temperature_field_comparison.png` | `thesis/figures/exp02_temperature_field_comparison.png` | exact 与 FFT9 数值解使用相同色标，展示解析温度场、误差分布最大值与 `y=0.5` 截线 |
| True Helmholtz denominator heatmap | `code/python/experiments/figures/exp05_denominator_heatmap.png` | `thesis/figures/exp05_denominator_heatmap.png` | 展示 selected modal pair `(2,3)/(3,2)` 下的 resonance strength `-log10|d_{p,q}|`，数值越大表示分母越小、近共振越强；图内包含 `p,q<=10` zoomed inset |
| Complex manufactured fields | `code/python/experiments/figures/exp06_complex_manufactured_fields.png` | `thesis/figures/exp06_complex_manufactured_fields.png` | 展示多模态 manufactured solution 的 exact field、FA/FFT9 数值解、两幅相同对数色标的误差图和 `y=0.5` 截线 |
| Complex manufactured convergence | `code/python/experiments/figures/exp06_complex_manufactured_convergence.png` | `thesis/figures/exp06_complex_manufactured_convergence.png` | 展示同一多模态补充问题中 FA 二阶与 FFT9 四阶收敛 |

生成命令：

```powershell
cd code/python
python -m experiments.generate_visualization_supplements
```

## 3. Thesis changes

已修改 `thesis/chapters/6_experiments.tex`：

- 在实验三内部加入“可视化补充”段落和 `fig:exp02_temperature_field_comparison`。
- 在实验六内部加入“可视化补充”段落和 `fig:exp05_denominator_heatmap`。
- 在实验三内部加入“多模态 manufactured solution 可视化补充”段落，以及 `fig:exp06_complex_fields` 和 `fig:exp06_complex_convergence`。
- 图注均说明其用途是空间结构或频域机制可视化，不作为新增收敛阶、性能比较或预处理器证据。
- denominator heatmap 的正文说明和图注已同步为 `-\log_{10}|d_{p,q}|`，并明确“数值越大表示分母越小、近共振越强”。

未新增 `实验七`，也未恢复旧 heatmap、restart、FACR 参数图或旧 timing 图。

同步更新：

- `README.md`：记录核心实验产物与可视化补充产物，并同步最终 PDF 为 56 pages。
- `FINAL_POLISH_REPORT.md`：说明补充图不改变 6 个核心实验结构，并同步最终 PDF 为 56 pages。

## 4. Verification

| Check | Command | Result |
|---|---|---|
| Figure generation | `cd code/python; python -m experiments.generate_visualization_supplements` | PASS |
| Complex manufactured script | `cd code/python; python -m experiments.exp06_complex_manufactured_visualization` | PASS: FA order 1.991/2.004, FFT9 order 3.988/3.989 |
| pytest | `cd code/python; python -m pytest -q` | PASS: 103 passed, 2 skipped, 2 warnings |
| XeLaTeX/BibTeX | `cd thesis; xelatex; bibtex; xelatex; xelatex` | PASS: `main.pdf` generated, 56 pages |
| LaTeX log scan | `Select-String thesis\main.log ...` | PASS: no fatal error, no undefined citation/reference, no multiply-defined label |

Non-blocking warnings:

- pytest 仍有既有 Neumann Poisson compatibility warnings，属于预期测试行为。
- XeLaTeX 仍有 SimSun 字体替代 warning 和 hyperref bookmark token warnings，不影响 PDF 生成、交叉引用或正文渲染。

## 5. Six-experiment structure

第 6 章仍只包含以下 6 个核心实验：

1. FFT 解与 sparse direct 解一致性；
2. 二阶与四阶收敛阶验证；
3. 非齐次 Dirichlet 边界验证；
4. Neumann 与 mixed 边界验证；
5. modified/true Helmholtz 谱指标与 Gaussian RHS 下 GMRES 行为；
6. true Helmholtz near-resonance 扫描。

所有补充图均归入实验三或实验六的“可视化补充”，不构成新增实验。

## 6. Final status

READY TO SUBMIT
