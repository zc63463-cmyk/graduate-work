# Recent Experiment Update Report

## 1. Scope

本报告记录最近一轮实验图表与正文修复，覆盖：

- 实验二 / exp03：Neumann 与 mixed 边界处理图表去冗余；
- 实验四 / exp07 + exp04：谱分母小分母风险图与 GMRES 行为展示优化；
- 实验五 / exp05：true Helmholtz near-resonance 模态放大与 GMRES 停滞图表优化。

本轮不新增核心实验、不改变数值算法、不改变第 6 章五个核心实验结构。

## 2. Experiment 3: Neumann/Mixed Boundary Visualization

### Files

- `code/python/experiments/exp03_neumann_mixed.py`
- `code/python/experiments/figures/exp03_neumann_mixed_summary.png`
- `code/python/experiments/figures/exp03_neumann_mixed_fields.png`
- `thesis/figures/exp03_neumann_mixed_summary.png`
- `thesis/figures/exp03_neumann_mixed_fields.png`
- `thesis/chapters/6_experiments.tex`

### Update

原 2x3 六子图 summary 图被压缩为 1x3 横版三子图：

1. `(a) Linf convergence, FA solver`：仅保留 FA solver 的四类 case 收敛曲线与 `O(h^2)` 参考线，避免 FA/CR/FACR-like 完全重叠造成视觉冗余。
2. `(b) Boundary constraint residuals`：只保留一条代表性 Neumann flux 曲线，并在 legend 中说明 `Neumann flux (all cases coincide)`；同时保留 Mixed `(N,D)` 的 Dirichlet boundary residual。
3. `(c) Neumann Poisson zero-mode check`：集中展示 Neumann Poisson 的 weighted mean 与 removed mean offset。

### Rationale

该图现在分别回答三个问题：

- 五点边界求解器是否保持二阶收敛；
- Neumann flux 与 mixed Dirichlet 边界约束是否闭合；
- Neumann Poisson 零模态归一化是否合理。

Fields 图保持不变，仍作为代表性解场和误差分布的可视化补充。

## 3. Experiment 4: Spectral Denominator Risk And GMRES Behavior

### Files

- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.png`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `thesis/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/chapters/6_experiments.tex`

### Update

图 10 不再使用全域 `d_{p,q}` 或全域 risk heatmap，而改为：

1. 上排：低频窗口 `p,q <= 12` 的小分母风险图
   \[
   R_{p,q}=-\log_{10}|d_{p,q}|.
   \]
2. 下排：三种 case 的前 20 个最小 `|d_{p,q}|` 排序曲线。

near-resonance case 中使用绿色 marker 标出目标简并模态 `(2,3)` 与 `(3,2)`；排序图显示 near-resonance 的前两个最小分母正是该模态对，且降至 `10^{-2}` 量级。

### Text Fixes

正文补充说明：

- 当 `|d_{p,q}| > 1` 时，`R_{p,q}` 为负值；负值只表示分母较大、小分母风险较低，并不代表“负风险”。
- true-away 情形虽然已出现分母正负分裂，但最近模态距离仍为 `O(10^1)`，因此没有形成显著小分母热点。
- risk 图只是 small-denominator risk visualization，不是条件数图。

### Rationale

全域 `127 x 127` 图会被大范围高频低风险背景淹没；低频 zoom + sorted denominator plot 更直接展示：

- modified Helmholtz 分母远离零；
- true-away 已 sign-changing 但仍远离小分母；
- true-near 在 `(2,3)/(3,2)` 处产生局部小分母热点。

## 4. Experiment 5: Near-Resonance Modal Amplification

### Files

- `code/python/experiments/exp05_true_helmholtz_resonance.py`
- `code/python/experiments/figures/exp05_multimode_resonance_summary.png`
- `code/python/experiments/figures/exp05_multimode_resonance_summary.pdf`
- `code/python/experiments/figures/exp05_resonance_gmres_history.png`
- `code/python/experiments/figures/exp05_resonance_gmres_history.pdf`
- `thesis/figures/exp05_multimode_resonance_summary.png`
- `thesis/figures/exp05_multimode_resonance_summary.pdf`
- `thesis/figures/exp05_resonance_gmres_history.png`
- `thesis/figures/exp05_resonance_gmres_history.pdf`
- `thesis/chapters/6_experiments.tex`

### Update

`exp05_multimode_resonance_summary` 被重画为聚焦 “目标模态放大 + GMRES 停滞” 的 2x2 图：

1. `(a) target pair amplification`：分别绘制 `|\hat u_{2,3}|` 与 `|\hat u_{3,2}|`，使用不同颜色，log scale，并给出 `C/|\delta|` 参考线和约 `1000x` 放大标注。
2. `(b) amplification ratio by target mode`：替代原先冗余的 target subspace energy 展示，比较不同目标模态集合相对于 `|\delta|=10^{-1}` 的放大比例。
3. `(c) dominant projection check`：用 full solution norm 与 dominant projection norm 的散点图压缩原 full/projection/difference 的重复信息。
4. `(d) GMRES(30) capped/stagnation`：横轴为 target modal amplification，纵轴为 final absolute residual；red x 表示 capped / not converged，并保留 absolute tolerance 参考。

GMRES residual history 图中的 capped/not-converged 末端点也统一使用红叉标记。

### Text Fixes

第 6 章对应图注和正文改为：

- 强调 `(2,3)` 与 `(3,2)` 的目标模态系数均按 `1/|\delta|` 放大；
- 删除对逐模态能量变化的重复解释；
- 将主要结论收束为“模态放大与无预处理 GMRES(30) 停滞相关”；
- 继续限定结论只适用于当前 Dirichlet 五点离散谱与无预处理 GMRES(30)。

## 5. Validation

### Commands Run

Figures were regenerated from existing CSV files only:

```powershell
$env:PYTHONPATH='code/python'
python - <<'PY'
from pathlib import Path
import pandas as pd
from experiments import exp05_true_helmholtz_resonance as exp05

results_dir = Path(exp05.get_results_dir())
figures_dir = exp05.get_figures_dir()

multimode = pd.read_csv(results_dir / exp05.PATCH5_MULTIMODE_CSV)
exp05.plot_multimode_summary(multimode, figures_dir)

history = pd.read_csv(results_dir / exp05.PATCH5_GMRES_HISTORY_CSV)
exp05.plot_resonance_gmres_history(history, figures_dir)
PY
```

Python tests:

```powershell
cd code/python
python -m pytest -q
```

LaTeX:

```powershell
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

### Results

- `pytest`: `103 passed, 2 skipped, 2 warnings`
- XeLaTeX/BibTeX/XeLaTeX/XeLaTeX: success
- Output PDF: `thesis/main.pdf`, 66 pages
- LaTeX log scan:
  - no fatal error
  - no undefined reference
  - no multiply-defined label
  - no overfull hbox
- `git diff --check`: no whitespace errors; only CRLF normalization warnings.

## 6. Final Status

READY TO USE

The latest experiment visual fixes are integrated at the source and figure level. They improve readability without adding experiments, changing numerical methods, or expanding the thesis claims.
