# PATCH3 Report

## 1. Scope

本轮完成 Patch3，并顺手完成 Patch2b 两个小修：

- 新增 `exp07_spectral_denominator_maps.py`，生成五点 Dirichlet 频域分母热图，用于可视化支撑第 4 章关于 modified Helmholtz 正定性、true Helmholtz 不定性与小分母近共振风险的理论分析。
- 修复 `exp06_accuracy_cost.py` 中 error-time 图的 GMRES 未收敛标记：由红色叉号改为黑色叉号，legend 明确写 `GMRES30 not converged`。
- 在第 6 章 exp06 精度--成本对比正文和 caption 中明确 timing scope，并说明该图是当前实现口径下的整体基准趋势，不是严格 kernel-to-kernel benchmark；由于 GMRES setup 被排除，该比较没有刻意偏向 FFT 方法。
- `exp06_accuracy_cost.py` 与 `exp07_spectral_denominator_maps.py` 均已补充 direct-file execution 兼容，支持从仓库根目录用 `PYTHONPATH=code/python python code/python/experiments/<script>.py` 形式运行。

## 2. Files Modified

- `code/python/experiments/exp06_accuracy_cost.py`
- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 3. Files Generated Or Updated

- `code/python/experiments/results/exp07_spectral_denominator_summary.csv`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.png`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `thesis/figures/exp07_spectral_denominator_heatmaps.pdf`
- `code/python/experiments/figures/exp06_accuracy_cost_error_time.png`
- `code/python/experiments/figures/exp06_accuracy_cost_error_time.pdf`
- `thesis/figures/exp06_accuracy_cost_error_time.png`
- `thesis/figures/exp06_accuracy_cost_error_time.pdf`
- `code/python/experiments/results/exp06_accuracy_cost.csv`
- `PATCH3_REPORT.md`

## 4. Patch2b Details

`exp06_accuracy_cost_error_time` 已将未收敛 GMRES30 点改为黑色 `x` 标记，并在 legend 中写明 `GMRES30 not converged`，避免与 FFT9 红色曲线混淆。

第 6 章已补充 timing scope 说明：

- direct solvers：计时为 solver call，包含内部 RHS/transform setup；
- GMRES30：计时为 gmres core solve，matrix/RHS setup excluded；
- 该图用于展示当前实现下的整体基准趋势，而非严格 kernel-to-kernel benchmark。

## 5. Exp07 Design

统一设置：

- Dirichlet 五点离散；
- `n = 129`, `N = 127`, `h = 1/128`；
- 频域分母 `d_{p,q} = lambda_{p,q} + sigma`；
- 热图使用 `d_{p,q}` 本身作为颜色值，true Helmholtz 情形用黑色等值线标出 `d=0`。

三种情形：

| case | sigma | purpose |
|---|---:|---|
| modified_sigma_plus_100 | 100.000000 | modified Helmholtz，全正分母 |
| true_away_from_resonance | -64.140291 | true Helmholtz，远离共振但不定 |
| true_near_resonance_23_32 | -128.266807 | true Helmholtz，靠近 `(2,3)/(3,2)` 简并模态对 |

论文正文明确避免将 `spectral_indicator` 称为矩阵条件数；它只作为频域分母指标或小分母风险指标。

## 6. Sanity Checks

`exp07_spectral_denominator_summary.csv` 已生成且非空。关键检查结果：

| case | d_min_abs | nearest mode | num_positive | num_negative | num_near_zero | spectral_indicator |
|---|---:|---|---:|---:|---:|---:|
| modified_sigma_plus_100 | 119.738218 | `(1,1)` | 16129 | 0 | 0 | 1.095325e+03 |
| true_away_from_resonance | 14.800691 | `(1,2)` | 16126 | 3 | 0 | 8.850136e+03 |
| true_near_resonance_23_32 | 0.010000 | `(2,3)` | 16121 | 8 | 2 | 1.309240e+07 |

验收点：

- modified case: `num_negative = 0`。
- true-away case: `num_positive > 0` 且 `num_negative > 0`，并且 `d_min_abs = 14.800691`，不接近 0。
- near-resonance case: 最近模态为 `(2,3)`，`d_min_abs = 1.0e-2`，符合 `delta = 1e-2`。

## 7. Validation

- Root command equivalent to `PYTHONPATH=code/python python code/python/experiments/exp07_spectral_denominator_maps.py`: PASS。
- Root command equivalent to `PYTHONPATH=code/python python code/python/experiments/exp06_accuracy_cost.py`: PASS。
- `python -m pytest -q`: PASS, `103 passed, 2 skipped, 2 warnings`。
- `xelatex -> bibtex -> xelatex -> xelatex`: PASS。
- PDF output: `main.pdf (59 pages)`。
- LaTeX log scan: no fatal error, no undefined citation/reference, no multiply-defined label, no overfull hbox.
- Non-blocking warnings: existing SimSun bold-italic font substitution and existing hyperref PDF-string warnings.
- `git diff --check`: PASS; only CRLF normalization warnings were reported.

## 8. Final Status

READY TO USE
