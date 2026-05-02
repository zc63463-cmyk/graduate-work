# EXP04 Experiment Data Audit

## 1. Scope

本报告只审计 `exp04_modified_vs_true` 的数据口径，不新增算法、不重跑实验、不修改论文正文或图像。

审计对象：

- CSV: `code/python/experiments/results/exp04_modified_vs_true.csv`
- Script: `code/python/experiments/exp04_modified_vs_true.py`
- Thesis figures:
  - `thesis/figures/exp04_min_denom_vs_sigma.png`
  - `thesis/figures/exp04_spectral_indicator_vs_sigma.png`
  - `thesis/figures/exp04_gmres_iters_vs_sigma.png`
- Thesis text: `thesis/chapters/6_experiments.tex`

## 2. Direct Answers

| Item | Audit result |
|---|---|
| Figure 8 n | `n=129`. The script plots the largest available `n` for each equation type. |
| Figure 9 n | `n=129`. The GMRES iteration plot also selects the largest available `n`. |
| Table 9 n | `n=33`. The caption explicitly says `n=33`, and the numbers match CSV rows with `n=33`. |
| GMRES tolerance | `tol=1e-10`, absolute residual tolerance. |
| GMRES restart | `restart=30`. |
| GMRES max_iter | Configured as `max_iter=500`. CSV capped rows report `gmres_iters=501`, due to iteration counting/callback convention. This is not `1000/1001`. |
| RHS | Gaussian non-eigenfunction RHS. |
| Figure/table n consistency | Not consistent: Figures 8-9 use `n=129`, Table 9 uses `n=33`. This is acceptable only if captions/text explicitly state the different representative grids. |
| Capped marker consistency | Table 9 status is consistent with `n=33` CSV. Figure 9 caption says capped marker, but the current script uses a plain `o-` plot and does not visibly mark capped/not-converged points. |

## 3. Artifact Inventory

| Thesis object | Data source CSV | Generation/source | n used | sigma values | GMRES settings | Notes |
|---|---|---|---:|---|---|---|
| Figure 8: `fig:exp04_spectral` | `exp04_modified_vs_true.csv` | `exp04_modified_vs_true.py`, plots `exp04_min_denom_vs_sigma.png` and `exp04_spectral_indicator_vs_sigma.png` | 129 | `-50,-10,-5,-1,1,10,100,1000` | Spectral indicators only; same CSV includes GMRES fields | Script selects `sub['n'].max()`. Caption currently does not state `n=129`. |
| Figure 9: `fig:exp04_gmres` | `exp04_modified_vs_true.csv` | `exp04_modified_vs_true.py`, plot `exp04_gmres_iters_vs_sigma.png` | 129 | `-50,-10,-5,-1,1,10,100,1000` | `tol=1e-10`, `restart=30`, `max_iter=500` | Script selects `n=129`. Many rows are capped/not converged, but plot code does not use a special capped marker. |
| Table 9: `tab:exp04_mod_true` | `exp04_modified_vs_true.csv` | Manual LaTeX table in `6_experiments.tex` | 33 | `-50,-10,-5,-1,1,10,100,1000` | `tol=1e-10`, `restart=30`, `max_iter=500` | Table values match CSV `n=33`. Only `sigma=-50` is marked near-resonance / not converged. |

## 4. CSV Audit

Only one exp04 CSV was found:

`code/python/experiments/results/exp04_modified_vs_true.csv`

CSV properties:

- Rows: 24
- Grid sizes: `n = 33, 65, 129`
- Modified Helmholtz sigma values: `1, 10, 100, 1000`
- True Helmholtz sigma values: `-1, -5, -10, -50`
- RHS type: `gaussian`
- Key columns: `sigma`, `equation_type`, `rhs_type`, `n`, `min_denominator`, `spectral_condition_indicator`, `gmres_iters`, `gmres_success`, `final_abs_residual`, `final_relative_residual`, `status`

Representative `n=33` rows, used by Table 9:

| sigma | type | min denominator | spectral indicator | GMRES iters | success | status |
|---:|---|---:|---:|---:|---|---|
| 1 | modified | 20.7234 | 394.399 | 256 | True | ok |
| 10 | modified | 29.7234 | 275.281 | 203 | True | ok |
| 100 | modified | 119.7234 | 69.0949 | 91 | True | ok |
| 1000 | modified | 1019.7234 | 8.99487 | 27 | True | ok |
| -1 | true | 18.7234 | 436.421 | 270 | True | ok |
| -5 | true | 14.7234 | 554.716 | 319 | True | ok |
| -10 | true | 9.72336 | 839.450 | 428 | True | ok |
| -50 | true | 0.786574 | 10326.1 | 501 | False | near_resonance |

Representative `n=129` rows, used by Figures 8-9:

| sigma | type | min denominator | spectral indicator | GMRES iters | success | status |
|---:|---|---:|---:|---:|---|---|
| 1 | modified | 20.7382 | 6319.41 | 501 | False | failed |
| 10 | modified | 29.7382 | 4407.20 | 501 | False | failed |
| 100 | modified | 119.7382 | 1095.32 | 501 | False | failed |
| 1000 | modified | 1019.7382 | 129.496 | 130 | True | ok |
| -1 | true | 18.7382 | 6993.80 | 501 | False | failed |
| -5 | true | 14.7382 | 8891.66 | 501 | False | failed |
| -10 | true | 9.73822 | 13456.5 | 501 | False | failed |
| -50 | true | 0.660400 | 198368 | 501 | False | near_resonance |

## 5. Script Audit

Relevant script facts from `code/python/experiments/exp04_modified_vs_true.py`:

- `n_list = [33, 65, 129]`
- `sigma_modified = [1, 10, 100, 1000]`
- `sigma_true_safe = [-1, -5, -10]`
- `sigma_true_near = [-50]`
- GMRES call uses `restart=30`, `max_iter=500`, `tol=1e-10`
- Figure 8 spectral plots use:
  - `n_max = sub['n'].max()`
  - therefore `n=129`
- Figure 9 GMRES plot uses:
  - `df_plot = df[df['gmres_iters'] > 0]`
  - `n_max = sub['n'].max()`
  - therefore `n=129`
- Figure 9 plot command is a plain `plt.plot(..., 'o-')`; there is no special marker for failed/capped points.

## 6. Figure and Caption Audit

### Figure 8

Current thesis caption:

`Modified/true Helmholtz 的最小频域分母与谱条件数指标`

Audit result:

- Data source is valid.
- Figure uses `n=129`.
- Caption does not say `n=129`.
- Since Table 9 uses `n=33`, the caption should explicitly state that Figure 8 is the largest-grid plot.

### Figure 9

Current thesis caption:

`Gaussian RHS 下无预处理 GMRES 迭代次数对比；达到迭代上限的点表示未在给定容差内收敛的 capped marker`

Audit result:

- Data source is valid.
- Figure uses `n=129`.
- Caption does not say `n=129`.
- Caption says capped marker, but the current script does not visibly distinguish capped/not-converged points.
- At `n=129`, most rows are not converged within `max_iter=500` and appear as `gmres_iters=501`.

### Table 9

Current thesis caption:

`Modified/true Helmholtz 谱指标与 GMRES 迭代次数（n=33）`

Audit result:

- Table values match the CSV rows with `n=33`.
- The `>500` entry for `sigma=-50` is consistent with CSV `gmres_iters=501`, `gmres_success=False`, `status=near_resonance`.
- Table status is internally consistent for `n=33`.

## 7. Consistency Issues

| Issue | Severity | Evidence | Risk | Minimal fix |
|---|---|---|---|---|
| Figure 8/9 use `n=129`, Table 9 uses `n=33`, but only Table 9 states its `n` explicitly. | B | Script selects `n_max=129`; Table caption says `n=33`. | Reader may compare figure and table statuses as if they use the same grid. | Use option B: keep current artifacts, but state in captions/text that Figures 8-9 use `n=129` and Table 9 uses `n=33`. |
| Figure 9 caption mentions capped marker, but plot script uses plain markers. | B | Script uses `plt.plot(..., 'o-')` with no failed/capped marker branch. | Caption overclaims visual encoding. | Either update plot to mark capped rows, or downgrade caption to say capped rows appear at `501` iterations. |
| GMRES configured `max_iter=500`, while CSV reports capped `gmres_iters=501`. | C | CSV capped rows have `501`; script uses `max_iter=500`. | Could be mistaken for `max_iter=1000/1001` or an off-by-one algorithm setting. | State that `501` is the iteration counter convention for capped runs under configured `max_iter=500`. |

## 8. Recommended Resolution

Recommended choice: **Option B**.

Keep Figures 8-9 at `n=129` and Table 9 at `n=33`, because this preserves the current artifacts and avoids regenerating exp04. The paper should explicitly describe the mixed grid usage:

- Figures 8-9: largest-grid visual behavior at `n=129`.
- Table 9: compact representative numeric table at `n=33`.

Suggested caption/text patch for a later edit:

1. Figure 8 caption:
   `Modified/true Helmholtz 的最小频域分母与谱条件数指标（Gaussian RHS 实验数据，图中取最大网格 n=129）。`

2. Figure 9 caption:
   `Gaussian RHS 下无预处理 GMRES 行为（n=129, tol=1e-10, restart=30, max_iter=500）；迭代数为 501 的点表示达到迭代上限后仍未达到绝对残差容差。`

3. Table 9 caption:
   `Modified/true Helmholtz 谱指标与 GMRES 迭代次数的代表性数值（n=33；GMRES 设置同图 9）。`

4. Add one explanatory sentence near Table 9:
   `图 8--图 9 使用最大网格 n=129 展示较大规模下的谱指标和 GMRES 压力；表 9 使用 n=33 给出同一参数组的紧凑代表性数值，因此图表中的收敛状态不应脱离网格规模直接逐项比较。`

If a future patch chooses **Option A**, then `exp04_modified_vs_true.py` should be changed to plot `n=33` for Figures 8-9 and the PNGs regenerated. That was not done in this audit.

## 9. Final Audit Status

`EXP04 DATA PROVENANCE AUDIT COMPLETE`

Main finding: exp04 data are traceable to one CSV and one script, but Figure 8/9 and Table 9 currently use different grid sizes. The safest minimal documentation fix is Option B: explicitly state `n=129` for the figures and `n=33` for the table, and clarify that capped `501` rows correspond to `max_iter=500`.
