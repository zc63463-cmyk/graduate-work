# PATCH4 Report

## 1. Scope

本轮升级 `exp04_modified_vs_true.py`，新增 Dirichlet 五点矩阵 condition check 与 GMRES residual history，不重构求解器，不覆盖历史 `exp04_modified_vs_true.csv`。

核心目标：

- 验证 Dirichlet 五点离散下 `max|d|/min|d|` 与真实稠密矩阵 `cond_2(A)` 一致。
- 为 modified/true Helmholtz 的代表性参数补充 GMRES(30) 残差历史，展示收敛、未收敛和 near-resonance 的过程差异。
- 在论文第 6 章中说明该等价性依赖当前 Dirichlet 五点正交 DST-I 对角化结构，不外推到 Neumann ghost-point 原始非对称矩阵、mixed 边界或一般 non-normal 系统。

## 2. Modified Files

- `code/python/experiments/exp04_modified_vs_true.py`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`
- `PATCH4_REPORT.md`

## 3. New CSV And Figures

CSV:

- `code/python/experiments/results/exp04_condition_check.csv`
- `code/python/experiments/results/exp04_gmres_history.csv`

Figures:

- `code/python/experiments/figures/exp04_condition_check.png`
- `code/python/experiments/figures/exp04_condition_check.pdf`
- `thesis/figures/exp04_condition_check.png`
- `thesis/figures/exp04_condition_check.pdf`
- `code/python/experiments/figures/exp04_gmres_history.png`
- `code/python/experiments/figures/exp04_gmres_history.pdf`
- `thesis/figures/exp04_gmres_history.png`
- `thesis/figures/exp04_gmres_history.pdf`

## 4. Condition Check

Settings:

- `n_list_dense = [17, 33]`
- `sigma_list = [-50, -10, -5, -1, 1, 10, 100, 1000]`
- Matrix: Dirichlet five-point sparse matrix converted to dense.
- Dense condition number: `np.linalg.cond(A_dense, 2)`.
- Spectral denominator indicator: `max(abs(lambda_pq + sigma)) / min(abs(lambda_pq + sigma))`.

Result:

- Maximum `relative_difference`: `2.9118238778189183e-13`.
- Threshold requested: `< 1e-8`.
- Status: PASS.

Interpretation:

For the current Dirichlet five-point system, the matrix is orthogonally diagonalized by DST-I, so the spectral denominator indicator equals `cond_2(A)` for nonsingular cases. The thesis text now states this equivalence and also states that it should not be directly extrapolated to non-orthogonal or non-normal systems.

## 5. GMRES Residual History

Settings:

- `n = 129`
- `restart = 30`
- `tol = 1e-10`
- `tol_type = absolute_residual`
- `max_iter = 500`
- RHS: Gaussian, same as exp04 main experiment.
- Solver: unpreconditioned GMRES(30).

Final status table:

| sigma | iterations | final_abs_residual | final_rel_residual | success | status |
|---:|---:|---:|---:|---|---|
| `+10` | 501 | `9.2189147326005e-03` | `5.139902531074e-04` | False | failed |
| `+1000` | 130 | `8.950862949220597e-11` | `4.990453265101621e-12` | True | ok |
| `-10` | 501 | `7.699054483454947e-01` | `4.2925214924108e-02` | False | failed |
| `-50` | 501 | `8.572596107135034e+00` | `4.779554829065356e-01` | False | near_resonance |

The value `501` follows the existing exp04 capped-count convention for a run that reaches `max_iter=500` without satisfying the absolute residual tolerance.

## 6. Thesis Integration

`thesis/chapters/6_experiments.tex` was updated to include:

- A condition-check paragraph and `exp04_condition_check` figure.
- A residual-history paragraph and `exp04_gmres_history` figure.
- A statement that GMRES stopping uses absolute residual tolerance; relative residual is only for normalized comparison.
- A scope limitation: the condition-number equivalence is only for the current Dirichlet five-point orthogonal diagonalization.

Original Figure 8/9/Table 9 were retained; Patch4 only adds supporting evidence.

## 7. Validation

- `PYTHONPATH=code/python python code/python/experiments/exp04_modified_vs_true.py`: PASS.
- `python -m pytest -q`: PASS, `103 passed, 2 skipped, 2 warnings`.
- `xelatex -> bibtex -> xelatex -> xelatex`: PASS.
- PDF output: `main.pdf (60 pages)`.
- LaTeX log scan: no fatal error, no undefined citation/reference, no multiply-defined label, no overfull hbox.
- Non-blocking warnings: existing SimSun bold-italic font substitution and existing hyperref PDF-string warnings.
- `git diff --check`: PASS; only CRLF normalization warnings were reported.

## 8. Known Limitations

- GMRES history is for unpreconditioned GMRES(30), not a preconditioned Krylov method study.
- GMRES stopping uses absolute residual tolerance `1e-10`; relative residual is plotted only for comparison.
- The condition-number equivalence applies only to the current Dirichlet five-point matrix with orthogonal DST-I diagonalization.
- The equivalence should not be extrapolated to Neumann ghost-point original nonsymmetric matrices, mixed boundary systems, non-orthogonal diagonalizations, or general non-normal operators.

## 9. Final Status

READY TO USE
