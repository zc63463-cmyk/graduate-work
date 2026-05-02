# Files To Patch

This file is an optional patch plan only. The audit did not modify thesis text, code, tests, experiments, CSV files, or figures.

| Priority | File | Location | Why | Minimal patch | Need rerun |
|---|---|---|---|---|---|
| A | `thesis/chapters/2_math_preliminary.tex` | `:315`, `:322`, `:507` | Chapter 2 FFT9 derivation/effective eigenvalue signs do not fully match the final Chapter 3 scheme and Python implementation. | Rewrite the affected FFT9 derivation so `L_h` approximates `Delta`, the scheme is `(-L_h + sigma R_h)u = R_h f`, the raw denominator is `-lambda_L + sigma lambda_R`, and the effective expansion uses `-lambda_L/lambda_R`. | no |
| B | `code/python/experiments/exp04_modified_vs_true.py` | `:331-349` | Exp04 GMRES plot includes capped/failed rows as ordinary iteration-count points. | Filter `gmres_success == True` for the main curve, or add explicit failure markers/caption text for capped runs. | yes, regenerate Exp04 figure |
| B | `code/python/experiments/exp05_true_helmholtz_resonance.py` | `:176-190` | Baseline rows hard-code `gmres_success=True`. | Store `info['success']` when GMRES returns, and `False` in the exception path. | yes, regenerate Exp05 CSV; regenerate figure if script rewrites it |
| C | `code/python/tests/test_07_krylov_baselines.py` | `:21-23` | Test header still claims GMRES iteration-ordering validation, while the implemented test only checks convergence/applicability sanity. | Update docstring/header to match the implemented scope. | no |

## Files Checked With No Patch Recommended

- `code/python/helmholtz_solver.py`: final five-point, FFT9 raw denominator, FFT9 boundary correction, and Dirichlet-only FFT9 support match the intended scheme.
- `code/python/gmres_solver.py`: GMRES stopping criterion is absolute residual tolerance, matching Chapter 6.
- `code/python/cyclic_reduction.py`: current comments distinguish FACR-like `O(N^2 log N)` implementation from classical FACR `O(N^2 log log N)` background.
- `code/python/experiments/exp02_nonhom_bc.py`: manufactured nonhomogeneous Dirichlet problem matches the paper and `utils.py`.
- `code/python/experiments/exp04_modified_vs_true.py`: Gaussian RHS choice is correct; only the failure-aware plotting presentation needs attention.
- `README.md` and `FINAL_POLISH_REPORT.md`: page count, experiment count, CSV list, and core PNG list are consistent.
