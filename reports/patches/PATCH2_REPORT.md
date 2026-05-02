# Patch-2 Report: exp06 Accuracy-Cost Comparison

## 1. Summary

Added `exp06_accuracy_cost.py` to generate an accuracy-cost comparison for one smooth homogeneous Dirichlet manufactured solution:

`u(x,y)=x^2(1-x)^2 y^2(1-y)^2`.

The comparison covers `sigma=0` and `sigma=10`:

- FA, CR, FACR-like: five-point second-order discretization.
- FFT9: Dirichlet compact fourth-order discretization.
- GMRES30: unpreconditioned restarted GMRES(30) applied to the five-point system.

This is a current-implementation baseline and does not represent all preconditioned Krylov methods.

## 2. Generated Files

CSV:

- `code/python/experiments/results/exp06_accuracy_cost.csv`

Figures:

- `code/python/experiments/figures/exp06_accuracy_cost_error_time.png`
- `code/python/experiments/figures/exp06_accuracy_cost_error_time.pdf`
- `code/python/experiments/figures/exp06_time_scaling.png`
- `code/python/experiments/figures/exp06_time_scaling.pdf`
- `thesis/figures/exp06_accuracy_cost_error_time.png`
- `thesis/figures/exp06_accuracy_cost_error_time.pdf`
- `thesis/figures/exp06_time_scaling.png`
- `thesis/figures/exp06_time_scaling.pdf`

Thesis integration:

- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 3. Experiment Settings

- `sigma_list = [0, 10]`
- `n_list_direct = [33, 65, 129, 257, 513]`
- `n_list_gmres = [33, 65, 129]`
- `warmup = 2`
- `repeat = 7`
- Time statistic: median and IQR
- Direct timing scope: `solver_call_including_internal_rhs_and_transform_setup`
- GMRES timing scope: `gmres_core_solve_only_matrix_and_rhs_setup_excluded`

GMRES settings:

- `restart = 30`
- `gmres_tol_type = absolute_residual`
- `gmres_tol = 1e-10`
- `maxiter = 2000`

## 4. GMRES Results

| sigma | n | iterations | converged | final abs residual | final rel residual |
|---:|---:|---:|---|---:|---:|
| 0 | 33 | 149 | True | `8.068e-11` | `4.687e-11` |
| 0 | 65 | 784 | True | `9.797e-11` | `2.761e-11` |
| 0 | 129 | 2000 | False | `1.583e-07` | `2.198e-08` |
| 10 | 33 | 130 | True | `9.387e-11` | `4.445e-11` |
| 10 | 65 | 613 | True | `9.739e-11` | `2.260e-11` |
| 10 | 129 | 1992 | True | `9.968e-11` | `1.145e-11` |

The non-converged GMRES30 point is marked with a red cross in the error-time figure.

## 5. Thesis Wording

The new subsection `精度--成本对比与复杂度实证` states:

- FA/CR/FACR-like solve five-point second-order systems.
- FFT9 is a Dirichlet compact fourth-order method.
- GMRES30 solves the five-point system, so its error reflects both discretization error and algebraic residual.
- The experiment is a current implementation baseline with unpreconditioned GMRES30.
- Direct-method and GMRES30 timing scopes are not a strict kernel-to-kernel benchmark; the CSV `timing_scope` field records the exact timing scope, and the thesis text now explicitly states this limitation.
- The time-scaling figure uses an `O(N^2 log N)` reference and does not claim FACR-like achieves `O(N^2 log log N)`.

Follow-up visual cleanup:

- The non-converged GMRES30 marker in `exp06_accuracy_cost_error_time` was changed from a red cross to a black cross to avoid visual confusion with the red FFT9 curve.

## 6. Validation

Commands run:

```powershell
cd code/python
python -m experiments.exp06_accuracy_cost
python -m pytest -q

cd ../../thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Results:

- exp06 script: PASS
- pytest: `103 passed, 2 skipped, 2 warnings`
- XeLaTeX/BibTeX: PASS
- Final PDF: `main.pdf`, 58 pages
- Final log scan: no fatal error, no undefined citation/reference, no multiply-defined label
- Follow-up XeLaTeX after marker/caption cleanup: PASS, `main.pdf`, 58 pages

Non-blocking warnings:

- Existing SimSun bold-italic font substitution warning.
- Existing hyperref bookmark warnings for math tokens.
- MiKTeX update notice.

## 7. Final Status

READY TO USE
