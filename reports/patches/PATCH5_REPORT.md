# PATCH5 Report

## 1. Summary

Patch5 upgraded `exp05_true_helmholtz_resonance.py` while preserving the original near-resonance experiment outputs. The patch adds a multimode near-resonance scan, DST-I modal projection checks, dominant-mode projection visualization, and near-resonance GMRES residual histories.

The Chapter 6 experiment structure remains unchanged: this is still part of Experiment 6, not a new core experiment.

## 2. Modified Files

- `code/python/experiments/exp05_true_helmholtz_resonance.py`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 3. New CSV Outputs

- `code/python/experiments/results/exp05_multimode_resonance.csv`
- `code/python/experiments/results/exp05_resonance_gmres_history.csv`

Existing outputs were preserved and regenerated:

- `code/python/experiments/results/exp05_resonance.csv`

## 4. New Figures

Code-side figures:

- `code/python/experiments/figures/exp05_multimode_resonance_summary.png`
- `code/python/experiments/figures/exp05_multimode_resonance_summary.pdf`
- `code/python/experiments/figures/exp05_dominant_mode_projection.png`
- `code/python/experiments/figures/exp05_dominant_mode_projection.pdf`
- `code/python/experiments/figures/exp05_resonance_gmres_history.png`
- `code/python/experiments/figures/exp05_resonance_gmres_history.pdf`

Thesis-side copies:

- `thesis/figures/exp05_multimode_resonance_summary.png`
- `thesis/figures/exp05_multimode_resonance_summary.pdf`
- `thesis/figures/exp05_dominant_mode_projection.png`
- `thesis/figures/exp05_dominant_mode_projection.pdf`
- `thesis/figures/exp05_resonance_gmres_history.png`
- `thesis/figures/exp05_resonance_gmres_history.pdf`

Existing Experiment 6 figures were preserved and regenerated:

- `code/python/experiments/figures/exp05_resonance.png`
- `code/python/experiments/figures/exp05_near_resonance_summary.png`
- `thesis/figures/exp05_resonance.png`
- `thesis/figures/exp05_near_resonance_summary.png`

## 5. Multimode Resonance Sanity Check

The multimode scan uses the off-center Gaussian RHS

`f(x,y)=exp(-40*((x-0.37)^2+(y-0.41)^2))`

normalized to unit RMS norm, and orthonormal DST-I modal coefficients. The checked target sets are `(1,1)`, `(2,3)/(3,2)`, and `(3,3)`.

Global sanity checks:

- CSV rows: `12`
- Maximum modal amplitude relative error: `3.319662256415239e-11`
- `(2,3)/(3,2)` target lambda spread: `0.0` for all tested detuning values

At `delta=1e-4`:

| mode_label | target_energy_fraction | shape_correlation | dominant_relative_error | gmres_iterations | gmres_status |
|---|---:|---:|---:|---:|---|
| `mode_11` | `0.9999999999916498` | `0.9999999999958252` | `2.889671157473289e-06` | `1001` | `near_resonance` |
| `mode_23_32` | `0.9999999999168744` | `0.9999999999584371` | `9.117325534474181e-06` | `1001` | `near_resonance` |
| `mode_33` | `0.9999999959215596` | `0.9999999979607794` | `6.386266853331281e-05` | `1001` | `near_resonance` |

These values confirm the DST formula

`u_hat[p,q] = f_hat[p,q] / (lambda[p,q] - kappa^2)`

and the expected `1/|delta|` target-mode amplification.

## 6. GMRES History Summary

The residual history uses target modes `(2,3)/(3,2)`, `n=65`, restart `30`, absolute residual tolerance `1e-10`, and the existing Experiment 6 capped-count convention.

| delta | final iteration | success | status | final_abs_residual | final_rel_residual |
|---:|---:|---|---|---:|---:|
| `1e-1` | `1001` | `False` | `near_resonance` | `0.18295758741233145` | `0.18295758741233142` |
| `1e-3` | `1001` | `False` | `near_resonance` | `0.18328608743115676` | `0.18328608743115674` |
| `1e-4` | `1001` | `False` | `near_resonance` | `0.18328769529899563` | `0.1832876952989956` |

## 7. Validation

Script execution:

- `PYTHONPATH=code/python python code/python/experiments/exp05_true_helmholtz_resonance.py`: PASS
- `python -m experiments.exp05_true_helmholtz_resonance` from `code/python`: PASS

Python tests:

- `python -m pytest -q`: PASS
- Result: `103 passed, 2 skipped, 2 warnings`
- Warnings are the existing Neumann Poisson compatibility warnings in tests.

LaTeX:

- `xelatex -interaction=nonstopmode -halt-on-error main.tex`: PASS
- `bibtex main`: PASS
- Final PDF: `thesis/main.pdf`, `62 pages`
- Log scan found no fatal errors, undefined citations/references, or multiply-defined labels.
- Non-blocking warnings remain: SimSun bold italic substitution and hyperref PDF bookmark token warnings.

Whitespace:

- `git diff --check`: PASS
- Only line-ending notices were printed (`LF will be replaced by CRLF`), with no whitespace errors.

## 8. Known Limitations

- The Gaussian RHS is fixed and off-center; modal projection magnitudes depend on this RHS.
- GMRES results are for unpreconditioned GMRES(30) only.
- The near-resonance scan is based on the finite-dimensional Dirichlet five-point discrete spectrum and is not a direct statement about the continuous spectrum limit.
- In this Gaussian setup, target energy fraction and shape correlation are already close to `1` across the scanned detuning range; the clearest quantitative evidence is the target modal amplitude following `1/|delta|` and the residual-history stagnation.

## 9. Final Status

READY TO USE
