# Experiment 6 Rollback Report

## 1. Goal

This round rolls the Experiment 6 main summary figure back from target mode `(4,4)` to the original degenerate modal pair `(2,3)/(3,2)`, while keeping the clearer 2x2 summary layout introduced during the visual redesign.

No new core experiment was added. Chapter 6 still keeps the six-core-experiment structure, and this remains part of Experiment 6: true Helmholtz near-resonance.

## 2. Why rollback to `(2,3)/(3,2)`

The original near-resonance case used the square-domain degeneracy

`lambda_(2,3)^h = lambda_(3,2)^h`.

Restoring this target pair keeps the experiment aligned with the earlier near-resonance design and table values. The figure uses `|u_hat_(2,3)|` as a representative modal coefficient for the degenerate pair, and the caption now states this explicitly.

## 3. Current version vs older versions

The current version differs from earlier visual variants in three ways:

- It does not use the denominator heatmap as the main figure.
- It does not use normalized fields in panels (c)--(d).
- It does not use a far-from-spectrum auto-selected control case.

Instead, panels (c)--(d) compare raw solution fields on the same near-resonance branch:

- `delta = 1e-1`: `kappa^2 = lambda_(2,3)^h + 1e-1`.
- `delta = 1e-4`: `kappa^2 = lambda_(2,3)^h + 1e-4`.

This makes the bottom-row comparison a fair detuning comparison for the same degenerate modal pair.

## 4. Why raw fields instead of normalized fields

Normalized fields made it easy to see spatial structure, but they hid the main near-resonance effect: amplitude growth caused by small denominators.

The current raw fields use independent colorbars. Their spatial structures may look similar because the same modal pair is being amplified. The important difference is the amplitude, which is shown quantitatively by panel (a) and by each panel title's `max |u|` and `min |d_pq|`.

## 5. Numerical comparison

Grid: `n = 65`, `N = 63`, `h = 1/64`.

Target modal pair: `(2,3)/(3,2)`.

Discrete eigenvalue:

`lambda_(2,3)^h = 128.11274946987737`.

| delta | min `|d_pq|` | `max |u|` | RMS `||u||_2` | `|u_hat_(2,3)|` | GMRES iterations | GMRES success | final abs residual |
|---:|---:|---:|---:|---:|---:|---|---:|
| `1e-1` | `1.000e-1` | `5.618e-1` | `2.223e-1` | `1.249e1` | `1001` | `False` | `1.445e0` |
| `1e-4` | `1.000e-4` | `5.639e2` | `2.223e2` | `1.249e4` | `1001` | `False` | `1.447e0` |

The representative modal coefficient satisfies the expected small-denominator pattern:

`|u_hat_(2,3)| * |delta| ≈ 1.249`.

## 6. Files updated

- `code/python/experiments/exp05_true_helmholtz_resonance.py`
- `code/python/experiments/results/exp05_resonance.csv`
- `code/python/experiments/figures/exp05_resonance.png`
- `code/python/experiments/figures/exp05_near_resonance_summary.png`
- `thesis/figures/exp05_resonance.png`
- `thesis/figures/exp05_near_resonance_summary.png`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 7. Validation

- Figure/CSV regeneration: `python -m experiments.exp05_true_helmholtz_resonance` - PASS.
- Pytest: `103 passed, 2 skipped, 2 warnings` - PASS.
- XeLaTeX/BibTeX chain: PASS.
- PDF output: `main.pdf (56 pages)`.
- Log scan: no fatal error, no undefined citation/reference, no multiply-defined label. Remaining warnings are non-blocking font substitution and hyperref bookmark-token warnings.
- Content scan: no `4,4`, `normalized`, `归一化`, `far-from-spectrum`, `实验七`, or `Experiment 7` remains in `thesis/chapters/6_experiments.tex`.
- `git diff --check`: no whitespace errors; Git only reported line-ending normalization warnings.

## 8. Final status

READY TO USE
