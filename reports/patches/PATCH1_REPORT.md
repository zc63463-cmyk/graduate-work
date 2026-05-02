# Patch-1 Report: exp04 Figure/Table Data Scope Fix

## 1. Summary

本轮只修复 `exp04_modified_vs_true` 的图表口径说明与 Figure 9 的 capped/not-converged 可视化标记。
未新增实验，未重构求解器，未改变 exp04 参数，也未重跑 GMRES 计算。

## 2. Files Modified

- `thesis/chapters/6_experiments.tex`
- `code/python/experiments/exp04_modified_vs_true.py`
- `code/python/experiments/figures/exp04_gmres_iters_vs_sigma.png`
- `thesis/figures/exp04_gmres_iters_vs_sigma.png`
- `thesis/main.pdf`

No change was made to:

- `code/python/experiments/results/exp04_modified_vs_true.csv`

## 3. Thesis Caption / Text Changes

Final scope:

- Figure 8 uses the largest grid: `n=129`.
- Figure 9 uses the largest grid: `n=129`.
- Table 9 uses representative small-grid values: `n=33`.

Updated captions:

- Figure 8 now states that the spectral indicators use Gaussian RHS experiment data and the largest grid `n=129`.
- Figure 9 now states `n=129`, `tol=1e-10`, `restart=30`, `max_iter=500`, and clarifies that iteration count `501` means the run reached the iteration cap without satisfying the absolute residual tolerance.
- Table 9 now states that it contains representative values at `n=33` and uses the same GMRES settings as Figure 9.

Added explanatory sentence near Table 9:

Figures 8--9 use `n=129` to show larger-grid spectral indicators and GMRES pressure, while Table 9 uses `n=33` for a compact representative numeric table; convergence statuses should not be compared entry-by-entry without accounting for grid size.

## 4. Figure 9 Plot Change

`code/python/experiments/exp04_modified_vs_true.py` now marks capped/not-converged rows specially:

- `gmres_success=False` or `status in ["failed", "near_resonance"]` is plotted with a red `x`.
- Converged points remain circular markers.
- A dashed reference line at `501` is labeled `501 = capped under max_iter=500`.

The existing CSV was loaded and the figure was regenerated from it. The script was not run in full, so experiment parameters and timing fields were not recomputed.

Generated / updated:

- `code/python/experiments/figures/exp04_gmres_iters_vs_sigma.png`
- `thesis/figures/exp04_gmres_iters_vs_sigma.png`

The exp04 script does not generate a PDF version of this figure, so no figure PDF was produced.

## 5. GMRES Settings

Exp04 GMRES settings remain:

- RHS: Gaussian non-eigenfunction RHS
- tolerance: `tol=1e-10`
- tolerance type: absolute residual tolerance
- restart: `30`
- max_iter: `500`

The CSV reports capped runs as `gmres_iters=501` due to the current iteration counter convention. This does not mean that `max_iter` is `1000` or `1001`.

## 6. Validation

Commands run:

```powershell
cd code/python
@'
import os
import shutil
import pandas as pd
from experiments.exp04_modified_vs_true import plot

csv_path = os.path.join('experiments', 'results', 'exp04_modified_vs_true.csv')
df = pd.read_csv(csv_path)
plot(df)
shutil.copy2(
    os.path.join('experiments', 'figures', 'exp04_gmres_iters_vs_sigma.png'),
    os.path.join('..', '..', 'thesis', 'figures', 'exp04_gmres_iters_vs_sigma.png'),
)
'@ | python -

cd ../../thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Validation result:

- Figure regeneration: PASS
- XeLaTeX/BibTeX compile: PASS
- PDF output: `main.pdf`, 57 pages
- Fatal LaTeX error: none
- Undefined citation/reference detected by final log scan: none
- Multiply-defined label detected by final log scan: none

Non-blocking warnings:

- Existing SimSun bold-italic font substitution warning.
- Existing hyperref bookmark warnings for math tokens.
- MiKTeX update notice.

## 7. Final Status

READY TO USE
