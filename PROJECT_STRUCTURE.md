# Project Structure

This repository contains a numerical PDE thesis project for FFT-based direct solvers and unpreconditioned GMRES baselines for
\[
(-\Delta+\sigma)u=f
\]
on rectangular domains.

## Stable Runtime Layout

The following paths are treated as stable because Python scripts, tests, LaTeX figures, and reports refer to them directly.
Do not move them without updating imports, figure references, tests, and report links together.

| Path | Purpose | Move policy |
|---|---|---|
| `code/python/` | Python solvers, tests, and experiment scripts | Do not move solver modules or experiment scripts in cleanup-only patches |
| `code/python/experiments/` | Experiment entrypoints, `results/`, and generated `figures/` | Keep current relative paths; scripts assume these locations |
| `code/lean4_formalization/` | Lean 4 auxiliary algebra formalization | Keep `.lake/` ignored; do not move theorem files without rerunning `lake build` |
| `thesis/` | LaTeX thesis source, figures, bibliography, and generated PDF | Keep `thesis/figures/` stable because chapter references use these filenames |
| `reports/ch6_experiment_report/` | Standalone Chapter 6 HTML/PPT explanation materials | Generated report package; rebuild through scripts in `reports/` |

## Documentation Layout

| Path | Purpose |
|---|---|
| `README.md` | Short delivery overview and validation commands |
| `PROJECT_STRUCTURE.md` | This project map and file-movement safety policy |
| `reports/REPORT_INDEX.md` | Index of audit, patch, and review reports |
| `PROJECT_REVIEW_AND_STRUCTURE_REPORT.md` | Latest project review and structure-cleanup record |
| `codex_prompts/` | Prompt drafts and handoff text for future Codex/GPT review |
| `notes/` | Long-form reading notes and derivation notes |

## Report Organization Policy

Many historical audit and patch reports are currently kept at the repository root to preserve links from prior review rounds.
They are indexed in `reports/REPORT_INDEX.md` instead of being moved in this pass.

For new work:

- Put project-wide audit reports under `reports/audit/`.
- Put generated presentation/HTML outputs under `reports/ch6_experiment_report/`.
- Keep root-level files for final delivery reports only when they are referenced by README or external review prompts.

## Runtime Safety Rules

- Do not move `code/python/experiments/results/` or `code/python/experiments/figures/` without rerunning all affected experiment scripts and updating thesis figure references.
- Do not move `thesis/figures/` assets without recompiling LaTeX and checking all figure references.
- Do not rename solver modules such as `helmholtz_solver.py`, `gmres_solver.py`, or `cyclic_reduction.py` without updating imports and tests.
- Treat `code/python/_archive/` as historical prototype code, not final experiment evidence.

## Validation Commands

From the repository root:

```powershell
python -m pytest -q code/python/tests
cd code/lean4_formalization
lake build
cd ../../thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Expected current status:

- Python tests: `105 passed, 2 skipped`
- Lean: build succeeds; existing unused simp warnings may remain
- LaTeX: `main.pdf` generated, currently 71 pages
