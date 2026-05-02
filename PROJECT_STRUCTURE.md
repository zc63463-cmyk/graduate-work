# Project Structure

This repository contains a numerical PDE thesis project for FFT-based direct solvers and unpreconditioned GMRES baselines for
\[
(-\Delta+\sigma)u=f
\]
on rectangular domains.

## Root Layout

The repository root is intentionally small:

| Path | Purpose |
|---|---|
| `README.md` | Delivery overview, main entrypoints, and validation commands |
| `PROJECT_STRUCTURE.md` | This file: file layout and movement safety policy |
| `code/` | Python solver code, experiment scripts, tests, and Lean 4 formalization |
| `thesis/` | LaTeX source, bibliography, figures, and compiled thesis PDF |
| `reports/` | Review reports, patch reports, generated HTML/PPT materials, and notes |
| `codex_prompts/` | Prompt drafts and handoff material for future Codex/GPT review |

## Stable Runtime Layout

The following paths are treated as stable because Python scripts, tests, LaTeX figures, and generated reports refer to them directly.
Do not move them without updating imports, figure references, tests, and report links together.

| Path | Purpose | Move policy |
|---|---|---|
| `code/python/` | Python solvers, tests, and experiment scripts | Do not move solver modules or experiment scripts in cleanup-only patches |
| `code/python/experiments/` | Experiment entrypoints, `results/`, and generated `figures/` | Keep current relative paths; scripts assume these locations |
| `code/lean4_formalization/` | Lean 4 auxiliary algebra formalization | Keep `.lake/` ignored; do not move theorem files without rerunning `lake build` |
| `thesis/` | LaTeX thesis source, figures, bibliography, and generated PDF | Keep `thesis/figures/` stable because chapter references use these filenames |
| `reports/ch6_experiment_report/` | Standalone Chapter 6 HTML/PPT explanation materials | Rebuild through scripts in `reports/` |

## Report Layout

Root-level reports have been moved under `reports/` so the repository homepage stays readable.

| Path | Purpose |
|---|---|
| `reports/REPORT_INDEX.md` | Human-readable index of review, audit, patch, project, defense, and experiment reports |
| `reports/project/` | Project onboarding, structure review, delivery polish, and high-level project reports |
| `reports/planning/` | Roadmaps, execution plans, optional patch plans, and future-work notes |
| `reports/patches/` | Patch reports, Q-series reports, and final audit patch records |
| `reports/experiments/` | Chapter 6 experiment review, figure-patch reports, asset indexes, and experiment evidence summaries |
| `reports/audit/` | Theory, formula, Lean, paper-code, and experiment audit reports |
| `reports/defense/` | Defense Q&A, PPT outlines, and GPT review question packs |
| `reports/literature/` | Literature and arXiv search notes |
| `reports/notes/` | Long-form reading notes and explanation notes |
| `reports/archive/` | Local historical drafts not used as active delivery material |

## Runtime Safety Rules

- Do not move `code/python/experiments/results/` or `code/python/experiments/figures/` without rerunning all affected experiment scripts and updating thesis figure references.
- Do not move `thesis/figures/` assets without recompiling LaTeX and checking all figure references.
- Do not rename solver modules such as `helmholtz_solver.py`, `gmres_solver.py`, or `cyclic_reduction.py` without updating imports and tests.
- Treat `code/python/_archive/` as historical prototype code, not final experiment evidence.
- Documentation-only organization patches should prefer moving Markdown reports under `reports/` and leave runtime paths unchanged.

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
