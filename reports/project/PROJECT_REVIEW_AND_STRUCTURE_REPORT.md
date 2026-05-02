# Project Review and Structure Report

## 1. Review Goal

This pass reviews the repository structure and adds navigation files so the project is easier to understand without moving runtime-sensitive files.

The main safety rule was:

> Do not move code, experiment CSVs, generated figures, or thesis assets in a way that could make scripts, imports, or LaTeX references fail.

## 2. Actions Taken

- Added `PROJECT_STRUCTURE.md` as the top-level project map and file-movement safety policy.
- Added `reports/REPORT_INDEX.md` to index root-level audit and patch reports without moving historical files.
- Added `reports/audit/README.md` as the intended home for future audit reports.
- Added `code/python/experiments/README.md` to document experiment scripts, CSV outputs, figures, and runtime safety.
- Updated `README.md` with a project navigation section and corrected the current PDF page count to 71 pages.

## 3. Files Deliberately Not Moved

No runtime files were moved.

In particular, the following locations remain stable:

- `code/python/`
- `code/python/experiments/`
- `code/python/experiments/results/`
- `code/python/experiments/figures/`
- `thesis/chapters/`
- `thesis/figures/`
- `code/lean4_formalization/`

Historical root-level reports were not moved because prior Codex/GPT prompts and review summaries may reference their current paths. They are now discoverable via `reports/REPORT_INDEX.md`.

## 4. Current Repository Organization

| Area | Current role |
|---|---|
| `code/python/` | Solver implementations, experiments, and Python tests |
| `code/lean4_formalization/` | Lean 4 auxiliary algebra formalization |
| `thesis/` | Thesis source, figures, bibliography, and compiled PDF |
| `reports/` | Chapter 6 HTML/PPT materials and report index |
| `codex_prompts/` | Prompt and patch text for future model review |
| `notes/` | Reading notes and derivation notes |

## 5. Validation Plan

The safe validation commands after this structure pass are:

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

## 6. Push Scope

The push should include:

- existing thesis/code/report changes already present in the working tree;
- the new structure/navigation files from this pass;
- no deletion or migration of runtime-sensitive files.

Because the repository is on branch `revise-sigma-fft9-krylov`, the safest publishing route is to commit on this branch and push it to `origin/revise-sigma-fft9-krylov`.
