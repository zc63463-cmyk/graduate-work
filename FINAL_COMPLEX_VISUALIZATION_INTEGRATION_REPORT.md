# Final Complex Visualization Integration Report

生成时间：2026-04-29

## 1. Summary

本轮只整理“多模态 manufactured solution 可视化补充”，不新增核心实验、不恢复旧图、不改变第 6 章 6 个核心实验结构。

当前第 6 章仍只包含六个核心实验：

1. FFT 解与 sparse direct 解一致性；
2. 二阶与四阶收敛阶验证；
3. 非齐次 Dirichlet 边界验证；
4. Neumann 与 mixed 边界验证；
5. modified/true Helmholtz 谱指标与 Gaussian RHS 下 GMRES 行为；
6. true Helmholtz near-resonance 扫描。

`exp06_complex_manufactured_visualization.py` 虽沿用 `exp06` 文件名，但已在 README 和报告中明确为 supplementary visualization script，不计入核心实验编号。

## 2. Mathematical setup

PDE:

```text
(-Delta + sigma) u = f, sigma = 10
```

Boundary:

```text
u = 0 on boundary
```

Exact solution:

```text
u(x,y) =
sin(pi*x)sin(pi*y)
+ 0.5 sin(2*pi*x)sin(3*pi*y)
+ 0.3 sin(5*pi*x)sin(4*pi*y)
```

RHS is constructed mode-by-mode:

```text
f(x,y) =
(2*pi^2 + sigma) * sin(pi*x)sin(pi*y)
+ 0.5*(13*pi^2 + sigma) * sin(2*pi*x)sin(3*pi*y)
+ 0.3*(41*pi^2 + sigma) * sin(5*pi*x)sin(4*pi*y)
```

Confirmed:

- no `f = C * u_exact` unified-coefficient construction;
- mode `(1,1)` coefficient is `2*pi^2 + sigma`;
- mode `(2,3)` coefficient is `13*pi^2 + sigma`;
- mode `(5,4)` coefficient is `41*pi^2 + sigma`.

## 3. Numerical verification

Sanity check at `n=65`:

| Method | Linf error | L2 error | Boundary max |
|---|---:|---:|---:|
| FA | `1.932e-03` | `7.254e-04` | `0.0` |
| FFT9 | `1.320e-06` | `6.537e-07` | `0.0` |

Fitted convergence orders:

| Method | Linf order | L2 order |
|---|---:|---:|
| FA | `1.991` | `2.004` |
| FFT9 | `3.988` | `3.989` |

The result excludes the old suspicious `max error ~= 1.0e+00` behavior.

## 4. Figure integration

Code-side outputs:

- `code/python/experiments/results/exp06_complex_manufactured_convergence.csv`
- `code/python/experiments/figures/exp06_complex_manufactured_fields.png`
- `code/python/experiments/figures/exp06_complex_manufactured_convergence.png`

Thesis-side figures:

- `thesis/figures/exp06_complex_manufactured_fields.png`
- `thesis/figures/exp06_complex_manufactured_convergence.png`

Thesis integration:

- Modified `thesis/chapters/6_experiments.tex`.
- Inserted under experiment three as `多模态 manufactured solution 可视化补充`.
- No `实验七` or seventh core experiment was added.
- The thesis wording is limited to the current smooth Dirichlet multimode manufactured solution and does not claim coverage of all complex physical problems.

Figure checks:

- exact / FA / FFT9 field panels use common `vmin/vmax`;
- FA / FFT9 error panels use common log-error color scale;
- error titles contain max error;
- line cut includes exact / FA / FFT9;
- new figures do not show `max error ~= 1.0e+00`.

## 5. Documentation synchronization

Synchronized documents:

- `README.md`
- `FINAL_POLISH_REPORT.md`
- `VISUALIZATION_SUPPLEMENT_REPORT.md`
- `COMPLEX_MANUFACTURED_EXPERIMENT_REPORT.md`

Current synchronized status:

- final PDF page count: `56 pages`;
- core experiment count: `6`;
- supplementary visualization status: explicit;
- no wording that treats the supplement as an additional core experiment;
- `exp06_complex_manufactured_visualization.py` documented as supplementary and not parallel to exp00-exp05 core experiments.

Startup/git hygiene notes:

- current branch: `revise-sigma-fft9-krylov`;
- latest commit before this pass: `7de4157 fix: consolidate experiment chapter and final package consistency`;
- LaTeX intermediates exist under `thesis/` but are not staged and should not be committed;
- `__pycache__` directories exist but are not staged and should not be committed;
- old audit reports remain untracked process artifacts and should be committed only if process history is desired.

## 6. Validation

Experiment script:

```text
python -m experiments.exp06_complex_manufactured_visualization
PASS
Sanity n=65: FA Linf=1.932e-03, FA L2=7.254e-04, FFT9 Linf=1.320e-06, FFT9 L2=6.537e-07
Fitted orders: FA Linf=1.991, FA L2=2.004, FFT9 Linf=3.988, FFT9 L2=3.989
```

Pytest:

```text
python -m pytest -q
103 passed, 2 skipped, 2 warnings
```

XeLaTeX/BibTeX:

```text
xelatex; bibtex; xelatex; xelatex
PASS: Output written on main.pdf (56 pages).
```

Log scan:

```text
no fatal error
no undefined citation/reference
no multiply-defined label
```

Non-blocking warnings:

- existing Neumann Poisson compatibility warnings in pytest;
- SimSun bold-italic font fallback;
- hyperref bookmark token warnings;
- MiKTeX update notice.

## 7. Final status

READY TO SUBMIT
