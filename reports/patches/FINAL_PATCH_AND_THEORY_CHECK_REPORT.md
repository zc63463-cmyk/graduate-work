# Final Patch and Theory Check Report

生成时间：2026-04-29

## 1. Summary

本轮以实际源文件为准完成最终确认，并做了 3 处安全小修：

- 将 Ch4 总览表中 true Helmholtz 的“正定性/适用迭代法”从过宽的“不定/MINRES-GMRES”改为“依谱而定”。
- 处理 C-001：将 Ch2 中 “对任意非零模式 h^4 项系数不为零” 改为 “该系数不恒为零，误差一般为 O(h^4)”。
- 处理 C-002：将 Ch3 中 FACR “第三项趋于常数” 改为 “第三项降为 O(N^2) 量级”，并再次明确当前 FACR-like 实现整体仍为 O(N^2 log N)。

同步更新：

- README PDF 页数：56 pages。
- FINAL_POLISH_REPORT PDF 页数：56 pages。
- PATCH_COMPLETION_REPORT：C 类风险已处理，Final status 改为 READY TO SUBMIT。

启动检查记录：

- 当前分支：`revise-sigma-fft9-krylov`。
- 最新 commit：`7de4157 fix: consolidate experiment chapter and final package consistency`。
- deleted 文件：无。
- `.lake/packages/proofwidgets`：存在 local changes，具体为 `lakefile.lean`。
- 存在 LaTeX 中间文件和 `__pycache__`，均不应提交。

## 2. Source-of-truth verification

| Item | Source verification | Result |
|---|---|---|
| A-001 Ch2 FFT9 `L_h` 符号 | `thesis/chapters/2_math_preliminary.tex` 中保留 `L_h \approx \Delta = \nabla^2`；主格式为 `(-L_h + \sigma R_h)u = R_h f`；分解为 `-L_h = L_h^{(5)} - h^2 C_h/6`。未发现 `L_h \approx -\Delta` 残留。 | PASS |
| A-002 三参数 `lambda_mod` | Ch2 有效特征值为 `(-\hat{\lambda}_L - \alpha h^2\hat{\lambda}_5 + \sigma(\hat{\lambda}_R+\beta h^2))/(\hat{\lambda}_R+\gamma h^2)`；未发现旧式正负号反转公式。 | PASS |
| A-003 Lean 附录公式和范围 | `appendix_lean4.tex` 与 Ch2 corrected fraction 一致；明确不覆盖完整 PDE 全局误差估计、DST/DCT 谱理论、Python 求解器和实验脚本正确性；六阶不可达性限定为三参数九点修正族。自有 Lean 源文件未发现 `sorry/admit/axiom/unsafe`。 | PASS |
| A-004 true Helmholtz 谱分类 | Ch4 给出完整分类：`kappa^2 < lambda_min` 正定，等于离散特征值奇异，谱区间内非特征值不定，大于 `lambda_max` 负定；摘要/引言/结论使用“可能不定/近共振风险”而非必然不定。 | PASS |
| B-001 near-resonance `delta` 符号 | Ch4 使用 `kappa^2=lambda_{p*,q*}+delta`，并写明 `d_{p*,q*}=-delta`、`\hat u=-\hat f/delta`；未发现 `d=delta` 残留。 | PASS |

Lean theorem count note：

- `SixthOrderImpossibility.lean` 可执行整数域 theorem 声明为 14 个；文件末尾另有注释块中的形式化 outline，不计入源码 theorem 数。
- `SixthOrderImpossibility/MathlibTest.lean` 实数域 theorem 声明为 10 个。

## 3. C-class refinements

| Item | Change | Result |
|---|---|---|
| C-001 h^4 系数全称措辞 | 改为 “该 h^4 项系数不恒为零；例如取 xi=eta 非零时该项非零，因此误差一般为 O(h^4)” | DONE |
| C-002 FACR 复杂度第三项 | 改为 “第三项降为 O(N^2) 量级；经典 FACR 理论最优可达 O(N^2 log log N)；当前 FACR-like 实现包含完整方向 DST，整体仍为 O(N^2 log N)” | DONE |

## 4. Tests

| Check | Command | Result |
|---|---|---|
| pytest | `cd code/python; python -m pytest -q` | PASS: 103 passed, 2 skipped, 2 warnings |
| FFT9 expansion | `cd code/python; python verify_fft9_expansion.py` | PASS: h2 coefficient is zero; effective eigenvalue is `xi^2+eta^2+O(h^4)` |
| Lean build | `cd code/lean4_formalization; lake build` | PASS: Build completed successfully |
| XeLaTeX/BibTeX | `cd thesis; xelatex; bibtex; xelatex; xelatex` | PASS: `main.pdf` generated, 56 pages |
| LaTeX log scan | `Select-String thesis\main.log ...` | PASS: no fatal error, no undefined citation/reference, no multiply-defined label, no final Overfull warning |
| Git whitespace | `git diff --check` | PASS: no whitespace error |

Non-blocking warnings:

- pytest emits expected Neumann Poisson compatibility warnings in compatibility tests.
- Lean emits unused `simp` argument warnings and reports `.lake/packages/proofwidgets` local changes.
- XeLaTeX emits SimSun bold-italic font fallback and hyperref bookmark token warnings.
- Git reports LF to CRLF normalization notices on Windows.

## 5. Files changed

Source/theory/report files changed in this final pass:

- `thesis/chapters/2_math_preliminary.tex`
- `thesis/chapters/3_fft_direct.tex`
- `thesis/chapters/4_helmholtz.tex`
- `README.md`
- `FINAL_POLISH_REPORT.md`
- `PATCH_COMPLETION_REPORT.md`
- `FINAL_PATCH_AND_THEORY_CHECK_REPORT.md`
- `thesis/main.pdf`

Existing modified files from previous patch rounds remain in the working tree:

- `code/python/experiments/exp05_true_helmholtz_resonance.py`
- `code/python/experiments/results/exp05_resonance.csv`
- `code/python/tests/test_07_krylov_baselines.py`
- `thesis/main.tex`
- `thesis/chapters/1_introduction.tex`
- `thesis/chapters/6_experiments.tex`
- `thesis/chapters/7_conclusion.tex`
- `thesis/chapters/appendix_lean4.tex`

Untracked audit/process reports remain present and should be committed only if process history is desired.

## 6. Remaining risks

No blocking mathematical or verification risks remain.

Commit hygiene notes:

- Do not submit LaTeX intermediates such as `.aux`, `.log`, `.blg`, `.out`, `.toc`.
- Do not submit `__pycache__`.
- Do not submit `.lake/packages/proofwidgets`.
- Old figure files still exist under `thesis/figures/`, but final thesis source no longer references the old heatmap/restart/FACR/timing-style figures.
- `rg` was unavailable in this shell due `Access is denied`; equivalent `Select-String` source checks were used.

## 7. Final status

READY TO SUBMIT
