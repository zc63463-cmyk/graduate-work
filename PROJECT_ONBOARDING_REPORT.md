# Project Onboarding Report

**Date**: 2026-04-29
**Project**: 矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究
**Phase**: 第一轮，只读审计

---

## 1. Repository Snapshot

| Item | Status |
|------|--------|
| Current branch | `revise-sigma-fft9-krylov` |
| Latest commit | `ce29f04 chore: integrate verified multimode visualization supplement` |
| Remote status | Ahead of `origin/revise-sigma-fft9-krylov` by 8 commits |
| Git status | **13 modified + 19 untracked** (uncommitted) |

### Modified (not staged)

- `README.md`, `FINAL_POLISH_REPORT.md`, `VISUALIZATION_SUPPLEMENT_REPORT.md`
- `thesis/main.pdf` (recompiled, now 57 pages)
- `thesis/chapters/6_experiments.tex`
- `thesis/figures/exp03_neumann_mixed.png`, `thesis/figures/exp05_resonance.png`
- `code/python/experiments/exp03_neumann_mixed.py`, `code/python/experiments/exp05_true_helmholtz_resonance.py`
- Corresponding CSV + PNG in `experiments/results/` and `experiments/figures/`

### Untracked files

- Audit/review reports: `EXPERIMENT_AUDIT_REPORT.md`, `PAPER_CODE_CONSISTENCY_REPORT.md`, `FORMULA_AUDIT_REPORT.md`, `FORMULA_AUDIT_MATRIX.md`, `LEAN_PROOF_MAP.md`, `LEAN_SCOPE_AUDIT.md`, `LEAN_APPENDIX_TABLE_SUGGESTION.md`
- New supplementary figures: `exp03_neumann_mixed_fields.png`, `exp03_neumann_mixed_summary.png`, `exp05_near_resonance_summary.png`
- `EXPERIMENT6_VISUAL_REDESIGN_REPORT.md`, `EXPERIMENT6_ROLLBACK_REPORT.md`, `EXPERIMENT_EVIDENCE_MATRIX.md`
- `CLAIM_DOWNGRADE_SUGGESTIONS.md`, `FILES_TO_PATCH.md`, `OPTIONAL_PATCH_PLAN.md`, `PATCH_COMPLETION_REPORT.md`, `QUESTIONS_FOR_GPT55_PRO.md`
- `codex_prompts/` directory

### LaTeX intermediates present

`thesis/*.aux`, `thesis/*.log`, `thesis/*.bbl`, `thesis/*.blg`, `thesis/*.toc`, `thesis/*.out`, `thesis/compile_log*.txt`, `thesis/compile_new.log`

### Assessment

Repository has **significant uncommitted changes** since the "READY TO SUBMIT" V2.3 state declared in `FINAL_POLISH_REPORT.md`. The supplementary visualization integrations are in the working tree but not committed. Multiple untracked audit reports and `codex_prompts/` suggest ongoing process work.

---

## 2. Project Structure

```
├── README.md                                    # V2.3, 声称 56 页 6 核心实验
├── FINAL_POLISH_REPORT.md                       # V2.3 最终清理报告
├── VISUALIZATION_SUPPLEMENT_REPORT.md            # 可视化补充报告
├── COMPLEX_MANUFACTURED_EXPERIMENT_REPORT.md     # 多模态 manufactured solution 报告
├── FINAL_COMPLEX_VISUALIZATION_INTEGRATION_REPORT.md
├── EXPERIMENT6_VISUAL_REDESIGN_REPORT.md         # (untracked)
├── EXPERIMENT6_ROLLBACK_REPORT.md                # (untracked)
├── various audit/consistency reports...          # (most untracked)
│
├── thesis/
│   ├── main.tex         # 主文件，含摘要、目录、7章+附录+参考文献
│   ├── main.pdf         # 编译输出 (57 pages after full compile)
│   ├── references.bib
│   └── chapters/
│       ├── 1_introduction.tex
│       ├── 2_math_preliminary.tex
│       ├── 3_fft_direct.tex
│       ├── 4_helmholtz.tex
│       ├── 5_gmres.tex
│       ├── 6_experiments.tex      # 6 核心实验 + 可视化补充
│       ├── 7_conclusion.tex
│       └── appendix_lean4.tex
│
├── code/python/
│   ├── helmholtz_solver.py        # 核心: FA/CR/FACR/FFT9/5pt (1463 lines)
│   ├── gmres_solver.py            # GMRES: restart + matrix-free (1065 lines)
│   ├── cyclic_reduction.py        # [deprecated] Poisson-only 原型
│   ├── fft9_complete.py           # [deprecated] Poisson-only FFT9 原型
│   ├── fft9_eigenvalue_proof.py   # SymPy 六阶不可达性证明
│   ├── verify_fft9_expansion.py   # SymPy FFT9 展开验证
│   ├── experiments/
│   │   ├── exp00_fft_vs_sparse.py
│   │   ├── exp01_convergence.py
│   │   ├── exp02_nonhom_bc.py
│   │   ├── exp03_neumann_mixed.py
│   │   ├── exp04_modified_vs_true.py
│   │   ├── exp05_true_helmholtz_resonance.py
│   │   ├── exp06_complex_manufactured_visualization.py  # 补充可视化
│   │   ├── generate_visualization_supplements.py
│   │   ├── utils.py
│   │   ├── results/               # CSV
│   │   └── figures/               # PNG
│   └── tests/
│       ├── test_01_eigenvalues.py
│       ├── test_02_nonhom_dirichlet.py
│       ├── test_03_fft_vs_spsolve.py
│       ├── test_04_fft9_vs_spsolve.py
│       ├── test_05_neumann_compatibility.py
│       ├── test_06_modified_true_helmholtz.py
│       └── test_07_krylov_baselines.py
│
└── code/lean4_formalization/
    ├── SixthOrderImpossibility.lean     # 六阶不可达性形式化证明
    ├── lakefile.toml
    └── lean-toolchain
```

---

## 3. Mathematical Mainline

### 3.1 Unified Sigma Framework

论文和代码统一使用：

$$(-\nabla^2 + \sigma) u = f$$

| $\sigma$ | 方程类型 | 特征 |
|----------|---------|------|
| $0$ | Poisson | 稳态扩散 |
| $+\kappa^2 > 0$ | Modified Helmholtz (screened Poisson) | 矩阵正定 |
| $-\kappa^2 < 0$ | True Helmholtz | 矩阵可能不定，近共振风险 |

代码实现：`helmholtz_solver.py:_resolve_sigma()` (line 55-84)，优先使用 `sigma` 参数，向后兼容 `k2`（映射为 `sigma = +k^2`）。

### 3.2 Five-Point Eigenvalues

**Dirichlet (DST-I)**：$N = n-2$ interior points

$$\mu_p = 2 - 2\cos\left(\frac{p\pi}{N+1}\right), \quad p = 1,\ldots,N$$

二维频域分母：

$$d_{p,q} = \frac{\mu_p + \mu_q}{h^2} + \sigma$$

**Neumann (DCT-I, ghost-point)**：$N = n$ full nodes

$$\lambda_k = \frac{2}{h^2}\left(1 - \cos\left(\frac{k\pi}{n-1}\right)\right), \quad k = 0,\ldots,n-1$$

实现：`helmholtz_solver.py:110-149`（`_eigenvalues_dirichlet`, `_eigenvalues_neumann`）

**无 double-counting 问题**：频域分母公式 `d_{p,q} = (mu_p + mu_q)/h^2 + sigma` 在代码中直接实现（line 648），无需额外 +2 修正。

### 3.3 FFT9 Compact Scheme

九点算子 $L_h$ 近似 $\nabla^2$（非 $-\nabla^2$）：

$$L_h = \frac{1}{6h^2}\begin{bmatrix}1 & 4 & 1 \\ 4 & -20 & 4 \\ 1 & 4 & 1\end{bmatrix}$$

等价系统（论文/代码自洽，`helmholtz_solver.py:988-1050`）：

$$(-L_h + \sigma R_h) u = R_h f \quad\iff\quad (L_h - \sigma R_h) u = -R_h f$$

代码使用后者（乘以 -1）。

Fourier symbol 四阶验证（`verify_fft9_expansion.py`）：

$$\frac{-\lambda_L}{\lambda_R} = (\xi^2 + \eta^2) + O(h^4)$$

$h^2$ 系数 = 0（`Is zero? True`），四阶精度确认。

**当前限制**：FFT9 仅 Dirichlet BC；Neumann/mixed 自动 fallback 到 FA（`helmholtz_solver.py:1007-1013`）。

### 3.4 Neumann DCT-I

Ghost-point: $u_{-1}=u_1$。非对称矩阵 G 对称化为 $S = D^{-1}GD$，由 DCT-I 对角化。

$D = \text{diag}(\sqrt{2}, 1, \ldots, 1, \sqrt{2})$

Neumann Poisson 零模态（$\lambda_0=0$）需要兼容条件 $\int f = 0$。`check_neumann_compatibility()` 检查离散积分。

### 3.5 Modified vs True Helmholtz

| 特性 | Modified ($\sigma > 0$) | True ($\sigma < 0$) |
|------|------------------------|---------------------|
| 频域分母 | $\lambda_{p,q} + \sigma > 0$ | $\lambda_{p,q} - \|\sigma\|$ 可为零 |
| 正定性 | 正定 | 可能不定 |
| 条件数 | 随 $\sigma$ 增大改善 | 近共振时恶化 |
| GMRES | 迭代下降 | 可能停滞/未收敛 |

Near-resonance: $|d_{p,q}| = |\lambda_{p,q} - \kappa^2| \approx 0$，解范数 $\propto 1/|d_{p,q}|$。

### 3.6 GMRES

当前为无预处理 Krylov 基线：
- `gmres_solver.py` 自实现 restart GMRES（Givens rotation）
- `build_helmholtz_matrix()` 构建 5-pt 稀疏矩阵
- `gmres_helmholtz_matfree()` 矩阵无关版本
- **不声称** shift-Laplacian、ILU 或任何预处理

---

## 4. Core Experiments (Chapter 6)

### 4.1 Six Core Experiments

| 编号 | 名称 | 脚本 | CSV | 主图 | 结论边界 |
|------|------|------|-----|---------|---------|
| 一 | FFT vs sparse 一致性 | `exp00_fft_vs_sparse.py` | `exp00_fft_vs_sparse.csv` | 表格 | 离散系统一致，非 PDE 收敛 |
| 二 | 二阶与四阶收敛阶 | `exp01_convergence.py` | `exp01_convergence.csv` | `exp01_convergence.png` | 五点二阶/FFT9 Dirichlet 四阶 |
| 三 | 非齐次 Dirichlet | `exp02_nonhom_bc.py` | `exp02_nonhom_bc.csv` | `exp02_nonhom_bc.png` | RHS 边界符号正确 |
| 四 | Neumann + mixed | `exp03_neumann_mixed.py` | `exp03_neumann_mixed.csv` | `exp03_neumann_mixed_summary.png` | 五点 ghost-point 二阶；非 FFT9 |
| 五 | Modified/True Helmholtz 谱 + GMRES | `exp04_modified_vs_true.py` | `exp04_modified_vs_true.csv` | 3 张图 | Gaussian RHS 下无预处理 GMRES 参数敏感 |
| 六 | Near-resonance 扫描 | `exp05_true_helmholtz_resonance.py` | `exp05_resonance.csv` | `exp05_near_resonance_summary.png` | 近共振放大/收敛困难 |

### 4.2 Supplementary Visualizations (不计入核心实验)

| 图 | 脚本 | 用途 |
|----|------|------|
| `exp02_temperature_field_comparison.png` | `generate_visualization_supplements.py` | 温度场空间结构 |
| `exp05_denominator_heatmap.png` | `exp05_true_helmholtz_resonance.py` | 分母热图备份 |
| `exp06_complex_manufactured_fields.png` | `exp06_complex_manufactured_visualization.py` | 多模态解析/数值解/误差 |
| `exp06_complex_manufactured_convergence.png` | `exp06_complex_manufactured_visualization.py` | 多模态收敛补充 |
| `exp03_neumann_mixed_fields.png` | `exp03_neumann_mixed.py` | Neumann/Mixed 解场 |

**命名歧义警告**：`exp06_complex_manufactured_visualization.py` 文件名含 `exp06`，但文档已注明为 supplementary，不计入核心实验编号。仍有混淆风险。

### 4.3 旧版实验

旧版 16 个历史实验（heatmap、restart、FACR 参数图、复杂度趋势图）已在 FINAL_POLISH_REPORT 声明废弃。

---

## 5. Code Architecture

### 5.1 Core Solver (`helmholtz_solver.py`)

- `fa_helmholtz()` — 2D DST/DCT 对角化直接解
- `cr_helmholtz()` — 1D 变换 + Thomas 求解 y 三对角
- `facr_helmholtz()` — CR + FA 混合，vectorized 奇偶约化
- `fft9_helmholtz()` — 九点紧致格式 (Dirichlet only)
- `point5_helmholtz()` — 五点基准
- `solve_helmholtz()` — 统一接口
- `check_resonance()` — 频域共振检测
- `check_neumann_compatibility()` — Neumann 兼容条件检测

### 5.2 GMRES Solver (`gmres_solver.py`)

- `gmres()` — 自实现 restart GMRES (Givens)
- `build_helmholtz_matrix()` — 稀疏矩阵构建 (Dirichlet/Neumann/mixed)
- `gmres_helmholtz()` — 完整 GMRES + Helmholtz 接口
- `gmres_helmholtz_matfree()` — LinearOperator 矩阵无关版本
- 测试问题：eigenfunction、polynomial、Gaussian、multimode

### 5.3 Deprecated Files

- `cyclic_reduction.py` — Poisson-only FA/CR/FACR 原型，已同步修正 RHS 边界符号
- `fft9_complete.py` — Poisson-only FFT9 原型
- 论文不应引用这些文件

### 5.4 Experiment Dependencies

```
utils.py (共用: test problems, equation_type, path tools)
  ├── exp00 → helmholtz_solver + gmres_solver.build_helmholtz_matrix
  ├── exp01 → helmholtz_solver
  ├── exp02 → helmholtz_solver
  ├── exp03 → helmholtz_solver (FA/CR/FACR)
  ├── exp04 → helmholtz_solver + gmres_solver
  ├── exp05 → helmholtz_solver + gmres_solver
  └── exp06 → helmholtz_solver (FA/FFT9)
```

---

## 6. Validation Status

### 6.1 Pytest

```text
103 passed, 2 skipped, 2 warnings
```

- 2 warnings: Neumann Poisson compatibility 预期 warning（故意不相容 RHS）
- Status: **PASS**

### 6.2 FFT9 Expansion Verification

```text
verify_fft9_expansion.py:
  h^2 coefficient = 0 → Is zero? True
  h^0 coefficient = eta^2 + xi^2 → Match: True
```

Status: **PASS**

### 6.3 FFT9 Eigenvalue Proof

确认 $c_2 = 0$（四阶），$c_4 \neq 0$（六阶不可达），Poisson 和 Helmholtz 均无通用解。

Status: **PASS**

### 6.4 Lean 4 Build

Per FINAL_POLISH_REPORT: `lake build → Build completed successfully`

Status: **PASS** (未重复验证，仍有 unused simp warnings)

### 6.5 LaTeX / PDF

**Full compile (xelatex → bibtex → xelatex → xelatex)**:

```text
Output written on main.pdf (57 pages).
0 LaTeX error
0 undefined citation
0 multiply-defined label
```

Non-blocking: SimSun font substitute warning, hyperref bookmark token warning.

### 6.6 README vs Reality Discrepancy

| 声明 | README/FINAL_POLISH | Actual |
|------|-------------------|--------|
| PDF 页数 | 56 pages | **57 pages** |
| 核心实验数 | 6 | 6 (consistent) |
| 测试通过 | 103 passed | 103 passed (consistent) |

**差异原因**：Uncommitted changes 向 `6_experiments.tex` 添加了 supplementary visualization 段落（温度场、多模态 manufactured solution），导致增加 1 页。

---

## 7. Risk Register

### A. Must Fix Before Submission

| # | 路径 | 问题 | 风险 | 建议 |
|---|------|------|------|------|
| A1 | `README.md:6` + `FINAL_POLISH_REPORT.md:7` | 声称 56 pages，实际 PDF 为 **57 pages** | README 与交付物不一致 | 更新页数为 57 |
| A2 | `thesis/chapters/6_experiments.tex` (modified but uncommitted) | 包含未经 Git 提交的补充可视化段落 | 论文定稿状态不确定 | 确认补充图段落是否全部保留，提交或回退 |
| A3 | All modified tracked files | 13 个文件未提交 | 当前 working tree 不代表任何明确版本 | 做一次完整的 `git commit` 收束当前状态，或 `git stash` 回退到已知状态 |

### B. Should Fix If Time Allows

| # | 路径 | 问题 | 风险 | 建议 |
|---|------|------|------|------|
| B1 | `exp06_complex_manufactured_visualization.py` | 文件名含 `exp06`，与六核心 `exp00`-`exp05` 命名重叠 | 读者/审核人混淆编号 | 重命名或在 README 醒目注明 |
| B2 | 19 untracked files | 包含十余个未决审计报告和 `codex_prompts/` | 混乱最终交付物 | 归档到 `notes/` 或 `.workbuddy/` |
| B3 | `thesis/main.log` 等编译中间文件 | tracked as modified | 不应随论文提交 | 确保 `.gitignore` 覆盖，最终提交排除 |
| B4 | `code/python/_archive/` | 未检视的存档目录 | 可能包含冗余数据 | 确认不干扰最终实验 |

### C. Non-blocking

| # | 路径 | 问题 | 建议 |
|---|------|------|------|
| C1 | Lean build | unused simp argument warnings | 不影响构建 |
| C2 | LaTeX | SimSun font + hyperref bookmark warnings | 不影响 PDF |
| C3 | `cyclic_reduction.py`, `fft9_complete.py` | deprecated 但未明确标注在何处 | 确认论文不引用 |
| C4 | `.lake/packages/proofwidgets` | local changes | 不提交 |

---

## 8. Recommended Next Actions

以下为建议的下一步行动（需确认后执行）：

1. **Resolve page-count mismatch**：README/FINAL_POLISH_REPORT 中的 56 改为 57；
2. **Commit or revert uncommitted changes**：所有 modified files 需要收束为一次干净提交；
3. **Clean untracked files**：将 audit 报告移到 `notes/` 或加入 `.gitignore`，删除 `codex_prompts/`（如不需保留）；
4. **Rename `exp06` script**：`exp06_complex_manufactured_visualization.py` → `supp_complex_manufactured_visualization.py`（可选，但建议避免编号混淆）；
5. **Final full compile check**：在最终提交后重新运行 `xelatex → bibtex → xelatex → xelatex` 确认 0 error/undefined；
6. **Verify Lean build**：在 final 提交状态下重新 `lake build`；
7. **Push to remote**：提交后推送（当前 ahead 8 commits）。

---

## 9. Questions for User

1. 当前 `6_experiments.tex` 的 uncommitted 修改（补充可视化段落和新增 figure includes）是否应全部保留进入最终论文？还是部分回退？
2. 十九个 untracked 文件（审计报告、patch 计划等）应如何处置：归档到 `notes/`、删除、还是提交为过程记录？
3. `exp06_complex_manufactured_visualization.py` 是否重命名为不含 `exp06` 的名字？
4. 是否需要在下一轮执行上述 Recommended Next Actions 中的修改，还是仅等待确认？
