# Pro Review Fix Completion Report

**Date**: 2026-04-28  
**Branch**: revise-sigma-fft9-krylov  
**Base commit**: 87259b2

---

## 1. Summary

本次修复针对 Pro 终稿裁决报告中指出的 A 类问题、额外硬错误和 B 类论文修复项，共完成 6 项 P0 修复、2 项 P1 降级、13 项 B 类修复以及 P4/P5 GMRES 残差问题处理。所有修复遵循"最小必要修复"原则，未推翻论文主线。

---

## 2. A-class Issues

### P2 (GMRES iterations 全为 1)
- **状态**: ✅ 已修复（方案 B：降级论文结论）
- **修复**: 在第 6.17 节中明确说明 iterations=1 是因为右端项 $f=\sin(2\pi x)\sin(3\pi y)$ 主要激发少数 Fourier 模态，不以该迭代次数判断一般 RHS 下的 GMRES 收敛难度，重点改为谱指标 $d_{\min}$ 和条件数对比。
- **文件**: `thesis/chapters/6_experiments.tex`

### P10 (FFT9 RHS/分母不一致)
- **状态**: ✅ 已修复
- **修复**: 统一采用 raw denominator $D_{p,q}^{\mathrm{raw}} = -\hat\lambda_L + \sigma\hat\lambda_R$，物理空间右端项 $\tilde{g}=(R_h f)_I + L_{IB}g_B - \sigma R_{IB}g_B$。Algorithm 3 已改写为 raw denominator 方案，Ch3 和 Ch4 中的频域公式已统一。
- **文件**: `thesis/chapters/3_fft_direct.tex`, `thesis/chapters/4_helmholtz.tex`

### P11 (五点 Dirichlet RHS 符号 + FFT9 边界修正符号)
- **状态**: ✅ 已修复
- **修复 1**: 五点非齐次 Dirichlet 右端项公式从减号改为加号：$\tilde{f}_{i,j} = f_{i,j} + \frac{1}{h^2}(\ldots)$，因为边界值在差分格式中以 $-g/h^2$ 出现在左端，移至右端变加号。
- **修复 2**: FFT9 非齐次 Dirichlet 边界修正公式符号修正：$+1/(6h^2)[\ldots]$ 为 $L_h$ 贡献（移项后加号），$-\sigma/12\,g_D$ 为 $\sigma R_h$ 贡献（移项后减号）。
- **文件**: `thesis/chapters/2_math_preliminary.tex`, `thesis/chapters/3_fft_direct.tex`

### P17 (near-resonance 实验缺失)
- **状态**: ✅ 已修复（方案 B：降级结论）
- **修复**: 在第 6.17 节和第 7 章中明确说明 near-resonance 作为理论分析，实验验证留作后续工作。不再声称"实验揭示近共振迭代困难"。
- **文件**: `thesis/chapters/6_experiments.tex`, `thesis/chapters/7_conclusion.tex`

---

## 3. Extra Hard Errors

### 第 6.16 非齐次 Dirichlet manufactured RHS
- **状态**: ✅ 已修复
- **修复**: 将 $f = ((2\pi)^2+(3\pi)^2+\sigma)u(x,y)$ 修正为 $f = ((2\pi)^2+(3\pi)^2+\sigma)\sin(2\pi x)\sin(3\pi y) + \sigma$，明确常数 1 的 Laplacian 为零。
- **文件**: `thesis/chapters/6_experiments.tex`

### 第 6.17 统一方程符号
- **状态**: ✅ 已修复
- **修复**: 将 $\nabla^2 u - \sigma u = f$ 和 $\nabla^2 u + k^2 u = f$ 改为 $(-\nabla^2 + \sigma)u = f$ 和 $(-\nabla^2 - \kappa^2)u = f$，与全文统一框架一致。
- **文件**: `thesis/chapters/6_experiments.tex`

### Rh_boundary_freq 引理
- **状态**: ✅ 已降级
- **修复**: 将未证明的 Lemma `Rh_boundary_freq` 降级为 Remark，删除 $\hat\lambda_R^{-1}\hat{b}$ 频域边界修正公式，改为强调物理空间移项方案。
- **文件**: `thesis/chapters/3_fft_direct.tex`

---

## 4. Experiments

### exp04 GMRES iterations
- **未重跑**。论文已降级相关结论，不以 iterations=1 支撑迭代行为对比。

### near-resonance 实验
- **未补**。论文已降级结论，near-resonance 作为理论分析，实验验证留作后续工作。

### exp02 非齐次 Dirichlet — ✅ 已完全修复
- **原问题**：代码使用 `test_problem_dirichlet_mode(sigma, m=2, n=3)` 生成 $u = \sin(2\pi x)\sin(3\pi y)$，但论文写 $u = \sin(2\pi x)\sin(3\pi y) + 1$。由于 $\sin(2\pi \cdot 0) = \sin(2\pi \cdot 1) = 0$，代码实际是齐次 Dirichlet（BC=0），而非声称的非齐次 Dirichlet。
- **修复**：统一为 $u = \sin(2\pi x)\sin(3\pi y) + 1$，$f = ((2\pi)^2+(3\pi)^2+\sigma)\sin(2\pi x)\sin(3\pi y) + \sigma$，$g_D = 1$ on $\partial\Omega$。
- **修改文件**：
  - `code/python/experiments/utils.py` — 新增 `test_problem_nonhom_dirichlet()` 函数
  - `code/python/experiments/exp02_nonhom_bc.py` — 使用新函数，更新文档字符串
  - `code/python/helmholtz_solver.py` — inline test 改为 u=sin+1, f 正确, bc=1
  - `code/python/tests/test_02_nonhom_dirichlet.py` — `_make_nonhom_problem` 改为 u=sin+1, bc=1
  - `code/python/tests/test_04_fft9_vs_spsolve.py` — `test_fft9_nonhom_bc` 改为 u=sin+1, bc=1
  - `thesis/chapters/6_experiments.tex` — 补充 $g_D=1$ 说明，$sin(2\pi x)sin(3\pi y)=0$ on $\partial\Omega$
- **重跑结果**：exp02 PASS，收敛阶 FA/CR/FACR ~2.00-2.03，FFT9 ~4.70
- **pytest**：103 passed, 2 skipped ✅
- **FFT9 vs sparse**：test_fft9_nonhom_bc PASS ✅（旧版 FAILED 因为声称非齐次但实际齐次）

---

## 5. B-class Fixes

| ID | 修复内容 | 状态 | 文件 |
|----|---------|------|------|
| B1 | manufactured solution 局限说明 | ✅ | 6_experiments.tex |
| B2 | FFT9 Neumann 支持范围 + 代码 warnings.warn | ✅ | 3_fft_direct.tex, helmholtz_solver.py |
| B3 | fft9_complete.py 弃用说明 | ✅ | fft9_complete.py |
| B4 | 非正方形域支持范围说明 | ✅ | 6_experiments.tex |
| B5 | 混合边界测试覆盖（无需修改，论文无"全面验证"措辞） | ✅ | — |
| B6 | 条件数渐近表述修正 | ✅ | 2_math_preliminary.tex |
| B7 | 六阶不可行补足说明 + Lean 4 验证范围 | ✅ | 2_math_preliminary.tex |
| B8 | "四个数量级"措辞修正 | ✅ | 6_experiments.tex |
| B9 | GMRES 收敛界适用前提 | ✅ | 5_gmres.tex |
| B10 | modified Helmholtz 结论限定 | ✅ | 6_experiments.tex, 7_conclusion.tex |
| B11 | Timing 无统计性说明 | ✅ | 6_experiments.tex |
| B12 | Lean 4 范围说明 | ✅ | 7_conclusion.tex |
| B13 | FACR-like 加速来源说明 | ✅ | 6_experiments.tex |

---

## 6. GMRES Residual Issues (P4/P5)

### P4: 残差类型
- **方案 B**: 不改代码，修改文档说明 `tol` 为 absolute residual tolerance。
- **文件**: `code/python/gmres_solver.py`（docstring 修改），`thesis/chapters/6_experiments.tex`（容差描述修改）

### P5: Givens residual vs 重启点 residual
- 无需大改，两种 residual 用途已在代码注释中说明。

---

## 7. Tests

| 测试 | 结果 |
|------|------|
| pytest (103 passed, 2 skipped) | ✅ PASS |
| verify_fft9_expansion.py | ✅ PASS |
| lake build (Lean 4) | ✅ PASS (Build completed successfully, 3 jobs) |
| XeLaTeX (xelatex + bibtex + xelatex × 2) | ✅ PASS (65 pages, 0 error, 0 undefined ref/cite, 0 multiply-defined label) |

---

## 8. Remaining Risks

### High Risk
（无）

### Medium Risk
- **near-resonance 未补实验**: 论文已降级结论为理论分析，实验验证留作后续工作。如评审要求实验证据，需补 exp05。**（已降级为 Low Risk，见下方）**

### Low Risk
- **near-resonance 未补实验**: 论文已降级结论为理论分析/后续工作，不影响论文正确性。
- **GMRES 仍使用 absolute residual**: 论文已明确说明，不声称 relative residual。
- **非正方形域未系统测试**: 论文已说明实验仅在 $[0,1]^2$ 下进行。
- **Timing 无统计性**: 论文已说明单次运行，绝对时间仅供参考。
- **SimSun bold italic font shape undefined**: LaTeX 警告，不影响内容，默认字体替代。

---

## 9. Modified Files

### 论文
- `thesis/chapters/2_math_preliminary.tex` — 五点 RHS 符号、条件数渐近、六阶补充说明
- `thesis/chapters/3_fft_direct.tex` — Algorithm 3 raw denominator、频域公式统一、Rh_boundary_freq 降级、边界修正符号、FFT9 Neumann 限定
- `thesis/chapters/4_helmholtz.tex` — D^(9) 改为 D^raw
- `thesis/chapters/5_gmres.tex` — GMRES 收敛界适用前提
- `thesis/chapters/6_experiments.tex` — 6.16 RHS 修正、6.17 符号统一、iterations=1 降级、near-resonance 降级、四数量级措辞、modified Helmholtz 限定、Timing 说明、FACR 说明、域限定、容差说明
- `thesis/chapters/7_conclusion.tex` — GMRES 结论限定、near-resonance 降级、Lean 4 范围说明
- `thesis/main.pdf` — 重新编译

### 代码
- `code/python/experiments/utils.py` — 新增 `test_problem_nonhom_dirichlet()` 函数
- `code/python/experiments/exp02_nonhom_bc.py` — 使用新的非齐次 Dirichlet 制造解
- `code/python/helmholtz_solver.py` — FFT9 Neumann 回退 warnings.warn + inline test 改为 u=sin+1
- `code/python/gmres_solver.py` — docstring 修改为 absolute residual
- `code/python/fft9_complete.py` — 弃用说明
- `code/python/tests/test_02_nonhom_dirichlet.py` — `_make_nonhom_problem` 改为 u=sin+1, bc=1
- `code/python/tests/test_04_fft9_vs_spsolve.py` — `test_fft9_nonhom_bc` 改为 u=sin+1, bc=1

---

## 10. Final Status

**READY TO SUBMIT**

- 所有 A 类问题和硬错误已修复或降级
- 所有 B 类修复已完成
- **exp02 非齐次 Dirichlet 一致性已完全修复**（论文/代码/u/f/bc 四统一）
- 全部测试通过（pytest 103 passed / Lean build PASS / XeLaTeX 65 pages 0 error）
- Low risk: near-resonance 未补实验（已降级结论）；GMRES absolute residual（已说明）
- 论文不再包含混搭 RHS/分母、错误符号、过度声称、论文/代码不一致等问题
