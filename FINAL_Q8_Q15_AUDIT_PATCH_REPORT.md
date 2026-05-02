# FINAL Q8--Q15 Audit Patch Report

## 1. 本轮目标

本轮按 `final_prompt_and_patch_text_for_codex_q8_q15.md` 做最后核查与最小文字修补。范围限定为论文文字、Lean 覆盖口径和最终核查报告：

- 不新增实验；
- 不改实验参数；
- 不重绘图；
- 不改 CSV/PNG/PDF 实验资产；
- 不修改 FFT9、GMRES 或其他数值算法。

## 2. FFT9 final code convention

```text
FFT9 final code convention: raw/theory
```

当前最终代码与论文理论系统一致，采用

\[
(-L_h+\sigma R_h)u=R_hf,
\]

\[
\tilde g=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\]

\[
D^{raw}_{p,q}=-\widehat\lambda_L(p,q)+\sigma\widehat\lambda_R(p,q).
\]

### 2.1 代码证据

`code/python/helmholtz_solver.py`：

- RHS 构造：
  - line 1039: `rhs_tilde = apply_Rh_full(F)`
  - line 1044: `rhs_tilde += bc_corr`
- 频域分母：
  - line 1064: `denom = -lam_L + sigma * lam_R4`
- 注释口径：
  - line 930: `This implementation uses the same raw sign convention as the thesis.`
  - line 1035--1038: 明确 `rhs_tilde = (R_h f)_I + L_IB g_B - sigma R_IB g_B` 与 `-lambda_L + sigma*lambda_R` 配套。

`code/python/tests/test_04_fft9_vs_spsolve.py`：

- line 10: sparse compact system 写为 `(-L_h + sigma * R_h) * u = R_h * f + bc_correction`
- line 54--60: 对照测试中同样使用 `rhs_tilde = apply_Rh_full(F)` 并加上 `bc_corr`

因此，本轮没有采用“当前 Python 实现为整体反号系统”的解释分支。

## 3. 第 3/4 章 implementation note 状态

本轮未修改第 3.4.3 / 3.4.5 / 4.1.2 附近的 implementation note，因为当前文本已经与仓库事实一致：

- 论文理论推导和当前 Python 主实现均采用 raw/theory 系统；
- 若其他实现整体乘以 \(-1\)，则 RHS 与 denominator 必须同步反号；
- 两种符号系统不能混用，尤其不能只反 RHS 或只反 denominator。

本轮仅在最终报告中明确记录该核查结论。

## 4. Lean theorem table 核查

### 4.1 当前 Lean theorem 事实

`code/lean4_formalization/SixthOrderImpossibility/MathlibTest.lean` 中：

- theorem `c2_helmholtz_forces_neg_gamma` 的结论是

\[
\alpha+\gamma=0,
\]

即 \(\alpha=-\gamma\)。

该 theorem 没有同时把 \(\beta=\gamma\) 作为结论。

### 4.2 论文附录处理

`thesis/chapters/appendix_lean4.tex` 中 theorem table 保持为：

```latex
Helmholtz 情形下 \(c_2=0\) 推出必要关系 \(\alpha=-\gamma\)
```

没有改成“Lean 证明 \(\alpha=-\gamma,\beta=\gamma\)”。

本轮新增说明：

- `c2_helmholtz_forces_neg_gamma` 直接验证的是 Helmholtz 情形下 \(c_2=0\) 所需参数关系中的 \(\alpha=-\gamma\) 代数核心；
- 若论文推导中使用 \(\beta=\gamma\)，该关系来自相应系数方程的纸面推导，而不是该 Lean 定理的直接结论。

这避免了“正文说两个关系，Lean table 只支持一个关系”的覆盖范围错配。

## 5. 摘要 / 引言 / 结论贡献边界

### 5.1 摘要最小修补

`thesis/main.tex` 中两处摘要句子已收紧：

1. FFT9 accuracy-cost 句子：

```latex
显示四阶 FFT9 在当前规则 Dirichlet 光滑问题与实现条件下，以相近复杂度获得显著更低误差。
```

2. near-resonance 句子：

```latex
near-resonance 扫描验证了当前 Dirichlet 五点离散中目标模态幅值遵循 \(1/|\delta|\) 放大规律，
并观察到该小分母结构与无预处理 GMRES(30) 的残差停滞相一致。
```

### 5.2 过强表述扫描

已扫描 `thesis/main.tex`、`thesis/chapters/1_introduction.tex`、`thesis/chapters/7_conclusion.tex` 和 `thesis/chapters/appendix_lean4.tex`，未发现以下风险表述：

- `PDE 理论完备`
- `PDE 理论完全`
- `FFT9 全局最优`
- `FFT9 最优`
- `所有九点格式不能六阶`
- `所有 Helmholtz 求解器`
- `所有 Krylov 方法`
- `严格证明 GMRES 停滞`
- `Lean 证明完整`
- `Lean 证明 FFT9`

当前摘要、引言和结论使用的是较稳妥的证据边界：

- 当前规则 Dirichlet 光滑问题与实现条件；
- 当前 Dirichlet 五点离散 near-resonance；
- 无预处理 GMRES(30)；
- 离散谱分析、实现校验和数值实验之间的相互印证或证据链；
- Lean 4 仅作为有限代数核心辅助验证。

## 6. 未修改范围

本轮未修改：

- `code/python/helmholtz_solver.py` 执行逻辑；
- `code/python/gmres_solver.py` 执行逻辑；
- Lean theorem 语句或证明代码；
- 实验脚本；
- CSV；
- PNG/PDF 实验图；
- HTML/PPT 讲解材料。

## 7. 验证结果

### 7.1 Python 编译

命令：

```powershell
python -m py_compile code/python/helmholtz_solver.py code/python/gmres_solver.py
```

结果：通过。

### 7.2 指定测试

命令：

```powershell
python -m pytest -q code/python/tests/test_04_fft9_vs_spsolve.py code/python/tests/test_02_nonhom_dirichlet.py code/python/tests/test_05_neumann_compatibility.py code/python/tests/test_01_eigenvalues.py
```

结果：

```text
59 passed
```

### 7.3 全量 tests

命令：

```powershell
python -m pytest -q code/python/tests
```

结果：

```text
105 passed, 2 skipped
```

### 7.4 Lean

命令：

```powershell
cd code/lean4_formalization
lake build
```

结果：

- build completed successfully；
- 仅保留既有 unused simp argument warnings；
- `proofwidgets` 依赖仓库有 local changes 提示；
- 无构建失败。

自有 Lean 文件扫描：

- 排除 `.lake` 后，未发现 `sorry`、`admit`、`axiom` 或 `unsafe`。

### 7.5 LaTeX

优先尝试：

```powershell
latexmk -xelatex -interaction=nonstopmode main.tex
```

结果：本机 MiKTeX 缺少 Perl script engine，`latexmk` 无法运行。这是本地工具链问题，不是论文错误。

随后按计划使用：

```powershell
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

结果：

- `thesis/main.pdf` 正常生成；
- 最终 PDF 页数：71 页；
- `main.log` 扫描未发现 fatal error；
- 未发现 undefined reference/citation；
- 未发现 multiply-defined label；
- 未发现 overfull hbox。

保留的非致命提示：

- SimSun font shape substitution warning；
- hyperref PDF string token warnings；
- MiKTeX update 提示；
- Lean unused simp argument warnings。

## 8. Diff / static checks

已确认：

- 没有把当前 FFT9 代码写成 global-sign-flipped 系统；
- 没有在 Lean theorem table 中补写 unsupported 的 \(\beta=\gamma\)；
- 摘要/引言/结论未发现指定的过强贡献边界表述。

## 9. 最终结论

本轮最终核查后，Q8--Q15 的口径已经进一步收紧：

- FFT9 final code convention 明确为 raw/theory；
- Lean theorem table 与实际 Lean theorem statement 保持一致；
- \(\beta=\gamma\) 不被归为该 Lean theorem 的直接验证结论；
- 摘要对 FFT9 accuracy-cost 和 near-resonance 的表述已限定到当前离散与实现条件；
- 所有验证命令均通过，除 `latexmk` 因本机缺 Perl 无法运行外，已通过手动 `xelatex/bibtex` 完整编译替代验证。
