# Q4 FFT9 非齐次 Dirichlet 符号系统说明报告

## 1. 修复目标

本轮处理 FFT9 非齐次 Dirichlet 边界修正项的符号说明问题。关键事实是：当前 `helmholtz_solver.py` 已经采用论文理论系统

\[
(-L_h+\sigma R_h)u=R_hf,
\]

而不是整体乘以 \(-1\) 的代码反号系统。因此本轮不修改数值算法，只补充论文与代码注释，明确当前实现与论文公式一致，并说明若其他实现采用整体反号系统，RHS 与 denominator 必须同步反号。

## 2. 当前实现采用的符号系统

论文和当前代码均采用 raw/theory 系统：

\[
\tilde g
=
(R_hf)_I
+L_{IB}g_B
-\sigma R_{IB}g_B,
\]

\[
D_{p,q}^{\mathrm{raw}}
=
-\hat{\lambda}_L(p,q)
+\sigma\hat{\lambda}_R(p,q).
\]

代码中的对应关系为：

- `compute_bc_correction_9pt` 返回 \(L_{IB}g_B-\sigma R_{IB}g_B\)；
- `fft9_helmholtz` 构造 `rhs_tilde = apply_Rh_full(F) + bc_corr`；
- `fft9_helmholtz` 使用 `denom = -lam_L + sigma * lam_R4`；
- `fft9_oer_helmholtz` 使用同一 RHS，并在块系数中采用 \(-L_h+\sigma R_h\)。

因此当前代码与论文公式一致，不存在 RHS 与 denominator 只反一个的混用。

## 3. 等价反号系统的说明

若某个实现将整个线性系统乘以 \(-1\)，则会得到等价系统

\[
(L_h-\sigma R_h)u=-R_hf.
\]

此时必须同时使用

\[
\tilde g_{\mathrm{code}}
=
-(R_hf)_I
-L_{IB}g_B
+\sigma R_{IB}g_B,
\]

以及

\[
D_{p,q}^{\mathrm{code}}
=
\hat{\lambda}_L(p,q)-\sigma\hat{\lambda}_R(p,q)
=
-D_{p,q}^{\mathrm{raw}}.
\]

由于分子和分母同时反号，两种写法给出的 \(\hat u_{p,q}\) 完全相同。错误情况是只反 RHS 或只反 denominator。本轮在论文和代码注释中加入该 remark，目的是防止后续审查时把等价反号写法误判为当前实现事实。

## 4. 修改内容

- `thesis/chapters/3_fft_direct.tex`
  - 在物理空间边界修正策略后补充：当前代码采用与论文一致的 raw 符号系统。
  - 增加等价反号系统公式，并强调 RHS 与 denominator 必须同步反号。
  - 在 FFT9 算法后增加实现符号系统 remark。

- `thesis/chapters/4_helmholtz.tex`
  - 补充当前 FFT9 实现同样使用 \(D_{p,q}^{\mathrm{raw}}\)。
  - 说明整体反号分母只作为代数等价形式，不能与 raw RHS 混用。

- `code/python/helmholtz_solver.py`
  - 增强 `compute_bc_correction_9pt` docstring。
  - 在 `fft9_helmholtz` 与 `fft9_oer_helmholtz` 的 RHS 构造处补充注释。
  - 未修改任何计算语句。

## 5. 未修改范围

- 未修改 FFT9 求解器的数值算法；
- 未改变 RHS/denominator 的实际符号；
- 未修改 sparse 对照测试矩阵；
- 未修改实验脚本、CSV、PNG/PDF 图像；
- 未新增实验。

## 6. 验证状态

已完成以下检查：

- `python -m py_compile code/python/helmholtz_solver.py`：通过。
- `cd code/python; python -m pytest -q tests/test_04_fft9_vs_spsolve.py tests/test_02_nonhom_dirichlet.py`：`21 passed in 1.30s`。
- 轻量非齐次 Dirichlet FFT9 manufactured-solution 收敛检查：通过，误差随网格加密呈四阶趋势。
  - \(n=17\)：\(1.085\times 10^{-5}\)
  - \(n=33\)：\(4.477\times 10^{-7}\)，阶数约 \(4.60\)
  - \(n=65\)：\(2.445\times 10^{-8}\)，阶数约 \(4.19\)
  - \(n=129\)：\(1.474\times 10^{-9}\)，阶数约 \(4.05\)
- LaTeX 完整编译链 `xelatex -> bibtex -> xelatex -> xelatex`：通过，`main.pdf` 正常生成。
- `main.log` 扫描：无 fatal error、undefined reference、multiply-defined label、overfull hbox；最终 PDF 为 68 页。
- `git diff --check`：通过；PowerShell 输出中仅有既有 CRLF 行尾提示。

## 7. 结论

Q4 不需要改变 FFT9 的实际计算逻辑。当前最终实现已经与论文 raw/theory 系统一致：

\[
\tilde g=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\qquad
D_{p,q}^{\mathrm{raw}}=-\hat{\lambda}_L+\sigma\hat{\lambda}_R.
\]

本轮修复的重点是把该事实在论文、代码注释和审查报告中写清楚，并明确等价反号系统只能作为整体同步反号的代数说明，不能作为本文当前实现口径。
