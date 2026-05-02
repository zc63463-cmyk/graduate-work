# Q2 FFT9 Raw Denominator Fix Report

## 1. 修复目标

本轮修复统一 FFT9 九点紧致格式的符号口径。论文和最终求解器均采用
\[
(-L_h+\sigma R_h)u=R_h f,
\qquad
L_h\approx \Delta,
\]
并以原始频域分母
\[
D^{\mathrm{raw}}_{p,q}
=
-\hat{\lambda}_L(p,q)+\sigma\hat{\lambda}_R(p,q)
\]
作为主表述。

## 2. raw denominator 与 effective denominator

若物理空间中已经构造
\[
\tilde g=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\]
并对 \(\tilde g\) 做 DST-I，则频域求解必须写为
\[
\hat u_{p,q}
=
\frac{\hat{\tilde g}_{p,q}}
{D^{\mathrm{raw}}_{p,q}}.
\]

等价的有效分母
\[
D^{(9)}_{p,q}
=
-\hat{\lambda}_L(p,q)/\hat{\lambda}_R(p,q)+\sigma
\]
只直接匹配原始右端项的频域系数 \(\hat f_{p,q}\)。若右端项已经包含 \(R_h f\)
或非齐次 Dirichlet 边界修正，则不能再用 \(\hat{\tilde g}/D^{(9)}\)。

## 3. 论文修改

已增强以下位置的 FFT9 符号说明：

- 第 2 章：补充 \(L_h\approx\Delta\)、\(\hat{\lambda}_L<0\)、\(-\hat{\lambda}_L>0\)、\(\hat{\lambda}_R>0\) 的说明；
- 第 3 章：将频域推导改为 raw denominator 优先，区分 \(\hat f\) 与 \(\hat{\tilde g}\)，并删除“齐次 Dirichlet 时 \(R_hf\) 仅涉及内节点值”的过强表述；
- 第 4 章：并列区分五点分母 \(d_{p,q}=\lambda^{(5)}_{p,q}+\sigma\) 与 FFT9 原始分母 \(D^{\mathrm{raw}}_{p,q}\)；
- 第 4 章：补充 FFT9 true Helmholtz 小分母条件
  \[
  -\hat{\lambda}_L(p,q)\approx \kappa^2\hat{\lambda}_R(p,q),
  \]
  等价于 \(-\hat{\lambda}_L/\hat{\lambda}_R\approx\kappa^2\)。

## 4. 代码修改

`code/python/helmholtz_solver.py` 中 FFT9 主实现已改为 raw convention：

- `compute_bc_correction_9pt` 返回论文物理空间修正
  \(L_{IB}g_B-\sigma R_{IB}g_B\)；
- `fft9_helmholtz` 构造 `rhs_tilde = apply_Rh_full(F) + bc_corr`；
- 频域分母使用 `denom = -lam_L + sigma * lam_R4`；
- `fft9_oer_helmholtz` 的块系数同步改为 \((-L_h+\sigma R_h)\) 对应符号。

`code/python/tests/test_04_fft9_vs_spsolve.py` 已同步改为构造 sparse 系统
\[
(-L_h+\sigma R_h)u=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B.
\]

## 5. 未修改范围

本轮未修改：

- 五点 FA/CR/FACR-like 求解器；
- GMRES 求解器；
- Neumann/mixed 边界处理；
- `code/python/_archive/*` 历史归档代码；
- 第 6 章 exp05 near-resonance 的核心结论。该实验仍限定为 Dirichlet 五点离散谱实验，不混入 FFT9 共振判据。

## 6. 验证状态

已完成：

- `python -m py_compile code\python\helmholtz_solver.py code\python\verify_fft9_expansion.py`；
- `python code/python/verify_fft9_expansion.py`，确认 \(-\lambda_L/\lambda_R\) 的 \(h^2\) 项消去；
- `python -m pytest -q tests/test_04_fft9_vs_spsolve.py tests/test_02_nonhom_dirichlet.py`，结果为 `21 passed`；
- `python -m pytest -q`，结果为 `103 passed, 2 skipped`，仅保留既有 Neumann compatibility warnings。
- LaTeX 完整编译完成，`main.pdf` 正常生成，最终页数为 68 页；
- `main.log` 扫描未发现 fatal error、undefined reference、multiply-defined label 或 overfull hbox，仅保留既有 SimSun 字体替代 warning。
