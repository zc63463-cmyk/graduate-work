# Q1--Q7 终审最小补丁报告

## 1. 本轮目标

本轮执行终审最小补丁，只处理 Q1--Q7 复核后仍需要收紧的表述和记号问题：

- 按当前仓库事实确认 FFT9 代码系统；
- 修复 Neumann 二维矩阵右乘非对称矩阵的转置记号；
- 将 FFT9 展开验证脚本的输出措辞收紧为 interior Fourier symbol consistency；
- 不新增实验，不改算法，不改实验参数，不重绘图，不修改 CSV。

## 2. FFT9 符号系统检查结论

当前 `code/python/helmholtz_solver.py` 中 `fft9_helmholtz` 主路径采用 raw/theory 系统：

```python
rhs_tilde = apply_Rh_full(F)
rhs_tilde += bc_corr
denom = -lam_L + sigma * lam_R4
U_hat[mask] = rhs_hat[mask] / denom[mask]
```

当前 `code/python/tests/test_04_fft9_vs_spsolve.py` 也按同一 raw 系统构造 sparse compact system：

```text
(-L_h + sigma * R_h) * u_int = (R_h * f) + bc_correction
```

因此本轮选择“情形 B：代码已经是 raw 系统”。未采用“理论 raw / 代码整体反号等价”的论文替换方案；论文中“当前代码采用 raw 系统”的事实性表述保留。

仍保留的安全说明是：若其他实现整体乘以 \(-1\)，则 RHS 与 denominator 必须同步反号，不能混用。

## 3. Neumann 二维矩阵记号修复

已修改：

- `thesis/chapters/2_math_preliminary.tex`
- `thesis/chapters/4_helmholtz.tex`

原先二维 ghost-point 写法中右乘项写作

\[
G_xU+UG_y+\sigma U=F.
\]

由于 ghost-point \(G_y\) 原始矩阵是非对称的，右乘方向需要显式写成转置：

\[
G_xU+UG_y^T+\sigma U=F.
\]

变量缩放仍为

\[
\tilde U=D_x^{-1}UD_y^{-1},
\qquad
\tilde F=D_x^{-1}FD_y^{-1}.
\]

若

\[
S_x=D_x^{-1}G_xD_x,\qquad
S_y=D_y^{-1}G_yD_y,
\]

则

\[
D_yG_y^TD_y^{-1}
=
(D_y^{-1}G_yD_y)^T
=
S_y^T
=
S_y.
\]

因此对称化系统保持为

\[
S_x\tilde U+\tilde U S_y+\sigma\tilde U=\tilde F.
\]

这一修复只改变论文记号，不改变 Neumann/DCT-I 求解算法。

## 4. FFT9 展开脚本措辞修复

已修改：

- `code/python/verify_fft9_expansion.py`

结论输出从

```text
FFT9 has 4th-order interior Fourier symbol consistency.
```

收紧为

```text
FFT9 has fourth-order interior Fourier symbol consistency.
```

该修改只影响输出文字，不改变任何符号展开或计算逻辑。

## 5. GPT-5.5 总报告 addendum

已更新：

- `Q1_Q7_GPT55_REVIEW_REPORT.md`

新增终审 addendum，说明当前仓库检查结果为 FFT9 raw 系统，因此 Q2/Q4 的审阅口径应以“理论与当前代码均为 raw；整体反号只作为其他实现的等价说明”为准。

## 6. 未修改范围

- 未修改 `fft9_helmholtz` 数值算法。
- 未修改 `test_04_fft9_vs_spsolve.py`。
- 未新增实验。
- 未运行 exp04/exp07 等实验脚本。
- 未重绘 PNG/PDF 图。
- 未修改 CSV 数据。
- 未把 FFT9 扩展到 Neumann/mixed。
- 未改变 Q1--Q7 已完成的理论主线。

## 7. 验证记录

执行的验证命令和结果：

```text
python -m py_compile code/python/helmholtz_solver.py
python -m py_compile code/python/verify_fft9_expansion.py
python code/python/verify_fft9_expansion.py
```

结果：通过。

```text
python -m pytest -q code/python/tests/test_04_fft9_vs_spsolve.py
python -m pytest -q code/python/tests/test_02_nonhom_dirichlet.py
python -m pytest -q code/python/tests/test_05_neumann_compatibility.py
python -m pytest -q code/python/tests/test_01_eigenvalues.py
```

结果：通过。

LaTeX 完整编译链：

```text
xelatex -> bibtex -> xelatex -> xelatex
```

结果：通过，`thesis/main.pdf` 正常生成。

LaTeX log 扫描结果：

- no fatal error；
- no undefined reference/citation；
- no multiply-defined label；
- no overfull hbox；
- 保留既有 SimSun 字体形状 warning。

最终 PDF 页数：70 页。

## 8. 本轮结论

终审最小补丁完成后，Q1--Q7 的剩余 blocking issue 已按当前仓库事实处理：FFT9 当前代码与 sparse 对照测试均采用 raw/theory 系统，不需要将论文改成“代码反号系统”；Neumann 二维矩阵记号已消除右乘非对称矩阵的转置歧义；FFT9 展开脚本不再使用可能被误读为完整精度证明的输出措辞。
