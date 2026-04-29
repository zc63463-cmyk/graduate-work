# Complex Manufactured Solution Experiment Report

生成时间：2026-04-29

## 1. Goal

本轮重做旧版 `Complex Test Problem`，定位为“多模态 manufactured solution 的可视化补充验证”，不是新的第 6 章核心实验编号。

旧图中出现 FA Error / FFT9 Error 的 `max=1.0e+00`，这与二阶/四阶离散求解器在光滑 Dirichlet manufactured solution 上的预期误差量级不符。此类现象通常提示 RHS 构造、full-grid/interior-grid 比较、边界补齐或误差计算存在问题。本轮先完成数学和代码一致性检查，再生成图像。

## 2. Mathematical setup

区域：
\[
\Omega=(0,1)^2.
\]

统一 PDE：
\[
(-\Delta+\sigma)u=f,\qquad \sigma=10.
\]

边界条件：
\[
u=0\quad\text{on }\partial\Omega.
\]

解析解：
\[
u(x,y)=
\sin(\pi x)\sin(\pi y)
+0.5\sin(2\pi x)\sin(3\pi y)
+0.3\sin(5\pi x)\sin(4\pi y).
\]

右端项逐模态构造：
\[
\begin{aligned}
f(x,y)=&
(2\pi^2+\sigma)\sin(\pi x)\sin(\pi y)\\
&+0.5(13\pi^2+\sigma)\sin(2\pi x)\sin(3\pi y)\\
&+0.3(41\pi^2+\sigma)\sin(5\pi x)\sin(4\pi y).
\end{aligned}
\]

关键检查：

- mode `(1,1)` 系数为 `2*pi^2 + sigma`。
- mode `(2,3)` 系数为 `13*pi^2 + sigma`。
- mode `(5,4)` 系数为 `41*pi^2 + sigma`。
- 未使用“整个 `u_exact` 乘统一系数”的错误 RHS。

## 3. Implementation

新增脚本：

- `code/python/experiments/exp06_complex_manufactured_visualization.py`

复用求解器：

- `fa_helmholtz` from `code/python/helmholtz_solver.py`
- `fft9_helmholtz` from `code/python/helmholtz_solver.py`

网格一致性检查：

- `fa_helmholtz` 和 `fft9_helmholtz` 均返回 full `n x n` grid。
- 脚本显式检查 exact / FA / FFT9 shape 均为 `(n,n)`。
- 脚本显式检查 interior shape 为 `(n-2,n-2)`。
- 脚本显式检查 exact / FA / FFT9 boundary max 均为零。
- 脚本检查 `U[i,j] = u(x_i,y_j)`，并在绘图时统一使用 `data.T` 与 `origin="lower"`，避免 imshow 方向与数组索引混淆。

脚本执行顺序：

1. `n=65` sanity check；
2. `[33,65,129,257]` 收敛阶验证；
3. 只有检查通过后才写 CSV 和 PNG。

## 4. Sanity check

`n=65, sigma=10` 结果：

| Method | Linf error | L2 error | Boundary max |
|---|---:|---:|---:|
| FA | `1.932e-03` | `7.254e-04` | `0.0` |
| FFT9 | `1.320e-06` | `6.537e-07` | `0.0` |

判断：

- 未出现 `max error ~ 1`。
- 边界误差为零。
- FFT9 明显优于 FA。
- 通过 sanity check，允许生成最终图像。

## 5. Convergence study

CSV：

- `code/python/experiments/results/exp06_complex_manufactured_convergence.csv`

结果表：

| n | h | FA Linf | FA L2 | FFT9 Linf | FFT9 L2 | FA order Linf | FA order L2 | FFT9 order Linf | FFT9 order L2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 33 | 3.125e-02 | 7.566e-03 | 2.919e-03 | 2.070e-05 | 1.027e-05 | - | - | - | - |
| 65 | 1.5625e-02 | 1.932e-03 | 7.254e-04 | 1.320e-06 | 6.537e-07 | 1.970 | 2.009 | 3.972 | 3.974 |
| 129 | 7.8125e-03 | 4.828e-04 | 1.811e-04 | 8.280e-08 | 4.104e-08 | 2.000 | 2.002 | 3.994 | 3.994 |
| 257 | 3.90625e-03 | 1.207e-04 | 4.525e-05 | 5.186e-09 | 2.568e-09 | 2.000 | 2.001 | 3.997 | 3.998 |

Log-log fitted orders:

- FA Linf: `1.991`
- FA L2: `2.004`
- FFT9 Linf: `3.988`
- FFT9 L2: `3.989`

结论：

- FA shows second-order convergence.
- FFT9 shows fourth-order convergence.

## 6. Figures generated

Code-side figures:

- `code/python/experiments/figures/exp06_complex_manufactured_fields.png`
- `code/python/experiments/figures/exp06_complex_manufactured_convergence.png`

Thesis figures:

- `thesis/figures/exp06_complex_manufactured_fields.png`
- `thesis/figures/exp06_complex_manufactured_convergence.png`

Figure design:

- Field plot uses a `2x3` layout: exact field, FA solution, FFT9 solution, FA error, FFT9 error, and `y=0.5` line cut.
- exact / FA / FFT9 fields share the same `vmin/vmax`.
- FA / FFT9 error panels share the same log-error color scale.
- Error titles use `abs_err.max()`.
- Convergence plot shows Linf and L2 error curves with `O(h^2)` and `O(h^4)` references.

## 7. Thesis integration

Modified:

- `thesis/chapters/6_experiments.tex`

Integration style:

- Added under experiment three as `多模态 manufactured solution 可视化补充`.
- Did not add `实验七`.
- Did not modify the six core experiment structure.
- Did not restore old unreliable figures.

The thesis text explains that old `O(1)` error complex-test figures are not trusted because they likely reflect RHS construction, full/interior grid comparison, or boundary-completion mistakes. The new figures are used only after sanity and convergence checks pass.

## 8. Validation results

Python tests:

```text
103 passed, 2 skipped, 2 warnings
```

The two warnings are the existing Neumann Poisson compatibility warnings and are non-blocking.

Experiment script:

```text
Sanity n=65: FA Linf=1.932e-03, FA L2=7.254e-04, FFT9 Linf=1.320e-06, FFT9 L2=6.537e-07
Fitted orders: FA Linf=1.991, FA L2=2.004, FFT9 Linf=3.988, FFT9 L2=3.989
```

LaTeX:

```text
Output written on main.pdf (56 pages).
```

Log scan:

- no fatal error;
- no undefined citation/reference;
- no multiply-defined label.

Non-blocking LaTeX warnings:

- SimSun bold-italic font fallback;
- hyperref bookmark token warnings;
- MiKTeX update notice.

## 9. Final assessment

SUCCESS
