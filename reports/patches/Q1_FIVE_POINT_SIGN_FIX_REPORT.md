# Q1 Five-Point Sign Fix Report

## 1. 修复目标

本轮修复针对五点正 Laplacian 离散格式与非齐次 Dirichlet 边界右端修正的符号一致性。统一连续模型为
\[
(-\Delta+\sigma)u=f,
\]
五点离散算子定义为
\[
L_h^{(5)}u_{i,j}
=
\frac{4u_{i,j}-u_{i-1,j}-u_{i+1,j}-u_{i,j-1}-u_{i,j+1}}{h^2},
\]
该算子近似的是 \(-\Delta u\)，不是 \(\Delta u\)。

## 2. residual / truncation error 方向

论文第 2.2 节已统一采用
\[
\tau_h(u)=\bigl(L_h^{(5)}+\sigma I\bigr)u-f,
\]
即将连续精确解代入离散算子后减去连续方程右端。

在该约定下，五点格式的局部截断误差主项为
\[
\tau_{i,j}
=
-\frac{h^2}{12}
\left(u_{xxxx}+u_{yyyy}\right)_{i,j}
+O(h^4).
\]

## 3. 非齐次 Dirichlet RHS 修正

第 2.5.1 节已经补充靠近左边界内点 \((1,j)\) 的推导。若 \(u_{0,j}=g_{0,j}\)，则
\[
\frac{4u_{1,j}-g_{0,j}-u_{2,j}-u_{1,j-1}-u_{1,j+1}}{h^2}
+\sigma u_{1,j}
=
f_{1,j}
\]
移项后得到
\[
\frac{4u_{1,j}-u_{2,j}-u_{1,j-1}-u_{1,j+1}}{h^2}
+\sigma u_{1,j}
=
f_{1,j}+\frac{g_{0,j}}{h^2}.
\]
因此，五点正 Laplacian 离散下非齐次 Dirichlet 边界值移至右端项后取正号。

## 4. 代码检查与修复范围

已确认活跃五点求解路径中的 Dirichlet RHS 边界修正使用正号：

- `code/python/cyclic_reduction.py`
- `code/python/helmholtz_solver.py`
- `code/python/gmres_solver.py`
- `code/python/experiments/exp00_fft_vs_sparse.py`
- `code/python/tests/test_03_fft_vs_spsolve.py`

本轮只清理了 `code/python/fft9_complete.py` 中 deprecated prototype 的两个五点 RHS 旧负号：

- `spectral` 分支中五点 Dirichlet RHS 修正由 `-=` 改为 `+=`；
- `poisson_5point_fft` helper 中五点 Dirichlet RHS 修正由 `-=` 改为 `+=`。

未修改以下内容：

- FFT9 紧致格式边界修正，因为其符号需要按 \((-L_h+\sigma R_h)u=R_h f\) 单独判断；
- Neumann ghost-point 相关代码；
- `code/python/fft9_eigenvalue_analysis.py` 中的 FFT9 边界修正；
- `code/python/_archive/*` 历史归档文件。

## 5. 验收口径

验收时应确认：

- 第 2 章中五点 residual 方向为 \(\tau_h=(L_h^{(5)}+\sigma I)u-f\)；
- 在该方向下，五点截断误差主项为负号；
- 非齐次 Dirichlet 五点 RHS 边界贡献为 \(+g/h^2\)；
- 活跃五点代码路径不再含有错误的 Dirichlet RHS `-=` 修正。
