# Q8/Q9 GMRES-Krylov 终审最小修复报告

## 1. 修复目标

本轮只修正 GMRES/Krylov 相关理论表述、复杂度记号和代码注释，不新增实验、不修改 GMRES 数值算法、不重绘图、不修改 CSV。

修复点：

- Q8：区分 restarted GMRES 的 restart length 与 total iterations；
- Q9：将 GMRES 正定型收敛率估计降级为背景 remark，并明确其不解释 true Helmholtz near-resonance 下 restarted GMRES(30) 停滞。

## 2. Q8 复杂度与内存口径

统一记号：

\[
N=\text{单方向内点数},\qquad M=N^2=\text{总未知量数}.
\]

对 restarted GMRES：

\[
r=\text{restart length},\qquad k=\text{total Arnoldi iterations}.
\]

第 3.5 节复杂度表已将 GMRES 行改为：

\[
\text{cost}=O(krN^2),\qquad r\text{ fixed}:O(kN^2),
\]

\[
\text{memory}=O(rN^2).
\]

表后说明同步补充：

- 五点稀疏矩阵一次 matrix-vector product 成本为 \(O(M)=O(N^2)\)；
- matvec-only 总成本为 \(O(kN^2)\)；
- 计入 Arnoldi 正交化后总成本为 \(O(krN^2)\)；
- Krylov basis 仅存储当前 restart cycle，因此内存是 \(O(rN^2)\)，不是 \(O(kN^2)\)；
- 本文实验使用 GMRES(30)，即 \(r=30\) 固定。

启发式估计已从旧的 \(m=O(h^{-1})\) 改为：

\[
k=O(N).
\]

于是 matvec-only 口径为

\[
O(N^3)=O(M^{3/2}),
\]

计入正交化为

\[
O(rN^3),
\]

固定 \(r\) 后仍为同一阶。

## 3. Q9 正定型 GMRES 收敛率背景

第 5.4 节中原“GMRES 收敛率估计” theorem 已降级为 remark：

\[
\text{正定型 GMRES 收敛率背景}.
\]

新的表述使用对称部分：

\[
H=\frac{A+A^T}{2},\qquad \mu=\lambda_{\min}(H)>0.
\]

full GMRES 的背景残差估计写为：

\[
\frac{\|\mathbf r_m\|_2}{\|\mathbf r_0\|_2}
\le
\left(
1-\frac{\mu^2}{\|A\|_2^2}
\right)^{m/2}.
\]

对于 SPD 矩阵：

\[
H=A,\qquad
\mu=\lambda_{\min}(A),\qquad
\|A\|_2=\lambda_{\max}(A),
\]

因此该估计化为：

\[
\frac{\|\mathbf r_m\|_2}{\|\mathbf r_0\|_2}
\le
\left(
1-\frac{1}{\kappa_2(A)^2}
\right)^{m/2}.
\]

论文同时说明该界通常较保守；对于 SPD 系统，CG 的经典条件数界更尖锐。本文保留该估计仅作为 GMRES 在正定型问题上的背景说明。

## 4. 适用范围收紧

第 5.4 节已将适用范围改为：

- Dirichlet Poisson：五点矩阵为 SPD，背景界只说明网格加密会恶化 Krylov 收敛保证；
- Dirichlet modified Helmholtz：\(A=A_0+\kappa^2I\) 为 SPD，固定网格和 \(\kappa>0\) 下条件数相对 Poisson 改善；
- Neumann modified Helmholtz：\(\sigma=+\kappa^2\) 将零模态分母移至正数，在对称化变量中可按正定问题理解；
- Neumann Poisson：需要 compatibility condition，并固定零均值代表元后才能在去零模态空间讨论正定型估计；
- True Helmholtz：除 \(\kappa^2<\lambda_{\min}\) 的正定区间外，矩阵可能奇异、不定或负定，正定型背景界一般不适用。

第 5.5 节和第 6 章实验解释已补充：

> true Helmholtz GMRES 停滞不由上述正定型 bound 解释，而应从离散谱小分母、Krylov 多项式逼近困难和 restarted GMRES 残差历史共同理解。

## 5. Restart 口径修改

第 5.5 节“重启策略”已改用 \(r\)：

- GMRES(\(r\)) 每 \(r\) 步后重启；
- \(r\) 越大，单个 restart cycle 保留的 Krylov 信息越多，但内存和正交化开销也更大；
- \(r\) 太小可能导致残差下降缓慢或停滞；
- 本文实验使用 \(r=30\)。

保留 `GMRES(30)` 和 `GMRES30` 作为实验标签，含义为 restart length \(r=30\)。

## 6. 代码注释修改

`code/python/gmres_solver.py` 仅修改 docstring/comment：

- `Restart parameter m (GMRES(m))` 改为 `Restart length r (GMRES(r))`；
- `GMRES restart parameter m` 改为 `GMRES restart length r`。

未修改内部变量名和执行逻辑。代码中局部变量 `m = min(restart, N, max_iter)` 仍作为当前 restart cycle 的有效长度使用。

## 7. 未修改范围

- 未修改 GMRES 算法；
- 未引入预处理 GMRES、MINRES 或新 Krylov 理论；
- 未新增实验；
- 未运行实验脚本；
- 未重绘图；
- 未修改 CSV；
- 未改变第 6 章中 GMRES(30)/GMRES30 的实验标签。

## 8. 验证状态

静态检查：

- 已确认论文中不再保留旧的 GMRES 复杂度/内存口径：
  - `O(mN^2)`；
  - `O(m N^2)`；
  - `O(mM)`；
  - `m 为 GMRES 的总迭代次数`；
  - `GMRES(m)`；
  - `重启参数 m`；
  - `存储 m 个 Krylov`。
- 已确认第 5.4 节不再使用 \(A+A^T\) 的 \(\lambda_{\min}/\lambda_{\max}\) 混合公式；只保留 \(H=(A+A^T)/2\) 的定义。

执行验证：

- `python -m py_compile code/python/gmres_solver.py`：通过。
- `python -m pytest -q code/python/tests`：通过，`105 passed, 2 skipped`。
- LaTeX 完整编译 `xelatex -> bibtex -> xelatex -> xelatex`：通过。
- `main.log` 扫描：
  - no fatal error；
  - no undefined reference/citation；
  - no multiply-defined label；
  - no overfull hbox。
- 最终 PDF 页数：71 页。

## 9. 本轮结论

Q8/Q9 后，GMRES 复杂度、内存、restart 记号和正定型收敛率背景的适用范围已经统一。论文不再把 restart length 和 total iterations 混用，也不再将正定型 GMRES bound 用作 true Helmholtz near-resonance 停滞的解释工具。
