# Q8--Q15 最近一轮问题修复总报告

## 1. 总体结论

本轮 Q8--Q15 主要是“终审口径收紧”而不是数值实验重做。修改集中在 GMRES/Krylov 理论表述、实验章证据边界、Lean 4 形式化覆盖范围、以及 FFT9 raw/theory 符号系统与代码实现的一致性说明。

本轮没有新增实验，没有修改实验参数，没有重绘图，没有修改 CSV/PNG/PDF 实验资产，也没有改变 FFT9、GMRES、Neumann 或 Helmholtz 求解算法。

当前总体状态：

- GMRES restart 复杂度符号已经统一；
- 正定型 GMRES 收敛率估计已经从 theorem 降级为背景 remark；
- exp06、exp07、exp05 的证据边界已经收紧；
- Lean 4 已限定为有限代数核心验证；
- FFT9 当前实现明确为 raw/theory 系统，不误写成代码反号系统；
- `thesis/main.pdf` 已正常生成，最终 71 页；
- LaTeX log 扫描未发现 fatal error、undefined reference/citation、multiply-defined label、overfull hbox。

## 2. 分项修复概览

| 编号 | 主题 | 原风险 | 修复结果 | 是否改算法/实验 |
|---|---|---|---|---|
| Q8 | restarted GMRES 复杂度与内存 | restart length 与 total iterations 混用，内存口径可能写成错误的 \(O(mN)\) 或旧 \(m\) 记号 | 统一为 restart length \(r\)、total Arnoldi iterations \(k\)、单方向内点 \(N\)、未知量 \(M=N^2\)；成本 \(O(krN^2)\)，固定 \(r\) 时 \(O(kN^2)\)，内存 \(O(rN^2)\) | 否 |
| Q9 | GMRES 收敛率理论边界 | 将正定型 GMRES 背景界写得像可解释 true Helmholtz near-resonance 停滞的定理 | 第 5.4 节降级为正定型背景 remark；使用 \(H=(A+A^T)/2\)、\(\mu=\lambda_{\min}(H)>0\)；明确不解释 true Helmholtz near-resonance 下 GMRES(30) 停滞 | 否 |
| Q10 | exp06 accuracy-cost 证据边界 | 容易被读成 FFT9 普遍优于 GMRES 或严格 kernel benchmark | 改为当前规则 Dirichlet 光滑 manufactured solution、当前实现和当前 timing scope 下的 error-time placement | 否 |
| Q11 | 图 10 risk map 口径 | risk map 容易被误读为条件数图 | 明确 \(R_{p,q}=-\log_{10}|d_{p,q}|\) 只是逐模态 small-denominator risk visualization；只有 Dirichlet 五点正交 DST-I 系统中 \(\max|d|/\min|d|\) 才与 \(\mathrm{cond}_2(A)\) 一致 | 否 |
| Q12 | near-resonance 与 GMRES 停滞因果 | “小分母导致 GMRES 停滞”表述过强 | 保留模态放大结论；GMRES 停滞改为“与小分母病态结构相一致，并有助于解释当前无预处理 GMRES(30) 观察” | 否 |
| Q13 | Lean 4 覆盖范围 | Lean 可能被理解为证明完整 PDE/全局误差/程序正确性 | 附录明确 Lean 只辅助验证三参数九点修正族六阶 local consistency 反证中的有限多项式代数核心 | 否 |
| Q14 | 摘要、引言、结论贡献边界 | “理论--实验闭环”可能暗示完整理论已由实验或 Lean 证明 | 改为“离散谱分析--实现校验--数值实验”的证据链或相互印证 | 否 |
| Q15 | FFT9 符号系统与代码事实 | 附件曾假设当前代码可能为整体反号系统，若照写会与仓库事实冲突 | 保持当前事实：代码和论文均采用 raw/theory 系统；只补充“若其他实现整体反号，RHS 与 denominator 必须同步反号，不能混用” | 否 |

## 3. Q8--Q9：GMRES/Krylov 口径修复

### 3.1 修复目标

Q8/Q9 修复的是 GMRES 理论和复杂度说明，不改变 GMRES 执行逻辑。

核心统一：

\[
N=\text{单方向内点数},\qquad M=N^2=\text{总未知量},
\]

\[
r=\text{restart length},\qquad k=\text{total Arnoldi iterations}.
\]

### 3.2 复杂度与内存

第 3.5 节复杂度表中，GMRES 行已改为：

\[
\text{cost}=O(krN^2),\qquad r \text{ fixed}: O(kN^2),
\]

\[
\text{memory}=O(rN^2).
\]

解释：

- 五点稀疏矩阵一次 matvec 成本是 \(O(N^2)=O(M)\)；
- restarted GMRES 每个 restart cycle 只保留当前 \(r\) 个 Krylov 向量；
- 因此内存是 \(O(rM)=O(rN^2)\)，不是 \(O(kN^2)\)，也不是 \(O(mN)\)；
- 若启发式假设总迭代数 \(k=O(N)\)，则 matvec-only 口径为 \(O(N^3)=O(M^{3/2})\)，计入正交化为 \(O(rN^3)\)，固定 \(r\) 后同阶。

### 3.3 正定型 GMRES 收敛率背景

第 5.4 节原“GMRES 收敛率估计”已经降级为背景 remark。采用的稳妥形式是：

\[
H=\frac{A+A^T}{2},\qquad \mu=\lambda_{\min}(H)>0,
\]

\[
\frac{\|\mathbf r_m\|_2}{\|\mathbf r_0\|_2}
\le
\left(1-\frac{\mu^2}{\|A\|_2^2}\right)^{m/2}.
\]

SPD 特例下：

\[
H=A,\qquad \mu=\lambda_{\min}(A),\qquad \|A\|_2=\lambda_{\max}(A),
\]

因此得到保守背景界：

\[
\frac{\|\mathbf r_m\|_2}{\|\mathbf r_0\|_2}
\le
\left(1-\frac{1}{\kappa_2(A)^2}\right)^{m/2}.
\]

论文同时说明：对 SPD 系统，CG 的经典条件数界通常更尖锐；该 GMRES 背景界不用于解释 true Helmholtz near-resonance 下 restarted GMRES(30) 的残差停滞。

### 3.4 代码层面

`code/python/gmres_solver.py` 只修改 docstring/comment，将 `GMRES(m)` 的 restart parameter 说明改为 restart length \(r\)。内部变量名和执行逻辑未改。

## 4. Q10--Q12：实验章证据边界收紧

### 4.1 Q10：exp06 accuracy-cost

第 6.4 节已经明确 exp06 是：

- 当前规则 Dirichlet 光滑 manufactured solution 上的核心基准证据；
- 当前实现和当前 timing scope 下的 error-time placement；
- 不是同一线性系统上的严格 kernel-to-kernel solver benchmark。

当前口径：

- FA/CR/FACR-like 和 GMRES30 求解五点二阶离散系统；
- FFT9 求解 Dirichlet 九点紧致四阶离散系统；
- FFT9 的更低误差来自高阶离散格式，不是把同一个五点系统“解得更精确”；
- GMRES30 是无预处理 restarted GMRES(30) 基线，不代表预处理 GMRES、MINRES 或所有 Krylov 方法。

### 4.2 Q11：图 10 small-denominator risk map

图 10 的解释已经收紧为：

\[
R_{p,q}=-\log_{10}|d_{p,q}|
\]

仅用于逐模态 small-denominator risk visualization。

明确限制：

- \(R_{p,q}<0\) 只表示 \(|d_{p,q}|>1\)，即分母较大、小分母风险较低；
- true-away 的 sign-changing 不等于 near-resonance；
- 图 10 不是解图、误差图，也不是条件数图；
- 只有当前 Dirichlet 五点正交 DST-I 系统中，\(\max|d|/\min|d|\) 才与 \(\mathrm{cond}_2(A)\) 一致；
- 条件数等价关系应结合后续 condition check 图理解，不能外推到 Neumann/mixed 或 non-normal 系统。

### 4.3 Q12：near-resonance 与 GMRES 停滞

实验五保留强支撑结论：

\[
\kappa^2=\lambda_{\mathrm{target}}^h+\delta
\quad\Longrightarrow\quad
d_{\mathrm{target}}=-\delta,
\]

若 RHS 在目标模态上的投影非零，则

\[
|\hat u_{\mathrm{target}}|
=
\frac{|\hat f_{\mathrm{target}}|}{|\delta|},
\]

所以目标模态按 \(1/|\delta|\) 放大。

GMRES 相关表述已经软化为：

- 当前无预处理 GMRES(30) 的残差停滞与该小分母病态结构相一致；
- 该结构有助于解释当前实验中观察到的停滞；
- 不说“小分母单独证明所有 GMRES/所有 Helmholtz 求解器都会停滞”；
- 不外推到预处理 Krylov、MINRES、shifted Laplacian 或连续谱极限。

## 5. Q13--Q15：Lean、贡献边界与 FFT9 符号系统

### 5.1 Q13：Lean 4 范围

Lean 附录现在明确：

- Lean 4 是辅助验证；
- 验证对象是三参数九点修正族六阶 local consistency 不可达性反证中的有限多项式系数条件；
- 具体覆盖 \(c_2=0\) 的参数必要关系、Poisson/Helmholtz 情形的代数路径、以及两个 Fourier 测试模式下 \(c_4=0\) 的矛盾；
- \((1,0)\)、\((1,1)\) 是 local Fourier symbol 的连续波数测试点，不是 Dirichlet DST-I 离散模态编号。

明确未覆盖：

- Fourier symbol 从模板导出的全过程；
- Taylor 展开建模；
- DST/DCT 谱理论；
- 边界闭合；
- 完整 PDE 全局误差；
- Python 求解器、实验脚本或数值实验结果的程序正确性。

### 5.2 Q14：贡献边界

摘要、引言和第六章总结中，较强的“理论--实验闭环”已经收紧为：

- 当前离散谱分析、实现校验和数值实验之间相互印证；
- “离散谱分析--实现校验--数值实验”的证据链。

这样保留了论文结构主线，但避免暗示实验或 Lean 已经完成完整 PDE 理论证明。

### 5.3 Q15：FFT9 raw/theory 系统

当前仓库事实：

- `fft9_helmholtz` 使用
  \[
  \tilde g=(R_h f)_I+L_{IB}g_B-\sigma R_{IB}g_B,
  \]
  并配合
  \[
  D^{raw}_{p,q}=-\hat\lambda_L+\sigma\hat\lambda_R.
  \]
- `test_04_fft9_vs_spsolve.py` 的 sparse compact system 对照也使用同一 raw 口径。

因此本轮没有把论文改成“当前代码整体反号系统”。第 3/4 章补充的是防误读说明：

- 整体乘以 \(-1\) 是等价代数写法；
- 若其他实现采用反号系统，RHS 与 denominator 必须同步反号；
- 不能只反右端项或只反频域分母；
- 本文当前理论推导和 Python 实现均采用 raw 系统。

## 6. 修改文件汇总

Q8--Q15 涉及的主要文件：

- `thesis/main.tex`
- `thesis/chapters/1_introduction.tex`
- `thesis/chapters/3_fft_direct.tex`
- `thesis/chapters/5_gmres.tex`
- `thesis/chapters/6_experiments.tex`
- `thesis/chapters/7_conclusion.tex`
- `thesis/chapters/appendix_lean4.tex`
- `code/python/gmres_solver.py`
- `Q8_Q9_GMRES_KRYLOV_PATCH_REPORT.md`
- `Q10_Q12_EXPERIMENT_EVIDENCE_PATCH_REPORT.md`
- `Q13_Q15_LEAN_CONCLUSION_CODE_CONSISTENCY_PATCH_REPORT.md`
- `Q8_Q15_RECENT_FIX_SUMMARY_REPORT.md`

注意：

- `code/python/helmholtz_solver.py` 当前已有 raw 系统注释与实现，但 Q13--Q15 这轮未改变其执行逻辑；
- 未修改实验脚本、CSV、PNG/PDF 实验图或 HTML/PPT 讲解材料。

## 7. 验证状态汇总

### 7.1 Q8--Q9 分项验证

分项报告记录：

- `python -m py_compile code/python/gmres_solver.py`：通过；
- `python -m pytest -q code/python/tests`：通过，`105 passed, 2 skipped`；
- LaTeX 完整编译通过；
- log 扫描无 fatal error、undefined reference/citation、multiply-defined label、overfull hbox；
- `main.pdf` 页数为 71 页。

### 7.2 Q10--Q12 分项验证

分项报告记录：

- 静态搜索确认旧风险表述已处理；
- LaTeX 完整编译通过；
- log 扫描无 fatal error、undefined reference/citation、multiply-defined label、overfull hbox；
- `main.pdf` 页数为 71 页。

### 7.3 Q13--Q15 分项验证

本轮已执行：

```powershell
cd code/lean4_formalization
lake build
```

结果：

- build completed successfully；
- 仅有既有 unused simp argument warnings；
- 有 `proofwidgets` 依赖仓库 local changes 提示；
- 无构建失败。

自有 Lean 文件检查：

- 未发现 `sorry`、`admit`、`axiom`、`unsafe`；
- `.lake` 依赖目录不纳入自有文件逃逸检查。

LaTeX：

- `xelatex -> bibtex -> xelatex -> xelatex` 通过；
- `thesis/main.pdf` 正常生成，最终 71 页；
- 当前 `main.log` 扫描仅见 `Output written on main.pdf (71 pages)`；
- 未发现 fatal error、undefined reference/citation、multiply-defined label、overfull hbox。

### 7.4 当前保留的非致命提示

- Lean：unused simp argument warnings；
- Lean：`proofwidgets` 依赖仓库 local changes 提示；
- LaTeX：SimSun font shape substitution warning；
- LaTeX：hyperref PDF string token warnings；
- MiKTeX：更新提示。

这些提示不影响本轮修复目标。

## 8. 给 GPT5.5 / 外部审阅者的重点审查问题

建议把以下问题交给 GPT5.5 或导师复审：

1. **GMRES 复杂度口径是否严谨？**
   在 \(N\) 为单方向内点数、\(M=N^2\)、restart length 为 \(r\)、total Arnoldi iterations 为 \(k\) 的约定下，论文写 \(O(krN^2)\)、固定 \(r\) 时 \(O(kN^2)\)、内存 \(O(rN^2)\) 是否合理？

2. **正定型 GMRES 背景界是否放在合适位置？**
   第 5.4 节将该界写成 remark，并使用 \(H=(A+A^T)/2\)、\(\mu=\lambda_{\min}(H)>0\)。该表述是否足够保守？是否清楚说明它不解释 true Helmholtz near-resonance 下 restarted GMRES(30) 停滞？

3. **exp06 的证据边界是否恰当？**
   现在论文只声称当前规则 Dirichlet 光滑 manufactured solution、当前实现和当前 timing scope 下 FFT9 有更好的 error-time placement。这个表述是否既不外推，又不削弱论文标题中的“对比研究”？

4. **图 10 risk map 与 condition check 是否分清？**
   图 10 是否已经充分说明自己是 small-denominator risk visualization，而不是 condition number heatmap？图 12/condition check 是否只限定在 Dirichlet 五点正交 DST-I 系统？

5. **near-resonance 与 GMRES 停滞的因果口径是否安全？**
   当前表述为“小分母模态放大和病态性与无预处理 GMRES(30) 残差停滞相一致，并有助于解释该现象”。这是否比“导致/造成 GMRES 停滞”更符合证据强度？

6. **Lean 4 覆盖范围是否仍有过强暗示？**
   附录是否已经明确 Lean 只验证有限多项式代数核心，不证明完整 PDE 全局误差、DST/DCT 谱理论、边界闭合或 Python 程序正确性？

7. **FFT9 raw/theory 系统是否与代码一致？**
   当前论文写
   \[
   \tilde g=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B,\qquad
   D^{raw}=-\hat\lambda_L+\sigma\hat\lambda_R.
   \]
   这是否与 `helmholtz_solver.py` 和 `test_04_fft9_vs_spsolve.py` 的实现一致？反号系统说明是否足够防止混用？

## 9. 仍需注意的边界

即使 Q8--Q15 已修复，论文最终口径仍应持续保持以下限制：

- GMRES 结论只针对无预处理 restarted GMRES(30)；
- FFT9 当前实现和实验验证只针对 Dirichlet 边界；
- FACR-like 当前实现仍按 \(O(N^2\log N)\) 处理；
- timing scope 不是严格 kernel benchmark；
- near-resonance 是离散谱实验，不是连续谱极限结论；
- Lean 4 是代数辅助验证，不是完整 PDE/程序正确性形式化证明；
- FFT9 raw 和整体反号系统只能整体切换，不能混用 RHS 与 denominator。

## 10. 结论

Q8--Q15 修复后，论文在理论表述、实验结论、形式化证明范围和代码符号系统四个方面都更稳健。最重要的变化不是新增结论，而是把已有结论放回其真实证据范围内：

- GMRES 复杂度与 restart 记号清楚；
- 正定型 GMRES 背景界不再解释不定 near-resonance 问题；
- exp06/exp07/exp05 的图表解释不再超出当前数据；
- Lean 4 不再被误读为完整 PDE 形式化；
- FFT9 符号系统与当前 Python 实现保持一致。

这份报告可作为提交 GPT5.5 或导师复核 Q8--Q15 修复的总索引。
