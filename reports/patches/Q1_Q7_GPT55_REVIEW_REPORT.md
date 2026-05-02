# Q1--Q7 理论与实验正确性修复总报告（给 GPT-5.5 审阅）

本文档用于交给 GPT-5.5 / GPT-5.5 Pro 做外部复核。它汇总了 Q1--Q7 这一轮针对论文理论、代码符号、实验解释和结论边界的修复内容，并给出逐项审阅 prompt。

> 终审 addendum：后续最小补丁再次检查当前仓库后，确认 `fft9_helmholtz` 与 `test_04_fft9_vs_spsolve.py`
> 均采用 FFT9 raw/theory 系统，即 `rhs_tilde=(R_hf)_I+bc_corr` 且
> `denom=-lambda_L+sigma*lambda_R`。因此本文档中 Q2/Q4 关于“当前代码 raw”的口径仍然成立；
> 整体反号系统只作为其他实现可采用的代数等价说明，不能理解为当前仓库实现。

建议审阅重点不是语言润色，而是：

1. 公式符号是否一致；
2. 推导结论是否成立；
3. 论文是否仍有过度声称；
4. 代码实现与论文公式是否同向；
5. 实验解释是否超过已有证据；
6. Lean 4 覆盖范围是否表述过强。

---

## 0. 统一背景

论文研究矩形区域上的 Poisson、modified Helmholtz 与 true Helmholtz 方程，统一模型为

\[
(-\Delta+\sigma)u=f.
\]

其中：

- Poisson: \(\sigma=0\)；
- modified Helmholtz: \(\sigma=+\kappa^2\)；
- true Helmholtz: \(\sigma=-\kappa^2\)。

论文涉及的主要数值对象：

- Dirichlet 五点二阶差分与 DST-I 对角化；
- Neumann ghost-point + DCT-I 对称化；
- Dirichlet 九点紧致 FFT9；
- FA / CR / FACR-like FFT 直接法；
- 无预处理 restarted GMRES(30)；
- true Helmholtz near-resonance 离散谱实验。

本轮 Q1--Q7 的共同原则：

- 不新增实验；
- 不扩大算法范围；
- 不把局部 Fourier symbol 结果写成完整全局误差定理；
- 不把 Dirichlet 五点正交对角化系统的条件数结论外推到 Neumann/mixed/non-normal 系统；
- 不把无预处理 GMRES(30) 的现象外推到所有 Krylov 方法。

---

## 1. 修改文件总览

### 论文 LaTeX

- `thesis/chapters/1_introduction.tex`
  - 收紧三参数九点修正族六阶不可达性的贡献表述。

- `thesis/chapters/2_math_preliminary.tex`
  - 统一五点格式 residual/truncation error 方向；
  - 补强非齐次 Dirichlet RHS 正号；
  - 收紧 FFT9 四阶 consistency 定理；
  - 收紧三参数九点修正族六阶不可达性定理；
  - 补强 Neumann ghost-point、DCT-I 缩放、zero mode；
  - 收紧 modified/true Helmholtz 条件数与谱解释。

- `thesis/chapters/3_fft_direct.tex`
  - 统一 FFT9 raw denominator；
  - 明确 \((R_h f)_I+L_{IB}g_B-\sigma R_{IB}g_B\) 与 \(D^{raw}\) 配对；
  - 收紧 FFT9 四阶为 interior Fourier symbol consistency + 实验验证整体收敛；
  - 说明整体反号系统只作为等价代数形式。

- `thesis/chapters/4_helmholtz.tex`
  - 区分五点分母与 FFT9 raw denominator；
  - 收紧 modified Helmholtz SPD 与条件数改善范围；
  - true Helmholtz 谱分类补全正定/奇异/不定/负定；
  - 补强 Neumann compatibility 与 zero mode；
  - 防止把 modified Helmholtz 与 true Helmholtz 的谱机制混淆。

- `thesis/chapters/6_experiments.tex`
  - 收紧 FFT9 四阶实验结论口径；
  - 收紧 risk map、condition check、near-resonance、GMRES 解释；
  - 明确 spectral denominator indicator 与 \(\mathrm{cond}_2(A)\) 的等价条件。

- `thesis/chapters/7_conclusion.tex`
  - 同步收紧 FFT9、Lean、modified/true Helmholtz、GMRES 与实验边界。

- `thesis/chapters/appendix_lean4.tex`
  - 将 Lean 4 角色限定为三参数九点修正族六阶不可达性反证中的代数核心验证。

### Python / Lean

- `code/python/fft9_complete.py`
  - Q1 中清理 deprecated prototype 五点 RHS 边界修正旧负号。

- `code/python/helmholtz_solver.py`
  - Q2 将 FFT9 主实现统一为 raw/theory 系统；
  - Q4 补强 raw/theory 系统与整体反号系统的注释；
  - Q6 修正/说明 Neumann compatibility 与 DCT-I 缩放口径。

- `code/python/tests/test_04_fft9_vs_spsolve.py`
  - Q2 同步 sparse compact system 到 \((-L_h+\sigma R_h)u=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B\)。

- `code/python/tests/test_05_neumann_compatibility.py`
  - Q6 增加/调整 pure Neumann weighted compatibility 与 mixed boundary 不触发 pure Neumann compatibility 的测试。

- `code/python/verify_fft9_expansion.py`
  - Q3 将注释/输出从宽泛 “fourth-order accuracy” 收紧为 “fourth-order interior Fourier symbol consistency”。

- `code/python/experiments/exp04_modified_vs_true.py`
  - Q7 增加注释，限定 spectral denominator indicator 作为 \(\mathrm{cond}_2(A)\) 的条件。

- `code/python/experiments/exp07_spectral_denominator_maps.py`
  - Q7 增加注释，说明 risk map 不是条件数图。

- `code/lean4_formalization/SixthOrderImpossibility.lean`
- `code/lean4_formalization/SixthOrderImpossibility/MathlibTest.lean`
  - Q5 仅收紧注释，不改 theorem 语句或证明逻辑。

### 分报告

- `Q1_FIVE_POINT_SIGN_FIX_REPORT.md`
- `Q2_FFT9_RAW_DENOMINATOR_FIX_REPORT.md`
- `Q3_FFT9_FOURTH_ORDER_SYMBOL_CONSISTENCY_REPORT.md`
- `Q4_FFT9_NONHOM_DIRICHLET_SIGN_SYSTEM_REPORT.md`
- `Q5_THREE_PARAMETER_FFT9_SIXTH_ORDER_SCOPE_REPORT.md`
- `Q6_NEUMANN_GHOST_DCT_ZERO_MODE_REPORT.md`
- `Q7_MODIFIED_TRUE_HELMHOLTZ_SPECTRUM_CONDITIONING_REPORT.md`

本文件是上述分报告的总审阅包。

---

## 2. Q1：五点格式截断误差符号与非齐次 Dirichlet RHS

### 修复目标

统一五点正 Laplacian 离散算子、residual / truncation error 方向，并明确非齐次 Dirichlet 边界值移至右端项后的正号。

### 当前论文约定

五点正 Laplacian 离散算子定义为

\[
L_h^{(5)}u_{i,j}
=
\frac{4u_{i,j}-u_{i-1,j}-u_{i+1,j}-u_{i,j-1}-u_{i,j+1}}{h^2},
\]

它近似 \(-\Delta u\)。

residual / truncation error 方向定义为

\[
\tau_h(u)=\bigl(L_h^{(5)}+\sigma I\bigr)u-f.
\]

在该方向下，连续精确解代入五点格式后有

\[
\tau_{i,j}
=
-\frac{h^2}{12}
\left(u_{xxxx}+u_{yyyy}\right)_{i,j}
+O(h^4).
\]

对于靠近左边界的内点 \((1,j)\)，若 \(u_{0,j}=g_{0,j}\)，则

\[
\frac{4u_{1,j}-g_{0,j}-u_{2,j}-u_{1,j-1}-u_{1,j+1}}{h^2}
+\sigma u_{1,j}
=f_{1,j}.
\]

移项得到

\[
\frac{4u_{1,j}-u_{2,j}-u_{1,j-1}-u_{1,j+1}}{h^2}
+\sigma u_{1,j}
=f_{1,j}+\frac{g_{0,j}}{h^2}.
\]

因此五点 Dirichlet 边界贡献移至 RHS 后取正号。

### 代码修复

- 活跃五点求解路径保持正号不变；
- 仅清理 `code/python/fft9_complete.py` 中 deprecated prototype 的五点 RHS 旧负号；
- 未机械修改 FFT9 紧致边界项、Neumann ghost-point 或 archive 文件。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q1 的五点差分符号是否一致。

统一 PDE 是 (-Delta + sigma)u=f。论文定义

    L_h^{(5)}u_ij=(4u_ij-u_{i-1,j}-u_{i+1,j}-u_{i,j-1}-u_{i,j+1})/h^2,

它近似 -Delta u。论文定义 residual/truncation error 方向为

    tau_h(u)=(L_h^{(5)}+sigma I)u-f.

论文写五点截断误差主项为

    tau_ij= -h^2/12 (u_xxxx+u_yyyy)_ij + O(h^4).

同时非齐次 Dirichlet 边界值移至 RHS 后写成 +g/h^2。

请判断：
1. 该 truncation error 主项负号是否正确？
2. 非齐次 Dirichlet RHS 的 +g/h^2 是否正确？
3. 当前 residual 方向是否足够清楚，是否还需要补充一句防误读说明？

请输出 verdict / corrected formula if needed / one-paragraph justification。
```

---

## 3. Q2：FFT9 raw denominator、符号约定与 Helmholtz 解释

### 修复目标

统一 FFT9 九点紧致格式的理论与代码符号，采用

\[
(-L_h+\sigma R_h)u=R_hf,\qquad L_h\approx\Delta.
\]

当物理空间 RHS 已经构造为

\[
\tilde g=(R_h f)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\]

频域求解必须配 raw denominator：

\[
D^{\mathrm{raw}}_{p,q}
=-\hat{\lambda}_L(p,q)+\sigma\hat{\lambda}_R(p,q).
\]

不能把已经包含 \(R_h f\) 和边界修正的 \(\hat{\tilde g}\) 再除以 effective denominator

\[
-\hat{\lambda}_L/\hat{\lambda}_R+\sigma.
\]

### Helmholtz 解释

- modified Helmholtz: \(\sigma=+\kappa^2\)，因 \(-\hat{\lambda}_L>0\)、\(\hat{\lambda}_R>0\)，raw denominator 为正；
- true Helmholtz: \(\sigma=-\kappa^2\)，FFT9 小分母条件为

\[
-\hat{\lambda}_L(p,q)\approx\kappa^2\hat{\lambda}_R(p,q).
\]

这不同于五点条件 \(\lambda_{p,q}^{(5)}\approx\kappa^2\)。

### 代码修复

- `helmholtz_solver.py` 中 FFT9 主实现改为 raw convention；
- `compute_bc_correction_9pt` 返回 \(L_{IB}g_B-\sigma R_{IB}g_B\)；
- `fft9_helmholtz` 使用 `rhs_tilde = apply_Rh_full(F) + bc_corr`；
- denominator 使用 `denom = -lam_L + sigma * lam_R4`；
- `test_04_fft9_vs_spsolve.py` sparse compact system 同步为 raw/theory 系统。

### 验证记录

Q2 分报告记录：

- `py_compile` 通过；
- `verify_fft9_expansion.py` 通过；
- `test_04_fft9_vs_spsolve.py` 与 `test_02_nonhom_dirichlet.py` 通过；
- 当时全量 pytest：`103 passed, 2 skipped`；
- LaTeX 编译通过。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q2 的 FFT9 符号系统是否正确。

论文设 L_h 是九点 Laplacian，约定 L_h approximates Delta。紧致格式写为

    (-L_h + sigma R_h)u = R_h f.

若物理空间已构造

    tilde_g = (R_h f)_I + L_IB g_B - sigma R_IB g_B,

则频域求解写为

    u_hat = tilde_g_hat / D_raw,
    D_raw = -lambda_L + sigma lambda_R.

modified Helmholtz 下 sigma=+kappa^2，因此 D_raw>0。
true Helmholtz 下 sigma=-kappa^2，因此小分母条件为

    -lambda_L ≈ kappa^2 lambda_R.

请判断：
1. `(-L_h+sigma R_h)u=R_h f` 是否与 `(-Delta+sigma)u=f` 一致？
2. raw denominator 是否必须与 `tilde_g=(R_hf)_I+bc` 配对？
3. modified/true Helmholtz 的分母解释是否准确？
4. 是否还存在把 effective denominator 与 raw RHS 混用的风险？

请输出 verdict / any correction / concise reasoning。
```

---

## 4. Q3：FFT9 四阶展开与 consistency 表述收紧

### 修复目标

保留 FFT9 有效特征值展开公式，但将理论结论限定为 interior / local Fourier symbol consistency，不把该推导单独写成完整 Dirichlet 全局误差定理。

### 当前定理口径

设

\[
\alpha=\xi h,\qquad \beta=\eta h.
\]

FFT9 有效离散特征值为

\[
\lambda_{\mathrm{discrete}}(\alpha,\beta)
=-\frac{\widehat{\lambda}_L(\alpha,\beta)}
{\widehat{\lambda}_R(\alpha,\beta)}
+\sigma.
\]

论文保留展开

\[
\lambda_{\mathrm{discrete}}
=
\xi^2+\eta^2+\sigma
-
\frac{(\xi^2+\eta^2)(3\xi^4-8\xi^2\eta^2+3\eta^4)}{720}h^4
+O(h^6).
\]

其中 \(h^2\) 项消去：

\[
A\cdot\frac{A}{12}-\frac{A^2}{12}=0,
\qquad A=\xi^2+\eta^2.
\]

### 表述降级

论文现在写为：

- 第 2 章：FFT9 在 interior Fourier symbol 意义下具有 \(O(h^4)\) consistency；
- 第 3 章：FFT9 Dirichlet 实现的整体四阶收敛由第 6 章 manufactured solution 数值实验验证；
- 第 4 章：true Helmholtz near-resonance 小分母可能放大整体误差，一致性不等于 near-resonance 下稳定四阶全局误差；
- 第 6 章：只说当前规则 Dirichlet 光滑 manufactured solution 上呈现四阶收敛。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q3 的 FFT9 Fourier symbol 展开和结论范围。

论文定义九点 Laplacian symbol:

    lambda_L = [-20 + 8(cos alpha + cos beta) + 4 cos alpha cos beta] / (6h^2),

compact averaging symbol:

    lambda_R = 2/3 + (1/6)(cos alpha + cos beta).

其中 L_h approximates Delta，所以 lambda_L leading term 是 -(xi^2+eta^2)。
论文写有效特征值

    lambda_discrete = -lambda_L/lambda_R + sigma
      = xi^2 + eta^2 + sigma
        - ((xi^2+eta^2)(3xi^4-8xi^2 eta^2+3eta^4))/720 h^4
        + O(h^6).

请判断：
1. h^2 项是否确实消去？
2. h^4 系数是否正确？
3. 该结果是否只能证明 interior/local Fourier symbol consistency？
4. 将完整 Dirichlet 四阶收敛放到第 6 章实验验证，是否是更稳妥的论文口径？

请输出 verdict / corrected expansion if needed / scope judgment。
```

---

## 5. Q4：FFT9 非齐次 Dirichlet RHS 符号系统说明

### 修复目标

说明当前代码实际采用的就是论文 raw/theory 系统，而不是附件中假设的“当前代码整体反号系统”。本轮不改变计算逻辑，只补充防混用说明。

### 当前代码与论文一致的系统

当前论文与代码均采用

\[
(-L_h+\sigma R_h)u=R_hf,
\]

\[
\tilde g
=(R_hf)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\]

\[
D^{\mathrm{raw}}_{p,q}
=-\hat{\lambda}_L(p,q)+\sigma\hat{\lambda}_R(p,q).
\]

对应代码语义：

- `compute_bc_correction_9pt` 返回 \(L_{IB}g_B-\sigma R_{IB}g_B\)；
- `fft9_helmholtz` 构造 `rhs_tilde = apply_Rh_full(F) + bc_corr`；
- `fft9_helmholtz` 使用 `denom = -lam_L + sigma * lam_R4`。

### 等价反号系统的 remark

若某实现整体乘以 \(-1\)，则必须同时使用

\[
\tilde g_{\mathrm{code}}
=
-(R_hf)_I-L_{IB}g_B+\sigma R_{IB}g_B,
\]

\[
D^{\mathrm{code}}_{p,q}
=\hat{\lambda}_L-\sigma\hat{\lambda}_R
=-D^{\mathrm{raw}}_{p,q}.
\]

不能只反 RHS 或只反 denominator。

### 验证记录

Q4 分报告记录：

- `py_compile helmholtz_solver.py` 通过；
- `test_04_fft9_vs_spsolve.py` 与 `test_02_nonhom_dirichlet.py` 通过；
- 轻量非齐次 Dirichlet FFT9 manufactured-solution 收敛检查呈四阶趋势；
- LaTeX 编译通过。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q4 的 FFT9 非齐次 Dirichlet RHS 符号说明是否稳妥。

当前论文和代码都采用 raw/theory 系统：

    (-L_h + sigma R_h)u = R_h f,
    tilde_g = (R_h f)_I + L_IB g_B - sigma R_IB g_B,
    D_raw = -lambda_L + sigma lambda_R.

论文另加 remark：若某实现整体乘以 -1，则必须同时使用

    tilde_g_code = -(R_h f)_I - L_IB g_B + sigma R_IB g_B,
    D_code = lambda_L - sigma lambda_R = -D_raw.

请判断：
1. raw/theory 系统中的 boundary correction 符号是否正确？
2. 整体反号系统的 remark 是否数学等价？
3. 论文是否应保留该 remark，还是会让读者觉得复杂？
4. 当前写法是否足以防止 RHS 与 denominator 只反一个的错误？

请输出 verdict / minimal wording suggestion。
```

---

## 6. Q5：三参数 FFT9 六阶不可达性与 Lean 覆盖范围

### 修复目标

将“三参数九点修正族六阶不可达性”限定为：

- 三参数九点修正族；
- \(\alpha,\beta,\gamma\in\mathbb R\) 是 mode-independent constants；
- local Fourier symbol / interior consistency 意义下；
- Lean 4 只验证反证中的代数核心。

### 当前定理口径

论文讨论三参数族：

\[
(-L_h-\alpha h^2L_h^{(5)})u
+\sigma(R_h+\beta h^2 I)u
=(R_h+\gamma h^2 I)f.
\]

结论：不存在一组与 Fourier mode 无关的统一常数 \((\alpha,\beta,\gamma)\)，使该族在 interior Fourier symbol 意义下达到六阶一致性。

反证逻辑：

1. 要求 \(c_2=0\) 给出必要参数关系；
2. 代入后要求 \(c_4=0\)；
3. 使用两个连续波数测试模式得到矛盾；
4. \((\xi,\eta)=(1,0)\) 是 local Fourier symbol 的连续波数测试，不是 Dirichlet DST-I 离散模态编号，因此允许某方向为 0。

### Lean 4 覆盖范围

Lean 4 只覆盖：

- \(c_2=0\) 的必要关系；
- 代入后两个 Fourier 测试模式的 \(c_4=0\) 条件矛盾；
- Poisson/Helmholtz 三参数族在统一常数参数下的 local symbol 六阶条件不相容。

Lean 4 不覆盖：

- Fourier symbol 从模板导出的全过程；
- Taylor 展开建模；
- DST-I 谱理论；
- 边界闭合；
- 完整 PDE 全局误差估计；
- Python 程序正确性；
- 数值实验正确性。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q5 的三参数 FFT9 六阶不可达性定理范围是否正确。

论文考虑三参数族

    (-L_h - alpha h^2 L_h^{(5)})u
      + sigma(R_h + beta h^2 I)u
      = (R_h + gamma h^2 I)f,

其中 alpha,beta,gamma 是与 Fourier mode 无关的统一常数。
论文通过 local Fourier symbol 展开，要求 c2=0，再要求 c4=0；代入 c2 条件后，用两个连续波数测试模式得到矛盾。

论文现在把结论限定为：

    该三参数九点修正族 + mode-independent constants + local Fourier symbol / interior consistency
    意义下无法达到六阶一致性。

Lean 4 只说验证该反证中的 c2/c4 代数矛盾。

请判断：
1. 这种反证是否足以证明 restricted family 的六阶 local consistency 不可达？
2. 是否必须保留“三参数族、统一常数参数、local symbol”三个限定？
3. “四阶是该族在统一常数参数约束下的自然精度上限”是否仍然过强？
4. Lean 4 覆盖范围这样写是否足够谨慎？

请输出 verdict / overclaim risks / suggested safer wording。
```

---

## 7. Q6：Neumann ghost-point、DCT-I 缩放与 zero mode

### 修复目标

统一 Neumann ghost-point 矩阵、DCT-I 对称化缩放、pure Neumann Poisson compatibility 和 zero-mode normalization 的论文/代码口径。

### 当前论文口径

一维三节点 ghost-point 矩阵示例为

\[
G=\frac1{h^2}
\begin{pmatrix}
2&-2&0\\
-1&2&-1\\
0&-2&2
\end{pmatrix},
\qquad
G\mathbf 1=0.
\]

使用对角缩放

\[
\tilde u=D^{-1}u,\qquad
\tilde f=D^{-1}f,\qquad
S=D^{-1}GD,
\]

得到

\[
S\tilde u=\tilde f,\qquad u=D\tilde u.
\]

二维写法为

\[
G_xU+UG_y+\sigma U=F,
\]

\[
\tilde U=D_x^{-1}UD_y^{-1},\qquad
\tilde F=D_x^{-1}FD_y^{-1},
\]

\[
S_x\tilde U+\tilde U S_y+\sigma\tilde U=\tilde F.
\]

频域求解：

\[
\widehat{\tilde U}_{k,l}
=
\frac{\widehat{\tilde F}_{k,l}}
{\lambda_k+\lambda_l+\sigma},
\qquad
U=D_x\tilde U D_y.
\]

pure Neumann Poisson compatibility 写为 trapezoid-weighted condition：

\[
\sum_{i,j}w_iw_j f_{ij}=0,
\qquad
w_0=w_{n-1}=1/2,\quad w_i=1.
\]

zero mode normalization 写为

\[
\widehat{\tilde U}_{0,0}=0,
\]

等价于恢复后的 \(U\) 满足 trapezoid weighted mean zero。

### modified / true Helmholtz 区分

- pure Neumann Poisson: \(d_{0,0}=0\)，需要 compatibility；
- modified Helmholtz: \(d_{0,0}=\kappa^2>0\)，不需要 compatibility；
- true Helmholtz: zero mode 分母为 \(-\kappa^2\)，通常不是 zero-mode compatibility 问题，但其他模态可能 near-resonance。

### 代码与测试修复

- `helmholtz_solver.py` 中 Neumann 注释统一为 `F_tilde = D_x^{-1} F D_y^{-1}`；
- `check_neumann_compatibility` 仅用于 pure Neumann Poisson；
- compatibility 诊断按 DCT-I / trapezoid endpoint weights 判断；
- mixed Dirichlet/Neumann Poisson 不触发 pure Neumann compatibility；
- `test_05_neumann_compatibility.py` 增加非对称 RHS 与 mixed boundary 测试。

### 验证记录

Q6 分报告记录：

- `py_compile helmholtz_solver.py` 通过；
- `tests/test_01_eigenvalues.py` 与 `tests/test_05_neumann_compatibility.py`: `38 passed`；
- 当时全量 pytest: `105 passed, 2 skipped`；
- LaTeX 编译通过；
- PDF 当时为 69 页，后续 Q7 最新编译为 70 页。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q6 的 Neumann ghost-point / DCT-I / zero mode 口径。

论文使用一维 ghost-point 矩阵

    G = (1/h^2) [[2,-2,0],[-1,2,-1],[0,-2,2]],
    G 1 = 0.

通过对角缩放

    tilde u = D^{-1}u, tilde f=D^{-1}f, S=D^{-1}GD,
    S tilde u = tilde f, u=D tilde u.

二维写法为

    G_x U + U G_y + sigma U = F,
    tilde U = D_x^{-1} U D_y^{-1},
    tilde F = D_x^{-1} F D_y^{-1},
    S_x tilde U + tilde U S_y + sigma tilde U = tilde F.

pure Neumann Poisson compatibility 写为 trapezoid weighted sum:

    sum_i,j w_i w_j f_ij = 0, w endpoints=1/2.

zero mode normalization 写为 hat{tilde U}_{0,0}=0，等价于恢复后 weighted mean zero。

请判断：
1. 该缩放方向是否正确？
2. compatibility 是否应使用 trapezoid endpoint weights，而不是普通 mean？
3. modified Helmholtz 是否不需要 pure Neumann compatibility？
4. mixed Dirichlet/Neumann Poisson 是否不应触发 pure Neumann compatibility？

请输出 verdict / corrected equations if needed / implementation-risk notes。
```

---

## 8. Q7：Modified/True Helmholtz 谱分类、条件数与谱分母指标

### 修复目标

将 modified Helmholtz 条件数改善、spectral denominator indicator 与 \(\mathrm{cond}_2(A)\) 的等价性严格限定在当前可证明范围内。

### modified Helmholtz 条件数范围

在固定网格和同一 Dirichlet 五点 SPD 系统中：

\[
A_\kappa=A_0+\kappa^2I,\qquad \kappa>0.
\]

若 \(0<\lambda_{\min}<\lambda_{\max}\)，则

\[
\kappa_2(A_\kappa)
=
\frac{\lambda_{\max}+\kappa^2}{\lambda_{\min}+\kappa^2}
<
\frac{\lambda_{\max}}{\lambda_{\min}}
=\kappa_2(A_0).
\]

当 \(\kappa=0\) 时退化为 Poisson，条件数相等。

该结论不自动外推到：

- 原始 Neumann ghost-point 非对称矩阵；
- mixed 边界；
- 非正交对角化系统；
- 一般 non-normal 系统。

### true Helmholtz 谱分类

true Helmholtz:

\[
A_h=A_0-\kappa^2I.
\]

按 \(\kappa^2\) 与离散谱关系可分为：

- \(\kappa^2<\lambda_{\min}\)：正定；
- \(\kappa^2=\lambda_{p,q}\)：奇异；
- \(\lambda_{\min}<\kappa^2<\lambda_{\max}\) 且不等于特征值：不定且非奇异；
- \(\kappa^2>\lambda_{\max}\)：\(A_h\) 负定，此时 \(-A_h\) 为 SPD，但 \(A_h\) 本身不是 SPD。

near-resonance 条件为

\[
|d_{p,q}|=|\lambda_{p,q}-\kappa^2|\ll1,
\]

即 \(\kappa^2\approx\lambda_{p,q}\)，不是 modified Helmholtz 的 \(\lambda_{p,q}+\kappa^2\approx0\)。

### spectral denominator indicator 与 cond₂

仅在当前 Dirichlet 五点正交 DST-I 系统中：

\[
A=Q^T\operatorname{diag}(d_{p,q})Q,\qquad Q^TQ=I.
\]

若 \(A\) 非奇异且 symmetric/normal，则奇异值为 \(|d_{p,q}|\)，所以

\[
\mathrm{cond}_2(A)
=
\frac{\max |d_{p,q}|}{\min |d_{p,q}|}.
\]

因此第 6 章 exp04 condition check 只是对当前 Dirichlet 五点系统的数值校验，不是一般条件数理论。

图 10 risk map 只可视化 small-denominator risk，不是条件数图。

### 验证记录

Q7 最新验证：

- `py_compile exp04_modified_vs_true.py exp07_spectral_denominator_maps.py` 通过；
- LaTeX 完整编译通过；
- `main.pdf` 最新为 70 页；
- log 扫描无 fatal error、undefined reference/citation、multiply-defined label、overfull hbox；
- 仅保留 SimSun 字体形状 warning；
- 未运行 exp04/exp07 实验脚本，未重绘图，未修改 CSV。

### 需要 GPT-5.5 审阅的问题

```text
请审查 Q7 的 modified/true Helmholtz 谱分类与条件数口径。

论文现在只在固定网格、同一 Dirichlet 五点 SPD 系统中声明：

    A_kappa = A_0 + kappa^2 I, kappa>0,
    cond_2(A_kappa) = (lambda_max+kappa^2)/(lambda_min+kappa^2)
                    < lambda_max/lambda_min = cond_2(A_0).

并明确 kappa=0 时相等。

true Helmholtz 写为 A_h=A_0-kappa^2 I，并分为：
positive definite / singular / indefinite nonsingular / negative definite。
当 kappa^2>lambda_max 时，-A_h 是 SPD，但 A_h 本身不是 SPD。

near-resonance 条件写为

    |d_pq| = |lambda_pq-kappa^2| << 1.

论文还说明 spectral denominator indicator = max|d|/min|d| 只有在 Dirichlet 五点正交 DST-I 对角化系统

    A=Q^T diag(d_pq) Q, Q^T Q=I

中才等于 cond_2(A)，不能外推到 Neumann ghost-point 原始非对称矩阵、mixed 边界或 non-normal 系统。

请判断：
1. modified Helmholtz 条件数改善证明和限制是否正确？
2. true Helmholtz 四类谱分类是否完整？
3. negative definite case 的表述是否足够防止误称 SPD？
4. spectral denominator indicator 与 cond_2(A) 的等价条件是否写得足够严格？
5. risk map 不是条件数图这一限制是否必要且正确？

请输出 verdict / any overclaim remaining / safer wording if needed。
```

---

## 9. 总体验证状态

### 最新已确认

- Q7 后 LaTeX 完整编译通过；
- `thesis/main.pdf` 正常生成；
- 最新 PDF 页数：70；
- LaTeX log 扫描：
  - no fatal error；
  - no undefined reference/citation；
  - no multiply-defined label；
  - no overfull hbox；
  - 剩余为 SimSun 字体形状 warning。
- Q7 涉及的两个实验脚本只做注释修改，`py_compile` 通过。
- `git diff --check` 对本轮相关文件通过，仅有 LF/CRLF 行尾提示。

### Q2/Q4/Q6 分报告中记录的测试

- Q2：
  - `test_04_fft9_vs_spsolve.py`、`test_02_nonhom_dirichlet.py` 通过；
  - 当时全量 pytest `103 passed, 2 skipped`。

- Q4：
  - `test_04_fft9_vs_spsolve.py`、`test_02_nonhom_dirichlet.py` 通过；
  - 轻量 FFT9 非齐次 Dirichlet manufactured-solution 收敛检查呈四阶趋势。

- Q6：
  - `test_01_eigenvalues.py`、`test_05_neumann_compatibility.py`: `38 passed`；
  - 当时全量 pytest `105 passed, 2 skipped`。

### 未做或不应做

- Q7 未重跑 exp04/exp07 实验脚本；
- 未重绘图；
- 未修改 CSV；
- 未新增实验；
- 未把 GMRES 结论扩展到预处理 GMRES、MINRES 或所有 Krylov 方法；
- 未把 FFT9 扩展到 Neumann/mixed；
- 未把 Lean 4 扩展成完整 PDE 形式化证明。

---

## 10. 建议给 GPT-5.5 的总审阅 prompt

可以直接复制下面这段作为总任务开头，然后附上本报告全文或逐项附上对应章节片段。

```text
你是一个数值 PDE / 有限差分 / 谱方法 / Krylov 方法审稿专家。请审阅我论文中 Q1--Q7 这轮理论与实验表述修复是否数学上稳妥。

论文统一模型为

    (-Delta + sigma)u=f,

其中 sigma=0 是 Poisson，sigma=+kappa^2 是 modified Helmholtz，sigma=-kappa^2 是 true Helmholtz。论文包含 Dirichlet 五点 DST-I、Neumann ghost-point + DCT-I、Dirichlet 九点紧致 FFT9、FA/CR/FACR-like 直接法、无预处理 restarted GMRES(30)，以及 true Helmholtz near-resonance 离散谱实验。

请按 Q1--Q7 逐项判断：

1. 公式是否正确；
2. 符号约定是否自洽；
3. 代码实现口径是否与论文公式一致；
4. 定理结论是否被限定到正确范围；
5. 实验解释是否超过已有数据；
6. 是否还存在应降级为 remark、implementation note 或 experiment observation 的表述。

请不要泛泛评价语言风格。请输出：

- 每个 Q 的 verdict: OK / minor wording risk / mathematical issue / needs more evidence；
- 若有问题，给出 corrected formula 或 safer wording；
- 最后列出 blocking issues 和 non-blocking suggestions。
```

---

## 11. GPT-5.5 应重点判定的潜在风险

下面这些不是我认为一定有错，而是最值得外部模型复核的点：

1. Q3 的 \(h^4\) 系数是否完全正确。
2. Q5 的三参数六阶不可达性反证是否需要更多假设，例如对 \(\sigma\)、mode 选择或参数独立性的说明。
3. Q6 中 DCT-I scaling 与 trapezoid weighted compatibility 是否与代码中 `F_adj` 的缩放方向完全一致。
4. Q7 中 modified Helmholtz 条件数改善是否应该进一步强调“固定 \(h\)”和“同一离散矩阵族”。
5. Q2/Q4 中 raw denominator 与 boundary correction 是否在所有 Dirichlet FFT9 路径，包括 `fft9_oer_helmholtz`，都完全配对。
6. 第 6 章实验解释是否仍有任何地方把 risk map、spectral denominator indicator 或 near-resonance residual history 说得过强。

---

## 12. 审阅后建议的输出格式

请 GPT-5.5 按如下结构返回：

```markdown
# Q1--Q7 External Review

## Executive Verdict
- Blocking issues:
- Non-blocking wording risks:
- Safe to proceed? yes/no/conditional

## Q1 Five-point Sign
- Verdict:
- Formula check:
- Suggested wording:

## Q2 FFT9 Raw Denominator
- Verdict:
- Formula check:
- Paper-code consistency:

## Q3 FFT9 Fourier Symbol Consistency
- Verdict:
- Expansion check:
- Scope check:

## Q4 FFT9 Nonhomogeneous Dirichlet Sign System
- Verdict:
- Boundary correction check:
- Remark usefulness:

## Q5 Three-Parameter Sixth-Order Scope
- Verdict:
- Proof-scope check:
- Lean-scope check:

## Q6 Neumann DCT-I / Zero Mode
- Verdict:
- Scaling check:
- Compatibility check:

## Q7 Helmholtz Spectrum / Conditioning
- Verdict:
- Modified Helmholtz:
- True Helmholtz:
- Condition-number indicator:

## Minimal Patch Suggestions
1.
2.
3.
```

---

## 13. 总结

Q1--Q7 的整体方向是把论文从“结果大体正确但容易被误读”收紧到“每个结论都有明确系统、边界条件、离散格式和实验口径”。最关键的降级包括：

- 五点 truncation error 明确 residual 方向；
- FFT9 raw denominator 与 RHS 配对；
- FFT9 四阶只在 interior Fourier symbol 层面作为理论 consistency，完整 Dirichlet 收敛由实验验证；
- 三参数六阶不可达性只限 restricted family；
- Lean 4 只验证代数核心；
- Neumann compatibility 使用 DCT-I/trapezoid weighted condition；
- modified Helmholtz 条件数改善只限同一 Dirichlet 五点 SPD 系统；
- spectral denominator indicator 与 \(\mathrm{cond}_2(A)\) 的等价只限正交 DST-I 对角化系统；
- true Helmholtz near-resonance 只作为离散谱小分母机制，不外推到连续谱极限或所有 Helmholtz 求解器。
