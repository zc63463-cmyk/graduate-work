# Questions for GPT-5.5 Pro

用途：把论文中最值得外部复核的理论、推导、形式化验证和实验结论边界拆成可直接复制给 GPT-5.5 Pro 的问题。
建议用法：不要一次性全部发送。每次发送 1--3 个问题，并附上对应论文片段或截图。优先发送 Q1--Q5，因为它们最接近“定理/推导正确性”。

统一上下文可在每次提问前附上：

```text
我在写一篇数值 PDE 硕士论文，研究矩形区域上 Poisson / modified Helmholtz / true Helmholtz 方程的 FFT 快速直接求解器与无预处理 GMRES(30) 对比。统一模型为

    (-Delta + sigma) u = f

其中 sigma=0 是 Poisson，sigma=+kappa^2 是 modified Helmholtz，sigma=-kappa^2 是 true Helmholtz。论文使用五点二阶差分、Dirichlet DST-I、Neumann ghost-point + DCT-I、Dirichlet 九点紧致四阶 FFT9、以及无预处理 restarted GMRES(30)。请只判断下面这段公式/推导/结论是否数学上稳妥，不需要评价语言风格。
```

---

## Formula Derivation

## Q1. 五点格式截断误差符号与 residual 方向是否一致？

**需要验证的内容**

- 论文第 2.2 节定义五点正 Laplacian 离散：

  ```text
  L_h^{(5)} u + sigma u = f
  ```

- residual / truncation error 方向定义为：

  ```text
  tau_h(u)=L_h^{(5)}u + sigma u - f
  ```

- 在该约定下，论文写五点格式截断误差主项：

  ```text
  tau_{i,j} = - h^2/12 (u_xxxx + u_yyyy) + O(h^4)
  ```

- 非齐次 Dirichlet 边界贡献移至 RHS 后写成：

  ```text
  f_tilde = f + boundary_values / h^2
  ```

**Prompt**

```text
请审查以下五点差分符号是否一致。

统一 PDE 是 (-Delta + sigma)u=f。论文把五点正 Laplacian 离散算子定义为

    L_h^{(5)}u = (4u_ij-u_{i-1,j}-u_{i+1,j}-u_{i,j-1}-u_{i,j+1})/h^2,

并定义 residual/truncation error 方向为

    tau_h(u)=L_h^{(5)}u + sigma u - f.

论文声称在该约定下

    tau_h(u)= - h^2/12 (u_xxxx+u_yyyy) + O(h^4).

同时，对于非齐次 Dirichlet 边界，边界值在左端以 -g/h^2 出现，移至右端后为 +g/h^2。

请判断：
1. 该 truncation error 主项符号是否正确？
2. `f_tilde = f + boundary_values/h^2` 的符号是否正确？
3. 如果要写得最稳妥，论文中应如何一句话定义 residual 方向？

请给出 verdict / corrected formula if needed / brief derivation。
```

---

## Q2. FFT9 紧致格式的基本方程与原始分母是否正确？

**需要验证的内容**

- 论文定义九点 Laplacian `L_h ≈ Delta`，不是 `-Delta`。
- 紧致格式写为：

  ```text
  (-L_h + sigma R_h) u = R_h f
  ```

- 频域原始分母：

  ```text
  D_raw(p,q) = -lambda_L(p,q) + sigma lambda_R(p,q)
  ```

- modified Helmholtz 下 `sigma=+kappa^2`，因为 `-lambda_L>0` 且 `lambda_R>0`，所以分母为正。
- true Helmholtz 下 `sigma=-kappa^2`，小分母条件为 `-lambda_L ≈ kappa^2 lambda_R`。

**Prompt**

```text
请审查我论文中 FFT9 九点紧致格式的符号。

设 L_h 是九点 Laplacian 模板，约定 L_h ≈ Delta，而不是 -Delta。R_h 是 compact averaging operator。论文写紧致格式为

    (-L_h + sigma R_h) u = R_h f.

Fourier/DST-I 空间中，L_h 的 symbol 为 lambda_L，R_h 的 symbol 为 lambda_R，因此频域原始分母写为

    D_raw(p,q) = -lambda_L(p,q) + sigma lambda_R(p,q).

对 modified Helmholtz, sigma=+kappa^2，论文说因 -lambda_L>0 且 lambda_R>0，D_raw>0。
对 true Helmholtz, sigma=-kappa^2，小分母条件是

    -lambda_L(p,q) ≈ kappa^2 lambda_R(p,q).

请判断：
1. `(-L_h + sigma R_h)u=R_h f` 是否与 `(-Delta+sigma)u=f` 一致？
2. `D_raw=-lambda_L+sigma lambda_R` 是否正确？
3. modified / true Helmholtz 的分母解释是否正确？
4. 是否需要额外说明 lambda_L 本身是负的，因为 L_h approximates Delta？

请给 verdict / corrected statement if needed / concise mathematical justification。
```

---

## Q3. FFT9 四阶有效特征值展开是否正确？

**需要验证的内容**

论文定理声称：

```text
lambda_discrete(alpha,beta)
  = -lambda_L/lambda_R + sigma
  = xi^2 + eta^2 + sigma
    - ((xi^2+eta^2)(3xi^4 - 8xi^2 eta^2 + 3eta^4))/720 * h^4
    + O(h^6).
```

其中 `alpha=xi h`, `beta=eta h`。

**Prompt**

```text
请审查这个 FFT9 四阶 Fourier symbol 展开是否正确。

论文定义：

    R_h symbol:
    lambda_R = 2/3 + (1/6)(cos alpha + cos beta)

    nine-point Laplacian symbol:
    lambda_L = [-20 + 8(cos alpha + cos beta) + 4 cos alpha cos beta] / (6h^2)

其中 L_h approximates Delta，所以 lambda_L 的 leading term 应为 -(xi^2+eta^2)，alpha=xi h, beta=eta h。

论文声称有效离散特征值为

    lambda_discrete = -lambda_L/lambda_R + sigma
    = xi^2 + eta^2 + sigma
      - ((xi^2+eta^2)(3xi^4 - 8xi^2 eta^2 + 3eta^4))/720 * h^4
      + O(h^6).

请判断：
1. h^2 项是否确实消去？
2. h^4 项系数是否正确？
3. 该结论是否足以说明该紧致格式对光滑 Dirichlet 问题具有四阶 consistency？
4. 这是否只说明 local Fourier symbol / interior consistency，而不是完整边界全局误差证明？

请给 verdict / corrected expansion if needed / brief derivation or reasoning。
```

---

## Q4. FFT9 非齐次 Dirichlet 边界修正符号是否稳妥？

**需要验证的内容**

论文算法写：

```text
tilde_g = (R_h f)_I + L_IB g_B - sigma R_IB g_B
```

代码实现中因为整体乘以 `-1` 解等价系统：

```text
(L_h - sigma R_h)u = -R_h f
```

所以代码边界修正使用相反号：

```text
-L_IB g_B + sigma R_IB g_B
```

**Prompt**

```text
请审查 FFT9 非齐次 Dirichlet 边界修正符号。

论文采用紧致系统

    (-L_h + sigma R_h)u = R_h f,

其中内点未知量与边界已知值分块后，论文算法写物理空间修正右端项为

    tilde_g = (R_h f)_I + L_IB g_B - sigma R_IB g_B.

代码实现为了方便，整体乘以 -1，求解等价系统

    (L_h - sigma R_h)u = -R_h f,

因此代码中的边界修正写成

    -L_IB g_B + sigma R_IB g_B.

请判断：
1. 论文算法中的 `+L_IB g_B - sigma R_IB g_B` 是否正确？
2. 代码中整体乘以 -1 后使用 `-L_IB g_B + sigma R_IB g_B` 是否等价？
3. 论文是否需要明确说明“代码实现使用等价符号系统，因此代码分母和 RHS 同时反号”以避免误解？

请输出 verdict / corrected formula if needed / minimal thesis wording。
```

---

## Q5. 三参数九点修正族六阶不可达性的定理范围是否正确？

**需要验证的内容**

论文考虑修正族：

```text
(-L_h - alpha h^2 L_h^{(5)})u
  + sigma(R_h + beta h^2 I)u
  = (R_h + gamma h^2 I)f.
```

论文声称不存在统一常数 `(alpha,beta,gamma)` 使该三参数族达到 `O(h^6)`。Lean 4 只验证代数核心：`c2=0` 推出参数关系，然后用两个 Fourier modes 的 `c4` 条件矛盾。

**Prompt**

```text
请审查以下“九点三参数修正族无法达到六阶”的定理表述是否范围足够准确。

论文考虑修正紧致格式：

    (-L_h - alpha h^2 L_h^{(5)})u
      + sigma(R_h + beta h^2 I)u
      = (R_h + gamma h^2 I)f,

其中 alpha,beta,gamma 是与 mode 无关的统一常数。通过 local Fourier symbol 展开，要求 h^2 项 c2=0，再看 h^4 项 c4 是否可同时为 0。论文用两个测试模式，例如 (xi,eta)=(1,0) 与 (1,1)，得到互相矛盾的 gamma 条件，从而说明不存在统一常数参数使所有模式同时达到 O(h^6)。

论文还说明这里的 (xi,eta) 是局部 Fourier symbol 的连续波数测试，不是 Dirichlet DST 离散模态编号，因此允许某一方向波数为 0。

请判断：
1. 这种反证逻辑是否足以证明“该三参数族不存在 mode-independent constants 达到六阶 local consistency”？
2. 是否必须把结论限定为“该三参数九点修正族 / 统一常数参数 / local Fourier symbol”？
3. 是否可以说“四阶是该族在统一常数参数约束下的自然精度上限”，还是这个措辞太强？
4. Lean 4 只验证 c2/c4 多项式矛盾时，论文应该如何描述 Lean 的证明范围？

请给 verdict / safest theorem wording / limitations that must be stated。
```

---

## Q6. Neumann ghost-point 矩阵与 DCT-I 对称化推导是否完整？

**需要验证的内容**

论文对齐次 Neumann 使用 ghost-point：

```text
u_{-1,j}=u_{1,j}
```

一维三个物理节点例子：

```text
G = (1/h^2) [[2,-2,0],[-1,2,-1],[0,-2,2]]
```

非对称矩阵用：

```text
D=diag(sqrt(2),1,...,1,sqrt(2)), S=D^{-1}GD
```

对称化后由 DCT-I 对角化。

**Prompt**

```text
请审查 Neumann ghost-point + DCT-I 对称化推导是否正确。

论文对齐次 Neumann 边界采用 ghost-point：

    (u_{-1,j}-u_{1,j})/(2h)=0 => u_{-1,j}=u_{1,j}.

一维三个物理节点时，正 Laplacian 矩阵写作

    G = (1/h^2) [[2,-2,0],
                 [-1,2,-1],
                 [0,-2,2]].

该矩阵非对称。论文引入

    D = diag(sqrt(2),1,...,1,sqrt(2)),
    S = D^{-1} G D,

并声称 S 可由 DCT-I 正交对角化。

请判断：
1. 这个 3x3 ghost-point Neumann 矩阵是否正确？
2. `D^{-1}GD` 的对称化方向是否正确？
3. DCT-I 对角化和 endpoint weights / weighted mean 的说明是否足够？
4. 对 pure Neumann Poisson, 是否必须明确 zero-mode compatibility 和 weighted mean normalization？

请给 verdict / any corrected matrix or scaling / minimal thesis wording。
```

---

## Q7. Modified 与 True Helmholtz 的谱分类和条件数表述是否稳妥？

**需要验证的内容**

论文使用五点 Dirichlet 频域分母：

```text
d_{p,q}=lambda_{p,q}+sigma
```

- modified: `sigma=+kappa^2`, SPD, all denominators positive。
- true: `sigma=-kappa^2`, denominators can change sign; near resonance if `kappa^2` close to a discrete eigenvalue。
- 条件数/谱分母指标等价只限 Dirichlet 五点正交 DST-I 系统。

**Prompt**

```text
请审查 modified Helmholtz 与 true Helmholtz 的谱分类表述。

对五点 Dirichlet 离散，论文写

    d_{p,q}=lambda_{p,q}+sigma,

其中 lambda_{p,q}>0 是 -Delta_h 的离散特征值。

modified Helmholtz: sigma=+kappa^2，所以 d_{p,q}>0，矩阵 SPD。
true Helmholtz: sigma=-kappa^2，所以 d_{p,q}=lambda_{p,q}-kappa^2，矩阵可能正定、不定或奇异；当 kappa^2 接近某个 lambda_{p,q} 时发生 near-resonance。

论文还说：在 Dirichlet 五点正交 DST-I 系统中，spectral denominator indicator 与 cond_2(A) 一致；但不外推到 Neumann ghost-point 原始非对称矩阵、mixed 边界或一般 non-normal 系统。

请判断：
1. modified / true Helmholtz 的谱分类是否正确？
2. “modified Helmholtz 条件数严格小于 Poisson 条件数”在 Dirichlet 五点 SPD 情况下是否成立？是否需要限定 kappa>0 和相同网格？
3. spectral denominator indicator 与 cond_2(A) 的等价性是否只应限于正交对角化的 symmetric normal 系统？
4. 论文目前这种限制是否足够保守？

请给 verdict / corrected or safer wording。
```

---

## GMRES / Krylov Theory

## Q8. GMRES 复杂度和内存表述是否还需更精确？

**需要验证的内容**

论文现在将 `N` 作为单方向内点数，总未知量 `M=N^2`。复杂度表写：

```text
GMRES cost: O(m N^2)
memory: O(m N^2)
```

其中 `m` 有时表示总迭代次数，有时也表示 restart length。需要 GPT 判断是否需要区分：

- restart length `r`
- total iterations `k`
- matrix-vector cost
- Arnoldi orthogonalization cost

**Prompt**

```text
请审查论文中 GMRES 复杂度表述是否足够准确。

论文中 N 表示单方向内点数，总未知量 M=N^2。当前复杂度表写 GMRES cost 为 O(mN^2)，memory 为 O(mN^2)，其中 m 解释为 GMRES 总迭代次数；重启 GMRES(m) 中 m 又表示 restart length。

更严格地说，对 sparse five-point matrix，单次 matvec 是 O(M)=O(N^2)。但 Arnoldi 正交化需要与当前 restart cycle 内已有 Krylov 向量正交，若 restart length 为 r，总迭代数为 k，则成本可能更接近 O(k r M) 或 O(kM + k r M)，memory 为 O(rM)，而不是 O(kM)。

请判断：
1. 若论文只是做高层复杂度比较，写 `O(mN^2)` 是否可以接受？
2. 是否应把符号改为：restart length 为 `r`，total iterations 为 `k`，memory `O(rN^2)`，cost roughly `O(k r N^2)` including orthogonalization, or `O(kN^2)` if ignoring orthogonalization / fixed r?
3. 为避免错误，表格和正文最安全的写法是什么？

请输出 verdict / recommended notation / minimal patch wording。
```

---

## Q9. GMRES 收敛率定理是否适合放在论文中？

**需要验证的内容**

论文第 5 章引用正定型 GMRES 收敛率估计，条件是 `A+A^T` positive definite，并明确该 bound 不适用于 true Helmholtz。

**Prompt**

```text
请审查论文第 5 章 GMRES 收敛率估计是否表述稳妥。

论文写：若 A+A^T 正定，则 GMRES residual 有一个正定型收敛率估计。随后说明该估计只适用于 Dirichlet Poisson、modified Helmholtz 等正定或强正定情形；对 Neumann Poisson 需处理 zero mode；对 true Helmholtz 矩阵不定，该正定型收敛界不适用，只能从谱分布/多项式逼近角度解释。

请判断：
1. 该定理作为背景是否数学上正确？
2. 对 SPD 矩阵，论文把 lambda_min/lambda_max 与 A+A^T 或 A 的特征值联系起来是否需要更精确？
3. 是否需要在第 5 章明确说“本文实验中的 true Helmholtz GMRES 现象不由该正定型 bound 解释”？
4. 是否应把该 theorem 改成 remark/background 而不是核心定理？

请给 verdict / safer theorem or remark wording。
```

---

## Experiment Evidence

## Q10. exp06 accuracy-cost 是否足以支撑“FFT 快速求解器与迭代法对比研究”？

**需要验证的内容**

exp06 比较：

- FA/CR/FACR-like：五点二阶 direct solvers。
- FFT9：Dirichlet 九点紧致四阶 direct solver。
- GMRES30：无预处理 restarted GMRES(30)，作用于五点二阶离散。
- 问题：光滑 homogeneous Dirichlet manufactured solution，`sigma=0` 与 `sigma=10`。
- timing scope：direct solver call 包含内部 RHS/transform setup；GMRES timing 排除 matrix/RHS setup。

**Prompt**

```text
请判断 exp06 accuracy-cost 实验是否足以支撑论文标题中的“FFT 快速求解器与迭代法对比研究”这一有限结论。

实验设置：
1. 光滑 Dirichlet manufactured solution；
2. sigma=0 Poisson 与 sigma=10 modified Helmholtz；
3. FA/CR/FACR-like 和 GMRES30 都作用于五点二阶离散；
4. FFT9 作用于 Dirichlet 九点紧致四阶离散；
5. 图中横轴是 median solve time，纵轴是 L_infty error，二者 log scale；
6. timing scope 不是严格 kernel benchmark：direct solvers 是 solver call，GMRES core solve 排除 matrix/RHS setup。

论文结论写得比较保守：
“在规则 Dirichlet 光滑问题和当前实现条件下，FFT9 以与 FFT 直接法相同量级的复杂度获得显著更低误差；GMRES30 是无预处理基准，结果不代表所有预处理 Krylov 方法。”

请判断：
1. 这个实验是否支持上述保守结论？
2. 是否可以说 FFT9 “更优”，还是必须写成“在当前设置下具有更好的 error-time placement”？
3. 由于 FFT9 和 GMRES30 解的是不同离散系统，是否需要避免“同一线性系统求解器性能比较”的说法？
4. timing scope 已经说明后，是否还存在明显公平性风险？

请输出 supported / partially supported / unsupported，并给最安全的 thesis wording。
```

---

## Q11. exp07 risk map 是否容易被误读为条件数图？

**需要验证的内容**

exp07 图 10 绘制：

```text
R_{p,q}=-log10 |d_{p,q}|
```

上排是低频窗口 risk map，下排是 sorted smallest `|d|`。论文反复强调：

- risk map 不是解图。
- risk map 不是误差图。
- risk map 不是 condition number plot。
- `R<0` 仅表示 `|d|>1`。

**Prompt**

```text
请审查 exp07 的 small-denominator risk map 解释是否数学上稳妥。

论文对五点 Dirichlet 频域分母

    d_{p,q}=lambda_{p,q}+sigma

绘制

    R_{p,q}=-log10 |d_{p,q}|.

图 10 上排只显示低频窗口 p,q<=12 的 risk map，下排显示最小 20 个 |d_{p,q}| 的排序曲线。论文解释：

modified: min|d|≈1.197e2，风险低；
true-away: 分母已正负分裂，但 min|d|≈1.480e1，sign-changing 不等于 near-resonance；
true-near: (2,3)/(3,2) 前两个最小分母降至 1e-2。

论文还强调：risk map 只是 small-denominator visualization，不等同于 cond_2(A)；若有 d=0 contour，也只是插值符号分界，不代表离散零分母模态。

请判断：
1. R=-log10|d| 作为 small-denominator risk visualization 是否合理？
2. `R<0` 的解释是否正确？
3. true-away “sign-changing but not near-resonance” 的说法是否正确？
4. 该图与 condition check 图的关系应如何表述，才能避免把 risk map 说成条件数图？

请给 verdict / safer caption or explanatory sentence。
```

---

## Q12. exp05 near-resonance 模态放大与 GMRES 停滞的因果表述是否过强？

**需要验证的内容**

exp05 支持：

```text
u_hat_{p,q}=f_hat_{p,q}/(lambda^h_{p,q}-kappa^2)
```

当：

```text
kappa^2=lambda_target^h+delta
```

目标模态按 `1/|delta|` 放大。GMRES(30) 在当前设置中 capped/not-converged。

需要判断是否可以写“小分母导致 GMRES 停滞”，还是应写“consistent with / helps explain”。

**Prompt**

```text
请审查 exp05 true Helmholtz near-resonance 的结论边界。

实验在 Dirichlet 五点离散下选择目标模态，例如 (2,3)/(3,2)，设

    kappa^2 = lambda_target^h + delta,
    sigma = -kappa^2.

频域公式为

    u_hat_{p,q}=f_hat_{p,q}/(lambda^h_{p,q}-kappa^2).

因此目标模态分母为 -delta，目标模态幅值应按 1/|delta| 放大。实验用 Gaussian RHS，显示 target modal amplitude 与 C/|delta| 一致，target subspace energy 和 shape correlation 接近 1，unpreconditioned restarted GMRES(30) 达到 capped/not-converged。

论文目前想表达：
“near-resonance 的小分母导致目标模态放大，并与无预处理 GMRES(30) 残差停滞一致；该结论限定于当前 Dirichlet 五点离散谱、当前 RHS 和无预处理 GMRES(30)，不外推到所有 Helmholtz 求解器或连续谱极限。”

请判断：
1. `1/|delta|` modal amplification 的理论解释是否正确？
2. target subspace energy / shape correlation 接近 1 是否可以解释为解场被目标模态主导？
3. “small denominator causes GMRES stagnation” 是否太强？是否应改为 “small denominator / ill-conditioning is consistent with and helps explain observed unpreconditioned GMRES stagnation”？
4. 当前限制性表述是否足够？

请给 verdict / safest wording / any missing caveat。
```

---

## Lean Proof Scope

## Q13. Lean 4 附录的证明范围是否表述准确？

**需要验证的内容**

Lean 4 当前只验证：

- `c4` 系数在少数模式上的显式表达式；
- 两个模式条件不可能同时成立；
- `c2=0` 推出参数关系；
- Poisson / Helmholtz 的代数矛盾路径。

Lean 不验证：

- PDE 全局误差估计；
- DST/DCT 谱理论；
- Python 程序正确性；
- 全文实验结论。

**Prompt**

```text
请审查我论文中 Lean 4 附录的范围表述是否准确。

Lean 4 formalization 只证明三参数九点修正族六阶不可达性的代数核心。具体包括：
1. c2=0 推出 alpha=-gamma，以及 Helmholtz 情形下 beta=gamma；
2. 代入后，c4 对两个模式给出互不相容的 gamma 条件；
3. 因此不存在统一实数 gamma，使所有模式 c4=0；
4. 源码没有 sorry/admit/axiom/unsafe，lake build 通过。

论文明确说 Lean 不覆盖：
- 完整 PDE 全局误差估计；
- DST/DCT 谱理论；
- Python 求解器实现正确性；
- 实验图表结论。

请判断：
1. 这样的 Lean scope 描述是否足够准确？
2. 论文可以说“Lean 4 辅助验证了六阶不可达性的代数核心”吗？
3. 是否应避免说“Lean 证明了 FFT9 四阶最优”？
4. 附录中 theorem table 的数学含义是否应更强调“local Fourier coefficient contradiction”？

请输出 verdict / safe wording / phrases to avoid。
```

---

## Thesis Conclusion Wording

## Q14. 摘要、引言、结论中“理论完备性”和“贡献点”是否过强？

**需要验证的内容**

论文贡献包括：

- 实现并验证 FFT direct solver family；
- Dirichlet FFT9 四阶；
- Neumann/mixed 五点边界处理；
- accuracy-cost comparison；
- true Helmholtz small-denominator / near-resonance；
- Lean 4 代数核心验证。

需要判断是否可以说“PDE 理论完备”，还是应说“离散化与数值验证闭环完备”。

**Prompt**

```text
请审查论文摘要/引言/结论中关于贡献和理论完备性的表述边界。

论文实际完成：
1. 在矩形区域均匀网格上研究 (-Delta+sigma)u=f；
2. 五点二阶 FA/CR/FACR-like；
3. Dirichlet FFT9 九点紧致四阶；
4. Neumann/mixed 五点 ghost-point + DCT-I/DST-I；
5. 与无预处理 GMRES(30) 的 accuracy-cost comparison；
6. true Helmholtz near-resonance 的离散谱小分母、模态放大和 GMRES residual history；
7. Lean 4 只验证三参数九点修正族六阶不可达性的代数核心。

论文不做：
- 一般区域或变系数 PDE；
- FFT9 Neumann/mixed 四阶；
- 预处理 Krylov 系统研究；
- 连续谱极限理论；
- 完整 PDE Sobolev 全局误差证明的形式化。

请判断：
1. 摘要/结论中是否可以说“理论分析、实现和实验形成闭环”？
2. 是否应避免说“PDE 相关理论完全完备”？
3. 最安全的总贡献表述应是什么？
4. 哪些词需要避免，例如 robust, optimal, proves, all Helmholtz, all Krylov？

请给 verdict / safe abstract-style wording / words or claims to avoid。
```

---

## Paper-Code Consistency

## Q15. 论文 FFT9 公式与代码反号实现是否需要在正文脚注/备注中进一步说明？

**需要验证的内容**

论文公式用：

```text
(-L_h + sigma R_h)u = R_h f
D_raw=-lambda_L+sigma lambda_R
```

代码为了求解等价系统使用：

```text
(L_h - sigma R_h)u = -R_h f
denom=lambda_L-sigma lambda_R
```

现在代码注释已说明，但论文正文主要采用原始系统。

**Prompt**

```text
请判断论文是否需要在 FFT9 方法章节额外加入 implementation note，说明代码使用整体乘以 -1 的等价系统。

论文数学公式统一采用：

    (-L_h + sigma R_h)u = R_h f,
    D_raw=-lambda_L+sigma lambda_R.

但 Python 代码 `fft9_helmholtz` 实现中使用等价系统：

    (L_h - sigma R_h)u = -R_h f,
    denom=lambda_L-sigma lambda_R.

同时 RHS 与边界修正也整体反号，所以数值解等价。

请判断：
1. 只在代码注释中说明是否足够？
2. 论文正文是否应加一条 remark：“implementation may multiply both sides by -1; this changes denominator and RHS signs simultaneously but not the solution”？
3. 如果加，应该放在第 3 章 FFT9 算法后，还是附录/代码说明中？

请输出 verdict / recommended placement / one-sentence implementation note。
```

---

## Recommended Asking Order

如果时间有限，建议按以下顺序问：

1. Q2 + Q3 + Q4：FFT9 公式、四阶展开、边界修正。这是最容易被评审盯住的核心推导。
2. Q5 + Q13：六阶不可达性与 Lean 4 范围。确认“证明了什么”和“没证明什么”。
3. Q8 + Q9：GMRES 理论和复杂度。避免第 5 章被抓符号/量纲问题。
4. Q10 + Q12：实验三和实验五的结论边界。确认实验支撑力度。
5. Q14：摘要/结论总口径。最后收束语言。
