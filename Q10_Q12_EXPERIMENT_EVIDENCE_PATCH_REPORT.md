# Q10--Q12 实验证据表述收紧补丁报告

## 1. 修改目标

本轮只收紧实验章和结论章的证据边界表述：

- Q10：将 exp06 的结论限定为当前规则 Dirichlet 光滑 manufactured solution、当前实现和当前计时口径下的 error-time placement。
- Q11：继续强调图 10 是逐模态 small-denominator risk map，而不是矩阵条件数图。
- Q12：保留 true Helmholtz near-resonance 的模态放大结论，但将 GMRES 停滞表述从强因果改为“相一致并有助于解释”。

本轮未新增实验、未重画图、未修改 CSV、未修改 Python 求解器或实验脚本。

## 2. 修改文件

- `thesis/chapters/6_experiments.tex`
- `thesis/chapters/7_conclusion.tex`
- `Q10_Q12_EXPERIMENT_EVIDENCE_PATCH_REPORT.md`

## 3. Q10：exp06 精度--成本证据边界

第 6.4 节已改为：

- exp06 是规则 Dirichlet 光滑问题上的核心基准证据；
- 比较对象是当前实现下不同离散/求解路径的 error-time placement；
- 该图不是同一线性系统上的严格 kernel-to-kernel solver benchmark；
- FA、CR、FACR-like 与 GMRES30 对应五点二阶离散系统；
- FFT9 对应 Dirichlet 九点紧致四阶离散系统；
- GMRES30 仅作为无预处理 restarted GMRES(30) 基线。

第 6.7 节和第 7 章同步收紧为：FFT9 在当前设置下具有更好的误差--时间位置，不代表预处理 GMRES、MINRES 或所有 Krylov 方法。

## 4. Q11：图 10 risk map 口径

图 10 附近已保留并强化：

- \(R_{p,q}=-\log_{10}|d_{p,q}|\) 是逐模态 small-denominator risk visualization；
- \(R_{p,q}<0\) 仅表示 \(|d_{p,q}|>1\)，即小分母风险较低；
- true-away 的 sign-changing 不等于 near-resonance；
- 图 10 不画解、不画误差，也不是条件数图；
- 只有当前 Dirichlet 五点正交 DST-I 系统中，\(\max |d|/\min |d|\) 才与 \(\mathrm{cond}_2(A)\) 一致；
- 条件数关系应结合图 12 的 condition check 理解。

## 5. Q12：near-resonance 与 GMRES 停滞因果软化

第 6.6 节已明确：

- 当 \(\kappa^2=\lambda_{\mathrm{target}}^h+\delta\) 时，目标模态分母为 \(-\delta\)；
- 只要右端项在目标模态上的投影非零，目标模态幅值按 \(1/|\delta|\) 放大；
- target subspace energy 与 shape correlation 接近 1 的含义是：当前 Gaussian RHS 与 detuning 区间内，解场主要由目标模态子空间主导；
- 无预处理 GMRES(30) 达到迭代上限这一现象，与离散小分母造成的模态放大和病态性相一致，并有助于解释当前实验中观察到的残差停滞；
- 该观察不外推为所有 Helmholtz 求解器、预处理 Krylov 方法或连续谱极限结论。

第 7 章同步改为：在本文无预处理 GMRES(30) 基准中，病态谱结构与观察到的残差停滞相一致。

## 6. 未修改范围

- 未新增实验；
- 未运行实验脚本；
- 未重绘任何图；
- 未修改 CSV；
- 未修改 Python 求解器；
- 未修改 HTML/PPT 讲解材料。

## 7. 验证结果

静态搜索已确认：

- 未残留裸写 `FFT9 更优`、`FFT9 优于 GMRES`；
- 未残留 `小分母导致 GMRES 停滞`、`造成当前基准迭代法停滞`、`谱原因`；
- 未残留 `第 5.1节的正定型背景界`；
- `所有 Helmholtz 求解器`、`预处理 Krylov` 均出现在限制或未来工作语境中；
- `条件数图` 仅出现在否定说明中。

LaTeX 验证：

- `xelatex -> bibtex -> xelatex -> xelatex`：通过；
- `main.log` 扫描：
  - no fatal error；
  - no undefined reference/citation；
  - no multiply-defined label；
  - no overfull hbox；
- 最终 PDF：`thesis/main.pdf`，71 页。

## 8. 结论

Q10--Q12 后，实验章对 exp06、exp07 和 exp05 的证据边界更加清晰：

- exp06 支持当前设置下 FFT9 更好的 error-time placement，而不是普遍 solver 优劣；
- exp07 的 risk map 只可视化逐模态小分母风险，不等同于条件数图；
- exp05 严格支持 near-resonance 模态放大，GMRES 停滞则表述为与小分母病态结构相一致并有助于解释当前无预处理 GMRES(30) 的观测结果。
