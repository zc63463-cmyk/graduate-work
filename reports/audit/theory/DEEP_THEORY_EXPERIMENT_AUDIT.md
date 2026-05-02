# 论文理论正确性与实验正确性深度审查报告

审查日期：2026-04-30
项目根目录：`C:\Users\20564\Desktop\Graduate\论文收集`

本报告只做审查，不修改论文正文、实验脚本、CSV 或图像。审查重点是理论符号、离散格式、代码实现、实验数据、图表证据和结论边界是否一致。

## 0. 本轮修复状态

根据本报告第 3 节列出的轻量问题，已完成以下修复：

- A-01：已将第 5 章 GMRES 重启策略中的 Arnoldi/Krylov 存储从 `O(mN)` 修正为 `O(mN^2)`，并补充以总未知量 `M=N^2` 表示时为 `O(mM)`。
- B-01：已将 `README.md`、`FINAL_POLISH_REPORT.md`、`FINAL_REVISION_REPORT.md`、`RECENT_EXPERIMENT_UPDATE_REPORT.md` 中的 PDF 页数同步为当前编译结果 `66 pages`。
- C-01：已修正 `verify_fft9_expansion.py` 最后结论打印中 `-\lambda_L` 主项符号的说明文字。
- C-02：已收紧 `helmholtz_solver.py` 顶层 docstring，明确 FFT9 compact fourth-order 路径仅实现并验证 Dirichlet，非 Dirichlet 请求会 warning 并回退到二阶 FA。

这些修复均为文本/说明层面，不改变数值算法、实验参数、CSV 或图像。

## 1. 总体结论

当前论文主线整体已经形成较稳的“理论--代码--实验”闭环：

- 统一方程 `(-\Delta+\sigma)u=f` 在引言、数学基础、FFT 直接法、Helmholtz、实验和结论中基本一致。
- 五点 Dirichlet、Neumann ghost-point + DCT-I、FFT9 Dirichlet 紧致格式、true Helmholtz 小分母、GMRES absolute residual 等关键口径与代码和实验大体匹配。
- 第 6 章的五个核心实验结构清楚，`exp00` 已正确降级为 implementation validation。
- `exp06` 兑现了题目中的 FFT direct solvers vs unpreconditioned GMRES(30) accuracy-cost comparison，但 timing scope 已正确收窄为当前实现趋势。
- `exp07+exp04+exp05` 对 true Helmholtz small denominator、modal amplification 和 GMRES stagnation 的证据链是连贯的，且没有明显外推到所有 Krylov 方法或连续谱极限。

需要优先处理的不是实验数据，而是少数文本/报告一致性问题。最重要的一处是第 5 章仍有旧的 GMRES memory 表述 `O(mN)`，与当前 `N` 表示单方向内点数、未知量为 `N^2` 的全文章法不一致。

## 2. 审查与验证命令

已执行：

```powershell
cd code/python
python -m pytest -q
```

结果：

- `103 passed, 2 skipped`
- Neumann compatibility 的 2 个 warning 是测试中预期触发的不兼容 RHS 诊断。

已执行：

```powershell
python code/python/verify_fft9_expansion.py
```

结果：

- 有效特征值 `-\lambda_L/\lambda_R` 的 `h^2` 项为 0。
- FFT9 四阶展开验证通过。
- 注意：该脚本最后的说明文字有一处符号打印陈旧，详见 C-01。

已执行：

```powershell
cd code/lean4_formalization
lake build
```

结果：

- build completed successfully。
- 无 `sorry` / `admit` / `axiom` / `unsafe`。
- 仅有 unused simp argument warnings，以及 `proofwidgets` package local changes warning。

已检查 LaTeX log：

- `Output written on main.pdf (66 pages).`
- 无 fatal error。
- 无 undefined reference / citation。
- 无 multiply-defined label。
- 无 overfull hbox。
- 仅有一个非致命 SimSun 字体形状 warning。

## 3. A/B/C/D 问题清单

| ID | 严重度 | 类型 | 位置 | 问题 | 最小修复 | 是否需要重跑实验 |
|---|---|---|---|---|---|---|
| A-01 | A | 理论复杂度表述 | `thesis/chapters/5_gmres.tex:179-180` | 仍写 Arnoldi 存储 `m` 个 `N` 维向量，内存 `O(mN)`。但全文当前 `N` 表示单方向内点数，总未知量为 `N^2`，应为 `m` 个 `N^2` 维向量，内存 `O(mN^2)`。 | 改为“由于未知量总数为 `N^2`，Arnoldi 过程需存储 `m` 个 `N^2` 维 Krylov 向量，内存为 `O(mN^2)`；若以总未知量 `M=N^2` 表示，则为 `O(mM)`。” | 否 |
| B-01 | B | 交付文档陈旧 | `README.md:6`, `FINAL_POLISH_REPORT.md:7,98`, `FINAL_REVISION_REPORT.md:156,177`, `RECENT_EXPERIMENT_UPDATE_REPORT.md:161` | 多处报告仍写 `62 pages`，但当前 PDF 为 66 pages。不会影响数学正确性，但会影响最终交付一致性。 | 统一改为当前页数 `66 pages`，或写“以当前编译结果为准”。 | 否 |
| B-02 | B | 审查报告陈旧 | `EXPERIMENT_AUDIT_REPORT.md`, `CLAIM_DOWNGRADE_SUGGESTIONS.md` | 旧审查报告中仍有早期 exp05/exp04 口径，例如旧 near-resonance CSV 行、旧图注建议等；它们不影响论文正文，但若作为最终材料会混淆当前版本。 | 标注为历史审查报告，或更新摘要说明已被后续 Patch/Phase 覆盖。 | 否 |
| C-01 | C | 脚本输出文字陈旧 | `code/python/verify_fft9_expansion.py:139-141` | 脚本计算结果正确，但结尾说明打印 `-lambda_L = -(xi^2+eta^2) ...`，符号应为 `-lambda_L = +(xi^2+eta^2) - ...`。这是输出说明文字问题，不是计算问题。 | 改 print 字符串，不改符号计算。 | 否 |
| C-02 | C | 代码顶层 docstring 可更精确 | `code/python/helmholtz_solver.py:16-22,36-37` | 顶层说明写 solvers support Dirichlet/Neumann/mixed，且 FFT9 标成 `DST/DCT`，容易让读者误以为 FFT9 也支持 Neumann/mixed 四阶；但函数级代码实际会对非 Dirichlet FFT9 warning + fallback。 | 顶层 docstring 增一句“FFT9 compact fourth-order path is Dirichlet-only; non-Dirichlet requests fall back to FA.” | 否 |

## 4. 理论正确性审查

### 4.1 统一 PDE 与符号

结论：一致。

- `thesis/chapters/1_introduction.tex:12-20` 定义 `(-\nabla^2+\sigma)u=f`，并区分 `sigma=0`, `+kappa^2`, `-kappa^2`。
- `thesis/chapters/2_math_preliminary.tex:11-19` 与引言一致。
- `code/python/helmholtz_solver.py:6-14` 与论文一致。

没有发现 modified Helmholtz 与 true Helmholtz 的符号互换问题。

### 4.2 五点 Dirichlet 离散与截断误差

结论：一致。

- `thesis/chapters/2_math_preliminary.tex:68-73` 明确定义 residual 方向 `tau_h(u)=L_h^{(5)}u+\sigma u-f`。
- 在该方向下，五点正 Laplacian 主项为 `-h^2/12(u_xxxx+u_yyyy)`，与当前正文一致。
- `thesis/chapters/2_math_preliminary.tex:613-623` 非齐次 Dirichlet 边界项移至 RHS 后为加号。
- `code/python/helmholtz_solver.py:388-413` Dirichlet RHS 使用 `F += bc/h^2`。
- `code/python/experiments/exp02_nonhom_bc.py:1-9` 使用 `u=sin(2*pi*x)sin(3*pi*y)+1`，`g_D=1`，RHS 包含常数项的 `sigma` 贡献。

### 4.3 FFT9 分母与边界修正

结论：论文公式与代码实现等价，但代码用乘以 `-1` 的等价系统实现，已在代码注释中说明。

- 论文采用原始系统 `(-L_h+\sigma R_h)u=R_h f`。
- `thesis/chapters/3_fft_direct.tex:257-280` 使用原始除数 `D_raw=-lambda_L+sigma lambda_R`。
- `thesis/chapters/3_fft_direct.tex:377-380` 算法步骤写 `+L_IB g_B - sigma R_IB g_B`，与论文系统一致。
- `code/python/helmholtz_solver.py:910-923` 明确说明实现求解等价的 `(L_h-\sigma R_h)u=-R_h f`，因此边界修正符号整体反向。
- `code/python/helmholtz_solver.py:1046-1054` 使用 `denom=lam_L-sigma*lam_R4`，与等价系统一致。
- `code/python/tests/test_04_fft9_vs_spsolve.py` 通过测试；`pytest` 全部通过。

风险仅在 `verify_fft9_expansion.py` 的最后说明文字符号陈旧，非公式/实现错误。

### 4.4 Neumann / mixed 边界

结论：核心理论与实现一致。

- `thesis/chapters/2_math_preliminary.tex:650-663` 说明 ghost-point Neumann 矩阵非对称，并以三个物理节点给出 3x3 例子，当前表述正确。
- `thesis/chapters/2_math_preliminary.tex:665-670` 说明通过 `D^{-1}GD` 对称化并用 DCT-I 对角化。
- `code/python/helmholtz_solver.py:176-184` 与 DCT-I 对称化尺度一致。
- `code/python/tests/test_01_eigenvalues.py`、`test_05_neumann_compatibility.py` 覆盖 DCT-I 特征值、正交性、compatibility 和 zero mode 行为。

论文没有声称 FFT9 已支持 Neumann/mixed 四阶，这是正确的边界。

### 4.5 GMRES 理论与停止准则

结论：停止准则一致；复杂度文字有一处 A 级旧表述。

- `code/python/gmres_solver.py:53-54` 明确 `tol` 是 absolute residual norm。
- `code/python/gmres_solver.py:181-182` 用 `residual < tol` 停止。
- `code/python/experiments/exp04_modified_vs_true.py` 与 `exp05_true_helmholtz_resonance.py` 的 CSV 中均写 `tol_type=absolute_residual`。
- `thesis/chapters/6_experiments.tex:480`、`618`、`622-628` 已清楚区分 absolute residual 和 relative residual。

原审查发现 `thesis/chapters/5_gmres.tex:179-180` 仍写 `O(mN)` memory；第 3 章表格已改为 `O(mN^2)`，所以第 5 章需要同步。该项已在本轮修复，当前正文改为 `O(mN^2)`，并补充总未知量记号 `M=N^2` 下为 `O(mM)`。

### 4.6 FACR-like 复杂度

结论：论文结论边界正确。

- `thesis/chapters/3_fft_direct.tex:419-420` 区分经典 FACR `O(N^2 log log N)` 与本文实现 `O(N^2 log N)`。
- `thesis/chapters/6_experiments.tex:39-40`, `314`, `667` 也保持谨慎。
- `code/python/cyclic_reduction.py:6-19` 已标为 historical prototype，并说明当前实现仍是 `O(N^2 log N)`。

没有发现当前论文声称 FACR-like 实测达到经典 `O(N^2 log log N)` 的问题。

### 4.7 Lean 4 附录

结论：范围表述正确。

- `thesis/chapters/appendix_lean4.tex:7-9` 明确 Lean 只验证代数核心，不验证完整 PDE 全局误差。
- `thesis/chapters/appendix_lean4.tex:75-80` 列出未覆盖 PDE 全局估计、DST/DCT 谱理论和 Python 程序正确性。
- `lake build` 成功，无 `sorry/admit/axiom/unsafe`。

Lean warning 是 unused simp arguments，不影响证明状态。

## 5. 实验正确性与证据链审查

| 章节 | 实验定位 | 脚本/CSV/图表 | 证据强度 | 审查结论 |
|---|---|---|---|---|
| 6.1 | implementation validation | `exp00_fft_vs_sparse.py`, `exp00_fft_vs_sparse.csv`, 表 `tab:exp00_sparse_consistency` | strong | 正确降级为离散系统一致性测试，没有作为连续 PDE 精度贡献。 |
| 6.2 | Dirichlet 收敛 | `exp01_convergence.py`, `exp02_nonhom_bc.py`, `exp01/02` CSV 与图，补充 multimode visual | strong | 五点二阶、FFT9 四阶、非齐次 Dirichlet 边界修正均有脚本/CSV/图支持；Fourier eigenfunction 已降级为 sanity check。 |
| 6.3 | Neumann/mixed 边界 | `exp03_neumann_mixed.py`, `exp03_neumann_mixed.csv`, 图 6/7, 表 `tab:exp03_neumann_mixed` | strong | 实验重点放在 boundary residual、flux residual、weighted mean、compatibility 上，未外推 FFT9。 |
| 6.4 | 精度--成本对比 | `exp06_accuracy_cost.py`, `exp06_accuracy_cost.csv`, 图 `exp06_accuracy_cost_error_time`, `exp06_time_scaling` | strong | 数据支持 FFT9 在当前光滑 Dirichlet 设置下以同量级 FFT cost 获得更低误差；timing scope 限制已写清。 |
| 6.5 | Helmholtz 谱结构与 GMRES | `exp07_spectral_denominator_maps.py`, `exp04_modified_vs_true.py`, 相关 CSV/图 | strong | risk map、sorted denominator、condition check、GMRES iterations/history 形成连贯证据链；condition equivalence 已限定 Dirichlet 五点正交系统。 |
| 6.6 | True Helmholtz near-resonance | `exp05_true_helmholtz_resonance.py`, `exp05_resonance.csv`, `exp05_multimode_resonance.csv`, history CSV/图 | strong | 数据支持 `1/|delta|` 模态放大、target subspace dominance 和无预处理 GMRES(30) capped/stagnation；结论未外推到所有 Helmholtz 求解器。 |

### 5.1 exp06 数据核对

`exp06_accuracy_cost.csv` 中 `sigma=10, n=129` 代表行与论文表格一致：

- FA: time `1.6689e-02`, error `5.5720e-07`
- CR: time `3.7654e-02`, error `5.5720e-07`
- FACR: time `1.5405e-02`, error `5.5720e-07`
- FFT9: time `3.7615e-02`, error `3.4408e-11`
- GMRES30: time `1.5724e+01`, error `5.5720e-07`

GMRES30 在 `sigma=0, n=129` 未收敛，在 `sigma=10, n=129` 接近上限但成功达到 absolute residual tolerance。论文图注“黑色叉号表示未达到绝对残差容差”与数据口径一致。

### 5.2 exp07 数据核对

`exp07_spectral_denominator_summary.csv` 支持图 10 的三个数量级判断：

- modified: `min |d| = 119.738218`，约 `1.197e2`
- true-away: `min |d| = 14.800691`，约 `1.480e1`
- true-near: `min |d| = 0.010000`，目标模态 `(2,3)/(3,2)`，`num_near_zero=2`

正文“sign-changing 不等于 near-resonance”和“risk map 不是条件数图”的限制正确。

### 5.3 exp05 数据核对

`exp05_resonance.csv` 支持表 `tab:exp05_resonance`：

- `delta=1e-1`: `solution_linf≈0.561769`, `gmres_iters=1001`, not converged
- `delta=1e-2`: `solution_linf≈5.637208`, `gmres_iters=1001`, not converged
- `delta=1e-3`: `solution_linf≈56.391614`, `gmres_iters=1001`, not converged
- `delta=1e-4`: `solution_linf≈563.939447`, `gmres_iters=1001`, not converged

`exp05_multimode_resonance.csv` 支持正文说法：

- target energy fraction 最小约 `0.996005`
- shape correlation 最小约 `0.998001`
- 三组目标集合均 `gmres_iterations=1001`, `gmres_success=False`

`exp05_resonance_gmres_history.csv` 支持“最终 absolute residual 约 `1.83e-1`，远高于 `1e-10`”的最新表述。

## 6. 图表与引用审查

第 6 章所有 `includegraphics` 指向的 `thesis/figures` 文件均存在：

- `exp01_convergence.png`
- `exp02_nonhom_bc.png`
- `exp02_temperature_field_comparison.png`
- `exp06_complex_manufactured_fields.png`
- `exp06_complex_manufactured_convergence.png`
- `exp03_neumann_mixed_summary.png`
- `exp03_neumann_mixed_fields.png`
- `exp06_accuracy_cost_error_time.png`
- `exp06_time_scaling.png`
- `exp07_spectral_denominator_heatmaps.png`
- `exp04_min_denom_vs_sigma.png`
- `exp04_spectral_indicator_vs_sigma.png`
- `exp04_condition_check.png`
- `exp04_gmres_iters_vs_sigma.png`
- `exp04_gmres_history.png`
- `exp05_near_resonance_summary.png`
- `exp05_multimode_resonance_summary.png`
- `exp05_dominant_mode_projection.png`
- `exp05_resonance_gmres_history.png`

LaTeX log 未发现 undefined reference、multiply-defined label 或 overfull hbox。

## 7. 高风险表述检查

当前论文正文中以下高风险点已基本处理得当：

- GMRES 已限定为 unpreconditioned restarted GMRES(30) 或对应实验设置，没有外推到所有 Krylov 方法。
- GMRES tolerance 已写为 absolute residual tolerance。
- `501` 与 `1001` capped convention 已区分。
- FFT9 已限定 Dirichlet。
- FACR-like 未声称当前实现达到 `O(N^2 log log N)`。
- exp07 risk map 未被称为条件数图。
- condition-number equivalence 已限定 Dirichlet 五点正交 DST-I 系统。
- near-resonance 已限定为离散谱实验，不写成连续谱极限。
- timing 已限定为当前实现趋势，不是严格 kernel benchmark。

## 8. 已执行修复与剩余建议

1. A-01 已修复：第 5 章 GMRES memory 已由 `O(mN)` 改为 `O(mN^2)`，并补充总未知量 `M=N^2` 口径。
2. B-01 已修复：README 和相关最终/更新报告的 PDF 页数已同步为当前编译结果 `66 pages`。
3. C-01 已修复：`verify_fft9_expansion.py` 的打印说明已更正为 `-\lambda_L = xi^2+eta^2 - ...`。
4. C-02 已修复：`helmholtz_solver.py` 顶层 docstring 已明确 FFT9 compact fourth-order 路径仅 Dirichlet。
5. 剩余建议：若后续把旧审查报告作为最终交付材料，可进一步在旧报告标题处标注“historical audit”，避免读者将早期建议误认为当前未处理问题。

## 9. GPT Pro 复核建议

若要让 GPT Pro 或导师做二次审查，建议把问题集中在以下 4 个点，而不是让对方泛泛重读全文：

1. 第 2 章 FFT9 推导中，`L_h≈\Delta`、`-L_h+\sigma R_h`、`D_raw=-lambda_L+sigma lambda_R` 与代码中乘以 `-1` 的实现是否完全等价？
2. 三参数九点修正族六阶不可行性的表述是否应限定为“统一常数参数、局部 Fourier symbol、该三参数族”？
3. 第 5 章 GMRES 收敛率估计是否需要进一步说明该正定型 bound 对本文 true Helmholtz 实验并不适用？
4. 第 6 章 exp06 的 accuracy-cost 表述是否已经充分避免把当前 timing 解释成严格公平 kernel benchmark？

## 10. 最终审查意见

论文理论与实验主体可以认为是“基本正确、证据链完整、边界表述谨慎”的状态。当前不需要新增实验，也不需要重跑大规模实验。本轮已修复第 5 章 GMRES memory 的 `O(mN)` 陈旧表述，并同步了页数、脚本输出说明和代码注释层面的轻量问题。
