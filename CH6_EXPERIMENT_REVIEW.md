# 第六章实验审查报告

## 0. 审查说明

本报告为 Phase 1 只读审查结果。审查对象包括：

- `thesis/chapters/6_experiments.tex`
- `thesis/chapters/2_math_preliminary.tex`
- `thesis/chapters/3_fft_direct.tex`
- `thesis/chapters/4_helmholtz.tex`
- `thesis/chapters/5_gmres.tex`
- `EXPERIMENT_ASSET_INDEX.md`
- `RECENT_EXPERIMENT_UPDATE_REPORT.md`
- `RISKMAP_PATCH_REPORT_V2.md`
- `PATCH4_REPORT.md`
- `PATCH5_REPORT.md`

说明：输入清单中写有 `thesis/chapters/4_helmholtz_fft.tex`，当前仓库实际文件为 `thesis/chapters/4_helmholtz.tex`，本报告按实际源文件审查。

本轮未修改论文正文、实验脚本、CSV、PNG 或 PDF；仅生成本报告。

## 1. 总体评价

第 6 章已经形成较完整的“理论--实现--实验”闭环。当前结构从早期补丁式扩展整理为五个核心实验，并将 exp00 降级为 implementation validation，整体排布是合理的。

### 1.1 五个核心实验结构

当前第 6 章结构为：

1. `6.1 实验设置、可复现性与实现校验`
2. `6.2 实验一：Dirichlet 离散格式收敛性验证`
3. `6.3 实验二：Neumann 与 mixed 边界处理`
4. `6.4 实验三：精度--成本对比与复杂度实证`
5. `6.5 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为`
6. `6.6 实验五：True Helmholtz 近共振模态放大`
7. `6.7 实验总结`

五个核心实验覆盖了论文主线所需的关键证据：Dirichlet 收敛、边界处理、精度--成本对比、Helmholtz 谱结构与 GMRES 行为、true Helmholtz 近共振机制。

### 1.2 exp00 降级是否合理

合理。第 6.1 节明确写明 exp00 只用于确认 FFT 类求解路径与 sparse direct 解求解同一个离散线性系统，不作为连续 PDE 数值精度贡献。该定位避免了把“实现一致性”误说成“收敛性证明”。

### 1.3 exp01/exp02 合并是否合理

合理。原 exp01 的 Dirichlet 收敛验证和原 exp02 的非齐次 Dirichlet 边界验证都服务于 Dirichlet 离散格式正确性。当前第 6.2 节把非齐次 Dirichlet 作为“边界子情形 / visual sanity check”，既保留了用户希望保留的两组图，又避免重复占用核心实验编号。

### 1.4 exp06 是否兑现“对比研究”

基本兑现。第 6.4 节将 FA、CR、FACR-like、FFT9 与无预处理 GMRES(30) 放在同一 manufactured Dirichlet 问题上做误差--时间对比，并明确记录 timing scope。该实验是论文标题“FFT 快速求解器与迭代法对比研究”的核心实验证据。

需要注意的是，正文已经正确说明这不是严格 kernel-to-kernel benchmark；direct solver timing 包含内部 RHS/transform setup，而 GMRES timing 排除了 matrix/RHS setup。因此该实验支持“当前实现下的整体基准趋势”，不应扩展为所有实现或所有预处理 Krylov 方法的性能结论。

### 1.5 exp07/exp04/exp05 是否支撑 true Helmholtz 小分母与 GMRES 困难

支撑较强。第 6.5 节用 exp07 的低频小分母风险图与排序图解释 modified/true Helmholtz 的分母结构；用 exp04 的谱分母指标、条件数校验和 GMRES residual history 展示 Gaussian RHS 下的迭代困难；第 6.6 节进一步用 exp05 的模态投影公式和多模态扫描验证近共振时的 `1/|delta|` 放大规律。

整体证据链为：

`d_{p,q} = lambda_{p,q} + sigma`
到
small-denominator risk
到
modal amplification
到
unpreconditioned GMRES stagnation。

这个逻辑闭环成立，但结论必须继续限定在当前 Dirichlet 五点离散谱、Gaussian RHS 设置与无预处理 GMRES(30) 基线范围内。

## 2. 理论正确性审查

| 审查项 | 当前状态 | 证据位置 | 判断 |
|---|---|---|---|
| 统一方程 `(-nabla^2 + sigma)u=f` | 已贯穿 | `2_math_preliminary.tex`, `4_helmholtz.tex`, `6_experiments.tex` | 正确 |
| modified Helmholtz `sigma=+kappa^2` | 明确解释为 SPD / 分母全正 | 第 2 章、第 4 章、第 6.5 节 | 正确 |
| true Helmholtz `sigma=-kappa^2` | 表述为依谱而定、可能不定、近共振 | 第 4 章谱分类，第 6.5/6.6 节 | 正确 |
| 五点 Dirichlet 分母 `d_{p,q}=lambda_{p,q}+sigma` | 与 exp07/exp04 图表一致 | 第 2、3、4、6 章 | 正确 |
| FFT9 原始分母 `D_raw=-lambda_L+sigma lambda_R` | 第 3 章采用 raw denominator，并说明与 `R_h f` 匹配 | `3_fft_direct.tex` | 正确 |
| Dirichlet 五点下 spectral denominator indicator 与 `cond_2(A)` | 正文明确限定为正交 DST-I 对角化系统 | 第 6.5 节、PATCH4 | 正确 |
| Neumann/mixed 是否外推 FFT9 | 第 6.3 节明确“不声称 FFT9 已支持 Neumann 或 mixed 四阶格式” | `6_experiments.tex` | 正确 |
| GMRES 是否限定为无预处理基线 | 第 6.1、6.4、6.5、6.6 节均有限定 | `6_experiments.tex`, `5_gmres.tex` | 正确 |
| GMRES tolerance 是否写成 absolute residual | 第 6.1、6.4、6.5/6.6 已说明绝对残差容差 | `6_experiments.tex` | 正确 |
| FACR-like 复杂度 | 第 3 章与第 6.4 节均说明当前实现为 `O(N^2 log N)`，经典 `O(N^2 log log N)` 只是理论背景 | `3_fft_direct.tex`, `6_experiments.tex` | 正确 |

### 2.1 需要留意的理论表述

1. `5_gmres.tex` 中关于 GMRES(30) 大量迭代的“根本原因”表述略强。谱分母结构确实是重要原因，但 restart、RHS 投影、停止准则和实现细节也会影响结果。建议后续成稿时改为“重要原因之一”或“本文实验中观察到的主要原因”。

2. `6_experiments.tex` 第 6.6 节有“Patch5 补充扫描”这样的过程性表述。数学上无误，但正式论文正文中建议改成“补充的多模态扫描”或“进一步的多模态扫描”，避免读者看到内部 patch 术语。

3. `RECENT_EXPERIMENT_UPDATE_REPORT.md` 中 exp05 图表说明仍有一段描述“amplification ratio / dominant projection check”的旧版本图形，而当前源码和第 6 章已回到 target energy fraction 与 shape correlation 版本。若该报告作为提交材料，需要同步；若仅作过程记录，则不影响论文。

## 3. 实验排布审查

### 3.1 实验一：Dirichlet 离散格式收敛性验证

作用清晰。该实验用多项式 manufactured solution 作为主要收敛性证据，说明五点 FA/CR/FACR-like 为二阶、FFT9 为四阶。Fourier eigenfunction track 被降级为 frequency-formula sanity check，这是合理的，因为单一特征模态过于特殊，不宜作为泛化收敛性主证据。

非齐次 Dirichlet 子情形放在 6.2 内也合理。它主要回答边界值贡献移至右端项时，符号和模板处理是否正确；不再作为独立核心实验，避免重复。

### 3.2 实验二：Neumann 与 mixed 边界处理

该实验有独立价值，不只是重复二阶收敛。它验证：

- ghost-point + DCT-I/DST-I 对称化路径；
- Neumann flux residual；
- mixed Dirichlet boundary residual；
- Neumann Poisson compatibility 与 weighted mean；
- 零模态处理。

当前图表已从六子图压缩为三子图，叙事更清楚：收敛、边界约束、零模态。正文也明确不外推 FFT9 到 Neumann/mixed。

### 3.3 实验三：精度--成本对比与复杂度实证

这是目前最能支撑论文题目“FFT 快速求解器与迭代法对比研究”的实验。它把直接 FFT 求解器、FFT9 和无预处理 GMRES(30) 放到同一误差--时间框架下比较，并保留了 timing scope 说明。

当前定位稳健：它比较的是当前实现下的整体基准趋势，不是严格内核级 benchmark，也不代表预处理 Krylov 方法。

### 3.4 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为

逻辑顺序合理：

1. 先用 exp07 展示频域分母风险结构；
2. 再用 exp04 展示谱分母指标随 sigma 的变化；
3. 再用 condition check 说明当前 Dirichlet 五点系统下 indicator 等于 `cond_2(A)`；
4. 最后用 GMRES iterations 和 residual history 说明迭代行为。

这使第 4 章理论的“modified 全正、true 可符号分裂、near resonance 小分母”得到了图像和数值双重支撑。

### 3.5 实验五：True Helmholtz 近共振模态放大

该实验收束得比较好。它不只报告 capped/not-capped，而是用离散谱展开

```text
u_hat[p,q] = f_hat[p,q] / (lambda[p,q] - kappa^2)
```

验证目标模态幅值按 `1/|delta|` 放大，并用 dominant projection 与 residual history 解释 near-resonance 解场和 GMRES 停滞。

正文已经加入稳健表述：在当前 Gaussian RHS 下，target subspace energy 和 shape correlation 在扫描区间内已接近 1，因此图 (b)(c) 主要说明目标模态已主导解场，而不是展示显著趋势。这是正确且必要的降级。

## 4. 图表解释充分性审查

| 图表 | 当前 caption 是否足够 | 正文读图说明 | 审查意见 |
|---|---|---|---|
| `fig:exp01_convergence` | 足够 | 足够 | 明确区分 Fourier track 与 polynomial track，定位稳健。 |
| `fig:exp02_nonhom_bc` | 基本足够 | 足够 | 已说明是边界符号和模板处理检查，不是独立核心实验。 |
| `fig:exp02_temperature_field_comparison` | 足够 | 足够 | 作为 physical-space boundary sanity check 合理。 |
| `fig:exp06_complex_fields` | 足够 | 足够 | 说明为 visual sanity check，避免把图像当作主要收敛证据。 |
| `fig:exp06_complex_convergence` | 足够 | 基本足够 | 作为补充收敛图可保留；不应压过主 exp01 证据。 |
| `fig:exp03_neumann_mixed` | 足够 | 足够 | 三子图叙事清楚：收敛、边界约束、零模态。 |
| `fig:exp03_neumann_mixed_fields` | 足够 | 足够 | 正确限定为解场与误差分布展示。 |
| `fig:exp06_accuracy_cost` | 足够 | 足够 | caption 明确 log scale、算法 marker、GMRES not converged；正文解释 timing scope。 |
| `fig:exp06_time_scaling` | 足够 | 足够 | FACR-like 未声称 `O(N^2 log log N)`，风险较低。 |
| `fig:exp07_spectral_denominator_maps` | 足够 | 很充分 | 已解释 `R<0`、低频 zoom、排序图、插值符号分界、非条件数图。 |
| `fig:exp04_spectral` | 基本足够 | 足够 | 图表与表格 n 不同的问题已说明：图 n=129，表 n=33。 |
| `fig:exp04_condition_check` | 足够 | 足够 | 可选在正文加最大 relative difference，但不是必要。 |
| `fig:exp04_gmres` | 足够 | 足够 | 501 capped convention 和绝对残差容差已说明。 |
| `fig:exp04_gmres_history` | 足够 | 足够 | 能补充最终迭代数看不到的停滞过程。 |
| `fig:true_helmholtz_near_resonance_summary` | 基本足够 | 基本足够 | 作为原 near-resonance summary 合格；更强解释由后续 multimode 和 history 图承担。 |
| `fig:exp05_multimode_resonance_summary` | 足够 | 很充分 | 图前后已解释四个 panel 的作用，尤其说明 (b)(c) near 1 的含义。 |
| `fig:exp05_dominant_mode_projection` | 基本足够 | 足够 | 图注简洁；正文已说明目标模态主导。 |
| `fig:exp05_resonance_gmres_history` | 足够 | 足够 | 明确无预处理 GMRES(30)、绝对残差容差和停滞过程。 |

### 4.1 重点图表判断

- 图 10 小分母风险图：当前版本是本章中解释最充分的图之一。低频 zoom + sorted denominator 比全域 heatmap 更适合表达 small-denominator risk。正文已经防止了“risk map = condition number”的误读。
- exp04 condition check：证据充分。`exp04_condition_check.csv` 中最大 relative difference 约为 `2.91e-13`，足以支持当前 Dirichlet 五点系统下的等价性。
- exp04 GMRES residual history：必要且有解释力。它弥补了仅看 iteration count 的不足。
- exp05 multimode resonance summary：当前 caption 与正文读图说明已经能解释为什么 (b)(c) 视觉变化小但仍有意义。
- exp05 resonance GMRES history：能支撑“near-resonance 下无预处理 GMRES(30) 停滞”的限定结论。
- exp06 error-time 图：对论文题目支撑强。timing scope 已说明，比较口径透明。
- exp03 Neumann/mixed summary：去冗余后清晰，不再用多条重合线分散读者注意力。

## 5. 高风险表述检查

| 风险点 | 当前处理 | 残余风险 |
|---|---|---|
| FACR-like 复杂度 | 已说明当前实现为 `O(N^2 log N)`，经典 FACR `O(N^2 log log N)` 只作理论背景 | 低 |
| GMRES 结论外推 | 第 6 章多处限定为无预处理 GMRES(30) / 当前设置 | 低 |
| FFT9 适用边界 | 明确限定为 Dirichlet，Neumann/mixed 不声称四阶 FFT9 | 低 |
| risk map 被误读为条件数图 | 图 10 正文和 caption 均说明不是条件数图 | 低 |
| near-resonance 被误写成连续谱极限 | 第 6.6 节明确限定为当前 Dirichlet 五点离散谱 | 低 |
| timing 被误写成公平 kernel benchmark | 第 6.4 节明确不是严格 kernel-to-kernel benchmark | 低 |
| condition-number equivalence 外推 | 第 6.5 节明确只适用于当前 Dirichlet 五点正交 DST-I 系统 | 低 |
| 过程性 patch 术语进入正文 | 第 6.6 节出现“Patch5 补充扫描” | 中低，建议改 |
| GMRES 停滞“根本原因”措辞 | 第 5 章有较强因果语气 | 中低，建议降为“重要原因之一” |
| 最近报告与当前图 16 口径不完全一致 | `RECENT_EXPERIMENT_UPDATE_REPORT.md` 对 exp05 图形描述可能是旧版 | 若该报告提交给导师/评审，建议同步；若仅过程记录，不影响正文 |

## 6. 建议修改列表

### A. 必须改

暂无 A 类问题。当前第 6 章实验结构、公式对应关系和主要图表证据链没有发现阻断提交的理论或证据错误。

### B. 建议改

1. 将第 6.6 节中的“Patch5 补充扫描”改为正式论文语气，例如“补充的多模态扫描”或“进一步的多模态扫描”。
   原因：数学无误，但 patch 术语属于内部过程记录，不适合出现在终稿正文。

2. 将第 5 章关于 GMRES(30) 大量迭代的“根本原因”改为“重要原因之一”或“本文实验中观察到的主要原因”。
   原因：near-resonance 谱结构确实关键，但 GMRES 行为还受 restart、RHS 投影、容差和实现细节影响。

3. 如果 `RECENT_EXPERIMENT_UPDATE_REPORT.md` 会作为最终交付材料之一，建议同步 exp05 图 16 的描述。
   原因：当前源码和第 6 章描述的是 target energy fraction / shape correlation 版本，而该报告中仍残留 amplification ratio / dominant projection check 版本的说明。

### C. 可不改

1. 可在 `fig:exp04_condition_check` 附近补一句“最大相对差约为 `2.91e-13`”，但当前图与正文已经足够说明一致性。

2. 可在 `fig:exp05_dominant_mode_projection` caption 中补一句“右图差值明显小于 full solution 主尺度”，但正文已有对应解释，不是必须。

3. 可统一中英文术语风格，例如 `implementation validation`、`visual sanity check`、`kernel-to-kernel benchmark` 是否保留英文。当前写法技术读者能理解，不影响数学正确性。

4. 可在文档索引或 README 中说明实际第 4 章文件名为 `4_helmholtz.tex`，避免与外部 prompt 中的 `4_helmholtz_fft.tex` 混淆。该问题不影响论文内容。

## 7. 审查结论

当前第 6 章已经具备正式论文实验章的结构，而不是补丁堆叠状态。五个核心实验之间分工清晰：

- exp00 仅作实现一致性校验；
- 实验一验证 Dirichlet 离散格式收敛和边界修正；
- 实验二验证 Neumann/mixed 边界处理与零模态；
- 实验三承担“FFT 快速求解器与迭代法对比研究”的核心实证；
- 实验四建立 modified/true Helmholtz 谱结构与 GMRES 行为之间的联系；
- 实验五用离散谱模态投影收束 near-resonance 小分母放大机制。

总体状态：**可提交级别，建议做少量 B/C 类文字清理**。
