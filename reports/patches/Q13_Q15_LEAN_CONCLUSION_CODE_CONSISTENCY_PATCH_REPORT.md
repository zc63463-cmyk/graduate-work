# Q13--Q15 Lean / 贡献边界 / FFT9 符号系统最小补丁报告

## 1. 本轮目标

本轮只做终审级文字与说明收紧，不新增实验、不改数值算法、不改实验参数、不重绘图、不改 CSV/PNG/PDF 实验资产。

核心目标有三点：

1. 将 Lean 4 的作用限定为三参数九点修正族六阶 local consistency 不可达性反证中的有限代数核心验证。
2. 将摘要、引言和第六章总结中的“理论--实验闭环”类表述收紧为“离散谱分析--实现校验--数值实验”的证据链或相互印证。
3. 根据当前仓库事实，继续采用 FFT9 raw/theory 符号系统口径，不把当前 Python 实现误写成整体反号系统。

## 2. 当前仓库事实：FFT9 使用 raw/theory 系统

经检查，当前 `code/python/helmholtz_solver.py` 中 FFT9 主实现采用：

\[
\tilde g=(R_h f)_I+L_{IB}g_B-\sigma R_{IB}g_B,
\qquad
D^{raw}_{p,q}=-\hat\lambda_L+\sigma\hat\lambda_R .
\]

具体表现为：

- `rhs_tilde = apply_Rh_full(F)` 后加上 `compute_bc_correction_9pt(...)` 返回的 raw 边界修正；
- 频域分母使用 `denom = -lam_L + sigma * lam_R4`；
- `tests/test_04_fft9_vs_spsolve.py` 中 sparse compact system 对照也使用 `rhs_tilde = R_h f + bc_corr` 的同一口径。

因此，本轮没有采用“当前代码反号系统”的论文替换方案。第 3 章和第 4 章继续写明：当前理论推导和当前 Python 实现均采用 raw 系统。反号系统只作为等价代数写法出现，并明确要求 RHS 与 denominator 必须同步反号，不能只反其中一项。

## 3. 修改文件与内容

### 3.1 `thesis/main.tex`

摘要中将 Helmholtz 实验的表述从较宽泛的“理论谱分析与数值现象相互印证”收紧为：

- 在当前离散谱分析、实现校验和数值实验之间相互印证；
- near-resonance 下的小分母结构与无预处理 GMRES(30) 残差停滞相一致。

这样避免将观察到的 GMRES 停滞写成对所有 Krylov/Helmholtz 求解器的强因果结论。

### 3.2 `thesis/chapters/1_introduction.tex`

将贡献点中的“理论--实验闭环”改为：

- “离散谱分析--实现校验--数值实验”的证据链。

这保留了论文实验章的组织逻辑，同时避免暗示 Lean 或实验已经完成完整 PDE 理论证明。

### 3.3 `thesis/chapters/3_fft_direct.tex`

保留当前 raw 系统口径，并补强反号系统防误读说明：

- 当前代码采用与理论推导一致的 raw 符号系统；
- 若某实现整体乘以 \(-1\)，则 RHS 与频域分母必须同步反号；
- 不能只反右端项或只反频域分母。

本轮未修改 FFT9 求解算法。

### 3.4 `thesis/chapters/4_helmholtz.tex`

在 FFT9 Helmholtz 分母说明处补充：

- 反号系统只能作为整体代数变换使用；
- 不能与本文 raw 右端项或 raw 分母交叉混用。

### 3.5 `thesis/chapters/6_experiments.tex`

将第六章总结中的“数值证据与前文理论形成闭环”收紧为：

- 数值证据与前文离散谱分析和实现校验形成证据链。

同时保留已有限制：FFT9 四阶表现仅在当前 Dirichlet 光滑 manufactured solution 上验证，GMRES 结论仅限无预处理 restarted GMRES(30)。

### 3.6 `thesis/chapters/appendix_lean4.tex`

Lean 附录进一步收紧为：

- Lean 4 是辅助验证；
- 验证对象是 local Fourier symbol 展开后的有限多项式系数条件；
- 表格只对应反证链条中的有限代数命题；
- \((1,0)\)、\((1,1)\) 是 local Fourier symbol 中的连续波数测试点，不是 Dirichlet DST-I 离散模态编号；
- Lean 不覆盖 Fourier symbol 从模板导出的全过程、Taylor 建模、DST/DCT 谱理论、边界闭合、完整 PDE 全局误差、Python 程序或实验正确性。

## 4. 未修改范围

本轮未修改：

- Lean theorem 语句或证明代码；
- `code/python/helmholtz_solver.py` 执行逻辑；
- FFT9 sparse 对照测试逻辑；
- 实验脚本、CSV、PNG、PDF 实验图；
- HTML/PPT 讲解材料；
- GMRES、Neumann、exp05/exp06/exp07 的实验数据。

## 5. 静态检查结果

已检查并确认：

- 未发现“Lean 证明完整 PDE 理论”“所有九点格式六阶不可达”“完整理论闭环”等过强表述；
- 未发现把当前 Python FFT9 实现写成“当前代码反号系统”的表述；
- 自有 Lean 文件中未发现 `sorry`、`admit`、`axiom` 或 `unsafe`。

说明：Lean 依赖目录 `.lake` 不纳入自有文件逃逸检查。

## 6. 验证命令与结果

### Lean

命令：

```powershell
cd code/lean4_formalization
lake build
```

结果：

- build completed successfully；
- 仅出现既有 unused simp argument warnings；
- 另有 `proofwidgets` 依赖仓库存在 local changes 的提示；
- 无构建失败。

### LaTeX

命令：

```powershell
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

结果：

- `main.pdf` 正常生成；
- 最终页数：71 页；
- log 扫描未发现 fatal error；
- 未发现 undefined reference/citation；
- 未发现 multiply-defined label；
- 未发现 overfull hbox。

保留的非致命信息：

- SimSun font shape substitution warning；
- hyperref PDF string token warnings；
- MiKTeX update 提示。

## 7. 结论

Q13--Q15 本轮补丁完成后，论文对 Lean 4 的覆盖范围、实验与理论的证据边界、以及 FFT9 raw/theory 符号系统与当前代码实现之间的关系更加明确。当前版本不会把 Lean 4 外推为完整 PDE 形式化证明，也不会把当前 FFT9 Python 实现误写成整体反号系统。
