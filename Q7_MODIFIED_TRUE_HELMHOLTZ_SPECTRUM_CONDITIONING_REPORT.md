# Q7 Modified/True Helmholtz 谱分类、条件数与谱分母指标口径收紧报告

## 1. 修复目标

本轮 Q7 只收紧理论表述、实验解释和代码注释，不新增实验、不改变数值算法、不改变实验参数，也不重绘图或改写 CSV。核心目标是避免把 Dirichlet 五点正交 DST-I 系统中的特殊结论外推到一般边界、非正交对角化或 non-normal 系统。

## 2. 论文理论章节修改

### 第 2 章

- 将 modified Helmholtz 条件数改善限定为固定网格、同一 Dirichlet 五点 SPD 系统、\(\kappa>0\) 下的结论。
- 明确 \(\kappa=0\) 时 modified Helmholtz 退化为 Poisson，条件数相等。
- 将 true Helmholtz 的谱分母指标解释为 near-resonance 小分母风险；当系统可能奇异、不定或负定时，不再按 SPD 条件数结论解释。

### 第 4 章

- 收紧 \(\sigma\) 分类表：true Helmholtz 的矩阵性质按谱可为正定、奇异、不定或负定；迭代法口径改为正定区间可用 CG，不定时考虑 MINRES/GMRES。
- 将 modified Helmholtz 定理重写为 Dirichlet 五点 \(A_\kappa=A_0+\kappa^2I\) 的 SPD 与 \(\kappa_2\) 改善定理，并说明 Neumann ghost-point 原始非对称矩阵、mixed 边界和一般 non-normal 系统不能直接套用该公式。
- 保留 true Helmholtz 四类谱分类，并补充 \(\kappa^2>\lambda_{\max}\) 时 \(-A_h\) 为 SPD，但 \(A_h\) 本身不是 SPD。
- near-resonance 条件统一为
  \[
  |d_{p,q}|=|\lambda_{p,q}-\kappa^2|\ll 1,
  \]
  即 \(\kappa^2\approx\lambda_{p,q}\)，不与 modified Helmholtz 的正分母情形混淆。

## 3. 第 6 章实验解释修改

- 图 10 周边补强：risk map 只可视化 small-denominator risk，不是一般条件数图。
- 图 11/12 周边显式写出当前 Dirichlet 五点系统的正交对角化结构
  \[
  A=Q^T\operatorname{diag}(d_{p,q})Q,\qquad Q^TQ=I.
  \]
  因此非奇异时奇异值才等于 \(|d_{p,q}|\)，并得到
  \[
  \mathrm{cond}_2(A)=\max |d_{p,q}|/\min |d_{p,q}|.
  \]
- 将 condition check 明确为当前 Dirichlet 五点系统的数值校验，不外推到 Neumann ghost-point 原始非对称矩阵、mixed 边界、非正交对角化或一般 non-normal 系统。

## 4. 第 7 章结论修改

- 同步限定 modified Helmholtz 的条件数改善只针对同一 Dirichlet 五点 SPD 系统。
- 同步限定 spectral denominator indicator 与 \(\mathrm{cond}_2(A)\) 的等价性依赖正交对角化与 symmetric/normal 结构。
- true Helmholtz 结论改为按离散谱分类：可能正定、奇异、不定或负定；near-resonance 的核心机制是小分母、模态放大和无预处理 GMRES(30) 停滞。

## 5. 代码注释修改

- `code/python/experiments/exp04_modified_vs_true.py`
  - 在 spectral denominator indicator 和 condition check 相关函数 docstring 中补充：该指标作为 \(\mathrm{cond}_2(A)\) 只适用于 Dirichlet 五点、正交 DST-I 对角化、非奇异的 symmetric/normal 系统。
  - 明确不要将其用于原始 Neumann ghost-point、mixed 边界、非正交对角化或一般 non-normal 系统。

- `code/python/experiments/exp07_spectral_denominator_maps.py`
  - 在脚本说明中补充：risk map 只展示 small-denominator risk，不是条件数图。
  - 说明只有 Dirichlet 五点 DST-I 正交对角化系统中，分母比值才与 \(\mathrm{cond}_2(A)\) 一致。

上述修改均为注释和解释层面，没有改变实验脚本参数、算法或输出文件名。

## 6. 未修改范围

- 未新增实验。
- 未运行 exp04 或 exp07 实验脚本。
- 未重绘任何 PNG/PDF 图。
- 未修改 CSV 数据。
- 未修改求解器算法、GMRES 参数、Helmholtz 参数或绘图参数。
- 未把 exp05 near-resonance 的 Dirichlet 五点离散谱结论扩展到 FFT9 或其他 Helmholtz 求解器。

## 7. 验证结果

- Python 语法检查：
  - `python -m py_compile code/python/experiments/exp04_modified_vs_true.py code/python/experiments/exp07_spectral_denominator_maps.py`
  - 结果：通过。

- LaTeX 编译：
  - 已执行 `xelatex -> bibtex -> xelatex -> xelatex`。
  - `thesis/main.pdf` 正常生成。
  - 最终 PDF 页数：70 页。

- LaTeX log 扫描：
  - 未发现 fatal error。
  - 未发现 undefined reference/citation。
  - 未发现 multiply-defined label。
  - 未发现 overfull hbox。
  - 剩余记录仅包含字体形状警告和 PDF 输出信息。

## 8. 本轮结论

Q7 后，modified Helmholtz 的条件数改善、true Helmholtz 的谱分类、risk map 的解释、spectral denominator indicator 与 \(\mathrm{cond}_2(A)\) 的等价性均已收紧到当前可证明和可验证的范围。全文不再把 Dirichlet 五点正交 DST-I 系统中的特殊条件数结论外推到 Neumann/mixed、非正交对角化或一般 non-normal 系统。
