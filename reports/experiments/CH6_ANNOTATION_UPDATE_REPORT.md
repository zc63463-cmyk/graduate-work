# CH6 Annotation Update Report

## 1. Scope

本轮执行 Phase 2：只为第 6 章关键图补充简短读图注解。

未新增实验，未修改图像，未修改 CSV，未改变数值算法，未改变理论结论。

## 2. Modified Files

- `thesis/chapters/6_experiments.tex`
- `CH6_ANNOTATION_UPDATE_REPORT.md`

## 3. Annotation Updates

### Figure `fig:exp03_neumann_mixed`

位置：Neumann/Mixed boundary verification summary 后。

新增说明：

- panel (a) 用于观察 FA 在不同边界/参数组合下是否保持二阶斜率；
- panel (b) 用于观察 Neumann flux residual 与 mixed Dirichlet residual 是否闭合；
- panel (c) 用于检查 Neumann Poisson 零模态归一化；
- CR/FACR-like 的重复结果保留在 CSV 中，主图不重复绘制。

### Figure `fig:exp06_accuracy_cost`

位置：精度--成本 error-time 图后。

新增说明：

- 读图原则为越靠左下越好；
- FFT9 的优势来自 Dirichlet 四阶离散精度，而不是迭代残差；
- GMRES30 是无预处理基线，其误差同时受离散误差和代数残差影响；
- timing scope 不是严格 kernel-to-kernel benchmark，且 GMRES setup 被排除。

### Figure `fig:exp07_spectral_denominator_maps`

位置：图 10 后。

新增说明：

- 上排低频 risk map 用于定位小分母热点；
- 下排 sorted `|d|` 图用于比较三类问题离零的数量级；
- `R_{p,q}<0` 仅表示 `|d_{p,q}|>1`，即分母较大、小分母风险较低，并不代表负风险；
- modified 约在 `10^2` 量级，true-away 约在 `10^1` 量级，near-resonance 的 `(2,3)/(3,2)` 降至 `10^{-2}` 量级。

### Figure `fig:exp04_condition_check`

位置：condition check 图后。

新增说明：

- 横轴为 spectral denominator indicator；
- 纵轴为 dense `cond_2(A)`；
- 点贴近 `y=x` 仅校验当前 Dirichlet 五点正交 DST-I 对角化系统；
- 不外推到 mixed/Neumann 原始矩阵或一般 non-normal 系统。

### Figure `fig:exp04_gmres_history`

位置：GMRES residual history 图后。

新增说明：

- absolute residual 是实际停止准则；
- relative residual 仅用于归一化比较；
- modified `sigma=1000` 作为条件较好的对照能降至容差；
- true Helmholtz 或 near-resonance 参数下出现下降缓慢和 capped；
- 最终迭代次数表不能替代残差历史。

### Figure `fig:exp05_multimode_resonance_summary`

位置：图 16 后。

新增说明：

- panel (a) 是核心证据，展示目标模态系数按 `1/|delta|` 放大；
- panel (b)(c) 已接近 `1`，主要说明目标模态子空间已经主导解场；
- panel (d) 只汇总无预处理 GMRES(30) 的 capped 计数；
- 该行为不外推到所有 Helmholtz 求解器。

### Figure `fig:exp05_resonance_gmres_history`

位置：near-resonance GMRES residual history 图后。

新增说明：

- 将表格中的“未收敛”展开为逐步残差过程；
- `1001` 是 `max_iter=1000` 下的 capped 计数约定；
- near-resonance 的问题体现为残差停滞过程，而不只是最终状态标签。

## 4. Validation

Commands run:

```powershell
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

During the chained run, one final XeLaTeX pass briefly failed with `Unable to open main.pdf`, consistent with a transient PDF file lock. A subsequent XeLaTeX pass succeeded without changing source files.

Final validation:

- XeLaTeX/BibTeX: PASS
- Output PDF: `thesis/main.pdf`
- PDF pages: 64
- Log scan:
  - no fatal error
  - no undefined reference
  - no multiply-defined label
  - no overfull hbox
- Non-blocking warning:
  - existing SimSun bold-italic font substitution warning
- `git diff --check`:
  - no whitespace errors
  - only existing LF/CRLF normalization notice

## 5. Final Status

READY TO USE
