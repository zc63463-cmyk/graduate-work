# EXP03 Visual Upgrade Report

生成时间：2026-04-29

## 1. Goal

本轮升级第 6 章“Neumann 与 mixed 边界条件验证”的指标和图像展示，并进一步对 summary figure 做去冗余优化：从原 2x3 六子图压缩为 1x3 三子图横版，不新增核心实验编号，不恢复旧版 16 个历史实验结构。

实验定位保持为：验证 ghost-point + DCT-I/DST-I 对称化框架下，FA、CR、FACR-like 五点求解器在 Neumann/mixed 边界下的正确性、二阶收敛、边界约束和 Neumann Poisson 零模态处理。

## 2. Metrics Added

`exp03_neumann_mixed.csv` 已扩展为 V2 指标表，新增或显式保留以下列：

- `linf_error`
- `l2_error`
- `observed_order_linf`
- `observed_order_l2`
- `boundary_flux_linf`
- `boundary_dirichlet_linf`
- `mean_offset`
- `weighted_mean_solution`
- `compatibility_integral`

同时保留旧兼容列 `error` 与 `rate`，分别对应 `linf_error` 和 `observed_order_linf`。

## 3. Figures Generated

已生成并同步复制到 `thesis/figures/`：

- `code/python/experiments/figures/exp03_neumann_mixed_summary.png`
- `thesis/figures/exp03_neumann_mixed_summary.png`
- `code/python/experiments/figures/exp03_neumann_mixed_fields.png`
- `thesis/figures/exp03_neumann_mixed_fields.png`

原图 `exp03_neumann_mixed.png` 保留并重新生成，但正文主图已替换为 1x3 summary 图。

当前 summary 图的三个 panel 为：

- `(a) Linf convergence, FA solver`：仅展示 FA 的四条 `Linf` 收敛曲线和 `O(h^2)` 参考线，避免 CR/FACR 与 L2 重叠线造成视觉冗余。
- `(b) Boundary constraint residuals`：仅绘制一条代表性 Neumann flux residual，并在图例中标注各 case 重合；同时绘制 Mixed `(N,D)` Dirichlet boundary residual。
- `(c) Neumann Poisson zero-mode check`：展示三种 solver 的 weighted mean solution 与 removed mean offset。

## 4. Numerical Results

FA 在 `n=65` 的代表性结果如下：

| Case | Linf | L2 | order Linf/L2 | flux residual |
|---|---:|---:|---:|---:|
| Pure Neumann, sigma=1 | 1.911e-04 | 9.704e-05 | 2.00 / 2.02 | 9.288e-05 |
| Pure Neumann, sigma=10 | 1.333e-04 | 6.767e-05 | 2.00 / 2.02 | 9.287e-05 |
| Mixed (N,D), sigma=10 | 1.333e-04 | 6.663e-05 | 2.00 / 2.00 | 9.287e-05 |
| Neumann Poisson, sigma=0 | 2.008e-04 | 1.020e-04 | 2.00 / 2.02 | 9.288e-05 |

在 `n=129`，FA 的 `observed_order_linf` 均约为 2.00，`observed_order_l2` 为约 2.00--2.01。Neumann boundary flux residual 从粗网格约 `4.69e-02` 降至 `1.16e-05`，呈稳定下降。

Neumann Poisson 的离散兼容条件积分为舍入误差量级或 0。误差比较已移除 trapezoid-weighted mean offset；FA 的 weighted mean solution 为约 `1e-18`，CR/FACR 的 mean offset 为可见但随网格加密下降的离散常数偏移，已在误差计算中扣除。

Pure Neumann `sigma=1` 若在粗网格 `Linf` 曲线上出现轻微折线，可由边界误差与高阶项影响解释；`L2` 曲线、boundary flux residual 和细网格拟合阶数仍支撑整体二阶结论。

## 5. Thesis Integration

已修改 `thesis/chapters/6_experiments.tex`：

- 实验四主图替换为 `figures/exp03_neumann_mixed_summary.png`。
- 新增 `figures/exp03_neumann_mixed_fields.png` 作为代表性解场与误差分布可视化补充。
- 表格增加 RMS 型 `L2` 误差和 `Linf/L2` 观测阶。
- 图注已同步为三子图叙事：收敛验证、边界约束、零模态处理；L2 保留在表格和 CSV 中作为补充证据，不再作为 summary 主图曲线。边界约束 panel 中仅保留一条代表性 Neumann flux 曲线，因为四类 case 的 flux residual 在当前同一 ghost-point 与二阶单边差分估算下数值重合。
- 正文明确本实验只验证受支持五点求解器的 Neumann/mixed 处理，不声称 FFT9 支持 Neumann/mixed 四阶验证。

README、`FINAL_POLISH_REPORT.md` 与 `VISUALIZATION_SUPPLEMENT_REPORT.md` 已同步为：第 6 章仍为 6 个核心实验，PDF 为 57 pages，可视化补充 PNG 数量更新为 5。

## 6. Validation

| Check | Command | Result |
|---|---|---|
| EXP03 script | `cd code/python; python -m experiments.exp03_neumann_mixed` | PASS |
| pytest | `cd code/python; python -m pytest -q` | PASS: 103 passed, 2 skipped, 2 warnings |
| XeLaTeX/BibTeX | `cd thesis; xelatex; bibtex; xelatex; xelatex` | PASS |
| PDF pages | `thesis/main.pdf` | 57 pages |
| LaTeX log scan | `Select-String thesis\main.log ...` | PASS: no fatal error, no undefined citation/reference, no multiply-defined label |

Non-blocking warnings:

- pytest 仍有既有 Neumann Poisson compatibility warnings，属于预期测试行为。
- XeLaTeX 仍有 SimSun bold italic font substitution 和 hyperref bookmark token warnings，不影响 PDF 生成或引用解析。

## 7. Final Status

READY TO USE
