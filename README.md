# 矩形区域 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究

当前分支：`revise-sigma-fft9-krylov`

最新论文 PDF：`thesis/main.pdf`（当前 71 页）

## 项目导航

根目录只保留两个入口文档：

- `README.md`：交付概览、运行入口和验收命令。
- `PROJECT_STRUCTURE.md`：项目结构地图、文件移动安全规则和报告归档策略。

主要内容位置如下：

| 内容 | 路径 |
|---|---|
| 论文 LaTeX 与 PDF | `thesis/` |
| Python 求解器、测试与实验脚本 | `code/python/` |
| Lean 4 辅助形式化验证 | `code/lean4_formalization/` |
| 第 6 章 HTML/PPT 讲解材料 | `reports/ch6_experiment_report/` |
| 审查、补丁、实验、答辩和项目报告索引 | `reports/REPORT_INDEX.md` |
| 实验脚本、CSV 与图像说明 | `code/python/experiments/README.md` |

## 文件整理原则

本轮项目整理只移动 Markdown 报告和说明材料，不移动运行资产：

- 不移动 `code/python/` 下的求解器、测试和实验入口。
- 不移动 `code/python/experiments/results/` 与 `code/python/experiments/figures/`。
- 不移动 `thesis/figures/` 或 LaTeX 章节文件。
- 不改动实验参数、CSV、PNG/PDF 图像生成逻辑。

因此现有 Python、Lean 和 LaTeX 入口仍按原路径运行。

## 第 6 章核心实验

第 6 章整理为 5 个核心实验；`exp00_fft_vs_sparse.py` 降级为 implementation validation，不作为核心实验编号。

| 实验 | 脚本 | 主要证据 |
|---|---|---|
| 实验一：Dirichlet 离散格式收敛性验证 | `exp01_convergence.py`, `exp02_nonhom_bc.py` | 二阶五点格式、FFT9 四阶 Dirichlet 实现、非齐次 Dirichlet 可视化检查 |
| 实验二：Neumann 与 mixed 边界处理 | `exp03_neumann_mixed.py` | ghost-point、DCT-I 对称化、flux residual、weighted mean |
| 实验三：精度--成本对比与复杂度实证 | `exp06_accuracy_cost.py` | error-time 与 time-scaling 图 |
| 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为 | `exp07_spectral_denominator_maps.py`, `exp04_modified_vs_true.py` | small-denominator risk map、condition check、GMRES residual history |
| 实验五：True Helmholtz 近共振模态放大 | `exp05_true_helmholtz_resonance.py` | near-resonance summary、dominant projection、GMRES history |

CSV 位于 `code/python/experiments/results/`；实验生成图位于 `code/python/experiments/figures/`；论文使用图位于 `thesis/figures/`。

## 关键限制

- FFT9 当前只针对规则矩形 Dirichlet 边界实现和验证。
- GMRES 结果限定为无预处理 restarted GMRES(30) 或对应实验设置下的无预处理重启 GMRES。
- FACR-like 当前实现仍按 `O(N^2 log N)` 处理，不声称达到经典 FACR 的 `O(N^2 log log N)`。
- near-resonance 是 Dirichlet 五点离散谱实验，不等同于连续谱极限。
- Lean 4 只验证三参数九点修正族六阶 local consistency 不可达反证中的有限代数核心，不证明完整 PDE 全局误差或 Python 程序正确性。

## 验收命令

从仓库根目录运行：

```powershell
python -m pytest -q code/python/tests

cd code/lean4_formalization
lake build

cd ../../thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

当前验证状态：

- Python tests：`105 passed, 2 skipped`
- Lean：`lake build` 通过；可能仍有既有 unused simp warning
- LaTeX：`main.pdf` 正常生成，当前 71 页
