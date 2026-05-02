# 矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究

当前版本：V2.3（实验章节删减收束与交付一致性清理版）
当前分支：`revise-sigma-fft9-krylov`
当前提交：随本次最终提交生成（见 Git 日志）
最新论文 PDF：`thesis/main.pdf`（71 pages）

## 项目导航

- 项目结构总图与运行安全约定：`PROJECT_STRUCTURE.md`
- 报告与审查文档索引：`reports/REPORT_INDEX.md`
- 第 6 章实验脚本、CSV 与图像说明：`code/python/experiments/README.md`
- 本轮结构整理审查记录：`PROJECT_REVIEW_AND_STRUCTURE_REPORT.md`

为避免破坏脚本和 LaTeX 引用，本轮结构整理不迁移 `code/`、`thesis/` 下的运行时文件、实验 CSV 或论文图像。历史根目录报告暂时保留原路径；后续新增审查报告优先放入 `reports/audit/`。

## 交付内容

- 论文正文与 PDF：`thesis/main.tex`、`thesis/main.pdf`
- Python 求解器与实验脚本：`code/python/`
- Lean 4 辅助形式化验证：`code/lean4_formalization/`
- 最终清理报告：`FINAL_POLISH_REPORT.md`
- 实验收束交接说明：`NEXT_DIALOG_EXPERIMENT_RESTRUCTURE.md`

## 第 6 章核心实验

第 6 章整理为 5 个核心实验；`exp00_fft_vs_sparse.py` 降级为 implementation validation，不再作为核心实验编号。论文不恢复旧版 16 个历史实验结构。

| 实验 | 脚本 | CSV | 图像 |
|------|------|-----|------|
| 实验一：Dirichlet 离散格式收敛性验证 | `exp01_convergence.py`, `exp02_nonhom_bc.py` | `exp01_convergence.csv`, `exp02_nonhom_bc.csv` | `exp01_convergence.png`, `exp02_nonhom_bc.png`, nonhom/multimode 可视化 |
| 实验二：Neumann 与 mixed 边界处理 | `exp03_neumann_mixed.py` | `exp03_neumann_mixed.csv` | `exp03_neumann_mixed_summary.png`, `exp03_neumann_mixed_fields.png` |
| 实验三：精度--成本对比与复杂度实证 | `exp06_accuracy_cost.py` | `exp06_accuracy_cost.csv` | `exp06_accuracy_cost_error_time.png`, `exp06_time_scaling.png` |
| 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为 | `exp07_spectral_denominator_maps.py`, `exp04_modified_vs_true.py` | `exp07_spectral_denominator_summary.csv`, `exp04_modified_vs_true.csv`, `exp04_condition_check.csv`, `exp04_gmres_history.csv` | `exp07_spectral_denominator_heatmaps.png`, `exp04_*` |
| 实验五：True Helmholtz 近共振模态放大 | `exp05_true_helmholtz_resonance.py` | `exp05_resonance.csv`, `exp05_multimode_resonance.csv`, `exp05_resonance_gmres_history.csv` | `exp05_near_resonance_summary.png`, `exp05_multimode_resonance_summary.png`, `exp05_dominant_mode_projection.png`, `exp05_resonance_gmres_history.png` |

Implementation validation：

- `exp00_fft_vs_sparse.py` / `exp00_fft_vs_sparse.csv`：只说明 FFT 求解路径与 sparse direct 解在离散系统层面一致，不作为连续 PDE 精度贡献。

CSV 位于 `code/python/experiments/results/`。
实验生成图位于 `code/python/experiments/figures/`，论文使用图位于 `thesis/figures/`。
当前交付包含 5 个核心实验对应的 CSV/PNG/PDF 证据；exp00 作为 implementation validation 保留在正文表格中。

此外提供若干 supplementary visualization scripts，用于生成非齐次 Dirichlet 温度场、near-resonance 机制备份图和多模态 manufactured solution 可视化：

- `exp02_temperature_field_comparison.png`：非齐次 Dirichlet 温度场、FFT9 数值解、误差分布与截线。
- `exp03_neumann_mixed_fields.png`：Neumann/mixed 边界代表性解场与误差分布；收敛阶结论仍由 `exp03_neumann_mixed_summary.png` 与 CSV 支撑。
- `exp05_denominator_heatmap.png`：true Helmholtz near-resonance 频域分母热图，作为备份可视化；正文实验五主图为 `exp05_near_resonance_summary.png`。
- `exp06_complex_manufactured_fields.png`：多模态 manufactured solution 的解析场、FA/FFT9 数值解与误差分布。
- `exp06_complex_manufactured_convergence.png`：多模态 manufactured solution 的 FA 二阶与 FFT9 四阶收敛补充验证。

其中 `exp06_complex_manufactured_visualization.py` 为 supplementary visualization script，
不计入 5 个核心实验编号。

## 范围限定

- FA、CR、FACR-like 与五点基准覆盖 Dirichlet、Neumann 和一种 mixed 边界条件。
- FFT9 四阶紧致求解器仅针对 Dirichlet 边界实现并验证。
- GMRES 仅作为无预处理 Krylov 基线，用于观察谱结构对迭代行为的影响。
- 当前 FACR-like 实现复杂度仍为 `O(N^2 log N)`；经典 FACR 的 `O(N^2 log log N)` 只作为理论背景。
- Lean 4 只验证三参数九点修正族六阶不可达性的代数核心步骤，不验证完整 PDE 全局误差估计、DST/DCT 谱理论或程序实现正确性。

## 验收命令

```powershell
cd code/python
python -m pytest -q
python verify_fft9_expansion.py

cd ../lean4_formalization
lake build

cd ../../thesis
xelatex -interaction=nonstopmode main.tex
bibtex main
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
```

预期状态：

- Python 活跃测试：`103 passed, 2 skipped`；Neumann compatibility 警告为预期测试行为。
- FFT9 展开验证：有效特征值的 h² 项消去。
- Lean 4：构建成功，可能仍有 unused simp argument warnings。
- LaTeX：0 error，0 undefined citation/reference，0 multiply defined label。
