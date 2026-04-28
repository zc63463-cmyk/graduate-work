# 矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究

当前版本：V2.3（实验章节删减收束与交付一致性清理版）  
当前分支：`revise-sigma-fft9-krylov`  
当前提交：随本次最终提交生成（见 Git 日志）  
最新论文 PDF：`thesis/main.pdf`（52 pages）

## 交付内容

- 论文正文与 PDF：`thesis/main.tex`、`thesis/main.pdf`
- Python 求解器与实验脚本：`code/python/`
- Lean 4 辅助形式化验证：`code/lean4_formalization/`
- 最终清理报告：`FINAL_POLISH_REPORT.md`
- 实验收束交接说明：`NEXT_DIALOG_EXPERIMENT_RESTRUCTURE.md`

## 第 6 章核心实验

第 6 章只保留 6 个核心实验，不再恢复旧版 16 个历史实验、旧 heatmap、restart 图、FACR 参数图或复杂度趋势图。

| 实验 | 脚本 | CSV | 图像 |
|------|------|-----|------|
| FFT 解与 sparse direct 解一致性 | `exp00_fft_vs_sparse.py` | `exp00_fft_vs_sparse.csv` | 表格 |
| 二阶与四阶收敛阶验证 | `exp01_convergence.py` | `exp01_convergence.csv` | `exp01_convergence.png` |
| 非齐次 Dirichlet 边界验证 | `exp02_nonhom_bc.py` | `exp02_nonhom_bc.csv` | `exp02_nonhom_bc.png` |
| Neumann 与 mixed 边界验证 | `exp03_neumann_mixed.py` | `exp03_neumann_mixed.csv` | `exp03_neumann_mixed.png` |
| modified/true Helmholtz 谱指标与 Gaussian RHS 下 GMRES 行为 | `exp04_modified_vs_true.py` | `exp04_modified_vs_true.csv` | `exp04_min_denom_vs_sigma.png`, `exp04_spectral_indicator_vs_sigma.png`, `exp04_gmres_iters_vs_sigma.png` |
| true Helmholtz near-resonance 扫描 | `exp05_true_helmholtz_resonance.py` | `exp05_resonance.csv` | `exp05_resonance.png` |

CSV 位于 `code/python/experiments/results/`。  
实验生成图位于 `code/python/experiments/figures/`，论文使用图位于 `thesis/figures/`。  
当前交付包含 exp00-exp05 共 6 个 CSV 和 7 张核心 PNG 图。

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
