# Final Polish Report

**项目**：矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究
**日期**：2026-04-29
**版本**：V2.4，五核心实验终局整理版
**当前分支**：`revise-sigma-fft9-krylov`
**最终 PDF**：`thesis/main.pdf`，66 pages
**最终状态**：READY TO SUBMIT

---

## 1. 总体结论

本轮清理的核心是把第 6 章从旧版 16 个历史实验和后续补丁式扩展整理为 5 个有脚本、CSV、图表和正文结论支撑的核心实验；`exp00` 降级为 implementation validation，并同步收紧摘要、引言、结论和实验总结的表述边界。

当前论文主线已经稳定为：在统一的 $(-\Delta+\sigma)u=f$ 框架下，验证规则矩形区域上 FFT 快速直接法的收敛阶、边界处理和精度--成本表现，并通过 modified/true Helmholtz 的谱结构、GMRES 残差历史与 near-resonance 模态投影观察 $\sigma$ 符号对求解难度的影响。

论文不再声称完成大范围系统性能比较、完整预处理迭代法研究、FFT9 的 Neumann/mixed 四阶推广，或 Lean 4 的完整 PDE 全局误差证明。

---

## 2. 验证结果

### Python 测试

运行目录：

```powershell
C:\Users\20564\Desktop\Graduate\论文收集\code\python
```

命令：

```powershell
python -m pytest -q
```

结果：

```text
103 passed, 2 skipped, 2 warnings
```

两条 warning 来自 Neumann Poisson compatibility 测试中故意触发的不相容 RHS，属于预期测试行为。

### FFT9 展开验证

命令：

```powershell
python verify_fft9_expansion.py
```

结果：有效特征值的 $h^2$ 项为 0，脚本确认 `Is zero? True`，FFT9 四阶展开验证通过。

### Lean 4

运行目录：

```powershell
C:\Users\20564\Desktop\Graduate\论文收集\code\lean4_formalization
```

命令：

```powershell
lake build
```

结果：

```text
Build completed successfully (3 jobs).
```

仍有 unused simp argument warnings；另有 `.lake/packages/proofwidgets` 依赖目录 local changes 提示。二者均不影响 Lean 源文件构建成功，不纳入本轮提交。

### LaTeX 与 PDF

运行目录：

```powershell
C:\Users\20564\Desktop\Graduate\论文收集\thesis
```

命令：

```powershell
xelatex -interaction=nonstopmode main.tex
bibtex main
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
```

结果：

```text
Output written on main.pdf (66 pages).
0 LaTeX error
0 undefined citation/reference
0 multiply defined label
0 Missing character
0 Overfull \hbox
```

日志中仍有 SimSun 字体替代 warning 和 hyperref bookmark 中数学符号移除 warning，不影响正文渲染、引用解析或 PDF 生成。

PDF 已用 `pdftoppm` 渲染第 6 章关键页和附录表 11 所在页目检：第 6 章表格/图像无裁切、无重叠；Lean 附录表 11 已固定在 `\textwidth` 内，状态列不再越界。

---

## 3. 第 6 章实验闭环

第 6 章当前整理为 5 个核心实验；`exp00` 作为 implementation validation 保留在 6.1，不再作为核心实验编号：

| 实验 | 支撑脚本 | CSV | 图像/表格 |
|------|----------|-----|-----------|
| 实验一：Dirichlet 离散格式收敛性验证 | `exp01_convergence.py`, `exp02_nonhom_bc.py` | `exp01_convergence.csv`, `exp02_nonhom_bc.csv` | `exp01_convergence.png`, `exp02_nonhom_bc.png`, nonhom/multimode 可视化 |
| 实验二：Neumann 与 mixed 边界处理 | `exp03_neumann_mixed.py` | `exp03_neumann_mixed.csv` | `exp03_neumann_mixed_summary.png`, `exp03_neumann_mixed_fields.png` |
| 实验三：精度--成本对比与复杂度实证 | `exp06_accuracy_cost.py` | `exp06_accuracy_cost.csv` | `exp06_accuracy_cost_error_time.png`, `exp06_time_scaling.png` |
| 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为 | `exp07_spectral_denominator_maps.py`, `exp04_modified_vs_true.py` | `exp07_spectral_denominator_summary.csv`, `exp04_modified_vs_true.csv`, `exp04_condition_check.csv`, `exp04_gmres_history.csv` | `exp07_spectral_denominator_heatmaps.png`, `exp04_*` |
| 实验五：True Helmholtz 近共振模态放大 | `exp05_true_helmholtz_resonance.py` | `exp05_resonance.csv`, `exp05_multimode_resonance.csv`, `exp05_resonance_gmres_history.csv` | `exp05_near_resonance_summary.png`, `exp05_multimode_resonance_summary.png`, `exp05_dominant_mode_projection.png`, `exp05_resonance_gmres_history.png` |

当前保留五张可视化补充图：

- `exp02_temperature_field_comparison.png`：用于展示非齐次 Dirichlet 温度场空间结构，不替代收敛阶证据。
- `exp03_neumann_mixed_fields.png`：用于展示 Neumann/mixed 边界代表性解场和误差分布，不替代收敛阶与边界残差证据。
- `exp05_denominator_heatmap.png`：用于展示 true Helmholtz 近共振分母模态分布，保留为备份可视化，不构成新增实验。
- `exp06_complex_manufactured_fields.png`：用于展示多模态 manufactured solution 的解析场、FA/FFT9 数值解与误差分布。
- `exp06_complex_manufactured_convergence.png`：用于展示同一多模态补充问题中 FA 二阶与 FFT9 四阶收敛。

这些补充图已并入对应核心实验的小节中，不新增“实验七”。

每个实验均补充或保留了结论边界：

- exp00 只证明离散 FFT-vs-sparse 一致性，不证明连续 PDE 收敛阶。
- 实验一验证 Dirichlet 收敛轨道；非齐次 Dirichlet 和多模态图作为 visual sanity check。
- 实验二验证受支持五点求解器的 Neumann/mixed 边界处理，不证明 FFT9 Neumann/mixed 四阶性。
- 实验三展示当前实现下的误差--时间趋势，不是严格 kernel-to-kernel benchmark。
- 实验四观察 Gaussian RHS 下无预处理 GMRES 的参数敏感性，不构成完整迭代法性能研究。
- 实验五只支持离散谱 near-resonance 模态放大结论，不外推为连续谱极限或一般波数性能比较。

旧版“实验七”至“实验十六”的结构、弱支撑图表和大范围性能比较表述未在正文中保留。

---

## 4. 已完成关键修复

### 论文正文

- 引言中明确：FA、CR、FACR-like 与五点基准覆盖 Dirichlet、Neumann 和 mixed；FFT9 仅针对 Dirichlet 边界实现并验证。
- FACR 复杂度表述已区分：当前 FACR-like 实现为 $O(N^2\log N)$；经典 FACR 的 $O(N^2\log\log N)$ 仅作为理论背景。
- GMRES 已限定为无预处理 Krylov 基线；true Helmholtz 近共振处改为“可能停滞或依赖参数/预处理”的保守表述。
- 第 6 章所有实验改为 5 个核心闭环，并加入“不证明什么”的限制句。
- Lean 附录表 11 改为固定在 `\textwidth` 内，移除 emoji glyph，避免 PDF 右侧越界。
- README 已同步为最终交付包说明，记录 62 页 PDF、5 个核心实验，并说明 exp00 为 implementation validation。

### 代码与数据

- `exp05_resonance.csv` 表头已将 `N` 改为 `n_interior`，避免 PowerShell `Import-Csv` 大小写重复列名问题。
- `helmholtz_solver.py` 与 `gmres_solver.py` 中 stale “subtract boundary contributions” 注释已改为 add/move to RHS。
- `cyclic_reduction.py` 保留为 deprecated/historical prototype，但同步修正非齐次 Dirichlet RHS 边界符号，避免误用。
- `.gitignore` 已增加 `PROJECT_REPORT_V1.0.md` 与 `Problem_To_Update.md`，二者保留本地但不进入最终交付提交。

---

## 5. Git 状态判断

建议纳入最终提交的文件：

```text
.gitignore
README.md
FINAL_POLISH_REPORT.md
NEXT_DIALOG_EXPERIMENT_RESTRUCTURE.md
code/python/cyclic_reduction.py
code/python/gmres_solver.py
code/python/helmholtz_solver.py
code/python/experiments/exp05_true_helmholtz_resonance.py
code/python/experiments/results/exp05_resonance.csv
thesis/main.tex
thesis/main.pdf
thesis/chapters/1_introduction.tex
thesis/chapters/2_math_preliminary.tex
thesis/chapters/3_fft_direct.tex
thesis/chapters/4_helmholtz.tex
thesis/chapters/5_gmres.tex
thesis/chapters/6_experiments.tex
thesis/chapters/7_conclusion.tex
thesis/chapters/appendix_lean4.tex
```

不建议提交：

```text
PROJECT_REPORT_V1.0.md
Problem_To_Update.md
thesis/*.aux
thesis/*.log
thesis/*.bbl
thesis/*.blg
__pycache__
.lake/packages/proofwidgets
temporary logs/rendered PNGs
```

---

## 6. Remaining Risks

High: none.

Medium: none.

Low:

- Lean 构建仍提示 unused simp arguments，可后续作为形式化代码清洁项处理。
- `.lake/packages/proofwidgets` 依赖目录存在 local changes 提示，本轮未纳入提交。
- LaTeX 存在字体替代与 hyperref bookmark warning，但最终 PDF 正文、表格和引用均正常。

---

**Final status: READY TO SUBMIT**
