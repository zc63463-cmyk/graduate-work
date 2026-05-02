# 下一对话快速启动任务书：GLM-5.1 论文实验收束后复核

## 0. 工作区

```text
C:\Users\20564\Desktop\Graduate\论文收集
```

当前分支：

```text
revise-sigma-fft9-krylov
```

最新已知提交：

```text
b1ae0c1 fix: correct table data inconsistency (n=65 si vs n=33 iters) and GMRES description
```

## 1. 当前状态

上一轮已经完成核心修正：第 6 章实验章节不再保留“实验一到实验十六”的历史堆叠，而是收束为 6 个核心实验。

最新 PDF：

```text
thesis/main.pdf
```

页数：

```text
52 pages
```

最新 Python 测试：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集\code\python
python -m pytest -q
# 103 passed, 2 skipped, 2 warnings
```

最新 LaTeX 日志检查：

```text
0 error
0 undefined citation/reference
0 multiply defined label
Output written on main.pdf (52 pages).
```

仍有原有 font/hyperref warnings，不影响编译成功和引用解析。

最新 Lean 验证：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集\code\lean4_formalization
lake build
# Build completed successfully (3 jobs)
```

仍有 unused simp argument warnings；另有 `.lake/packages/proofwidgets` 依赖目录 local changes 提示，不属于本轮正文修改。

## 2. 已完成改动

### 第 6 章实验结构

`thesis/chapters/6_experiments.tex` 已改为：

```text
6.1 实验设置与可复现性说明
6.2 实验一：FFT 解与 sparse direct 解的一致性
6.3 实验二：二阶与四阶收敛阶验证
6.4 实验三：非齐次 Dirichlet 边界条件验证
6.5 实验四：Neumann 与混合边界条件验证
6.6 实验五：Modified 与 True Helmholtz 的谱指标和 GMRES 行为
6.7 实验六：True Helmholtz 近共振扫描
6.8 实验总结
```

已从正文移除旧版弱支撑或重复实验引用：

```text
gmres_convergence.png
efficiency_comparison.png
wavenumber_effects.png
error_heatmaps.png
error_grid_refinement.png
solution_surface.png
condition_number.png
facr_parameter.png
gmres_restart.png
complexity_scaling.png
complex_test_heatmap.png
wavenumber_heatmap.png
```

图片文件本身未删除，只是不再被正文引用。

### 同步修改的论文文件

```text
thesis/main.tex
thesis/chapters/1_introduction.tex
thesis/chapters/2_math_preliminary.tex
thesis/chapters/3_fft_direct.tex
thesis/chapters/4_helmholtz.tex
thesis/chapters/6_experiments.tex
thesis/chapters/7_conclusion.tex
thesis/chapters/appendix_lean4.tex
thesis/main.pdf
```

重点：

- 摘要改为核心验证 + 谱指标/GMRES 行为观察。
- 引言不再承诺大范围系统性能比较。
- 结论将 GMRES 降为 Gaussian RHS 和 near-resonance 下的行为观察。
- 条件数公式残留已修。
- FA 边界贡献措辞已修。
- “无条件稳定”已降为“分母恒正/无小分母奇异”等更精确表述。
- Lean 附录表 11 已改为按 `\textwidth` 缩放并固定位置，修复 PDF 右侧越界裁切；字体不支持的 emoji checkmark 也已移除。

### 同步修改的代码/数据文件

```text
code/python/cyclic_reduction.py
code/python/experiments/exp05_true_helmholtz_resonance.py
code/python/experiments/results/exp05_resonance.csv
code/python/gmres_solver.py
code/python/helmholtz_solver.py
```

重点：

- `exp05_resonance.csv` 表头从 `N` 改为 `n_interior`，避免 PowerShell `Import-Csv` 重名列问题。
- `helmholtz_solver.py` / `gmres_solver.py` 中 stale “subtract boundary contributions” 注释已修为 add/move to RHS。
- `cyclic_reduction.py` 作为 deprecated 原型保留，但边界符号也已改正，避免误用。

## 3. 新对话优先复核清单

请从以下命令开始，不要回到旧版 67 页/16 实验状态：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集
git status --short
```

然后核验第 6 章是否已收束：

```powershell
Select-String -Path thesis\chapters\6_experiments.tex -Pattern '\\subsection|实验七|实验八|实验九|实验十|实验十一|实验十二|实验十三|实验十四|实验十五|实验十六|gmres_convergence|efficiency_comparison|wavenumber_effects|condition_number|facr_parameter|gmres_restart|complexity_scaling|complex_test_heatmap|wavenumber_heatmap'
```

期望结果：

- 只看到 6 个主实验 subsection。
- 不应出现“实验七”到“实验十六”。
- 不应出现旧弱支撑图表引用。

继续核验测试：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集\code\python
python -m pytest -q
```

继续核验 LaTeX：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集
Select-String -Path thesis\main.log -Pattern '^! |undefined|Undefined|multiply defined|Output written|Warning'
```

如修改了任何 `.tex` 文件，重新编译：

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集\thesis
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
```

## 4. 当前仍需人工决定

1. 是否把 `NEXT_DIALOG_EXPERIMENT_RESTRUCTURE.md` 提交进仓库。
2. 是否保留、删除或忽略未跟踪参考文件：

```text
PROJECT_REPORT_V1.0.md
Problem_To_Update.md
```

3. 是否处理 Lean 依赖目录 `.lake/packages/proofwidgets` 的 local changes 提示。当前建议不纳入本轮论文终稿提交。

4. 是否清理旧版未引用图片。当前建议：不要急着删除，先保留历史材料，只确保正文不再引用。

## 5. 建议提交范围

论文终稿提交建议包含：

```powershell
git add FINAL_POLISH_REPORT.md
git add code/python/cyclic_reduction.py code/python/gmres_solver.py code/python/helmholtz_solver.py
git add code/python/experiments/exp05_true_helmholtz_resonance.py code/python/experiments/results/exp05_resonance.csv
git add thesis/main.tex thesis/main.pdf
git add thesis/chapters/1_introduction.tex thesis/chapters/2_math_preliminary.tex thesis/chapters/3_fft_direct.tex thesis/chapters/4_helmholtz.tex thesis/chapters/6_experiments.tex thesis/chapters/7_conclusion.tex thesis/chapters/appendix_lean4.tex
```

交接文档可选：

```powershell
git add NEXT_DIALOG_EXPERIMENT_RESTRUCTURE.md
```

不建议默认提交：

```text
PROJECT_REPORT_V1.0.md
Problem_To_Update.md
```

## 6. 给下一位审查模型的关键提醒

不要只复核 P0/P1 小修。当前最重要的是判断：

- 第 6 章是否已经从实验堆叠收束为少而硬的闭环实验。
- 摘要、引言、结论是否都与 6 个核心实验一致。
- 是否仍有“系统性能比较”“大范围实验对比”“无条件稳定”等过强结论残留。
- exp04 是否明确是 Gaussian 非特征 RHS。
- exp05 是否明确是 true Helmholtz near-resonance，不与 modified Helmholtz 混写。

当前判断：实验删减收束已经完成，后续工作应以最终一致性复核、提交裁决和必要的小修为主，不应再扩展新实验。
