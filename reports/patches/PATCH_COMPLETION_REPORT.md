# Patch Completion Report

生成时间：2026-04-29

## 1. 修复列表

本轮根据 `OPTIONAL_PATCH_PLAN.md` 复核 A 类问题和明确安全的 B 类问题，并在最终 source-of-truth verification 中补齐两个非阻塞 C 类小修。结论：A/B 项已完成，C 类措辞风险已处理。

| Issue | Status | 修改原因 |
|---|---|---|
| A-001 Ch2 FFT9 `L_h` 符号约定 | 已完成 | 保持 `L_h \approx \Delta`，避免把九点 Laplacian 和正的 `-\Delta` 算子混写。当前写法为 `-L_h u = R_h(-\nabla^2 u)+O(h^4)`，并使用 `-L_h = L_h^{(5)} - h^2 C_h/6`。 |
| A-002 Ch2 三参数 `lambda_mod` 符号 | 已完成 | 与主格式 `(-L_h+\sigma R_h)u=R_h f` 保持一致，避免有效特征值分子漏掉负号或 `\beta` 项。 |
| A-003 Lean 附录三参数公式 | 已完成 | 附录公式已与 Ch2 corrected fraction 对齐，并保留 Lean 只覆盖代数核心步骤的范围说明。 |
| A-004 Ch4 true Helmholtz 谱分类 | 已完成 | 旧的“只要 `\kappa^2>\lambda_{\min}` 就不定”不完整；当前按正定、奇异、不定、负定四类描述。 |
| A-004 跨章节 true Helmholtz 范围措辞 | 已完成 | 将宽泛“矩阵不定”改为“可能不定/近共振风险”，避免对所有参数区间过度概括。 |
| B-001 Ch4 near-resonance `delta` 局部符号 | 已完成 | 保持 `\kappa^2=\lambda_{p^*,q^*}+\delta` 约定时，分母应为 `d_{p^*,q^*}=-\delta`，解分量应为 `-\hat f/\delta`。 |
| C-001 Ch2 `h^4` 系数措辞 | 已完成 | 将“对任意非零模式不为零”降级为“该系数不恒为零、误差一般为 `O(h^4)`”，避免未证明的全称断言。 |
| C-002 Ch3 FACR 第三项复杂度措辞 | 已完成 | 将“第三项趋于常数”改为“第三项降为 `O(N^2)` 量级”，并明确当前 FACR-like 实现整体仍为 `O(N^2 log N)`。 |

## 2. 修改文件

按 `OPTIONAL_PATCH_PLAN.md` 涉及并已满足的文件：

- `thesis/chapters/2_math_preliminary.tex`
- `thesis/chapters/appendix_lean4.tex`
- `thesis/chapters/4_helmholtz.tex`
- `thesis/chapters/3_fft_direct.tex`
- `thesis/main.tex`
- `thesis/chapters/1_introduction.tex`
- `thesis/chapters/7_conclusion.tex`
- `README.md`
- `FINAL_POLISH_REPORT.md`

本轮验证重新生成：

- `thesis/main.pdf`

本轮更新报告：

- `PATCH_COMPLETION_REPORT.md`

本轮未新增实验，未恢复旧图，未扩大 Lean 4 贡献范围。

## 3. 是否重跑实验

本轮未重跑实验。

原因：`OPTIONAL_PATCH_PLAN.md` 明确 A/B 修复不需要实验重跑，且本轮限制为不新增实验、不恢复旧图、不大改主线。

## 4. pytest 结果

命令：

```powershell
cd code/python
python -m pytest -q
```

结果：PASS

- 103 passed
- 2 skipped
- 2 warnings

说明：warnings 来自既有 Neumann Poisson compatibility 测试，属于预期兼容性警告，非阻塞。

## 5. verify_fft9 结果

命令：

```powershell
cd code/python
python verify_fft9_expansion.py
```

结果：PASS

- `h^2` coefficient verified as `0`
- 输出结论为 `-lambda_L / lambda_R = xi^2+eta^2 + O(h^4)`

## 6. lake build 结果

命令：

```powershell
cd code/lean4_formalization
lake build
```

结果：PASS

- Build completed successfully
- 3 jobs completed

非阻塞 warning：

- `SixthOrderImpossibility.lean` 有 unused `simp` argument warnings。
- `.lake/packages/proofwidgets` 有 local changes notice。

## 7. XeLaTeX 结果

命令：

```powershell
cd thesis
xelatex -interaction=nonstopmode main.tex
bibtex main
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
```

结果：PASS

- `main.pdf` 成功生成。
- PDF 页数：56 pages。
- `bibtex main` 成功完成。

非阻塞 warning：

- SimSun bold-italic font fallback。
- hyperref PDF bookmark token warnings。
- MiKTeX update notice。

未发现：

- Undefined citation。
- Undefined reference。
- final overfull box warning。

## 8. Remaining risks

- C-001/C-002 已在本轮作为安全 C 类小修处理。
- 未发现需要人工裁决的新数学矛盾，因此未新增 GPT-5.5 Pro 问题。
- 当前工作区仍包含此前审计输出文件和前一轮修复产生的代码/结果文件改动；本报告仅记录本轮基于 `OPTIONAL_PATCH_PLAN.md` 的补查与验证。

## 9. Final status

READY TO SUBMIT

确认：

- 未新增实验。
- 未恢复旧图。
- 未大改论文主线。
- 报告与 README 已同步最终 56 pages PDF。
- 未扩大 Lean 贡献范围。
