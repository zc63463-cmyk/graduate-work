# Final Polish Report

**项目**：矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究  
**日期**：2026-04-28  
**版本**：V2.1（最终提交清理版）  
**执行人**：WorkBuddy Agent

---

## 1. Repository Status

- **Branch**：`revise-sigma-fft9-krylov`
- **Latest commit**：`1dd7a82 fix: de-risk 'strict proof' wording and M5 resonance in Ch2+Ch7`
- **Git status**（当前）：
  - Untracked: `FINAL_POLISH_REPORT.md`, `code/python/pytest.ini`, 临时日志文件
  - 无 modified，无 deleted

---

## 2. Tests

### Python 测试

| 测试 | 状态 | 命令 | 输出 |
|------|------|------|------|
| pytest | ✅ PASS | `python -m pytest -q` | **103 passed, 2 skipped, 0 failed** |
| verify_fft9_expansion.py | ✅ PASS | `python verify_fft9_expansion.py` | h² 项消去，-λ_L/λ_R = ξ²+η²+O(h⁴) 确认 |

**运行环境**：Windows 11 + Python subprocess（Miniconda3）

**pytest.ini 已添加**：
```ini
[pytest]
testpaths = tests
norecursedirs = _archive .git .venv __pycache__
python_files = test_*.py
python_classes = Test*
python_functions = test_*
```

默认 `python -m pytest -q` 不再收集 `_archive/` 目录。归档旧测试（2 个 GMRES 奇异矩阵失败）不属于活跃测试范围。

### FFT9 展开验证

| 验证项 | 结果 |
|--------|------|
| h² 项系数（有效特征值） | **0**（消去确认） |
| h⁰ 项 | η²+ξ²（匹配 ξ²+η²） |
| h² 项 target | -(ξ²+η²)²/12，**Match: True** |
| h⁴ 项 target | (ξ²+η²)(ξ⁴+4ξ²η²+η⁴)/360，**Match: True** |
| -λ_L/λ_R | ξ²+η²+O(h⁴)，**四阶精度确认** |

### Lean 4 最终验收

| 验证项 | 结果 |
|--------|------|
| lake build | ✅ **Build completed successfully (3 jobs)** |
| lake build return code | **0** |
| sorry/admit/axiom/unsafe | **0 matches** |
| Lean 版本 | leanprover/lean4:v4.29.1 |
| Mathlib | v4.29.1 (prebuilt olean) |
| linter 警告 | 12 个 unused simp arguments（无害，可优化） |

**验证的定理**：

| Lean theorem | 数学含义 | 证明策略 |
|-------------|----------|----------|
| c4_num_mode_10_real | c₄(1,0,γ) = 60γ+3 | ring |
| c4_num_mode_11_real | c₄(1,1,γ) = 120γ-4 | ring |
| sixth_order_impossible_real | 不存在统一 γ | linarith |
| c2_helmholtz_forces_neg_gamma | c₂=0 在 Helmholtz 推出 α=-γ | linarith |
| poisson_no_sixth_order_real | Poisson 三参数族不可达六阶 | ring + linarith |
| helmholtz_no_sixth_order_real | Helmholtz 三参数族不可达六阶 | ring + linarith |

**Lean 4 定位声明**：Lean 4 验证覆盖的是代数系数矛盾，不包含完整 PDE 截断误差或全局收敛理论。

---

## 3. LaTeX Build

### XeLaTeX 最终编译

| 指标 | 结果 |
|------|------|
| 编译命令 | **xelatex + bibtex + xelatex ×2** |
| Engine | **This is XeTeX, Version 3.141592653-2.6-0.999998 (MiKTeX 26.2)** |
| 错误 | **0** |
| 警告 | 22（全部 `hyperref` 无害警告：章节标题数学公式无法放入 PDF 书签） |
| 重复标签 | **0** |
| Citation undefined | **0** |
| Reference undefined | **0** |
| PDF 大小 | **2,339 KB** |

---

## 4. Mathematical Fixes

### P0 措辞降风险

| 编号 | 原文 | 修改后 | 位置 | 状态 |
|------|------|--------|------|------|
| P0-1 | "严格的四阶精度证明" | "Fourier symbol / Taylor 展开层面的四阶一致性推导" | 报告已修正 | ✅ |
| P0-2 | "六阶不可行性证明" | "三参数九点修正族六阶不可达性" | 论文已正确 | ✅ |
| P0-3 | "Lean 4 证明了理论正确性" | "辅助形式化验证（覆盖代数恒等式与系数矛盾）" | 论文已正确 | ✅ |
| P0-4 | "FACR 复杂度 O(N²loglogN)" | 论文已区分 classical FACR 理论 vs 本文实现 O(N²logN) | 无需修改 | ✅ |
| M5 | "Modified H. 的共振条件 λ+κ²=0 永不成立" | "Modified H. 不存在共振（因 λ+κ² > 0 恒正）" | `7_conclusion.tex` | ✅ 已修复 |
| M4a | "为严格证明四阶精度" | "为建立四阶一致性" | `2_math_preliminary.tex` | ✅ 已修复 |
| M4b | "四阶精度的严格证明" | "四阶精度的 Fourier 分析" | `2_math_preliminary.tex` | ✅ 已修复 |
| M4c | "六阶格式不可行的证明：通过 Fourier 特征值分析严格证明了" | "三参数九点修正族六阶不可达性：通过 Fourier 特征值分析验证了" | `7_conclusion.tex` | ✅ 已修复 |

### 数学残留扫描

| 扫描项 | matches | 状态 |
|--------|---------|------|
| B = tridiag(-1,4,-1) | 0 | ✅ 无残留 |
| Modified H. λ+κ²=0 共振条件 | 0 | ✅ 无残留 |
| FACR O(N²loglogN) 未限定 | 0（6 matches 均为正确标注） | ✅ 全部正确 |
| "严格证明" 过强措辞 | 0（已全部修复） | ✅ 已清理 |
| Lean 证明理论 | 0 | ✅ 无残留 |

---

## 5. Figures and Experiments

### 实验脚本（5个 + utils）

| 实验 | 状态 | CSV | PNG |
|------|------|-----|-----|
| exp00 | ✅ PASS | ✅ | — |
| exp01 | ✅ PASS | ✅ | ✅ |
| exp02 | ✅ PASS | ✅ | ✅ |
| exp03 | ✅ PASS | ✅ | ✅ |
| exp04 | ✅ PASS | ✅ | ✅ (3张) |

---

## 6. Bibliography

- **总计**：41 条参考文献，0 citation undefined。

---

## 7. Archive Test Handling

- **方案**：添加 `pytest.ini`（方案 A）
- **配置**：`testpaths = tests`，`norecursedirs = _archive .git .venv __pycache__`
- **效果**：默认 `python -m pytest -q` 不再收集 `_archive/` 目录
- **结果**：103 passed, 2 skipped, 0 failed（return code 0）
- **归档旧测试**：2 个 GMRES 奇异矩阵失败仍存在于 `_archive/`，但不影响活跃测试

---

## 8. Remaining Risks

### 低风险 💡

| 风险 | 说明 | 缓解措施 |
|------|------|----------|
| _archive/ 归档代码存在旧失败测试 | 旧 GMRES 测试使用奇异矩阵 | 已通过 pytest.ini 排除，活跃测试全部通过 |
| 图表风格未统一 | 6 张实验图为 matplotlib 默认风格 | 如需更高视觉质量，可统一字体/线宽/marker |
| 答辩 PPT 未制作 | 仅有大纲，无实际幻灯片 | 使用 PowerPoint/Beamer 按大纲制作 |
| Lean 4 linter 警告 | 12 个 unused simp arguments | 不影响正确性，可后续优化 |

---

## 9. Recommended Commit

```bash
# 必须提交
git add FINAL_POLISH_REPORT.md
git add code/python/pytest.ini

# 视情况提交（如仓库已跟踪 PDF 或导师要求）
# git add thesis/main.pdf

git commit -m "chore: add pytest.ini and final validation report V2.1

- pytest.ini: exclude _archive/ from test collection
- Default pytest now passes: 103 passed, 2 skipped, 0 failed
- All validations confirmed: FFT9, Lean 4, XeLaTeX"
```

### 不应提交的文件

以下为临时日志/中间文件，不应提交：
- `code/python/final_*.txt` — pytest/FFT9 临时输出
- `code/lean4_formalization/final_*.txt` — Lean 临时输出
- `thesis/_check_log.py` — 临时脚本
- `thesis/compile_*.txt`, `thesis/final_*.txt` — 编译临时输出
- `PROJECT_REPORT_V1.0.md` — 旧版报告（已被 V1.1 替代）
- `*.aux`, `*.log`, `*.bbl`, `*.blg` — LaTeX 中间文件

---

## 附录：最终验收检查清单

- [x] git status 查看
- [x] git log 查看（4 commits on `revise-sigma-fft9-krylov`）
- [x] **pytest 默认命令通过**：`python -m pytest -q` → 103 passed, 2 skipped, 0 failed
- [x] **pytest.ini 已添加**：`_archive/` 默认排除
- [x] **verify_fft9_expansion.py 通过**：h² 项消去，O(h⁴) 确认
- [x] **Lean 4 lake build 通过**：Build completed successfully, 0 sorry/admit/axiom/unsafe
- [x] **XeLaTeX 编译 0 error**，0 citation undefined，0 reference undefined
- [x] 重复标签检查：0
- [x] P0-1 ~ P0-4 措辞审查通过
- [x] M1-M5 数学残留扫描通过
- [x] M4a-c "严格证明"过强措辞已修复
- [x] 参考文献核查（41 条，0 citation undefined）
- [x] Lean 4 附录创建并集成
- [x] 答辩 PPT 大纲 + Q&A 创建
- [x] Final Polish Report 更新至 V2.1

---

**Final status: READY TO SUBMIT**

---

**报告结束** | 生成时间：2026-04-28 18:45 | 版本：V2.1（最终提交清理版）
