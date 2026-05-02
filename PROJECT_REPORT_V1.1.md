# 项目阶段性汇报 V1.1

**项目名称**：矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究  
**汇报日期**：2026-04-28  
**阶段**：Phase 1~3 完成，Phase 4 待启动  
**版本变更**：V1.0 → V1.1（新增项目健康检查、文件清单、冒烟测试验证）

---

## 一、项目概述

本项目研究矩形区域上 Poisson 方程（σ=0）、Modified Helmholtz 方程（σ>0）与 True Helmholtz 方程（σ<0）的快速数值求解方法。核心贡献在于：

1. **σ 统一框架**：将三类方程统一为 `(-∇² + σ)u = f`，所有求解器通过 `sigma` 参数统一调用
2. **FFT9 四阶紧致差分**：实现九点紧致差分格式的频域快速求解，通过 Fourier symbol / Taylor 展开层面的四阶一致性推导，并结合 SymPy 验证与数值收敛实验确认四阶精度；对**三参数九点修正族**的六阶不可达性进行了 Fourier 特征值分析
3. **完整求解器套件**：FA / CR / FACR (l) / FFT9 / GMRES，覆盖 Dirichlet + Neumann + Mixed 边界条件
4. **Lean 4 辅助代数验证**：对三参数九点修正族六阶不可行性的关键代数步骤进行了机器可检查验证（覆盖代数恒等式与系数矛盾，不宣称完整形式化 PDE 截断误差或全局收敛理论）

---

## 二、已完成工作

### ✅ Phase 1：LaTeX 框架与公式修正（2026-04-27 完成）

| 编号 | 修正项 | 状态 |
|------|--------|------|
| A1 | 非齐次 Dirichlet BC 符号 `F-=bc/h²` → `F+=bc/h²`（20处） | ✅ 已修复，测试通过 |
| A2 | Kronecker B 矩阵零元素清理 | ✅ 已清理 |
| A3+A4 | σ 参数统一框架，`_resolve_sigma()`，k² 向后兼容 | ✅ 已实现 |
| A5 | `L_h ≈ Δ` 统一表述 | ✅ 已统一 |
| A6 | FACR 复杂度描述降级（不声称 O(N²loglogN)） | ✅ 已修正 |
| A7 | FFT9 方案B 边界处理统一 | ✅ 已统一 |
| A8 | GMRES λ_min 修正（Dirichlet 2π²+σ，Neumann σ） | ✅ 已修正 |
| A9 | DCT-I 完整公式补充 | ✅ 已补充 |
| A10 | 六阶证明路径修正（c₂=0 → α=-γ → c₄ 矛盾） | ✅ 已完成 |
| A11 | 条件数公式 4/π² | ✅ 已修正 |
| A12 | R_h 边界值说明 | ✅ 已补充 |
| A13 | Neumann 兼容性检查 + Poisson 零模式不误报共振 | ✅ 已修复 |

**产出**：`thesis/` 全文 59 页，MiKTeX 编译 0 错误。

---

### ✅ Phase 2：FFT9 数学推导深化（2026-04-28 完成）

#### Ch2 扩充（§2.3 九点紧致差分格式）

| 子节 | 内容 | 行数 |
|------|------|------|
| 2.3.1 | R_h 算子 Taylor 展开完整推导（h² 项消去机制） | ~40行 |
| 2.3.1 | L_h stencil 逐项推导 `[1,4,1;4,-20,4;1,4,1]/(6h²)` | ~25行 |
| 2.3.2 | 四阶精度定理证明（Fourier 特征值 h 幂级数展开到 O(h⁶)） | ~30行 |
| 2.3.3 | 六阶不可行性定理证明（c₄ 系数显式，对不同模式矛盾） | ~15行 |

#### Ch3 扩充（§3.4 FFT9 四阶紧致差分法）

| 子节 | 内容 | 行数 |
|------|------|------|
| 3.4.1 | 广义特征值问题 `(L_h - σR_h)u = -R_h f` 频域对角化完整推导 | ~40行 |
| 3.4.2 | 边界条件的 R_h 传播效应形式化描述 | ~20行 |
| 3.4.4 | 奇偶约化（OER）的块三对角视角 | ~40行 |

#### 数值验证

- **`verify_fft9_expansion.py`**（SymPy）：验证 h² 项消去，O(h⁴) 精度确认 ✅

---

### ✅ Phase 3：代码实现与实验验证（2026-04-28 完成）

#### 求解器代码

| 文件 | 大小 | 内容 | 状态 |
|------|------|------|------|
| `helmholtz_solver.py` | 63KB | FA / CR / FACR / FFT9，统一 σ 框架，Dirichlet+Neumann+Mixed BC | ✅ |
| `cyclic_reduction.py` | 19KB | CR, FA, FACR(l) 求解器，向量化实现 | ✅ (已修改) |
| `fft9_complete.py` | 20KB | FFT9 四阶/谱精度/OER 求解器 | ✅ |
| `gmres_solver.py` | 35KB | GMRES + Givens 旋转，多项式测试问题 | ✅（3个bug已修复）|

#### 实验脚本（5个 + utils）

| 实验 | 文件 | 大小 | 目的 | 状态 |
|------|------|------|------|------|
| exp00 | `exp00_fft_vs_sparse.py` | 10KB | FFT 直接法 vs SciPy sparse 一致性验证 | ✅ 28/28 PASS |
| exp01 | `exp01_convergence.py` | 6KB | Poisson + Modified H. + True H. 三轨收敛阶 | ✅ 二阶/四阶符合理论 |
| exp02 | `exp02_nonhom_bc.py` | 7KB | 非齐次 Dirichlet BC 验证（A1 bug fix） | ✅ 全部 PASS |
| exp03 | `exp03_neumann_mixed.py` | 9KB | Neumann & Mixed BC 验证 | ✅ 全部 PASS |
| exp04 | `exp04_modified_vs_true.py` | 11KB | Modified vs True Helmholtz 谱对比 | ✅ 运行成功 |
| utils | `utils.py` | 5KB | 6个测试问题 + 辅助函数 | ✅ |

#### 实验数据（CSV）

| 文件 | 大小 | 行数 |
|------|------|------|
| `exp00_fft_vs_sparse.csv` | 2.3 KB | 29 |
| `exp01_convergence.csv` | 4.6 KB | 61 |
| `exp02_nonhom_bc.csv` | 4.2 KB | 65 |
| `exp03_neumann_mixed.csv` | 5.8 KB | 61 |
| `exp04_modified_vs_true.csv` | 2.5 KB | 22 |

#### 实验图表（PNG）

| 文件 | 大小 | 内容 |
|------|------|------|
| `exp01_convergence.png` | 134 KB | 三轨收敛曲线 |
| `exp02_nonhom_bc.png` | 157 KB | 非齐次 BC 误差曲线 |
| `exp03_neumann_mixed.png` | 103 KB | Neumann/Mixed 收敛阶 |
| `exp04_gmres_iters_vs_sigma.png` | 43 KB | GMRES 迭代次数对比 |
| `exp04_min_denom_vs_sigma.png` | 61 KB | 频域分母最小值 |
| `exp04_spectral_indicator_vs_sigma.png` | 65 KB | 谱条件数指标 |

---

### ✅ Lean 4 形式化验证

| 文件 | 大小 | 内容 | 状态 |
|------|------|------|------|
| `SixthOrderImpossibility.lean` | 5KB | Int 版本，12 定理，`omega` 证明 | ✅ 全部验证通过 |
| `MathlibTest.lean` | 6KB | ℝ 版本，8 定理，`ring` + `linarith` 证明 | ✅ 全部验证通过 |

**关键定理**（ℝ 版本）：
1. `c4_num_mode_10_real`：c₄(1,0,γ) = 60γ + 3
2. `c4_num_mode_11_real`：c₄(1,1,γ) = 120γ - 4
3. `sixth_order_impossible_real`：★ 无实数 γ 使所有模式 c₄=0
4. `poisson_no_sixth_order_real`：★ Poisson 六阶不可行（ℝ）
5. `helmholtz_no_sixth_order_real`：★ Helmholtz 六阶不可行（ℝ，一般 k²）

---

## 三、项目健康检查（V1.1 新增）

### 冒烟测试

| 测试 | 结果 |
|------|------|
| `verify_fft9_expansion.py` | ✅ 通过（h² 项消去，O(h⁴) 确认） |
| `exp00_fft_vs_sparse` | ✅ 28/28 PASS |
| `exp01_convergence` | ✅ 运行成功，收敛阶符合理论 |

### TODO/FIXME 扫描

| 扫描范围 | 结果 |
|----------|------|
| 所有 `.py` / `.tex` / `.md` 文件 | 仅 2 个 TODO，均在 `code/python/_archive/fft9_true.py`（归档代码，非活跃） |
| 活跃代码 | ✅ 0 TODO / 0 FIXME / 0 XXX |

### LaTeX 编译状态

| 指标 | 结果 |
|------|------|
| 错误 | **0** |
| 重复标签 | **0**（`tab:modified_vs_true` 已修复） |
| 警告 | 21 个（全部 `hyperref` 无害警告：章节标题中数学公式无法放入 PDF 书签） |
| PDF 大小 | 3.59 MB |
| 最后编译 | 2026-04-28 17:18 |

### Git 状态

| 类型 | 文件 |
|------|------|
| **已修改** | `ADVANCED_ROADMAP.md`, `code/python/cyclic_reduction.py`, `notes/Cyclic_Reduction_FACR_精读.md`, `thesis/chapters/2_math_preliminary.tex`, `thesis/chapters/3_fft_direct.tex`, `thesis/chapters/6_experiments.tex`, `thesis/main.pdf` |
| **未跟踪** | `PROJECT_REPORT_V1.0.md`, `PROJECT_REPORT_V1.1.md`, `code/python/experiments/`, `code/python/verify_fft9_expansion.py`, `thesis/figures/exp*.png`, `thesis/chapters/figures/`, `thesis/compile_*.txt`, `code/python/tests/_test_output*.txt` |
| **当前分支** | `revise-sigma-fft9-krylov` |

---

## 四、论文当前状态

### 论文章节

| 章节 | 文件 | 大小 | 行数 | 内容 | 状态 |
|------|------|------|------|------|------|
| Ch1 | `1_introduction.tex` | 7 KB | 117 | 引言 | ✅ |
| Ch2 | `2_math_preliminary.tex` | 32 KB | 668 | 数学基础（已深化） | ✅ |
| Ch3 | `3_fft_direct.tex` | 19 KB | 400 | FFT 直接法（已深化） | ✅ |
| Ch4 | `4_helmholtz.tex` | 11 KB | 257 | Helmholtz 方程 | ✅ |
| Ch5 | `5_gmres.tex` | 11 KB | 241 | GMRES | ✅ |
| Ch6 | `6_experiments.tex` | 31 KB | 736 | 15 个实验（含新增十四、十五） | ✅ |
| Ch7 | `7_conclusion.tex` | 7 KB | 119 | 结论 | ✅ |

### 论文指标

- **总页数**：~60 页
- **编译状态**：`main.pdf` 3.59 MB，0 错误
- **图片总数**：21 张 PNG（14 原有 + 7 新增实验图）
- **表格总数**：~15 个
- **参考文献**：41 条

### 参考文献清单（41 条）

| 类别 | 条目数 | 代表文献 |
|------|--------|----------|
| 核心算法 | 8 | SaadSchultz1986, Houstis1979a/b, Swarztrauber1977, Strang1999, Hockney1965, Buzbee1970, Temperton1979, SwarztrauberSweet1979 |
| 紧致差分 | 3 | Gupta1984, Spotz1998, Sutmann2007, Lynch1977 |
| Helmholtz | 3 | Bayliss1983, Ernst2011, Ihlenburg1998 |
| 迭代法 | 3 | Saad2003, Walker1988, Niu2009 |
| 中文文献 | 13 | 闵涛2017, 邵文婷2017, 刘镓源2019 等 |
| 经典教材 | 6 | LeVeque2007, Trefethen2000, GolubVanLoan2013, Quarteroni2007, Samarskii2001, Briggs1995 |

---

## 五、待完成工作

### ⏳ Phase 4：最终打磨

| 任务 | 说明 | 优先级 | 预估时间 |
|------|------|--------|----------|
| 全文编译验证 | 交叉引用、参考文献标注全验证 | 高 | 0.5h |
| 语言润色 | 中英术语统一，语句通顺 | 中 | 2h |
| 图表风格统一 | 字体、颜色、线型统一 | 中 | 1h |
| Lean 4 附录整理 | 将形式化验证结果整理为附录 | 低 | 2h |

### ⏳ Phase 5：答辩准备

| 任务 | 说明 | 优先级 | 预估时间 |
|------|------|--------|----------|
| 答辩 PPT | 制作 15-20 页答辩幻灯片 | 高 | 4h |
| 演示代码 | 准备现场演示脚本 | 中 | 1h |
| 常见问题 | 准备答辩 Q&A | 中 | 2h |

---

## 六、关键技术决策记录

### 决策 1：σ 统一框架

**决策**：将所有方程统一为 `(-∇² + σ)u = f`  
**理由**：
- Modified Helmholtz：`σ = +κ²`（Yukawa 方程），矩阵正定
- True Helmholtz：`σ = -k²`，矩阵不定，存在共振风险
- Poisson：`σ = 0`，作为特例自然包含

**实现**：所有 solver 函数接受 `sigma` 参数（优先）或 `k2`（向后兼容）

### 决策 2：六阶不可行性证明策略

**决策**：采用 Fourier 特征值分析（数值验证 + Lean 4 形式化）  
**理由**：
- 纯代数方法需要处理大量符号运算，易出错
- Fourier 分析可以直接计算 c₄ 系数对不同模式 (p,q) 的依赖
- Lean 4 形式化提供了机器可检查的证明

### 决策 3：实验设计策略

**决策**：每个求解器 + 每个边界条件 + 每个方程类型 = 独立实验  
**理由**：
- 系统性验证 σ 统一框架的正确性
- 每个 bug fix 都有对应的验证实验（如 exp02 验证 A1 fix）
- 实验代码与论文章节一一对应，便于维护

---

## 七、已知问题清单

| 问题 | 影响 | 解决方案 | 状态 |
|------|------|----------|------|
| `exp04` 中 modified Helmholtz 部分误差计算无意义（用错测试函数）| 数据记录为 `nan`，不影响其他指标 | 改用 `np.nan`，或实现精确解 | ✅ 已修复 |
| `tab:modified_vs_true` 标签重复定义 | LaTeX 编译警告 | Ch6 中重命名为 `tab:exp15_mod_vs_true` | ✅ 已修复 |
| `thesis/chapters/figures/` 目录冗余 | 图片重复存储 | 已同步到 `thesis/figures/`；建议删除 `chapters/figures/` | ⚠️ 待清理 |

---

## 八、下一步计划

### 近期（1-2 天）

1. **清理冗余目录**：删除 `thesis/chapters/figures/`（图片已复制到正确位置）
2. **语言润色**：检查全文术语一致性，优化表达方式
3. **图表优化**：统一实验图表风格

### 中期（3-5 天）

1. **答辩 PPT 制作**：制作 15-20 页答辩幻灯片
2. **演示代码准备**：准备现场演示脚本
3. **Lean 4 附录**：将形式化验证结果整理为论文附录

### 长期（1-2 周）

1. **代码发布准备**：整理代码仓库，编写 README 和使用文档
2. **投稿准备**：根据目标期刊调整论文格式
3. **最终验收**：全文通读，确保逻辑连贯、无错别字

---

## 九、附录：完整文件清单

### Python 代码（20 个文件）

```
code/python/
├── correct_poisson_solver.py              7 KB
├── cyclic_reduction.py                   19 KB  [MODIFIED]
├── derivation_audit_v2.py                14 KB
├── derivation_audit_v3.py                14 KB
├── fft9_complete.py                      20 KB
├── fft9_eigenvalue_analysis.py           12 KB
├── fft9_eigenvalue_proof.py               5 KB
├── gmres_diagnosis.py                     4 KB
├── gmres_solver.py                       35 KB
├── gmres_validate.py                      4 KB
├── gmres_validate_fixed.py                5 KB
├── helmholtz_solver.py                   63 KB
├── poisson_fft_solver.py                  8 KB
├── test_exp3_gmres.py                     4 KB
├── test_neumann_dct.py                    4 KB
├── test_new_problems.py                   5 KB
├── thesis_derivation_audit.py            27 KB
├── thesis_final_verification.py           5 KB
├── verify_fft9_expansion.py               5 KB  [NEW]
├── working_poisson_solver.py              9 KB
└── experiments/
    ├── __init__.py                          0 KB
    ├── utils.py                             5 KB
    ├── exp00_fft_vs_sparse.py              10 KB
    ├── exp01_convergence.py                 6 KB
    ├── exp02_nonhom_bc.py                   7 KB
    ├── exp03_neumann_mixed.py               9 KB
    ├── exp04_modified_vs_true.py           11 KB
    ├── results/                          (5 CSV)
    └── figures/                          (6 PNG)
```

### Lean 4 形式化

```
code/lean4_formalization/
├── SixthOrderImpossibility.lean           5 KB  (Int, 12 定理)
├── SixthOrderImpossibility/
│   └── MathlibTest.lean                   6 KB  (ℝ, 8 定理)
└── .lake/                              (Mathlib 依赖)
```

### 论文 LaTeX

```
thesis/
├── main.tex                              14 KB
├── chapters/
│   ├── 1_introduction.tex                 7 KB  (117 lines)
│   ├── 2_math_preliminary.tex            32 KB  (668 lines) [MODIFIED]
│   ├── 3_fft_direct.tex                  19 KB  (400 lines) [MODIFIED]
│   ├── 4_helmholtz.tex                   11 KB  (257 lines)
│   ├── 5_gmres.tex                       11 KB  (241 lines)
│   ├── 6_experiments.tex                 31 KB  (736 lines) [MODIFIED]
│   ├── 7_conclusion.tex                   7 KB  (119 lines)
│   └── figures/                        (7 PNG, 冗余待清理)
├── figures/                            (21 PNG)
├── references.bib                      (41 条目)
└── main.pdf                           3.59 MB
```

### 笔记与文档

```
notes/
├── Cyclic_Reduction_FACR_精读.md          7 KB  (212 lines)
├── FFT9_精读总结报告.md                    18 KB  (599 lines)
├── FFT9_论文1_精读.md                     10 KB  (327 lines)
├── FFT9_论文2_算法实现_精读.md               12 KB  (441 lines)
├── GMRES_算法详解.md                       9 KB  (306 lines)
├── GMRES_论文精读.md                      13 KB  (445 lines)
└── lean4_formalization.md                 3 KB  (76 lines)
```

---

**报告结束** | 生成时间：2026-04-28 17:22 | 版本：V1.1 | 上次更新：V1.0 (2026-04-28 17:00)
