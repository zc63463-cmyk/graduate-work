# Q5 三参数 FFT9 六阶不可达性范围收紧报告

## 1. 修复目标

本轮修复的目标是收紧“三参数九点修正族六阶不可达性”的理论口径与 Lean 4 覆盖范围。核心结论保持不变，但明确限定为：

- 三参数九点修正族；
- \(\alpha,\beta,\gamma\) 为与 Fourier mode 无关的统一常数；
- local Fourier symbol / interior consistency 意义下的六阶不可达性。

该结论不表述为所有九点格式、完整 Dirichlet 边值问题全局误差，或完整 PDE 理论的六阶不可达性。

## 2. 论文修改内容

- `thesis/chapters/2_math_preliminary.tex`
  - 将定理标题改为“三参数九点修正族的六阶 local consistency 不可达性”。
  - 定理陈述中显式加入 mode-independent constants 与 interior Fourier symbol 限定。
  - 保留 Poisson 与 Helmholtz 的 \(c_2 \to c_4\) 反证逻辑。
  - 补充说明 \((\xi,\eta)=(1,0)\) 是 local Fourier symbol 的连续波数测试，不是 Dirichlet DST-I 离散模态编号。
  - 将“四阶是自然精度上限”收紧为“在该三参数九点修正族、mode-independent constants 和 local Fourier symbol 意义下”的自然上限。

- `thesis/chapters/1_introduction.tex`
  - 将贡献点从宽泛的“六阶不可行性/四阶最优选择”改为局部符号层面的 restricted-family 结论。

- `thesis/chapters/7_conclusion.tex`
  - 将结论同步限定为三参数九点修正族的六阶 local consistency 不可达性。
  - 明确不排除更大模板、更多参数或特殊边界闭合下的高阶格式。

- `thesis/chapters/appendix_lean4.tex`
  - 将 Lean 4 描述改为“验证反证中的代数核心”。
  - 表格中的 theorem 含义改为 local symbol 六阶条件不相容。
  - 明确未覆盖 Fourier symbol 推导、Taylor 展开建模、DST/DCT 谱理论、边界闭合、完整 PDE 全局误差估计和 Python 程序正确性。

## 3. Lean 注释修改

- `code/lean4_formalization/SixthOrderImpossibility.lean`
  - 收紧文件头注释，说明 Lean 文件只验证 restricted three-parameter family 的代数核心。
  - 将“9-point compact finite difference format cannot achieve 6th-order accuracy”类注释改为 restricted local Fourier-symbol obstruction。

- `code/lean4_formalization/SixthOrderImpossibility/MathlibTest.lean`
  - 将主定理注释改为 algebraic core / restricted-family local-symbol obstruction。
  - 未修改 theorem 名称、声明结构或证明代码。

## 4. 未修改范围

- 未修改 FFT9、五点、GMRES、Neumann/mixed 等求解算法；
- 未修改实验脚本、CSV、PNG/PDF 图像；
- 未新增实验；
- 未改变 Lean theorem 的证明逻辑；
- 未把 Lean 4 覆盖范围扩展到 PDE 全局误差、谱理论或程序正确性。

## 5. 静态检查

已搜索以下风险表述：

- `所有九点格式不可能`
- `九点格式不可能达到六阶`
- `Lean 4 证明了`
- `Lean 4 证明完整`
- `完整 PDE 六阶`
- `完整 PDE 收敛定理`
- `自然精度上限`
- `六阶不可行`
- `六阶不可达`

检查结果：未发现无条件声称“所有九点格式不可能六阶”或“Lean 4 证明完整 PDE 六阶不可达”的表述。保留的“自然精度上限”已经带有三参数族、mode-independent constants 与 local Fourier symbol 限定。

## 6. 验证状态

- `lake build`：通过。
  - Lean 输出保留若干既有 `unused simp argument` warning。
  - `proofwidgets` 依赖目录提示存在 local changes；不影响本项目自有 Lean 文件构建通过。
- 自有 Lean 文件逃逸词检查：未发现 `sorry`、`admit`、`axiom` 或 `unsafe`。
- LaTeX 完整编译链 `xelatex -> bibtex -> xelatex -> xelatex`：通过。
- `main.log` 扫描：无 fatal error、undefined reference、multiply-defined label、overfull hbox；仅保留既有 SimSun 字体警告。
- `main.pdf` 正常生成，共 68 页。

## 7. 结论

Q5 完成后，论文中的六阶不可达性结论已经限定为：

\[
\text{三参数九点修正族}
+
\text{mode-independent constants}
+
\text{local Fourier symbol / interior consistency}.
\]

Lean 4 的角色也已限定为验证该反证中的 \(c_2/c_4\) 多项式矛盾，不再被表述为完整 PDE、全局误差、边界闭合或程序正确性的形式化证明。
