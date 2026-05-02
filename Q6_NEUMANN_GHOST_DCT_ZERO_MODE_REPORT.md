# Q6 Neumann Ghost-Point / DCT-I / Zero-Mode 修复报告

## 1. 修复目标

本轮只处理 Neumann ghost-point 矩阵、DCT-I 对称化缩放、纯 Neumann Poisson compatibility 与 zero-mode normalization 的表述和诊断一致性。
未新增实验，未修改实验参数，未重绘图，未修改 CSV/PNG/PDF 实验资产。

## 2. 论文修改

- `thesis/chapters/2_math_preliminary.tex`
  - 保留三物理节点 ghost-point 矩阵
    \[
    G=\frac1{h^2}\begin{pmatrix}2&-2&0\\-1&2&-1\\0&-2&2\end{pmatrix}.
    \]
  - 补充 \(G\mathbf 1=0\)，说明常数向量是纯 Neumann Poisson 的零模态。
  - 补充 \(S=D^{-1}GD\) 的三节点对称矩阵形式。
  - 明确一维缩放变量
    \(\tilde u=D^{-1}u,\ \tilde f=D^{-1}f,\ S\tilde u=\tilde f,\ u=D\tilde u\)。
  - 将二维流程改为矩阵形式：
    \[
    G_xU+UG_y+\sigma U=F,\quad
    \tilde U=D_x^{-1}UD_y^{-1},\quad
    \tilde F=D_x^{-1}FD_y^{-1}.
    \]
  - 频域求解统一写成
    \[
    \widehat{\tilde U}_{k,l}
    =\widehat{\tilde F}_{k,l}/(\lambda_k+\lambda_l+\sigma),
    \quad U=D_x\tilde U D_y.
    \]

- `thesis/chapters/4_helmholtz.tex`
  - 明确 DCT-I 变换作用在对称化变量 \(\tilde U\) 上。
  - 纯 Neumann Poisson compatibility 写为 trapezoid-weighted condition。
  - zero mode normalization 改为 \(\widehat{\tilde U}_{0,0}=0\)，并说明其等价于恢复后解的 weighted mean zero。
  - 区分三种情形：
    - pure Neumann Poisson: \(d_{0,0}=0\)，需要 compatibility；
    - modified Helmholtz: \(d_{0,0}=\kappa^2>0\)，不需要 compatibility；
    - true Helmholtz: 零模态分母为 \(-\kappa^2\)，但其他模态仍可能 near-resonance。

- `thesis/chapters/6_experiments.tex`
  - 实验二补充说明：weighted mean diagnostic 对应 \(\widehat{\tilde U}_{0,0}=0\)，不是普通算术平均。

## 3. 代码与测试修改

- `code/python/helmholtz_solver.py`
  - Neumann 注释统一为当前实际流程：
    \(F_{\rm tilde}=D_x^{-1}FD_y^{-1}\)，不再把 FA/DCT-I 路径写成额外 \(h^2\) 缩放。
  - `check_neumann_compatibility` 仅对 pure Neumann Poisson 执行；mixed Dirichlet/Neumann Poisson 不再被误判为需要纯 Neumann compatibility。
  - 对已缩放 RHS `F_adj = F / d_x / d_y`，compatibility sum 改为
    `sum(F_adj / d_x / d_y)`，即
    `sum(F / d_x^2 / d_y^2)`，与 trapezoid endpoint weights 一致。
  - FA/CR/FACR-like 的 pure Neumann Poisson 路径均调用同一 compatibility 诊断。

- `code/python/tests/test_05_neumann_compatibility.py`
  - 将直接调用 `check_neumann_compatibility` 的测试改为传入已缩放 RHS。
  - 新增非对称 RHS 测试：trapezoid-weighted integral 为零但普通求和不为零，确认诊断使用 DCT-I/trapezoid 权重。
  - 新增 mixed Poisson 测试：含 Dirichlet 方向时不触发 pure Neumann compatibility warning。

## 4. 验证结果

- `python -m py_compile code/python/helmholtz_solver.py`: passed.
- `cd code/python; python -m pytest -q tests/test_01_eigenvalues.py tests/test_05_neumann_compatibility.py`: 38 passed.
- `cd code/python; python -m pytest -q`: 105 passed, 2 skipped.
- LaTeX 完整编译 `xelatex -> bibtex -> xelatex -> xelatex`: passed.
- `main.log` 扫描：
  - no fatal error;
  - no undefined reference/citation;
  - no multiply-defined label;
  - no overfull hbox.
- `thesis/main.pdf` 正常生成，页数为 69 pages.

## 5. 未修改范围

- 未改变 Neumann ghost-point 离散方向。
- 未重构 DCT-I 求解算法。
- 未修改 FFT9、Dirichlet、GMRES 或实验脚本。
- 未重跑大规模实验。
- 未修改任何 CSV、PNG 或 PDF 实验图资产。
