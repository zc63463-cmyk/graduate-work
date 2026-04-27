# 完备高可执行计划：PDE-FFT 论文修订

> 生成时间：2026-04-27  
> 整合依据：GPT 5.5 审稿意见 (`PDE_FFT_Thesis_Roadmap.md`) + 独立审稿报告 + Lean 4 校验 + 项目现状分析  
> 目标：将论文从"有系统性矛盾"修复为"理论-代码-实验-结论完全自洽"的毕业论文

---

## 零、根本问题诊断

**论文定义了 modified Helmholtz $(-\nabla^2+k^2)u=f$，但行文时按 true Helmholtz $(-\nabla^2-\kappa^2)u=f$ 的困难讨论。**

| 层面 | 现状 | 应改为 |
|------|------|--------|
| 理论 | 共振条件 $\lambda_{p,q}+k^2=0$ 对 modified H. 永不成立 | modified: $\lambda+k^2>0$ 恒成立；true: $\lambda-k^2=0$ 有共振 |
| 实验 | GMRES 迭代次数随 $k^2$ 增大而**减少** | modified H. 这是正确行为；true H. 才会增加 |
| 结论 | "GMRES 迭代次数随波数增大而显著增加" | modified: 迭代减少；true: 可能增加 |

---

## 一、8 阶段执行路线

```
阶段 0：备份与建立修订分支              [Day 0, 30分钟]
阶段 1：理论符号冻结                    [Day 1-2, 2天]
阶段 2：核心公式与论文理论重写           [Day 3-5, 3天]
阶段 3：代码 bug 修复与 σ 参数重构       [Day 5-7, 2.5天]
阶段 4：单元测试体系建设                [Day 8-9, 1.5天]
阶段 5：核心实验重构                    [Day 10-12, 3天]
阶段 6：图表与数据产物规范化             [Day 13, 1天]
阶段 7：Lean 4 附录整理 + 最终验收       [Day 14, 1天]
```

---

## 阶段 0：备份与建立修订分支

### 操作

```powershell
cd C:\Users\20564\Desktop\Graduate\论文收集
git init   # 如果还没有 git
git add .
git commit -m "backup: original thesis and solver before theory revision"
git checkout -b revise-sigma-fft9-krylov
```

### 验收

- [ ] `git log` 显示备份提交
- [ ] 当前分支为 `revise-sigma-fft9-krylov`
- [ ] 可以一键回退到修改前版本

---

## 阶段 1：理论符号冻结

### 1.1 目标

全文统一 PDE 形式为：

$$(-\nabla^2 + \sigma)u = f$$

其中：
- $\sigma = 0$：Poisson 方程
- $\sigma = +\kappa^2$：modified Helmholtz / screened Poisson
- $\sigma = -\kappa^2$：true Helmholtz

频域分母统一为：$d_{p,q} = \lambda_{p,q} + \sigma$

### 1.2 需要修改的 LaTeX 文件清单

| 文件 | 行数 | 修改内容 |
|------|------|----------|
| `main.tex` | L82-96 | 摘要中 $(-\nabla^2+k^2)$ 改为 $(-\nabla^2+\sigma)u=f$，解释 σ 三种取值 |
| `1_introduction.tex` | L10-11, L57-69 | 方程定义改用 σ；Helmholtz 求解部分区分 modified/true |
| `2_math_preliminary.tex` | L9-14 | 方程定义改为 $(-\nabla^2+\sigma)u=f$ |
| `2_math_preliminary.tex` | L82-101 | 矩阵表示改为 $B=\mathrm{tridiag}(-1,2,-1)$，$A=\frac{1}{h^2}(I\otimes B+B\otimes I)+\sigma I$ |
| `2_math_preliminary.tex` | L103-155 | 特征值改为 $\mu_p=2-2\cos\theta_p$，二维 $d_{p,q}=\frac{\mu_p+\mu_q}{h^2}+\sigma$，条件数改为 $4/(\pi^2 h^2)$ |
| `2_math_preliminary.tex` | L302-306 | 紧致格式改为 $(-L_h+\sigma R_h)u=R_h f$ |
| `3_fft_direct.tex` | L18-37 | FA 分母改为 $\lambda_{p,q}+\sigma$ |
| `3_fft_direct.tex` | L358 | FACR 复杂度改为 $O(N^2\log N)$ |
| `4_helmholtz.tex` | 全文(78行) | **大幅重写**：4.1 modified H.(SPD)，4.2 true H.(不定) |
| `5_gmres.tex` | L144-157 | λ_min 改为 $\approx 2\pi^2+k^2$（Dirichlet）或 $\approx k^2$（Neumann modified） |
| `6_experiments.tex` | L149-151 | "大波数困难"描述修正 |
| `7_conclusion.tex` | L60-63 | 结论第3条完全重写 |

### 1.3 搜索与替换命令

```powershell
# 搜索所有需要修改的位置
Select-String -Path thesis\chapters\*.tex -Pattern "Helmholtz" | Select-Object -Property LineNumber, Line
Select-String -Path thesis\chapters\*.tex -Pattern "共振" | Select-Object -Property LineNumber, Line
Select-String -Path thesis\chapters\*.tex -Pattern "大波数" | Select-Object -Property LineNumber, Line
Select-String -Path thesis\chapters\*.tex -Pattern "shifted" | Select-Object -Property LineNumber, Line
Select-String -Path thesis\chapters\*.tex -Pattern "k\^2" | Select-Object -Property LineNumber, Line
```

### 1.4 验收标准

- [ ] 全文使用 $(-\nabla^2+\sigma)u=f$ 统一形式
- [ ] modified H. 和 true H. 分开定义
- [ ] 共振只出现在 true H. 章节
- [ ] 摘要和结论不再把 $(-\nabla^2+k^2)$ 说成 indefinite Helmholtz
- [ ] 全文不存在 "$\lambda+k^2\approx 0$" 作为 modified H. 共振条件
- [ ] 结论章能解释：modified H. 中 $\kappa^2$ 增大通常改善条件数

---

## 阶段 2：核心公式与论文理论重写

### 2.1 五点差分矩阵 (消除 double counting)

**当前问题**：论文使用 $T=\mathrm{tridiag}(-1,4,-1)$，导致特征值 $\lambda_p=4-2\cos\theta_p$，二维公式 $(\lambda_p+\lambda_q)/h^2+k^2$ 看起来像 double counting。

**修正**：改用 $B=\mathrm{tridiag}(-1,2,-1)$：

$$A_h = \frac{1}{h^2}(I\otimes B + B\otimes I) + \sigma I$$

一维特征值：$\mu_p = 2-2\cos\frac{p\pi}{N+1}$

二维特征值：$d_{p,q} = \frac{\mu_p+\mu_q}{h^2}+\sigma$

**修改位置**：
- [ ] Ch2 §2.2.2 矩阵定义改用 $B$
- [ ] Ch2 §2.2.3 DST-I 特征值使用 $\mu_p$
- [ ] Ch3 §3.1.1 Kronecker 表达式改为 $I\otimes B+B\otimes I$
- [ ] 全文不再出现错误的"$(4-2\cos\theta_p)+(4-2\cos\theta_q)$"分母

### 2.2 条件数公式

**修正**：

$$\kappa(A_h) \sim \frac{4}{\pi^2 h^2}$$

（不是 $2/(\pi^2 h^2)$）

Modified H. 条件数：$\kappa = \frac{\lambda_{\max}+\kappa^2}{\lambda_{\min}+\kappa^2}$

- [ ] 条件数公式常数改为 $4/(\pi^2 h^2)$
- [ ] modified H. 条件数写为 $(\lambda_{\max}+\kappa^2)/(\lambda_{\min}+\kappa^2)$
- [ ] true H. 条件数不按 SPD 条件数解释

### 2.3 FFT9 符号与频域公式

**当前问题**：论文中 $L_h$ 的符号约定混乱，有时写 $L_h\approx -\nabla^2$，有时 $L_h\approx\Delta$。

**修正**（方案 B，与代码一致）：
- $L_h \approx \Delta$（不是 $-\Delta$）
- 紧致格式：$(-L_h+\sigma R_h)u = R_h f$
- 频域求解：$\hat{u}_{p,q} = \frac{\widehat{R_hf}_{p,q}}{-\hat{\lambda}_L+\sigma\hat{\lambda}_R}$

- [ ] Ch2 §2.3 改为 $L_h\approx\Delta$
- [ ] FFT9 正定形式写为 $(-L_h+\sigma R_h)u=R_h f$
- [ ] Algorithm FFT9 统一为方案 B
- [ ] 非齐次 Dirichlet FFT9 使用 $R_h f$ + boundary correction
- [ ] 不再出现 $g=-R_hf$ 却除以 $-\lambda_L/\lambda_R+\sigma$ 的混搭写法

### 2.4 Neumann DCT-I 与兼容条件

- [ ] DCT-I 公式补充 $\alpha_j\alpha_k$ 双端点权重
- [ ] Neumann Poisson 写明兼容条件：$\sum w_iw_jf_{i,j}=0$
- [ ] 求解流程写明 $\hat{u}_{0,0}=0$
- [ ] 解满足 weighted mean zero

### 2.5 FACR 复杂度

**决策：降级为 $O(N^2\log N)$**

- [ ] 当前实现复杂度改为 $O(N^2\log N)$
- [ ] 表格中 FACR-like 不再标 $O(N^2\log\log N)$
- [ ] classical FACR 作为理论说明或展望
- [ ] 实验图中不再声称已验证 $O(N^2\log\log N)$

### 2.6 六阶不可行证明

**修正**：保留为"受限三参数族不可行"，不泛化为所有六阶紧致格式不可行。

- [ ] 章节标题改为"三参数九点修正族的六阶不可行性"
- [ ] 删除"c₂ 无法消去"的错误论证
- [ ] 改为"先消 c₂（令 $\alpha=-\gamma$），再由 c₄ 模态依赖性得矛盾"
- [ ] 不再写"六阶格式普遍不可行"

### 2.7 Ch4 大幅重写 (78行 → 250+行)

当前 Ch4 仅 78 行，是论文最薄弱的章节。需要：

| 新增内容 | 预计行数 |
|----------|----------|
| σ 参数统一框架定义 | 40 |
| 4.1 modified H.：SPD 性质、条件数分析、无 Dirichlet 共振 | 50 |
| 4.2 true H.：不定矩阵、共振条件 $\lambda=\kappa^2$、近共振分析 | 50 |
| 频域除数统一公式 $D=\lambda_{p,q}+\sigma$ | 20 |
| Lean 4 验证结论引用 | 20 |
| σ 参数的代码实现细节 | 30 |
| 数值例子 | 40 |

- [ ] Ch4 扩充到 250+ 行
- [ ] 共振检测只用于 $\sigma=-\kappa^2$（true H.）
- [ ] 引用 Lean 4 `modified_no_resonance` 和 `true_helmholtz_resonance_iff`

---

## 阶段 3：代码修复与 σ 参数重构

### 3.1 修复非齐次 Dirichlet 符号 bug（最高优先级）

**搜索**：
```powershell
Select-String -Path .\code\python\*.py -Pattern "F.*-=" | Select-Object -Property LineNumber, Line
Select-String -Path .\code\python\*.py -Pattern "bc_l|bc_r|bc_b|bc_t" | Select-Object -Property LineNumber, Line
```

**修改**：所有 `F -= bc/h**2` 改为 `F += bc/h**2`

已确认位置（11处）：

| 文件 | 行号 | 修改 |
|------|------|------|
| `helmholtz_solver.py` | L269 | `F[0, :] -= bc_l_full[1:-1]/h**2` → `+=` |
| `helmholtz_solver.py` | L270 | `F[-1, :] -= bc_r_full[1:-1]/h**2` → `+=` |
| `helmholtz_solver.py` | L277 | `F[:, 0] -= bc_b_full[1:-1]/h**2` → `+=` |
| `helmholtz_solver.py` | L278 | `F[:, -1] -= bc_t_full[1:-1]/h**2` → `+=` |
| `helmholtz_solver.py` | L319 | `F[0, :] -= bc_l_full/h**2` → `+=` |
| `helmholtz_solver.py` | L320 | `F[-1, :] -= bc_r_full/h**2` → `+=` |
| `gmres_solver.py` | L391 | `F[0, :] -= bc_l[1:-1]/h2` → `+=` |
| `gmres_solver.py` | L392 | `F[-1, :] -= bc_r[1:-1]/h2` → `+=` |
| `gmres_solver.py` | L393 | `F[:, 0] -= bc_b[1:-1]/h2` → `+=` |
| `gmres_solver.py` | L394 | `F[:, -1] -= bc_t[1:-1]/h2` → `+=` |
| `gmres_solver.py` | L418-419 | Neumann BC 的 2 处同样修改 |

- [ ] 上述 11 处全部修改
- [ ] 代码中加注释说明为什么是加号

### 3.2 增加 PDEParams 类

新建 `code/python/pde_params.py`：

```python
from dataclasses import dataclass
from typing import Literal

@dataclass(frozen=True)
class PDEParams:
    kind: Literal["poisson", "modified_helmholtz", "true_helmholtz"]
    kappa: float = 0.0

    @property
    def sigma(self) -> float:
        if self.kind == "poisson":
            return 0.0
        if self.kind == "modified_helmholtz":
            return self.kappa ** 2
        if self.kind == "true_helmholtz":
            return -(self.kappa ** 2)
        raise ValueError(f"Unknown PDE kind: {self.kind}")
```

- [ ] 所有 solver 不再直接使用 `k2`，改用 `params.sigma`
- [ ] 所有实验 CSV 同时记录 `kind, kappa, sigma`

### 3.3 五点 FFT 分母修正

```python
theta = np.arange(1, N+1) * np.pi / (N+1)
mu = 2 - 2 * np.cos(theta)
denom = (mu[:, None] + mu[None, :]) / h**2 + sigma
```

- [ ] Dirichlet 5pt FFT 使用 `mu=2-2cos(theta)`
- [ ] `denom = (mu_x + mu_y)/h^2 + sigma`
- [ ] true H. 时 near denominator warning 可触发

### 3.4 FFT9 实现方案 B

```python
def lambda_L_9pt(alpha, beta, h):
    return (-20 + 8*np.cos(alpha) + 8*np.cos(beta) + 4*np.cos(alpha)*np.cos(beta)) / (6*h*h)

def lambda_R(alpha, beta):
    return 2/3 + (np.cos(alpha) + np.cos(beta)) / 6

denom = -lambda_L + sigma * lambda_R
u_hat = Rh_f_hat / denom
```

- [ ] FFT9 内部注释写明 $L_h\approx\Delta$
- [ ] FFT9 使用 `denom = -lambda_L + sigma*lambda_R`
- [ ] FFT9 与 sparse compact matrix 解一致

### 3.5 Neumann compatibility

```python
def trapezoid_weights(n):
    w = np.ones(n)
    w[0] = 0.5
    w[-1] = 0.5
    return w

def weighted_mean_neumann(F):
    nx, ny = F.shape
    wx = trapezoid_weights(nx)
    wy = trapezoid_weights(ny)
    W = wx[:, None] * wy[None, :]
    return np.sum(W * F) / np.sum(W)

def project_neumann_compatible(F):
    return F - weighted_mean_neumann(F)
```

- [ ] Neumann Poisson 检查 `weighted_mean(f)`
- [ ] 不兼容 RHS 报错或自动投影
- [ ] $\hat{u}_{0,0}=0$
- [ ] 解满足 `weighted_mean(u) < 1e-12`

---

## 阶段 4：单元测试体系建设

### 4.1 测试文件结构

```
code/python/tests/
  test_01_eigenvalues.py       # Bv = μv，‖Bv-μv‖∞ < 1e-12
  test_02_nonhom_dirichlet.py  # u=exp(x+y)，验证 F+=bc/h²
  test_03_fft_vs_spsolve.py    # FFT vs scipy.sparse.linalg.spsolve ≤ 1e-10
  test_04_fft9_vs_spsolve.py   # FFT9 vs compact sparse ≤ 1e-10
  test_05_neumann_compatibility.py  # weighted mean < 1e-12
  test_06_modified_true_helmholtz.py # modified 分母恒正；true 可接近零
  test_07_krylov_baselines.py  # SPD 用 CG；true H. 不用 CG
```

### 4.2 必须通过的测试指标

| 测试 | 验收 |
|------|------|
| 特征值 | $\|Bv-\mu v\|_\infty<10^{-12}$ |
| 非齐次 Dirichlet | 收敛阶 ≈ 2 |
| FFT vs spsolve | $\le 10^{-10}$ |
| FFT9 vs spsolve | $\le 10^{-10}$ |
| Neumann 兼容 | weighted mean $<10^{-12}$ |
| modified 分母 | $\lambda+\kappa^2>0$ 恒正 |
| true 共振 | $\min|\lambda-\kappa^2|$ 可接近 0 |
| Krylov 基线 | SPD 用 CG；true H. 不用 CG |

- [ ] `pytest -q` 全部通过
- [ ] 重点测试 `test_02_nonhom_dirichlet.py`

---

## 阶段 5：核心实验重构

### 5.1 实验目录

```
code/python/experiments/
  exp00_verify_fft_vs_sparse.py
  exp01_convergence_order.py
  exp02_nonhom_dirichlet.py
  exp03_neumann_mixed_bc.py
  exp04_modified_true_helmholtz.py
  exp05_true_helmholtz_resonance.py
  exp06_krylov_fair_comparison.py
  exp07_complexity_scaling.py
```

### 5.2 关键实验设计

#### exp00：FFT vs sparse direct（离散一致性）

- **目标**：证明 FFT solver 解的是同一个离散系统
- **对象**：5pt Dirichlet, 5pt Neumann modified H., FFT9 Dirichlet, mixed DST/DCT
- **验收**：$\|u_{\text{fft}}-u_{\text{spsolve}}\|_\infty \le 10^{-10}$

#### exp02：非齐次 Dirichlet（验证 A1 Bug 修复）

- **测试函数**：$u=e^{x+y}$（边界值不为零）
- **验收**：5pt 二阶收敛；FFT9 与 sparse compact 一致

#### exp04：modified vs true Helmholtz

- **modified**：$(-\nabla^2+\kappa^2)u=f$ → 分母恒正，条件数改善，CG/GMRES 迭代下降
- **true**：$(-\nabla^2-\kappa^2)u=f$ → 分母可能为零，矩阵不定，MINRES/GMRES 适用

#### exp05：true Helmholtz 近共振

- **固定模态**：$p=2, q=3$
- **扫描**：$\kappa^2 = \lambda_{p,q}^h + \delta$，$\delta \in \{-10^{-1}, -10^{-2}, -10^{-3}, 10^{-3}, 10^{-2}, 10^{-1}\}$
- **RHS**：使用固定 RHS $f=\sin(p\pi x)\sin(q\pi y)$，不用共振模态作为 manufactured solution
- **验收**：min_abs_denominator 随 $\delta\to 0$ 而下降；解范数/迭代数有可观察变化

#### exp06：Krylov 公平对比

- **modified H. SPD**：CG, GMRES, Jacobi-CG, ILU-GMRES
- **true H. 不定**：MINRES, GMRES, shifted-Laplacian GMRES
- **统一条件**：$x_0=0$，$\text{tol}=10^{-10}$，same A, same b

- [ ] 所有实验输出 CSV
- [ ] 所有图从 CSV 生成
- [ ] 5pt 收敛阶在 [1.8, 2.2]
- [ ] FFT9 收敛阶在 [3.7, 4.3]

---

## 阶段 6：图表与数据产物规范化

### 6.1 禁止项

- [ ] 空坐标轴
- [ ] 没有数据点的图
- [ ] 只有一行或一列色块的 heatmap
- [ ] 没有参考斜率的 convergence 图
- [ ] 手工填表后再画图

### 6.2 必须项

| 图类型 | 必须包含 |
|--------|----------|
| convergence 图 | log-log, 数据点, 参考 $O(h^2)$/$O(h^4)$, legend, grid |
| resonance 图 | min_denominator vs $\kappa^2$, solution_norm vs $\kappa^2$, krylov_iterations vs $\kappa^2$ |
| heatmap | 至少 5×5 参数网格, 色条 $\log_{10}$ |
| timing 图 | setup/solve/total, 参考 $O(M\log M)$ |

---

## 阶段 7：Lean 4 附录整理 + 最终验收

### 7.1 Lean 4 论文写法

**推荐写法**：

> 为降低离散谱推导中的符号错误风险，本文使用 Lean 4 + mathlib 对若干核心代数恒等式进行了辅助形式化校验。校验内容包括：一维五点差分 sine 模态特征值 $2-2\cos\theta$、二维 Kronecker 和特征值叠加、九点紧致格式 $L_h$ 的 Fourier symbol、$R_h$ 的 Fourier symbol、modified Helmholtz 与 true Helmholtz 的频域分母差异，以及 Neumann ghost-point 矩阵边界项经 $\sqrt{2}$ 缩放后的对称化。Lean 项目在无 `sorry`、`admit`、`axiom`、`unsafe` 的条件下通过 `lake build`。

**禁止写法**：~~"本文所有理论推导均已由 Lean 4 严格证明"~~

### 7.2 最终验收清单

#### 理论验收

- [ ] 全文统一 PDE：$(-\nabla^2+\sigma)u=f$
- [ ] $\sigma=0/+\kappa^2/-\kappa^2$ 定义清楚
- [ ] modified H. 无 Dirichlet 共振
- [ ] true H. 共振条件为 $\lambda=\kappa^2$
- [ ] 五点谱公式无 double counting
- [ ] FFT9 的 $L_h\approx\Delta$
- [ ] FFT9 RHS 与分母一致
- [ ] 非齐次 Dirichlet RHS 为加号
- [ ] Neumann DCT-I 端点权重完整
- [ ] Neumann Poisson 有兼容条件和零均值
- [ ] FACR 当前实现不声称 $O(N^2\log\log N)$
- [ ] CG/MINRES/GMRES 适用性区分清楚

#### 代码验收

- [ ] `pytest -q` 全部通过
- [ ] 非齐次 Dirichlet 测试通过
- [ ] 5pt FFT vs spsolve $\le 10^{-10}$
- [ ] FFT9 vs spsolve $\le 10^{-10}$
- [ ] Neumann compatibility 测试通过
- [ ] modified H. 分母恒正
- [ ] true H. near-resonance warning 可触发

#### 实验验收

- [ ] 5pt 收敛阶在 [1.8, 2.2]
- [ ] FFT9 收敛阶在 [3.7, 4.3]
- [ ] modified H. 下 $\kappa^2$ 增大迭代通常下降
- [ ] true H. near resonance 下 min denominator 下降
- [ ] timing 区分 setup/solve/total
- [ ] 图全部从 CSV 生成
- [ ] 无空图

#### Lean 4 验收

- [ ] `lake build` 通过
- [ ] `grep` 不存在 sorry/admit/axiom/unsafe
- [ ] 论文只声明"关键代数公式辅助校验"
- [ ] 不声明"完整 PDE 理论已形式化证明"

#### LaTeX 编译

- [ ] 0 错误
- [ ] 公式编号连贯
- [ ] 页数约 55-60 页

---

## 二、两周日历排期

| 日期 | 阶段 | 具体任务 | 交付物 |
|------|------|----------|--------|
| Day 0 | 0 | git 备份 + 建分支 | 备份提交 |
| Day 1 | 1 | Ch1+Ch2 方程定义改 σ | 2_tex_diff |
| Day 2 | 1 | Ch3+Ch4+Ch5+Ch7 符号统一 | 全文 σ 冻结 |
| Day 3 | 2 | Ch2 矩阵(B)+特征值(μ)+条件数(4/π²) | 2_tex_math |
| Day 4 | 2 | Ch2 FFT9(Lh≈Δ)+Ch3 频域公式 | 2_tex_fft9 |
| Day 5 | 2 | Ch4 大幅重写(78→250行)+Ch7 结论 | 2_tex_helmholtz |
| Day 5-6 | 3 | A1 Bug 修复(11处)+PDEParams | 代码修复 |
| Day 6-7 | 3 | σ 参数重构+FFT9 方案 B+Neumann | 代码重构 |
| Day 8-9 | 4 | 7 个单元测试编写+运行 | tests/ |
| Day 10-12 | 5 | 8 个实验重跑+CSV 输出 | results/ |
| Day 13 | 6 | 图表规范化+从 CSV 生成 | figures/ |
| Day 14 | 7 | Lean 4 附录+全文校对+LaTeX 编译 | 终稿 |

---

## 三、各章节行数预期

| 章节 | 当前行数 | 目标行数 | 变化 |
|------|----------|----------|------|
| Ch1 引言 | 103 | ~130 | +27 (σ 框架介绍) |
| Ch2 数学基础 | 601 | ~680 | +79 (B矩阵+条件数+FFT9符号修正) |
| Ch3 FFT直接法 | 369 | ~400 | +31 (FACR降级+FFT9方案B) |
| Ch4 Helmholtz | 78 | ~260 | +182 (大幅扩充) |
| Ch5 GMRES | 229 | ~280 | +51 (λ_min修正+CG基线+true H.) |
| Ch6 实验 | 626 | ~650 | +24 (新增非齐次+true H.实验) |
| Ch7 结论 | 108 | ~120 | +12 (结论修正) |
| **总计** | **2114** | **~2520** | **+406** |

---

## 四、关键决策（已确定）

| 决策 | 选择 | 理由 |
|------|------|------|
| FACR 复杂度 | 方案 B：改声明为 $O(N^2\log N)$ | 更务实，当前实现确实是 DST-CR hybrid |
| σ 参数 | 代码中也实现 true H. ($\sigma=-\kappa^2$) | 改动小，支撑 true H. 实验 |
| Ch4 扩充 | 250 行 | 毕业论文需要足够的理论深度 |
| 六阶证明 | "三参数族不可行"，不泛化 | 精确对应证明范围 |
| L_h 符号 | $L_h\approx\Delta$ | 与代码一致，与 Lean 4 一致 |
| Lean 4 定位 | 附录亮点，"辅助形式化校验" | 诚实准确，不夸大 |

---

## 五、完成定义

本项目完成的标准不是"图很多"，而是以下 6 点同时满足：

1. **理论符号统一**：$(-\nabla^2+\sigma)u=f$
2. **代码离散系统可信**：FFT vs spsolve 通过
3. **收敛阶可信**：5pt 二阶，FFT9 四阶
4. **Helmholtz 分类可信**：modified 和 true 分开
5. **Krylov 对比公平**：modified 有 CG，true 有 MINRES/GMRES
6. **Lean 4 附录克制准确**：只证明关键代数公式，不夸大

最终论文一句话定位：

> 本文在统一的 $(-\nabla^2+\sigma)u=f$ 框架下，研究矩形区域上 Poisson、modified Helmholtz 与 true Helmholtz 方程的 FFT 快速直接求解器和 Krylov 迭代法；通过 sparse direct 一致性、收敛阶、边界条件、近共振和复杂度实验，验证不同算法在规则区域上的精度、效率和适用边界。
