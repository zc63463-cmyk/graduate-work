---
title: "GMRES: 求解非对称线性系统的广义最小残差算法"
authors: "Youcef Saad, Martin H. Schultz"
year: 1986
journal: "SIAM Journal on Scientific and Statistical Computing"
tags:
  - 数值线性代数
  - Krylov子空间方法
  - 迭代法
  - 非对称矩阵
  - 大型稀疏系统
  - 算法
status: 精读完成
date: 2026-04-25
---

# GMRES: 求解非对称线性系统的广义最小残差算法

> [!info] 论文基本信息
> - **标题**: GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems
> - **作者**: Youcef Saad¹, Martin H. Schultz¹
> - **单位**: ¹耶鲁大学计算机科学系
> - **发表**: SIAM Journal on Scientific and Statistical Computing, Vol.7, No.3, 1986
> - **页数**: 14页
> - **被引**: >10,000次 (Google Scholar)
> - **PDF**: [[GMRES_translated.pdf]]

## 一、核心贡献

### 1.1 要解决什么问题？

求解大型稀疏**非对称**线性方程组：

$$Ax = b$$

其中 $A \in \mathbb{R}^{N \times N}$ 是非对称矩阵（可能不定）。

**现有方法的缺陷**：
- **GCR (广义共轭残差法)**: 当矩阵A不是正实数（对称部分不是正定）时，可能**崩溃**（除零错误）
- **Orthodir/Orthomin**: 数值稳定性差，误差累积
- **MINRES**: 仅适用于对称不定矩阵

### 1.2 GMRES的创新点

1. **永不崩溃**（除非已收敛到精确解）
2. **最优性保证**: 在Krylov子空间上最小化残差L2范数
3. **存储效率高**: 比GCR节省约一半存储空间
4. **计算效率高**: 比GCR少约1/3算术运算

### 1.3 理论基础

GMRES基于**Arnoldi过程**（非对称版本的Lanczos过程）：

$$K_k = \text{span}\{v_1, Av_1, A^2v_1, ..., A^{k-1}v_1\}$$

通过Arnoldi过程构建Krylov子空间的正交基 $\{v_1, v_2, ..., v_k\}$。

---

## 二、算法原理

### 2.1 Arnoldi过程（算法1）

**目标**: 构建Krylov子空间 $K_k$ 的一组正交基

**输入**: 矩阵 $A$, 初始向量 $v_1$ (满足 $\|v_1\|_2 = 1$)

**输出**: 正交基 $\{v_1, ..., v_k\}$ 和上Hessenberg矩阵 $H_k$

```
算法1: Arnoldi方法

1. 开始: 选择初始向量v₁，满足 ||v₁||₂ = 1

2. 迭代: 对于 j = 1, 2, ..., k
   wⱼ = Avⱼ
   hᵢⱼ = (wⱼ, vᵢ),  i = 1, 2, ..., j
   wⱼ = wⱼ - Σᵢ₌₁ʲ hᵢⱼvᵢ
   hⱼ₊₁ⱼ = ||wⱼ||₂
   vⱼ₊₁ = wⱼ / hⱼ₊₁ⱼ
```

**关键关系**（公式3）:

$$AV_k = V_{k+1} H_k$$

其中：
- $V_k = [v_1, v_2, ..., v_k]$ 是 $N \times k$ 矩阵
- $H_k$ 是 $(k+1) \times k$ 上Hessenberg矩阵
- $H_k$ 的元素由Arnoldi过程生成的 $h_{i,j}$ 构成

### 2.2 GMRES算法（算法3）

**目标**: 在Krylov子空间 $K_k$ 中找到近似解 $x_k$，使得残差范数 $\|b - Ax_k\|_2$ 最小

**数学形式**:

$$x_k = x_0 + V_k y_k$$

其中 $y_k$ 最小化：

$$J(y) = \|\beta e_1 - H_k y\|_2$$

这里 $\beta = \|r_0\|_2$, $e_1 = (1, 0, ..., 0)^T \in \mathbb{R}^{k+1}$

```
算法3: GMRES

1. 开始: 
   选择x₀，计算 r₀ = b - Ax₀
   v₁ = r₀ / ||r₀||₂
   β = ||r₀||₂

2. Arnoldi迭代: 对于 j = 1, 2, ..., k
   （执行算法1的第2步）

3. 形成近似解:
   x_k = x₀ + V_k y_k
   其中y_k通过最小化 ||βe₁ - H_k y||₂ 获得
```

### 2.3 最小二乘问题的求解

**使用Givens旋转进行QR分解**:

$$Q_k H_k = R_k$$

其中 $Q_k$ 是旋转矩阵的乘积。

**重要发现**（命题1）:

> 近似解 $x_k$ 的残差范数等于变换后右侧向量 $g_k = Q_k \beta e_1$ 的第 $(k+1)$ 个分量。
> 
> **这意味着无需显式计算 $x_k$，就能获得残差范数！**

### 2.4 重启GMRES（算法4）

**问题**: 当 $k$ 增大时，存储和运算成本呈**平方增长**。

**解决方案**: 每 $m$ 步重启一次

```
算法4: 重启GMRES(m)

1. 开始: 选择x₀，计算 r₀ = b - Ax₀，v₁ = r₀ / ||r₀||₂

2. Arnoldi迭代: 对于 j = 1, 2, ..., m
   （执行算法1的第2步）

3. 形成近似解:
   x_m = x₀ + V_m y_m
   其中y_m最小化 ||βe₁ - H_m y||₂

4. 重启:
   计算 r_m = b - Ax_m
   if ||r_m|| 足够小: 停止
   else: 
      x₀ := x_m
      v₁ := r_m / ||r_m||₂
      转到步骤2
```

---

## 三、理论基础

### 3.1 不会崩溃的定理（命题2）

GMRES在第 $j$ 步产生精确解 **当且仅当** 以下四个等价条件成立：

1. 算法在第 $j$ 步失效（崩溃）
2. $v_{j+1} = 0$
3. $h_{j+1,j} = 0$
4. 初始残差向量 $r_0$ 的极小多项式的次数等于 $j$

**推论3**: 对于 $N \times N$ 问题，GMRES最多 $N$ 步收敛。

### 3.2 收敛性分析（定理5）

**假设**: 矩阵 $A$ 有 $v$ 个特征值 $\lambda_1, \lambda_2, ..., \lambda_v$ 的实部非负，其余特征值位于以 $C$ 为中心、半径为 $R$ 的圆内（$C > 0$）。

**上界**:

$$\|r_m\| \leq \left(\frac{R}{C}\right)^m \cdot \kappa(X) \cdot \|r_0\|$$

其中 $\kappa(X) = \|X\|_2 \|X^{-1}\|_2$ 是特征向量矩阵的条件数。

**关键洞察**:
- 左半平面特征值少 → 收敛快
- $m$ 足够大时保证收敛
- 实际收敛所需 $m$ 可能远小于理论预测

### 3.3 与GCR的等价性

GMRES在数学上等价于**广义共轭残差法（GCR）**，但：
- GMRES数值更稳定
- GMRES存储需求仅为GCR的一半
- GMRES每步运算量更少

---

## 四、数值实验

### 4.1 测试问题

**PDE**: 二维椭圆偏微分方程

$$-\frac{\partial}{\partial x}(b\frac{\partial u}{\partial x}) - \frac{\partial}{\partial y}(c\frac{\partial u}{\partial y}) + d\frac{\partial u}{\partial x} + e\frac{\partial u}{\partial y} + fu = g$$

在单位正方形上离散化，使用五点差分格式。

**参数**:
- $n$: 每边内部节点数
- $N = n^2$: 矩阵维度
- $\gamma, \beta$: 调节矩阵对称性程度的参数

**预处理**: MILU (Modified Incomplete LU)

### 4.2 实验结果

**实验1** ($n=48, N=2304, \gamma=50, \beta=1$):

| 方法 | 收敛? | 乘法运算量 |
|------|--------|-----------|
| Orthomin(1) | ❌ 不收敛 | - |
| Orthomin(5) | ❌ 不收敛 | - |
| GCR(1) | ✅ | 中等 |
| **GMRES(5)** | ✅ | **与GCR(1)相当** |

**实验2** ($n=18, N=324, \gamma=50, \beta=-20$):

| 方法 | 收敛? |
|------|--------|
| Orthomin(k), k=1~10 | ❌ 全部发散 |
| GCR(1), GCR(2), GCR(3) | ❌ 发散 |
| GMRES(2), GMRES(3), GMRES(4) | ❌ 发散 |
| **GMRES(5)** | ✅ 开始收敛 |
| **GMRES(20)** | ✅ 性能最佳 |

**关键发现**:
- GMRES(k)在 $k \geq 5$ 时开始收敛，且随 $k$ 增大显著改善
- 对于困难问题，需要更大的重启参数 $m$

---

## 五、实现细节

### 5.1 高效实现技巧

1. **使用改进的Gram-Schmidt**
   - 比经典Gram-Schmidt数值更稳定
   - 参考文献[15]: Stewart (1973)

2. **逐步更新QR分解**
   - 在Arnoldi过程的每一步更新 $H_k$ 的QR分解
   - 无需显式计算 $x_k$ 就能获得残差范数

3. **避免计算 $v_{m+1}$**
   - 重启时，通过 $v_i$ 和 $Av_m$ 的线性组合计算残差
   - 节省 $(2m+1)N$ 次乘法运算

### 5.2 存储需求对比

| 方法 | 存储需求 |
|------|----------|
| GCR(m-1) | $(2m+1)N$ |
| **GMRES(m)** | **$(m+2)N$** |

**示例** ($m=20, N=1000$):
- GCR: 需要存储 41,000 个向量元素
- GMRES: 需要存储 22,000 个向量元素（节省约46%）

### 5.3 运算量对比

每步平均运算量（忽略 $NZ$ 项）：

| 方法 | 乘法运算量 |
|------|-----------|
| GCR(m-1) | $(3m+5)N$ |
| **GMRES(m)** | **$(m+3+1/m)N$** |

**示例** ($m=20$):
- GCR: $65N$ 次乘法/步
- GMRES: $23.05N$ 次乘法/步（节省约64%）

---

## 六、优缺点总结

### 6.1 优点

- ✅ **理论完备**: 不会崩溃（除非已收敛）
- ✅ **最优性保证**: 残差范数在Krylov子空间上最小
- ✅ **存储高效**: 比GCR节省约一半存储
- ✅ **计算高效**: 比GCR少约1/3算术运算
- ✅ **残差范数免费**: 无需显式计算解就能监控收敛
- ✅ **广泛适用**: 适用于各种非对称问题

### 6.2 缺点

- ❌ **重启可能停滞**: 对于某些不定问题，重启GMRES可能不收敛到零
- ❌ **m选择困难**: 重启参数 $m$ 需要经验或试错
- ❌ **理论界保守**: 收敛性上界可能过于悲观
- ❌ **大规模问题**: 对于超大问题($N > 10^6$)，仍需结合预处理

---

## 七、实践建议

### 7.1 参数选择

1. **重启参数 $m$**:
   - 简单问题: $m = 10 \sim 20$
   - 中等难度: $m = 20 \sim 50$
   - 困难问题: $m = 50 \sim 100$ 或更大

2. **停止准则**:
   - 基于残差范数: $\|r_k\|_2 / \|r_0\|_2 < \epsilon$
   - 利用命题1免费获得残差范数

3. **预处理**:
   - ILU (Incomplete LU): 简单有效
   - MILU (Modified ILU): 对某些问题更好
   - 多网格预处理: 对于PDE离散化系统非常有效

### 7.2 常见陷阱

1. ⚠️ **忘记归一化 $v_1$**
   ```python
   v1 = r0 / np.linalg.norm(r0)
   ```

2. ⚠️ **Arnoldi过程中 $h_{j+1,j} = 0$**
   - 这意味着算法已完成（解是精确的）
   - 不要继续计算 $v_{j+1}$

3. ⚠️ **重启时丢弃有用信息**
   - 考虑使用**灵活GMRES (FGMRES)**
   - 或使用** deflation **技术

4. ⚠️ **m太小导致停滞**
   - 监控残差范数
   - 如果残差下降太慢，增大 $m$

### 7.3 调试技巧

1. **检查正交性**:
   ```python
   # V_k^T V_k 应该接近单位矩阵
   orth_error = np.linalg.norm(V.T @ V - np.eye(k))
   ```

2. **验证Arnoldi关系**:
   ```python
   # ||AV_k - V_{k+1} H_k||_F 应该很小
   arnoldi_error = np.linalg.norm(A @ V - V_new @ H)
   ```

3. **监控残差历史**:
   ```python
   plt.semilogy(residuals)
   plt.xlabel('Iteration')
   plt.ylabel('Residual norm')
   ```

---

## 八、扩展阅读

### 8.1 Saad的其他工作

- Saad (1981): "Krylov subspace methods for solving large unsymmetric linear systems"
- Saad & Schultz (1983): "Conjugate gradient-like algorithms for solving nonsymmetric linear systems"
- Saad (1993): 《Iterative Methods for Sparse Linear Systems》教科书（第6章详细讲解GMRES）

### 8.2 现代发展

1. **Flexible GMRES (FGMRES)**:
   - 允许预处理矩阵在每次迭代中变化
   - 适用于非线性预处理

2. **Deflated GMRES**:
   - 移除Krylov子空间中的"坏方向"
   - 加速收敛

3. **Augmented GMRES**:
   - 在Krylov子空间中加入特定方向
   - 针对已知特征值或特征向量

### 8.3 相关算法对比

| 算法 | 适用矩阵 | 理论保证 | 存储 | 稳定性 |
|------|----------|----------|------|--------|
| **GMRES** | 非对称 | ✅ 残差非增 | 中等 | ✅ 好 |
| BiCGStab | 非对称 | ❌ 可能震荡 | 低 | ✅ 较好 |
| CG | 对称正定 | ✅ 最优 | 低 | ✅ 很好 |
| MINRES | 对称不定 | ✅ 残差非增 | 低 | ✅ 很好 |
| IDR(s) | 非对称 | ❌ 较弱 | 低 | ✅ 较好 |

---

## 九、代码实现

### 9.1 Python实现

参见 [[GMRES_Python实现]]

### 9.2 C++实现

参见 [[GMRES_C++实现]]

---

## 十、总结

> [!quote] 论文摘要（翻译）
> 我们提出一种求解线性方程组的迭代算法，该算法具有在每一步骤中最小化残差向量在Krylov子空间上L2范数的特性。该算法源自构建Krylov子空间正交基的Arnoldi过程。可视为Paige和Saunders MINRES算法的推广，其理论等价于广义共轭残差法（GCR）和orthodir方法。相较于GCR和orthodir，新算法展现出多项优势。

**核心要点**:
1. GMRES = Arnoldi过程 + 最小残差原理
2. 永不崩溃（除非已收敛到精确解）
3. 存储和计算效率高于GCR
4. 实际使用必须重启（GMRES(m)）
5. $m$ 的选择是艺术与科学的结合

**影响**:
- 数值线性代数领域的**里程碑式工作**
- 成为MATLAB `gmres` 函数、PETSc、Hypre等库的核心算法
- 后续发展出FGMRES、deflated GMRES等多种变体

---

## 参考文献

1. Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. *SIAM Journal on Scientific and Statistical Computing*, 7(3), 856-869.

2. Saad, Y. (1993). *Iterative Methods for Sparse Linear Systems*. PWS Publishing Company.

3. Greenbaum, A. (1997). *Iterative Methods for Solving Linear Systems*. SIAM.

4. Kelley, C. T. (1995). *Iterative Methods for Linear and Nonlinear Equations*. SIAM.

---

**最后更新**: 2026-04-25
