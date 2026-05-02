# Cyclic Reduction & FACR 算法精读笔记

> 论文: *A Fast Poisson Solver — Cyclic Reduction + Fourier Analysis (FACR)*  
> Swarztrauber, SIAM Review, 1977

---

## 1. 问题背景

求解泊松方程 $-\nabla^2 u = f$ 在 $[0,1]^2$ 上，Dirichlet 边界条件。  
五点差分离散后，得到块三对角系统：

$$
A \mathbf{u} = \mathbf{f}, \quad A = I \otimes T + T \otimes I
$$

其中 $B = \text{tridiag}(-1, 2, -1)$，系统规模 $N^2 \times N^2$。
（注：若保留块三对角视角，对角块可记为 $T=\mathrm{tridiag}(-1,4,-1)$，
但此时一维特征值为 $\mu_p=2-2\cos\frac{p\pi}{N+1}$，而非 $4-2\cos\theta_p$。）

---

## 2. 三种求解方法

### 2.1 Fourier Analysis (FA) — O(N² log N)

1. 对每行做 1D-DST (type I)，变换到 Fourier 空间
2. 每个模式 $(p,q)$ 的特征值 $\lambda_p + \lambda_q$，直接除以特征值得解
3. 逆 DST 还原

```python
# 2D-DST → 除以特征值 → 逆 2D-DST
Fh = dst2d(F_adj)
Uh = Fh / (lambda_x[:,None] + lambda_y[None,:])
U = idst2d(Uh)
```

### 2.2 Cyclic Reduction (CR) — O(N² log N)，与 FA 等价

1. 仅在 x 方向做 DST → 解耦为 N 个标量三对角系统
2. 每个三对角系统用 Thomas 算法求解

**CR 的核心思想**：不是直接 Thomas，而是通过消去奇数行递归简化系统。

### 2.3 FACR(l) — O(N² log log N)，最优

混合策略：
1. x 方向 DST 解耦
2. l 步 CR 消去奇数行，系统缩小到 $N/2^l$
3. 解小型简化系统（Thomas）
4. 逐步回代恢复所有值

最优 $l \approx \log_2(\log_2(N))$，复杂度 $O(N^2 \log\log N)$。

---

## 3. Cyclic Reduction 详解

### 3.1 单步 CR 消元

原始系统（对称三对角）：
$$
e_{j-1} u_{j-1} + d_j u_j + e_j u_{j+1} = b_j, \quad j = 0, 1, \ldots, M-1
$$

消去奇数下标 $u_{2k+1}$：从行 $2k+1$ 解出
$$
u_{2k+1} = \frac{b_{2k+1} - e_{2k} u_{2k} - e_{2k+1} u_{2k+2}}{d_{2k+1}}
$$

代入偶数行 $2k$ 的方程，得到简化系统：

| 位置 | 新对角线 $d'_k$ | 新次对角线 $e'_k$ | 新右端 $b'_k$ |
|------|-----------------|-------------------|---------------|
| 内部 | $d_{2k} - \frac{e_{2k-1}^2}{d_{2k-1}} - \frac{e_{2k}^2}{d_{2k+1}}$ | $-\frac{e_{2k} e_{2k+1}}{d_{2k+1}}$ | $b_{2k} - \frac{e_{2k-1}}{d_{2k-1}} b_{2k-1} - \frac{e_{2k}}{d_{2k+1}} b_{2k+1}$ |
| 边界 (k=0) | $d_0 - \frac{e_0^2}{d_1}$ | $-\frac{e_0 e_1}{d_1}$ | $b_0 - \frac{e_0}{d_1} b_1$ |
| 边界 (k=M'-1) | $d_{M-1} - \frac{e_{M-2}^2}{d_{M-2}}$ | $-\frac{e_{M-3}e_{M-2}}{d_{M-2}}$ | $b_{M-1} - \frac{e_{M-2}}{d_{M-2}} b_{M-2}$ |

**关键发现**：简化后系统仍为对称三对角，但对角线和次对角线都是变系数的！

### 3.2 多步 CR 的关键问题

第一步 CR 后：
- 次对角线 $e' = -a_0^2/d_0$（常数，因初始系数常数）
- 对角线 $d'$ 在边界和内部不同

**第二步 CR 必须使用变系数公式**，因为 $d'$ 不再处处相等。此时：
- 简化后的次对角线 $e''$ 也变成变系数
- 第三步及以后所有步都处理一般变系数对称三对角

**这就是之前 l≥2 实现失败的根本原因**：旧代码只用单一 `a_per_row` 跟踪次对角线，无法捕捉变系数次对角线。

### 3.3 正确做法：完整跟踪对称三对角结构

每步 CR 跟踪两个数组：
- `d[0..M-1]`：对角线
- `e[0..M-2]`：次对角线（$e[k]$ 连接行 $k$ 和行 $k+1$）

简化公式（统一处理边界）：
```python
# 对每个偶数行 j = 2k
d_new[k] = cur_d[j]
b_new[k] = cur_b[j]

if j > 0:  # 左奇邻居
    d_new[k] -= cur_e[j-1]**2 / cur_d[j-1]
    b_new[k] -= (cur_e[j-1] / cur_d[j-1]) * cur_b[j-1]

if j < M-1:  # 右奇邻居
    d_new[k] -= cur_e[j]**2 / cur_d[j+1]
    b_new[k] -= (cur_e[j] / cur_d[j+1]) * cur_b[j+1]

# 次对角线
e_new[k] = -cur_e[2k] * cur_e[2k+1] / cur_d[2k+1]
```

### 3.4 回代公式

从最深层到第 0 层，逐层恢复奇数下标值：
$$
u_{2k+1} = \frac{b_{2k+1} - e_{2k} u_{2k} - e_{2k+1} u_{2k+2}}{d_{2k+1}}
$$

---

## 4. 向量化实现

### 4.1 跨模式向量化

所有 N 个 Fourier 模式共享相同的消元模式，仅对角线 $d_\text{eig}[p]$ 不同。

状态用 2D 数组表示：
- `cur_d`: shape `(N_modes, M)`
- `cur_e`: shape `(N_modes, M-1)`  
- `cur_b`: shape `(N_modes, M)`

每步 CR 使用 numpy 花式索引，一次处理所有模式。

### 4.2 性能对比 (n=513, N=511)

| 方法 | 时间 | 误差 |
|------|------|------|
| FA (2D-DST) | 56ms | 3.14e-06 |
| CR (DST+Thomas) | 412ms | 3.14e-06 |
| FACR(l=3, opt) | 124ms | 3.14e-06 |
| FACR(l=7) | 82ms | 3.14e-06 |

> Python 中 FACR 仍慢于 FA，因为 scipy 的 DST 使用优化 FFT 库。  
> 在 Fortran/C 实现中，FACR 的大 N 优势才会体现。

---

## 5. 收敛性验证

五点差分格式，二阶收敛率 ~2.01：

| n | h | error | rate |
|---|---|-------|------|
| 9 | 0.125 | 1.01e-01 | — |
| 17 | 0.0625 | 2.43e-02 | 2.07 |
| 33 | 0.03125 | 6.01e-03 | 2.02 |
| 65 | 0.01562 | 1.50e-03 | 2.01 |
| 129 | 0.00781 | 3.75e-04 | 2.01 |
| 257 | 0.00391 | 9.36e-05 | 2.01 |
| 513 | 0.00195 | 2.34e-05 | 2.01 |

---

## 6. 踩坑记录

### 6.1 次对角线符号错误（关键 BUG）

**错误**：内部行消元后 $d' = d - 2a^2/d$，次对角线写成 $+a^2/d$  
**正确**：$e' = -a^2/d$（负号！）

这导致 FACR l=1 给出 O(1) 误差。修复后 l=1 立刻机器精度。

### 6.2 边界行 d' 不同

边界行只有一个奇邻居，$d' = d - a^2/d$（不是 $d - 2a^2/d$）。

### 6.3 l≥2 失败：单一 a_per_row 无法跟踪变系数次对角线

**根本原因**：第一步 CR 后对角线变系数 → 第二步 CR 产生的次对角线也变系数 → 用单一 `a_per_row` 无法表示。

**解决方案**：完整跟踪对称三对角结构 `(d[], e[])`，每步 CR 保存当前层的 `(d, e, b)` 用于回代。

### 6.4 回代不能用原始 a0, d0

旧代码回代时使用原始层系数 `a_bs=a0, d_bs=d0`，但多步 CR 后方程的系数已改变。

**正确做法**：每层保存自己的 `(d, e, b)`，回代时使用该层系数：
$$
u_j^{(\text{odd})} = \frac{b_j - e_{j-1} u_{j-1} - e_j u_{j+1}}{d_j}
$$

---

## 7. 复杂度分析

| 方法 | 操作量 | 说明 |
|------|--------|------|
| FA | $O(N^2 \log N)$ | 2D-DST 主导 |
| CR | $O(N^2 \log N)$ | 等价于 FA |
| FACR(l) | $O(N^2 \log N) + O(N^2 \cdot l) + O(N \cdot N/2^l)$ | DST + l步CR + Thomas |
| FACR(opt) | $O(N^2 \log\log N)$ | $l \approx \log_2\log_2 N$ |

FACR 的优势在 $N \gg 1$ 时才显现：Thomas 部分 $O(N^2/2^l)$ 随 l 指数下降。

---

*实现文件: `GMRES_Implementation/python/cyclic_reduction.py`*
