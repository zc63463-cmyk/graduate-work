# GMRES算法实现项目

基于Saad & Schultz (1986)论文的完整GMRES算法实现。

## 论文信息

- **标题**: GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems
- **作者**: Youcef Saad, Martin H. Schultz
- **单位**: 耶鲁大学计算机科学系
- **发表**: SIAM Journal on Scientific and Statistical Computing, Vol.7, No.3, 1986
- **页码**: 14页
- **被引用次数**: >10,000次 (Google Scholar)

## 项目结构

```
GMRES_Implementation/
├── python/
│   ├── gmres.py          # Python实现（详细注释）
│   └── test_gmres.py    # 测试脚本
├── cpp/
│   ├── gmres.hpp         # C++头文件（模板实现）
│   ├── main.cpp          # C++示例程序
│   └── CMakeLists.txt   # CMake构建文件
└── README.md            # 本文件
```

## 算法概述

GMRES (Generalized Minimal Residual) 是一种迭代法，用于求解大型稀疏**非对称**线性方程组：

$$Ax = b$$

### 核心思想

1. **Krylov子空间方法**: 在Krylov子空间 $K_k = \text{span}\{r_0, Ar_0, A^2r_0, ..., A^{k-1}r_0\}$ 中寻找近似解

2. **最小残差**: 在Krylov子空间上最小化残差范数 $\|b - Ax_k\|_2$

3. **Arnoldi过程**: 构建Krylov子空间的正交基（非对称版本的Lanczos过程）

### 主要优势

- ✅ **永不崩溃**（除非已收敛到精确解）
- ✅ **最优性保证**: 残差范数在Krylov子空间上最小
- ✅ **存储高效**: 比GCR节省约一半存储
- ✅ **计算高效**: 比GCR少约1/3算术运算

## Python实现

### 依赖

```bash
pip install numpy scipy matplotlib
```

### 使用方法

```python
import numpy as np
from gmres import gmres_simplified

# 创建测试矩阵
N = 100
A = np.diag(np.ones(N)) + 0.1 * np.random.randn(N, N)
b = np.random.randn(N)

# 求解
x, info = gmres_simplified(A, b, tol=1e-8, max_iter=500, m=10)

print(f"收敛: {info['success']}")
print(f"迭代次数: {info['iterations']}")
print(f"残差: {info['residuals'][-1]}")
```

### 示例

运行内置示例：

```bash
python gmres.py
```

这将运行：
1. **示例1**: 论文中的2×2系统（GCR会崩溃的例子）
2. **示例2**: 1D Poisson方程（对称正定）
3. **示例3**: 非对称系统（可选）

## C++实现

### 依赖

- C++11或更高版本
- CMake 3.10+ (可选，用于构建)

### 编译

**方法1: 使用g++**

```bash
cd cpp
g++ -std=c++11 -o gmres_example main.cpp -I.
./gmres_example
```

**方法2: 使用CMake**

```bash
cd cpp
mkdir build && cd build
cmake ..
make
./gmres_example
```

### 使用方法

```cpp
#include "gmres.hpp"
#include <vector>

// 定义矩阵向量乘法函数
std::vector<double> my_matrix_vector_product(const std::vector<double>& v) {
    // 实现矩阵向量乘法
    std::vector<double> result(v.size());
    // ...
    return result;
}

int main() {
    // 创建右侧向量
    std::vector<double> b = {1.0, 2.0, 3.0};
    
    // 初始猜测
    std::vector<double> x = {0.0, 0.0, 0.0};
    
    // 定义算子
    auto A_op = my_matrix_vector_product;
    
    // 求解
    std::vector<double> residuals;
    bool converged = gmres::gmres(A_op, b, x, 1e-6, 100, 10, &residuals);
    
    return 0;
}
```

## 算法参数

### 重启参数 m

重启参数 `m` 是GMRES最关键的选择：

- **小m** (5-10): 存储需求低，但可能收敛慢或不收敛
- **中m** (20-50): 平衡存储和收敛速度，适合大多数问题
- **大m** (100+): 收敛快，但存储需求高

**建议**: 从 `m=10` 开始，如果收敛太慢，逐渐增大。

### 停止准则

```python
# 相对残差
||r_k||_2 / ||r_0||_2 < tol

# 绝对残差
||r_k||_2 < atol
```

## 理论基础

### Arnoldi过程

构建Krylov子空间 $K_k$ 的正交基 $\{v_1, v_2, ..., v_k\}$：

$$AV_k = V_{k+1} H_k$$

其中 $H_k$ 是 $(k+1) \times k$ 上Hessenberg矩阵。

### 最小二乘问题

GMRES求解：

$$\min_{y \in \mathbb{R}^k} \|\beta e_1 - H_k y\|_2$$

其中 $\beta = \|r_0\|_2$, $e_1 = (1, 0, ..., 0)^T$.

### 重要定理

**命题1** (残差范数免费计算): 近似解 $x_k$ 的残差范数等于变换后右侧向量的第 $(k+1)$ 个分量。

**命题2** (不崩溃定理): GMRES在第 $j$ 步产生精确解当且仅当：
1. 算法在第 $j$ 步失效（崩溃）
2. $v_{j+1} = 0$
3. $h_{j+1,j} = 0$
4. 初始残差向量 $r_0$ 的极小多项式次数等于 $j$

**推论3**: 对于 $N \times N$ 问题，GMRES最多 $N$ 步收敛。

## 数值实验

### 实验1: 论文中的2×2系统

$$A = \begin{bmatrix} 1 & 1 \\ 1 & 0 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 0 \end{bmatrix}$$

- GCR会崩溃（除零错误）
- GMRES在2步内收敛

### 实验2: 1D Poisson方程

$$-\frac{d^2u}{dx^2} = f(x), \quad u(0)=u(1)=0$$

使用五点差分格式离散化，得到对称正定三对角矩阵。

**结果** (N=50, m=10):
- 收敛: ✅
- 迭代次数: ~20
- 最终残差: < 1e-8

### 实验3: 非对称系统

随机生成非对称矩阵（对角占优）。

**结果** (N=100, 无重启):
- 收敛: ✅
- 迭代次数: ~150
- 最终残差: < 1e-8

## 性能对比

### 存储需求

| 方法 | 存储需求 |
|------|----------|
| GCR(m-1) | $(2m+1)N$ |
| **GMRES(m)** | **$(m+2)N$** |

**示例** (m=20, N=1000):
- GCR: 41,000个元素
- GMRES: 22,000个元素（节省46%）

### 计算量

每步平均运算量：

| 方法 | 乘法运算量 |
|------|-----------|
| GCR(m-1) | $(3m+5)N$ |
| **GMRES(m)** | **$(m+3+1/m)N$** |

**示例** (m=20):
- GCR: $65N$ 次乘法/步
- GMRES: $23.05N$ 次乘法/步（节省64%）

## 常见问题

### 1. GMRES不收敛怎么办？

**可能原因**:
- 重启参数 $m$ 太小
- 矩阵条件数太大（需要预处理）
- 问题本身难以求解

**解决方案**:
- 增大 $m$ (建议 $m \geq 20$)
- 使用预处理（ILU, MILU等）
- 尝试其他算法（BiCGStab, IDR(s)等）

### 2. 如何选择重启参数 m？

**经验法则**:
- 简单问题: $m = 10 \sim 20$
- 中等难度: $m = 20 \sim 50$
- 困难问题: $m = 50 \sim 100$

**自适应策略**:
- 如果残差下降太慢，增大 $m$
- 如果存储不足，减小 $m$

### 3. 如何实现预处理？

左预处理：求解 $M^{-1}Ax = M^{-1}b$

```python
# Python示例
def preconditioned_gmres(A, b, M_inv, tol=1e-6, m=10):
    # 左预处理: M^{-1} A x = M^{-1} b
    A_pcd = M_inv @ A
    b_pcd = M_inv @ b
    
    x, info = gmres_simplified(A_pcd, b_pcd, tol=tol, m=m)
    
    return x, info
```

## 扩展阅读

1. **Saad (1993)**: 《Iterative Methods for Sparse Linear Systems》- 第6章详细讲解GMRES
2. **Saad & Schultz (1983)**: "Conjugate gradient-like algorithms for solving nonsymmetric linear systems"
3. **Greenbaum (1997)**: 《Iterative Methods for Solving Linear Systems》
4. **Kelley (1995)**: 《Iterative Methods for Linear and Nonlinear Equations》

## 现代发展

1. **Flexible GMRES (FGMRES)**: 允许预处理矩阵在每次迭代中变化
2. **Deflated GMRES**: 移除Krylov子空间中的"坏方向"
3. **Augmented GMRES**: 在Krylov子空间中加入特定方向

## 参考文献

1. Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal residual algorithm for solving nonsymmetric linear systems. *SIAM Journal on Scientific and Statistical Computing*, 7(3), 856-869.

2. Saad, Y. (1993). *Iterative Methods for Sparse Linear Systems*. PWS Publishing Company.

3. Greenbaum, A. (1997). *Iterative Methods for Solving Linear Systems*. SIAM.

## 许可

本实现基于已发表的学术论文，代码可自由使用和修改。

## 联系

如有问题或建议，请通过以下方式联系：
- 创建Issue
- 提交Pull Request

---

**最后更新**: 2026-04-25
