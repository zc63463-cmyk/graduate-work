---
title: "Algorithm 543: FFT9 - 算法实现论文精读"
authors: "E. N. Houstis, T. S. Papatheodorou"
journal: "ACM Transactions on Mathematical Software (TOMS)"
year: 1979
volume: "Vol.5, No.4"
pages: "490-493"
doi: "10.1145/355853.355865"
tags: ["FFT9", "Fortran", "算法实现", "ACM算法", "椭圆方程"]
date: 1979-12-01
---

# Algorithm 543: FFT9 - Fast Solution of Helmholtz-Type PDEs

## 📋 论文基本信息

- **标题**: Algorithm 543: FFT9, Fast Solution of Helmholtz-Type Partial Differential Equations [D3]
- **作者**: 
  - E. N. Houstis (Purdue University)
  - T. S. Papatheodorou (Clarkson College of Technology)
- **发表**: ACM Transactions on Mathematical Software (TOMS), Volume 5, Issue 4, December 1979, Pages 490-494
- **DOI**: https://doi.org/10.1145/355853.355865
- **ACM算法编号**: 543
- **引用次数**: 8次
- **下载量**: 539次
- **编程语言**: Fortran
- **分类**: CR Categories 5.17, D3

---

## 📖 论文与论文1的关系

这篇论文是**论文1的补充**，提供了FFT9算法的**完整Fortran实现**。

| 论文 | 内容 | 页数 |
|------|------|------|
| **论文1** | 方法描述、数值实验、理论分析 | 12页 |
| **论文2** | Fortran代码实现、使用说明 | 5页 |

**引用关系**：
> "The algorithm given here is a complement to [1] where the description, test results, and references are given."
> 参考文献[1]: Houstis, E.N., AND PAPATHEODOROU, T.S. High-order fast elliptic equation solver. ACM Trans. Math. Software 5, 4 (Dec. 1979), 431-441.

---

## 💻 算法实现结构

### Fortran程序模块

FFT9算法由以下Fortran子程序组成：

| 模块名 | 记录数 | 功能描述 |
|--------|--------|----------|
| `MAIN` | 299 | 主程序，控制算法流程 |
| `RGHTSD` | 208 | 计算右侧向量（Right-hand side） |
| `EQSOL` | 77 | 方程求解模块 |
| `DISCRT` | 95 | 离散化模块（Discretization） |
| `EVENRD` | 55 | 偶数行归约（Even reduction） |
| `ODDRD` | 54 | 奇数行求解（Odd reduction） |
| `CRED` | 76 | 循环归约求解（Cyclic reduction） |
| `STOREX` | 18 | 存储管理X |
| `FETCHX` | 18 | 提取存储X |
| `TFOLD` | 138 | 傅里叶变换折叠（Fourier fold） |
| `FOUR` | 51 | 傅里叶分析/合成（Fourier analysis/synthesis） |
| `NEG` | 12 | 处理负数 |
| `ZERO` | 13 | 清零操作 |
| `SETF` | 41 | 设置参数F |
| `SUMARY` | 76 | 汇总输出（Summary） |
| `PDERGHT` | 7 | PDE右侧函数定义（用户提供） |
| `TRUE T` | 5 | 精确解函数（用户提供，可选） |
| `PDE T` | 9 | PDE系数定义（用户提供） |
| `DATA D` | 3 | 测试数据 |
| `BCOND` | 7 | 边界条件定义（用户提供） |

**总计**: 约1500行Fortran代码

---

## 📐 使用方法

### 1. 用户提供的子程序

FFT9需要用户编写以下Fortran子程序：

#### (1) PDE系数定义 - `SUBROUTINE PDE(X,Y,CVALUS)`

```fortran
SUBROUTINE PDE(X,Y,CVALUS)
REAL CVALUS(7)
CVALUS(1) = CUXX   ! u_xx的系数
CVALUS(2) = CUYY   ! u_yy的系数 (论文中似乎是CUYY，但代码显示为索引2)
CVALUS(3) = ...    ! 其他系数
CVALUS(6) = CU     ! u的系数
RETURN
END
```

**注意**：根据论文1的公式(1.1): $\alpha u_{xx} + \beta u_{yy} + \gamma u = f$
- `CVALUS(1)` = $\alpha$
- `CVALUS(3)` = $\beta$ (代码中的索引可能需要验证)
- `CVALUS(6)` = $\gamma$

#### (2) 右侧函数 - `REAL FUNCTION PDERGHT(X,Y)`

```fortran
REAL FUNCTION PDERGHT(X,Y)
PDERGHT = f(X,Y)   ! 定义 f(x,y)
RETURN
END
```

#### (3) 边界条件 - `REAL FUNCTION BCOND(I,X,Y,BVALUS)`

```fortran
REAL FUNCTION BCOND(I,X,Y,BVALUS)
REAL BVALUS(4)
GO TO (101,102,103,104), I

101 BVALUS(4) = C   ! 边1: u = C
    BCOND = BVALUS(4)
    RETURN

102 BVALUS(4) = G   ! 边2: u = G
    BCOND = BVALUS(4)
    RETURN

103 BVALUS(4) = C   ! 边3: u = C
    BCOND = BVALUS(4)
    RETURN

104 BVALUS(4) = C   ! 边4: u = C
    BCOND = BVALUS(4)
    RETURN
END
```

**边界编号**：
- I=1: 边1
- I=2: 边2
- I=3: 边3
- I=4: 边4

#### (4) 精确解（可选）- `REAL FUNCTION TRUE(X,Y)`

```fortran
REAL FUNCTION TRUE(X,Y)
TRUE = u_exact(X,Y)   ! 精确解（如果已知）
RETURN
END
```

### 2. 输入参数

#### 区域和网格定义

```fortran
READ(5,1000) SX, SY, NGRIDX, NGRIDY
1000 FORMAT(2F10.4, 2I5)
```

- `SX, SY`: 矩形区域的边长
- `NGRIDX, NGRIDY`: 网格线数（满足 $N_x = 2^{IQX} + 1, N_y = 2^{IQY} + 1$）

#### 输出控制

```fortran
READ(5,1002) LEVEL, NRUNS, ORDER
1002 FORMAT(3I2)
```

- `LEVEL`: 输出级别
  - `= 0`: 仅打印最大误差
  - `= 1`: 打印最大误差和最大相对误差（如果已知精确解）
  - `= 2`: 还打印精确解和近似解在内部网格点的值
- `NRUNS`: 连续运行次数（每次网格尺寸减半）
- `ORDER`: 差分格式阶数
  - `= 4`: 选择四阶差分近似
  - `= 6`: 选择六阶差分近似（仅用于Poisson型算子）

### 3. 存储需求

假设 `3 ≤ IQX, IQY ≤ 7`（即 $5 \leq N_x, N_y \leq 257$）

**数组维度公式**：

- `WORK`: $7 + IQ \times \max(NGRIDX, NGRIDY) + 2 \times NGRIDX \times NGRIDY$
- `Z, Y`: $7 + IQ \times \max(NGRIDX, NGRIDY) + NGRIDX \times NGRIDY$
- `CORE(K)`: 近似解在节点K的值，维度 $NCORE = NX + 2 + (NX+1) \times (NY+1)$
- `PTINT()`: 半格点处的G值，维度 $NPTINT = NX \times NY$（仅用于六阶格式）
- `GRIDX, GRIDY`: 网格坐标，维度 $NGRIDX = NX+1, NGRIDY = NY+1$
- `AKX(I+1)`: 对角块的第I个特征值除以非对角块的第I个特征值的平方减2，维度 $NAKXD = \max(NX, NY) + 2$
- `ED(I+1)`: 非对角块的第I个特征值，维度 $NEDD = \max(NX, NY) + 2$
- `INDEX`: FFT分析中使用的索引向量，维度 $MAX(NX, NY)$
- `SINES`: FFT分析和合成用的正弦向量，维度 $MAX(NX, NY)$

---

## 🔧 算法核心模块详解

### 1. 主程序 `MAIN`

**功能**：
- 读取输入参数
- 调用离散化模块
- 控制求解流程
- 输出结果

**核心变量**：
- `ORDER`: 4（四阶）或6（六阶）
- `NX, NY`: 内部网格点数（$N_x - 1, N_y - 1$）
- `HX, HY`: 网格尺寸

### 2. 离散化模块 `DISCRT`

**功能**：
- 根据 `ORDER` 选择差分格式
- 计算9点模板的权重
- 构造差分方程

**四阶格式权重**（公式2.1 in 论文1）：
$$
a = \frac{12}{p^2 + q}, \quad b = \frac{12p^2}{p^2 + q} - 2
$$
其中 $p = h_x^2/\alpha, q = h_y^2/\beta$

### 3. 右侧计算 `RGHTSD`

**功能**：
- 计算差分方程的右侧向量 $b_j$
- 对于六阶格式，在半格点计算 $f$

### 4. 偶数行归约 `EVENRD`

**功能**：
- 实现奇数-偶数归约的偶数步
- 将偶数行方程解耦

**数学原理**（论文1第2.3节）：
$$
h^2 u_{j-2} + (2h^2 - F^2) u_j + h^2 u_{j+2} = b_j^*
$$

### 5. 奇数行求解 `ODDRD`

**功能**：
- 计算奇数行方程的右侧
- 调用 `CRED` 求解三对角系统

**数学原理**：
$$
F u_j = b_j - h(u_{j-1} + u_{j+1})
$$

### 6. 循环归约求解 `CRED`

**功能**：
- 使用递归循环归约求解三对角系统
- 基于Hockney的算法[4]

**参考文献**：
> Hockney, R.W. "The potential calculation and some applications." Meth. Comput. Phys. 9 (1970), 135-211.

### 7. 傅里叶变换 `FOUR`

**功能**：
- 执行快速傅里叶变换（FFT）
- 用于傅里叶分析和合成

**变换公式**（论文1公式2.2）：
$$
b_j^* = \sum_{k=1}^N \tilde{b}_k V_k, \quad u_j = \sum_{k=1}^N \tilde{u}_{jk} V_k
$$

其中 $V_k = \sqrt{\frac{2}{N+1}} (\sin k\theta, \sin 2k\theta, ..., \sin Nk\theta)^T$

---

## 📊 数值实验摘要

论文1中提供了详细的数值实验结果。关键发现：

### 1. 精度对比

| 问题类型 | FFT9(4阶) | FFT9(6阶) | NCAR(2阶) | GLFFTP(4阶) | BUNPT(4阶) |
|----------|------------|------------|-------------|--------------|--------------|
| 光滑Poisson | ✓✓✓ | ✓✓✓✓ | ✓ | ✓✓ | ✓✓ |
| 奇异导数 | ✓✓ | ✓✓✓✓ | ✓ | - | ✓✓ |
| 边界层 | ✓✓ | - | ✓ | - | ✓✓ |
| 波前 | ✓ | ✓✓✓ | - | - | ✓ |

### 2. 效率对比

- **FFT9 vs NCAR**: 快100-300倍
- **FFT9 vs GLFFTP**: 快28倍
- **FFT9 vs BUNPT**: 快3-17倍（取决于问题）

---

## 💡 使用建议

### 1. 何时使用FFT9

✅ **适合**：
- 矩形区域上的椭圆PDE
- Dirichlet边界条件
- 常数系数
- 需要高精度和快速求解

❌ **不适合**：
- 变系数问题
- 复杂几何区域
- Neumann或Robin边界条件

### 2. 阶数选择

- **四阶格式**：适用于Helmholtz型方程（$\gamma \neq 0$）
- **六阶格式**：仅适用于Poisson方程（$\alpha = \beta = 1, \gamma = 0$）

### 3. 网格选择

- 使用 $N = 2^k + 1$ 的网格
- 对于光滑解，较粗的网格即可达到高精度
- 对于奇异解，需要较细的网格

---

## 🔄 与论文1的对照阅读

| 论文1（理论） | 论文2（实现） |
|--------------|--------------|
| 第2.1节：网格模块 | `DISCRT` - 网格设置 |
| 第2.2节：离散化模块 | `DISCRT` - 差分格式 |
| 第2.3节：方程求解模块 | `EVENRD`, `ODDRD`, `CRED` |
| 公式(2.2)：傅里叶变换 | `FOUR` - FFT子程序 |
| 第3节：实验结果 | `SUMARY` - 输出汇总 |

**建议阅读顺序**：
1. 先读论文1，理解算法原理
2. 再读论文2，了解实现细节
3. 结合两篇论文，实现自己的版本（Python/C++）

---

## 📝 Fortran代码特点（1979年）

### 1. 编码风格

- 使用 `FORMAT` 语句进行格式化I/O
- 使用 `COMMON` 块进行全局变量传递
- 使用 `EQUIVALENCE` 共享存储
- 固定格式（Fixed-format）Fortran 66/77

### 2. 数值计算特点

- 单精度算术（Single precision）
- 在CDC 6500计算机上测试
- 使用 `COS` 函数计算特征值
- 使用改进的Gram-Schmidt正交化（虽然主要用于FFT）

### 3. 现代改进方向

- 改为自由格式Fortran 90/95
- 使用模块（Module）代替COMMON块
- 双精度/四精度算术
- 并行化（OpenMP/MPI）
- 使用现代FFT库（FFTW）

---

## 🔗 相关资源

### 获取完整代码

1. **ACM算法分发服务**：
   - 见期刊内页订单表
   - 或查阅 "Collected Algorithms from ACM"

2. **在线获取**：
   - DOI: https://doi.org/10.1145/355853.355865
   - ACM Digital Library

### 相关算法

- **ACM Algorithm 541**: 可能与本文相关
- **NCAR软件包**: 5点块循环归约方法
- **ELLAPCK**: 椭圆方程求解器系统

---

## 📚 精读心得

### 1. 论文结构评价

**优点**：
- ✅ 理论与实现分离，便于阅读
- ✅ 详细的子程序说明
- ✅ 完整的输入输出规范

**不足**：
- ❌ Fortran代码可读性较差（固定格式）
- ❌ 缺乏现代预处理技术
- ❌ 仅适用于常数系数

### 2. 历史意义

这篇论文是**1970年代快速椭圆求解器**的代表作：
- 展示了高阶差分的优势
- 推广了FFT在PDE求解中的应用
- 为现代快速求解器（如FFTW, PFFT）奠定基础

### 3. 现代复现建议

如果要复现该算法，建议：

1. **使用Python/NumPy**：
   - 利用NumPy的FFT（`np.fft`）
   - 使用稀疏矩阵存储
   - 利用向量化操作

2. **使用C++**：
   - 使用Eigen或Armadillo库
   - 调用FFTW进行FFT
   - 使用模板实现泛型编程

3. **并行化**：
   - 使用OpenMP并行FFT
   - 使用MPI进行域分解
   - GPU加速（CUDA/OpenCL）

---

## 🔗 相关笔记

- [[FFT9_论文1_精读|FFT9论文1 - 理论与实验]]
- [[GMRES_论文精读|GMRES算法论文精读]]
- [[椭圆方程数值解综述|椭圆方程数值方法综述]]

---

**精读日期**: 2026-04-25
**精读人**: WorkBuddy AI Assistant
