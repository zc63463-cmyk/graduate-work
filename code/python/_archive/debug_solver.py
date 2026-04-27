#!/usr/bin/env python3
"""
调试脚本 - 找出为什么谱方法只有2阶收敛
"""

import numpy as np
from scipy.fft import dst, idst

def test_1d_poisson():
    """
    测试1D泊松求解器 - 打印中间结果
    """
    print("="*60)
    print("1D泊松求解器调试")
    print("="*60)
    
    # 使用N=5
    n = 5
    L = 1.0
    x = np.linspace(0, L, n)
    h = L / (n - 1)
    
    print(f"\n网格: n={n}, h={h:.6f}")
    print(f"x = {x}")
    
    # 精确解: u = sin(πx)
    u_exact = np.sin(np.pi * x)
    print(f"u_exact = {u_exact}")
    
    # 右侧: f = π² sin(πx)
    f = np.pi**2 * np.sin(np.pi * x)
    print(f"f = {f}")
    
    # 内部点
    f_int = f[1:-1]
    N = n - 2  # = 3
    
    print(f"\n内部点数量: N = {N}")
    print(f"f_int = {f_int}")
    
    # 正弦变换
    f_hat = dst(f_int, type=1, norm='ortho')
    print(f"\nf_hat = {f_hat}")
    
    # 特征值
    k_vals = np.arange(1, N+1)
    lambda_k = 2.0 / h**2 * (np.cos(k_vals * np.pi / (N + 1)) - 1.0)
    print(f"\nlambda_k = {lambda_k}")
    print(f"-lambda_k = {-lambda_k}")
    
    # 求解
    u_hat = f_hat / (-lambda_k)
    print(f"\nu_hat = {u_hat}")
    
    # 逆正弦变换
    u_int = idst(u_hat, type=1, norm='ortho')
    print(f"\nu_int (数值解) = {u_int}")
    
    # 精确解 (内部点)
    u_exact_int = u_exact[1:-1]
    print(f"u_exact_int (精确解) = {u_exact_int}")
    
    # 误差
    error = np.max(np.abs(u_int - u_exact_int))
    print(f"\n误差 = {error:.6e}")
    
    # 理论分析
    print("\n" + "="*60)
    print("Theory Analysis:")
    print("="*60)
    print(f"For u = sin(pi*x), we have:")
    print(f"  u_xx = -pi^2 sin(pi*x)")
    print(f"  -u_xx = pi^2 sin(pi*x) = f")
    print(f"\nAfter discretization, eigenvalues are:")
    print(f"  lambda_k = 2/h^2 * (cos(k*pi/(N+1)) - 1)")
    print(f"\nFor k=1:")
    print(f"  lambda_1 = 2/h^2 * (cos(pi/(N+1)) - 1)")
    print(f"  lambda_1 = {lambda_k[0]:.6f}")
    print(f"  -lambda_1 = {-lambda_k[0]:.6f}")
    print(f"  Theoretical value: pi^2 = {np.pi**2:.6f}")
    print(f"  Relative error: {abs(-lambda_k[0] - np.pi**2)/np.pi**2:.6e}")


def test_2d_poisson_small():
    """
    测试2D泊松求解器 - 使用小网格
    """
    print("\n" + "="*60)
    print("2D泊松求解器调试 (小网格)")
    print("="*60)
    
    # 使用N=5
    n = 5
    L = 1.0
    x = np.linspace(0, L, n)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = L / (n - 1)
    
    # 精确解: u = sin(πx)sin(πy)
    u_exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # 右侧: f = 2π² sin(πx)sin(πy)
    f = 2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # 内部点
    f_int = f[1:-1, 1:-1]
    N = n - 2  # = 3
    
    # 2D正弦变换
    f_hat = np.zeros((N, N))
    for i in range(N):
        f_hat[i, :] = dst(f_int[i, :], type=1, norm='ortho')
    for j in range(N):
        f_hat[:, j] = dst(f_hat[:, j], type=1, norm='ortho')
    
    print(f"\nf_hat[0,0] = {f_hat[0,0]:.6f}")
    
    # 特征值
    k_vals = np.arange(1, N+1)
    lambda_x = 2.0 / h**2 * (np.cos(k_vals * np.pi / (N + 1)) - 1.0)
    lambda_y = 2.0 / h**2 * (np.cos(k_vals * np.pi / (N + 1)) - 1.0)
    
    # 求解
    u_hat = np.zeros((N, N))
    for kx in range(N):
        for ly in range(N):
            lambda_total = lambda_x[kx] + lambda_y[ly]
            u_hat[kx, ly] = f_hat[kx, ly] / (-lambda_total)
    
    print(f"u_hat[0,0] = {u_hat[0,0]:.6f}")
    
    # 逆2D正弦变换
    u_int = np.zeros((N, N))
    for kx in range(N):
        u_int[kx, :] = idst(u_hat[kx, :], type=1, norm='ortho')
    for ly in range(N):
        u_int[:, ly] = idst(u_int[:, ly], type=1, norm='ortho')
    
    # 完整解
    u_num = np.zeros((n, n))
    u_num[1:-1, 1:-1] = u_int
    
    # 误差
    error = np.max(np.abs(u_num - u_exact))
    print(f"\n最大误差 = {error:.6e}")
    
    # 检查u_num是否等于u_exact
    print("\n数值解 u_num:")
    print(u_num)
    print("\n精确解 u_exact:")
    print(u_exact)


if __name__ == "__main__":
    test_1d_poisson()
    test_2d_poisson_small()
