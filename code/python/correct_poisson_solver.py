#!/usr/bin/env python3
"""
正确的泊松方程求解器 - 用于验证和对比

这个脚本实现：
1. 谱方法求解器 (基于正弦变换) - 正确的实现
2. 5点差分 + FFT求解器 - 简化的FFT类算法
3. 9点差分格式 - FFT9算法的核心

先确保基础求解器正确，然后再实现复杂的FFT9算法
"""

import numpy as np
from scipy.fft import dst, idst
import time

def solve_poisson_spectral(nx, ny, f_func, bc_func, sx=1.0, sy=1.0):
    """
    使用谱方法求解泊松方程 -u_xx - u_yy = f
    
    参数:
        nx, ny: 网格点数
        f_func: 右侧函数 f(x,y)
        bc_func: 边界条件函数 u(x,y) on ∂Ω
        sx, sy: 区域尺寸 [0,sx] × [0,sy]
        
    返回:
        U: 数值解 (nx × ny)
    """
    # 坐标
    x = np.linspace(0, sx, nx)
    y = np.linspace(0, sy, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # 网格间距
    hx = sx / (nx - 1)
    hy = sy / (ny - 1)
    
    # 右侧函数
    F = f_func(X, Y)
    
    # 初始化解 (设置边界条件)
    U = np.zeros((nx, ny))
    U[0, :] = bc_func(x[0], y)
    U[-1, :] = bc_func(x[-1], y)
    U[:, 0] = bc_func(x, y[0])
    U[:, -1] = bc_func(x, y[-1])
    
    # 内部点
    N = nx - 2
    M = ny - 2
    
    # 调整右侧 (考虑边界条件)
    # 对于谱方法，我们需要解一个齐次问题
    # 令 u = u_hom + u_part，其中 u_part 满足边界条件
    
    # 方法: 直接求解齐次问题 (假设边界条件为0)
    # 如果边界条件非0，我们需要使用其他方法
    
    F_int = F[1:-1, 1:-1].copy()
    
    # 如果边界条件非0，我们需要调整右侧
    # 但对于许多测试问题，边界条件为0
    
    # 2D正弦变换
    F_hat = np.zeros((N, M))
    
    # 对x方向进行正弦变换
    for i in range(N):
        F_hat[i, :] = dst(F_int[i, :], type=1, norm='ortho')
    
    # 对y方向进行正弦变换
    for j in range(M):
        F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
    
    # 特征值
    # λ_k = -2/h² * (cos(kπ/(N+1)) - 1)
    U_hat = np.zeros((N, M))
    
    for k in range(N):
        lambda_x = 2.0 / hx**2 * (np.cos((k+1) * np.pi / (N + 1)) - 1.0)
        for l in range(M):
            lambda_y = 2.0 / hy**2 * (np.cos((l+1) * np.pi / (M + 1)) - 1.0)
            lambda_total = lambda_x + lambda_y
            if abs(lambda_total) > 1e-12:
                U_hat[k, l] = -F_hat[k, l] / lambda_total  # 注意负号
    
    # 逆2D正弦变换
    U_int = np.zeros((N, M))
    
    for k in range(N):
        U_int[k, :] = idst(U_hat[k, :], type=1, norm='ortho')
    
    for l in range(M):
        U_int[:, l] = idst(U_int[:, l], type=1, norm='ortho')
    
    U[1:-1, 1:-1] = U_int
    
    return U


def solve_poisson_5point_fft(nx, ny, f_func, bc_func, sx=1.0, sy=1.0):
    """
    使用5点差分格式 + FFT求解泊松方程
    
    这是对FFT9算法的简化 (使用5点格式而不是9点格式)
    """
    # 坐标
    x = np.linspace(0, sx, nx)
    y = np.linspace(0, sy, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # 网格间距
    hx = sx / (nx - 1)
    hy = sy / (ny - 1)
    
    # 右侧函数
    F = f_func(X, Y)
    
    # 初始化解
    U = np.zeros((nx, ny))
    U[0, :] = bc_func(x[0], y)
    U[-1, :] = bc_func(x[-1], y)
    U[:, 0] = bc_func(x, y[0])
    U[:, -1] = bc_func(x, y[-1])
    
    # 内部点
    N = nx - 2
    M = ny - 2
    
    # 构建调整后的右侧 (考虑边界条件)
    F_adj = np.zeros((N, M))
    
    for i in range(N):
        for j in range(M):
            x_idx = i + 1
            y_idx = j + 1
            
            F_adj[i, j] = F[x_idx, y_idx]
            
            # 边界条件贡献
            if x_idx == 1:
                F_adj[i, j] -= U[0, y_idx] / hx**2
            if x_idx == nx - 2:
                F_adj[i, j] -= U[-1, y_idx] / hx**2
            if y_idx == 1:
                F_adj[i, j] -= U[x_idx, 0] / hy**2
            if y_idx == ny - 2:
                F_adj[i, j] -= U[x_idx, -1] / hy**2
    
    # 傅里叶正弦变换
    F_hat = np.zeros((N, M))
    
    for i in range(N):
        F_hat[i, :] = dst(F_adj[i, :], type=1, norm='ortho')
    
    for j in range(M):
        F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
    
    # 求解
    U_hat = np.zeros((N, M))
    
    for k in range(N):
        lambda_x = 2.0 / hx**2 * (np.cos((k+1) * np.pi / (N + 1)) - 1.0)
        for l in range(M):
            lambda_y = 2.0 / hy**2 * (np.cos((l+1) * np.pi / (M + 1)) - 1.0)
            lambda_total = lambda_x + lambda_y
            if abs(lambda_total) > 1e-12:
                U_hat[k, l] = -F_hat[k, l] / lambda_total
    
    # 逆傅里叶正弦变换
    U_int = np.zeros((N, M))
    
    for k in range(N):
        U_int[k, :] = idst(U_hat[k, :], type=1, norm='ortho')
    
    for l in range(M):
        U_int[:, l] = idst(U_int[:, l], type=1, norm='ortho')
    
    U[1:-1, 1:-1] = U_int
    
    return U


def test_solver():
    """
    测试求解器
    """
    print("="*60)
    print("测试泊松方程求解器")
    print("="*60)
    
    # 测试问题1: u = sin(πx)sin(πy), 边界条件为0
    print("\n测试问题1: u = sin(πx)sin(πy) (零边界条件)")
    
    def u_exact1(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def f1(x, y):
        return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def bc1(x, y):
        return 0.0  # 零边界条件
    
    print("\n谱方法求解器:")
    for n in [5, 9, 17, 33, 65]:
        U = solve_poisson_spectral(n, n, f1, bc1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        U_exact = u_exact1(X, Y)
        
        error = np.max(np.abs(U - U_exact))
        print(f"  N={n}: 误差 = {error:.6e}")
    
    print("\n5点差分 + FFT求解器:")
    for n in [5, 9, 17, 33, 65]:
        U = solve_poisson_5point_fft(n, n, f1, bc1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        U_exact = u_exact1(X, Y)
        
        error = np.max(np.abs(U - U_exact))
        print(f"  N={n}: 误差 = {error:.6e}")
    
    # 测试问题2: u = x(1-x)y(1-y), 边界条件非0
    print("\n" + "="*60)
    print("测试问题2: u = x(1-x)y(1-y) (非零边界条件)")
    
    def u_exact2(x, y):
        return x * (1.0 - x) * y * (1.0 - y)
    
    def f2(x, y):
        # -∇²u = - (u_xx + u_yy)
        # u_xx = y(1-y) * (-2)
        # u_yy = x(1-x) * (-2)
        return -2.0 * y * (1.0 - y) - 2.0 * x * (1.0 - x)
    
    def bc2(x, y):
        return u_exact2(x, y)
    
    print("\n谱方法求解器 (注意: 当前实现假设零边界条件)")
    print("对于非零边界条件，需要使用其他方法")
    
    # 对于这个测试，我们暂时跳过，因为我的实现还不支持非零边界条件


if __name__ == "__main__":
    test_solver()
