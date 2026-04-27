#!/usr/bin/env python3
"""
诊断脚本 - 找出谱方法求解器的错误
"""

import numpy as np
from scipy.fft import dst, idst

def test_dst():
    """
    测试正弦变换能否精确表示 sin(πx)
    """
    print("="*60)
    print("测试正弦变换")
    print("="*60)
    
    for n in [5, 9, 17, 33]:
        # 网格
        L = 1.0
        x = np.linspace(0, L, n)
        h = L / (n - 1)
        
        # 函数 u = sin(πx)
        u = np.sin(np.pi * x)
        
        # 内部点
        u_int = u[1:-1]
        N = n - 2
        
        # 正弦变换
        u_hat = dst(u_int, type=1, norm='ortho')
        
        # 逆正弦变换
        u_int_recon = idst(u_hat, type=1, norm='ortho')
        
        # 误差
        error = np.max(np.abs(u_int - u_int_recon))
        print(f"N={n}: 重建误差 = {error:.6e}")
        
        # 检查 u_hat 是否只有一个非零系数
        # 对于 u = sin(πx)，应该只有 k=1 的系数非零
        print(f"  u_hat[0] = {u_hat[0]:.6f}, u_hat[1] = {u_hat[1]:.6f}")
        print(f"  u_hat 最大值位置: {np.argmax(np.abs(u_hat))}")


def test_poisson_solver():
    """
    测试泊松求解器 - 使用最简单的方法
    """
    print("\n" + "="*60)
    print("测试泊松求解器 (最简单版本)")
    print("="*60)
    
    # 问题: -u_xx = f, u(0)=u(1)=0
    # 精确解: u = sin(πx)
    # 右侧: f = π² sin(πx)
    
    for n in [5, 9, 17, 33, 65]:
        # 网格
        L = 1.0
        x = np.linspace(0, L, n)
        h = L / (n - 1)
        
        # 右侧
        f = np.pi**2 * np.sin(np.pi * x)
        
        # 内部点
        f_int = f[1:-1]
        N = n - 2
        
        # 正弦变换 (x方向)
        f_hat = dst(f_int, type=1, norm='ortho')
        
        # 特征值
        # 对于 -u_xx，特征值是 λ_k = 2/h² * (cos(kπ/(N+1)) - 1)
        # 注意: 这是负数
        lambda_k = 2.0 / h**2 * (np.cos(np.arange(1, N+1) * np.pi / (N + 1)) - 1.0)
        
        # 求解
        u_hat = f_hat / (-lambda_k)  # 注意负号
        
        # 逆正弦变换
        u_int = idst(u_hat, type=1, norm='ortho')
        
        # 完整解
        u_num = np.zeros(n)
        u_num[0] = 0.0
        u_num[-1] = 0.0
        u_num[1:-1] = u_int
        
        # 精确解
        u_exact = np.sin(np.pi * x)
        
        # 误差
        error = np.max(np.abs(u_num - u_exact))
        print(f"N={n}: 误差 = {error:.6e}")
        
        # 计算收敛阶数
        if n > 5:
            h_prev = 1.0 / (prev_n - 1)
            h_curr = 1.0 / (n - 1)
            rate = np.log(prev_error/error) / np.log(h_prev/h_curr)
            print(f"  收敛阶数: {rate:.2f}")
        
        prev_n = n
        prev_error = error


def test_poisson_2d():
    """
    测试2D泊松求解器
    """
    print("\n" + "="*60)
    print("测试2D泊松求解器")
    print("="*60)
    
    # 问题: -∇²u = f, u=0 on ∂Ω
    # 精确解: u = sin(πx)sin(πy)
    # 右侧: f = 2π² sin(πx)sin(πy)
    
    for n in [5, 9, 17, 33]:
        # 网格
        L = 1.0
        x = np.linspace(0, L, n)
        y = np.linspace(0, L, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        h = L / (n - 1)
        
        # 右侧
        f = 2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # 内部点
        f_int = f[1:-1, 1:-1]
        N = n - 2
        
        # 2D正弦变换
        f_hat = np.zeros((N, N))
        
        # 对x方向进行正弦变换
        for i in range(N):
            f_hat[i, :] = dst(f_int[i, :], type=1, norm='ortho')
        
        # 对y方向进行正弦变换
        for j in range(N):
            f_hat[:, j] = dst(f_hat[:, j], type=1, norm='ortho')
        
        # 特征值
        # λ_kl = -2/h² * (cos(kπ/(N+1)) - 1) - 2/h² * (cos(lπ/(N+1)) - 1)
        k_vals = np.arange(1, N+1)
        l_vals = np.arange(1, N+1)
        
        lambda_x = 2.0 / h**2 * (np.cos(k_vals * np.pi / (N + 1)) - 1.0)
        lambda_y = 2.0 / h**2 * (np.cos(l_vals * np.pi / (N + 1)) - 1.0)
        
        # 求解
        u_hat = np.zeros((N, N))
        for k_idx, kx in enumerate(k_vals):
            for l_idx, ly in enumerate(l_vals):
                lambda_total = lambda_x[k_idx] + lambda_y[l_idx]
                u_hat[k_idx, l_idx] = f_hat[k_idx, l_idx] / (-lambda_total)
        
        # 逆2D正弦变换
        u_int = np.zeros((N, N))
        
        for k_idx in range(N):
            u_int[k_idx, :] = idst(u_hat[k_idx, :], type=1, norm='ortho')
        
        for l_idx in range(N):
            u_int[:, l_idx] = idst(u_int[:, l_idx], type=1, norm='ortho')
        
        # 完整解
        u_num = np.zeros((n, n))
        u_num[1:-1, 1:-1] = u_int
        
        # 精确解
        u_exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # 误差
        error = np.max(np.abs(u_num - u_exact))
        print(f"N={n}: 误差 = {error:.6e}")


if __name__ == "__main__":
    test_dst()
    test_poisson_solver()
    test_poisson_2d()
