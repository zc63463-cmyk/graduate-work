#!/usr/bin/env python3
"""
Simple debug script - ASCII only
"""

import numpy as np
from scipy.fft import dst, idst

def test_1d():
    """Test 1D Poisson solver with small N"""
    print("="*60)
    print("1D Poisson Solver Test (N=5)")
    print("="*60)
    
    n = 5
    L = 1.0
    x = np.linspace(0, L, n)
    h = L / (n - 1)
    
    print(f"Grid: n={n}, h={h:.6f}")
    
    # Exact solution: u = sin(pi*x)
    u_exact = np.sin(np.pi * x)
    
    # RHS: f = pi^2 * sin(pi*x)
    f = np.pi**2 * np.sin(np.pi * x)
    
    print(f"\nExact solution at internal points:")
    print(f"  u_exact[1:4] = {u_exact[1:-1]}")
    
    print(f"\nRHS at internal points:")
    print(f"  f[1:4] = {f[1:-1]}")
    
    # Internal points
    f_int = f[1:-1]
    N = n - 2  # = 3
    
    # DST
    f_hat = dst(f_int, type=1, norm='ortho')
    print(f"\nDST of f_int:")
    print(f"  f_hat = {f_hat}")
    
    # Eigenvalues
    k_vals = np.arange(1, N+1)
    lambda_k = 2.0 / h**2 * (np.cos(k_vals * np.pi / (N + 1)) - 1.0)
    
    print(f"\nEigenvalues (my formula):")
    print(f"  lambda_k = {lambda_k}")
    print(f"  -lambda_k = {-lambda_k}")
    
    # Theoretical eigenvalues of A (for -u_xx)
    # lambda_A = (2/h^2) * (1 - cos(k*pi/(N+1)))
    lambda_A = (2.0 / h**2) * (1.0 - np.cos(k_vals * np.pi / (N + 1)))
    print(f"\nEigenvalues of A (theoretical):")
    print(f"  lambda_A = {lambda_A}")
    
    # Check relation
    print(f"\nCheck: lambda_k = -lambda_A?")
    print(f"  lambda_k + lambda_A = {lambda_k + lambda_A}")
    
    # Solve
    u_hat = f_hat / lambda_A  # Use correct eigenvalues
    print(f"\nSolution in frequency domain:")
    print(f"  u_hat = {u_hat}")
    
    # Inverse DST
    u_int = idst(u_hat, type=1, norm='ortho')
    print(f"\nNumerical solution at internal points:")
    print(f"  u_num = {u_int}")
    
    # Error
    error = np.max(np.abs(u_int - u_exact[1:-1]))
    print(f"\nMax error = {error:.6e}")


def test_2d_small():
    """Test 2D Poisson solver with small N"""
    print("\n" + "="*60)
    print("2D Poisson Solver Test (N=5)")
    print("="*60)
    
    n = 5
    L = 1.0
    x = np.linspace(0, L, n)
    y = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = L / (n - 1)
    
    # Exact solution: u = sin(pi*x)*sin(pi*y)
    u_exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # RHS: f = 2*pi^2 * sin(pi*x)*sin(pi*y)
    f = 2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # Internal points
    f_int = f[1:-1, 1:-1]
    N = n - 2  # = 3
    
    # 2D DST
    f_hat = np.zeros((N, N))
    for i in range(N):
        f_hat[i, :] = dst(f_int[i, :], type=1, norm='ortho')
    for j in range(N):
        f_hat[:, j] = dst(f_hat[:, j], type=1, norm='ortho')
    
    print(f"f_hat[0,0] = {f_hat[0,0]:.6f}")
    
    # Eigenvalues
    k_vals = np.arange(1, N+1)
    lambda_x = (2.0 / h**2) * (1.0 - np.cos(k_vals * np.pi / (N + 1)))
    lambda_y = (2.0 / h**2) * (1.0 - np.cos(k_vals * np.pi / (N + 1)))
    
    # Solve
    u_hat = np.zeros((N, N))
    for kx in range(N):
        for ly in range(N):
            lambda_total = lambda_x[kx] + lambda_y[ly]
            u_hat[kx, ly] = f_hat[kx, ly] / lambda_total
    
    print(f"u_hat[0,0] = {u_hat[0,0]:.6f}")
    
    # Inverse 2D DST
    u_int = np.zeros((N, N))
    for kx in range(N):
        u_int[kx, :] = idst(u_hat[kx, :], type=1, norm='ortho')
    for ly in range(N):
        u_int[:, ly] = idst(u_int[:, ly], type=1, norm='ortho')
    
    # Full solution
    u_num = np.zeros((n, n))
    u_num[1:-1, 1:-1] = u_int
    
    # Error
    error = np.max(np.abs(u_num - u_exact))
    print(f"\nMax error = {error:.6e}")
    
    # Check if solution is correct
    print("\nNumerical solution:")
    print(u_num)
    print("\nExact solution:")
    print(u_exact)


if __name__ == "__main__":
    test_1d()
    test_2d_small()
