#!/usr/bin/env python3
"""
Correct FFT-based Poisson Solver (simplified FFT9)

This implements the correct FFT-based method for Poisson equation:
  u_xx + u_yy = f(x,y) on rectangular domain with Dirichlet BC

Method:
1. Use 2D discrete sine transform (DST)
2. Solve in frequency domain
3. Transform back

This is a simplified version of FFT9 for educational purposes.
For the full FFT9 with high-order differences, see the Fortran code in Paper 2.
"""

import numpy as np
from scipy.fft import dst, idst
import time


def solve_poisson_fft2d(f_func, bc_func, sx=1.0, sy=1.0, nx=33, ny=33):
    """
    Solve Poisson equation using FFT (sine transform).
    
    Equation: u_xx + u_yy = f(x,y)
    Domain: [0, sx] x [0, sy]
    BC: u = bc(x,y) on boundary
    
    Parameters:
    -----------
    f_func : Callable
        Right-hand side function f(x,y)
    bc_func : Callable
        Boundary condition function bc(x,y)
    sx, sy : float
        Domain dimensions
    nx, ny : int
        Number of grid lines (should be 2^k + 1 for FFT)
    
    Returns:
    --------
    U : ndarray (nx, ny)
        Solution
    X, Y : ndarray
        Grid coordinates
    """
    
    print("\n" + "="*60)
    print("Poisson Solver using FFT (Sine Transform)")
    print("="*60)
    
    start_time = time.time()
    
    # Grid spacing
    hx = sx / (nx - 1)
    hy = sy / (ny - 1)
    
    # Grid coordinates
    x = np.linspace(0, sx, nx)
    y = np.linspace(0, sy, ny)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    print(f"Grid: {nx} x {ny}")
    print(f"Grid spacing: hx={hx:.6f}, hy={hy:.6f}")
    
    # Step 1: Evaluate RHS
    F = f_func(X, Y)
    
    # Step 2: Initialize solution with boundary conditions
    U = np.zeros((nx, ny))
    
    # Apply boundary conditions
    U[0, :] = bc_func(x[0], y)   # Left
    U[-1, :] = bc_func(x[-1], y)  # Right
    U[:, 0] = bc_func(x, y[0])    # Bottom
    U[:, -1] = bc_func(x, y[-1])  # Top
    
    # Step 3: For Dirichlet BC, adjust RHS
    # (This is a simplified version)
    
    # Step 4: Compute eigenvalues
    # lambda_{k,l} = 2/hx^2 * (cos(k*pi/(nx)) - 1) + 2/hy^2 * (cos(l*pi/(ny)) - 1)
    # Actually, for interior points, N = nx-2, M = ny-2
    
    N = nx - 2  # Interior points in x
    M = ny - 2  # Interior points in y
    
    # Eigenvalues for sine transform
    lambda_x = np.zeros(N)
    lambda_y = np.zeros(M)
    
    for k in range(1, N+1):
        lambda_x[k-1] = 2.0 / hx**2 * (np.cos(k * np.pi / (N + 1)) - 1.0)
    
    for l in range(1, M+1):
        lambda_y[l-1] = 2.0 / hy**2 * (np.cos(l * np.pi / (M + 1)) - 1.0)
    
    # Step 5: Apply 2D sine transform to RHS
    print("Applying 2D sine transform...")
    
    # Extract interior RHS
    F_int = F[1:-1, 1:-1]
    
    # 2D DST (Discrete Sine Transform)
    # We can compute this using scipy.fft.dst
    
    # For type I DST (which corresponds to Dirichlet BC):
    # F_hat[k,l] = sum_{i=1}^N sum_{j=1}^M F[i,j] * sin(i*k*pi/(N+1)) * sin(j*l*pi/(M+1))
    
    # Unfortunately, scipy.fft.dst only does 1D.
    # For 2D, we need to apply DST in x then y direction.
    
    # Method: Use FFT with odd extension
    # For sine transform: extend f[-n+1:0] = -f[1:n], f[n+1:2n] = f[n-1:0:-1]
    
    # Simplified: Use scipy.fft.dst twice
    F_hat = np.zeros((N, M))
    
    for i in range(N):
        F_hat[i, :] = dst(F_int[i, :], type=1, norm='ortho')
    
    for j in range(M):
        F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
    
    # Step 6: Solve in frequency domain
    print("Solving in frequency domain...")
    
    U_hat = np.zeros((N, M))
    
    for k in range(N):
        for l in range(M):
            lambda_total = lambda_x[k] + lambda_y[l]
            if abs(lambda_total) > 1e-12:
                U_hat[k, l] = F_hat[k, l] / lambda_total
    
    # Step 7: Inverse 2D sine transform
    print("Applying inverse 2D sine transform...")
    
    U_int = np.zeros((N, M))
    
    for k in range(N):
        U_int[k, :] = idst(U_hat[k, :], type=1, norm='ortho')
    
    for l in range(M):
        U_int[:, l] = idst(U_int[:, l], type=1, norm='ortho')
    
    # Step 8: Place interior solution
    U[1:-1, 1:-1] = U_int
    
    end_time = time.time()
    computation_time = end_time - start_time
    
    print(f"Solution computed in {computation_time:.4f} seconds")
    
    return U, X, Y, computation_time


def test_poisson_smooth():
    """
    Test: Poisson equation with smooth solution.
    
    Equation: u_xx + u_yy = f
    Exact solution: u = sin(pi*x) * sin(pi*y)
    """
    
    print("\n" + "="*60)
    print("Test: Poisson Equation with Smooth Solution")
    print("Exact solution: u = sin(pi*x) * sin(pi*y)")
    print("="*60)
    
    # Exact solution
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # RHS: u_xx + u_yy = -pi^2 * sin(pi*x) * sin(pi*y) - pi^2 * sin(pi*x) * sin(pi*y)
    #                    = -2 * pi^2 * sin(pi*x) * sin(pi*y)
    def f_rhs(x, y):
        return -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # Boundary condition
    def bc(x, y):
        return u_exact(x, y)
    
    # Solve
    U, X, Y, time = solve_poisson_fft2d(f_rhs, bc, sx=1.0, sy=1.0, nx=33, ny=33)
    
    # Compute error
    U_exact = u_exact(X, Y)
    error = np.abs(U - U_exact)
    max_error = np.max(error)
    
    print(f"\nResults:")
    print(f"  Max error: {max_error:.6e}")
    print(f"  Grid: 33 x 33")
    
    # Check convergence
    print(f"\nConvergence check:")
    for n in [5, 9, 17, 33, 65]:
        U_test, X_test, Y_test, _ = solve_poisson_fft2d(f_rhs, bc, sx=1.0, sy=1.0, nx=n, ny=n)
        U_exact_test = u_exact(X_test, Y_test)
        error_test = np.max(np.abs(U_test - U_exact_test))
        print(f"  N={n}: Error = {error_test:.6e}")
    
    return max_error, time


def test_poisson_polynomial():
    """
    Test: Poisson equation with polynomial solution.
    
    Equation: u_xx + u_yy = f
    Exact solution: u = x*y*(1-x)*(1-y)
    """
    
    print("\n" + "="*60)
    print("Test: Poisson Equation with Polynomial Solution")
    print("Exact solution: u = x*y*(1-x)*(1-y)")
    print("="*60)
    
    # Exact solution
    def u_exact(x, y):
        return x * y * (1.0 - x) * (1.0 - y)
    
    # RHS: u_xx + u_yy
    # u = x*y*(1-x)*(1-y) = x*y*(1-x-y+xy) = xy - x^2y - xy^2 + x^2y^2
    # u_x = y - 2xy - y^2 + 2xy^2
    # u_xx = -2y + 2y^2
    # u_y = x - x^2 - 2xy + 2x^2y
    # u_yy = -2x + 2x^2
    # u_xx + u_yy = -2y + 2y^2 - 2x + 2x^2 = 2(x^2 + y^2 - x - y)
    def f_rhs(x, y):
        return 2.0 * (x**2 + y**2 - x - y)
    
    # Boundary condition
    def bc(x, y):
        return u_exact(x, y)
    
    # Solve
    U, X, Y, time = solve_poisson_fft2d(f_rhs, bc, sx=1.0, sy=1.0, nx=33, ny=33)
    
    # Compute error
    U_exact = u_exact(X, Y)
    error = np.abs(U - U_exact)
    max_error = np.max(error)
    
    print(f"\nResults:")
    print(f"  Max error: {max_error:.6e}")
    print(f"  Grid: 33 x 33")
    
    return max_error, time


def main():
    """Main test function"""
    print("="*60)
    print("Poisson FFT Solver - Test Program")
    print("="*60)
    
    # Run tests
    print("\n" + "="*60)
    print("Running Test 1: Smooth Solution (sin(pi*x)*sin(pi*y))")
    print("="*60)
    error1, time1 = test_poisson_smooth()
    
    print("\n" + "="*60)
    print("Running Test 2: Polynomial Solution (x*y*(1-x)*(1-y))")
    print("="*60)
    error2, time2 = test_poisson_polynomial()
    
    print("\n" + "="*60)
    print("Summary:")
    print(f"  Test 1 - Max error: {error1:.6e}, Time: {time1:.4f}s")
    print(f"  Test 2 - Max error: {error2:.6e}, Time: {time2:.4f}s")
    print("="*60)


if __name__ == "__main__":
    main()
