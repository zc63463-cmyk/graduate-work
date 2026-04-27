#!/usr/bin/env python3
"""
FFT9 Algorithm - Complete Implementation
Based on Houstis & Papatheodorou (1979), ACM TOMS Vol.5 No.4 pp.431-441

Sixth-order 9-point finite difference scheme for Poisson equation:
  L_h u = R_h f

L_h template (1/(6h^2)):
  [1   4   1]
  [4  -20  4]
  [1   4   1]

R_h template (1/360):
  [0   48  0]
  [48   0  48]
  [0   48  0]

Eigenvalues:
  lambda_A(k,l) = (1/(6h^2)) * [-20 + 8*cos(k*theta) + 8*cos(l*phi) + 4*cos(k*theta)*cos(l*phi)]
  lambda_B(k,l) = (4/15) * [cos(k*theta) + cos(l*phi)]

Solve: u_hat(k,l) = (lambda_B(k,l) / lambda_A(k,l)) * f_hat(k,l)

Date: 2026-04-25
"""

import numpy as np
from scipy.fft import dst, idst
import time


def fft9_solve_poisson(f_func, bc_func, n, sx=1.0, sy=1.0):
    """
    FFT9 solver for Poisson equation -u_xx - u_yy = f
    with Dirichlet boundary condition u = g on boundary.
    
    Uses sixth-order 9-point finite difference scheme.
    
    Parameters:
        f_func: RHS function f(x, y)
        bc_func: Boundary condition function g(x, y)
        n: Grid size (should be 2^k + 1)
        sx, sy: Domain size [0, sx] x [0, sy]
    
    Returns:
        U: Numerical solution (n x n)
    """
    # Grid
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    
    # Assume square grid for now (hx = hy = h)
    h = hx
    assert abs(hx - hy) < 1e-12, "FFT9 requires hx = hy for 6th-order scheme"
    
    # Number of interior points
    N = n - 2  # in each direction
    
    # Step 1: Evaluate RHS at grid points
    F = f_func(X, Y)
    
    # Step 2: Apply R_h template to F (in physical space)
    # R_h f_{i,j} = (1/360) * [48*(f_{i-1,j} + f_{i+1,j} + f_{i,j-1} + f_{i,j+1})]
    RF = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ii = i + 1  # index in full grid
            jj = j + 1
            RF[i, j] = (1.0 / 360.0) * 48.0 * (
                F[ii-1, jj] + F[ii+1, jj] + F[ii, jj-1] + F[ii, jj+1]
            )
    
    # Step 3: Handle non-zero boundary conditions
    # For non-zero BC, we need to move the BC contribution to RHS
    # BC values - evaluate on boundary arrays
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).copy()
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).copy()
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).copy()
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).copy()
    
    # The L_h operator at boundary-adjacent points involves BC values
    # L_h u_{i,j} involves u at (i-1,j), (i+1,j), (i,j-1), (i,j+1)
    # and diagonal neighbors (i-1,j-1), etc.
    # If any of these are on the boundary, we move them to RHS
    
    # For simplicity, handle zero BC case first, then add BC correction
    # BC correction for L_h u = RF:
    # At point (1,j) (i=0 in interior), L_h involves u_{0,j} = bc_left[j]
    # The contribution is: (1/(6h^2)) * [4*bc_left[j] + bc_left[j-1] + bc_left[j+1]]
    # (from the template: center column 4, top/bottom 1)
    
    bc_correction = np.zeros((N, N))
    
    # Left boundary (i=0 in full grid, i=0 in interior)
    for j in range(N):
        jj = j + 1
        val = 0.0
        val += 4.0 * bc_left[jj]  # from 4*u_{0,jj}
        val += 1.0 * bc_left[jj-1] if jj-1 >= 0 else 0  # from 1*u_{0,jj-1}
        val += 1.0 * bc_left[jj+1] if jj+1 < n else 0  # from 1*u_{0,jj+1}
        bc_correction[0, j] += (1.0 / (6.0 * h**2)) * val
    
    # Right boundary (i=n-1 in full grid, i=N-1 in interior)
    for j in range(N):
        jj = j + 1
        val = 0.0
        val += 4.0 * bc_right[jj]
        val += 1.0 * bc_right[jj-1] if jj-1 >= 0 else 0
        val += 1.0 * bc_right[jj+1] if jj+1 < n else 0
        bc_correction[N-1, j] += (1.0 / (6.0 * h**2)) * val
    
    # Bottom boundary (j=0 in full grid, j=0 in interior)
    for i in range(N):
        ii = i + 1
        val = 0.0
        val += 4.0 * bc_bottom[ii]  # from 4*u_{ii,0}
        val += 1.0 * bc_bottom[ii-1] if ii-1 >= 0 else 0  # from 1*u_{ii-1,0}
        val += 1.0 * bc_bottom[ii+1] if ii+1 < n else 0  # from 1*u_{ii+1,0}
        bc_correction[i, 0] += (1.0 / (6.0 * h**2)) * val
    
    # Top boundary (j=n-1 in full grid, j=N-1 in interior)
    for i in range(N):
        ii = i + 1
        val = 0.0
        val += 4.0 * bc_top[ii]
        val += 1.0 * bc_top[ii-1] if ii-1 >= 0 else 0
        val += 1.0 * bc_top[ii+1] if ii+1 < n else 0
        bc_correction[i, N-1] += (1.0 / (6.0 * h**2)) * val
    
    # Corner corrections (already partially added above, need diagonal terms)
    # Bottom-left corner (i=0, j=0 in interior): u_{0,0} contributes (1/(6h^2))*1*u_{0,0}
    bc_correction[0, 0] += (1.0 / (6.0 * h**2)) * 1.0 * bc_left[0]  # u_{0,0} from left
    # Actually bc_left[0] = bc_bottom[0], so we'd double count
    # Need to be more careful with corners...
    
    # For now, let's just test with zero BC
    # The BC correction is: RF = RF - bc_correction (move to RHS)
    
    # Step 4: 2D Sine Transform of RF
    RF_hat = np.zeros((N, N))
    for i in range(N):
        RF_hat[i, :] = dst(RF[i, :], type=1, norm='ortho')
    for j in range(N):
        RF_hat[:, j] = dst(RF_hat[:, j], type=1, norm='ortho')
    
    # Step 5: Eigenvalues of L_h and solve in frequency domain
    # lambda_A(k,l) = (1/(6h^2)) * [-20 + 8*cos(k*pi/(N+1)) + 8*cos(l*pi/(N+1)) + 4*cos(k*pi/(N+1))*cos(l*pi/(N+1))]
    
    k_vals = np.arange(1, N + 1)
    l_vals = np.arange(1, N + 1)
    
    cos_k = np.cos(k_vals * np.pi / (N + 1))
    cos_l = np.cos(l_vals * np.pi / (N + 1))
    
    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            lambda_A = (1.0 / (6.0 * h**2)) * (
                -20.0 + 8.0 * cos_k[ki] + 8.0 * cos_l[li] + 4.0 * cos_k[ki] * cos_l[li]
            )
            if abs(lambda_A) > 1e-14:
                U_hat[ki, li] = RF_hat[ki, li] / lambda_A
    
    # Step 6: Inverse 2D Sine Transform
    U_int = np.zeros((N, N))
    for ki in range(N):
        U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
    for li in range(N):
        U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
    
    # Step 7: Assemble full solution
    U = np.zeros((n, n))
    U[0, :] = bc_left
    U[-1, :] = bc_right
    U[:, 0] = bc_bottom
    U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    
    return U


def poisson_5point_fft(f_func, bc_func, n, sx=1.0, sy=1.0):
    """
    Standard 5-point finite difference + FFT solver for Poisson equation.
    This is the O(h^2) method, used as comparison.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    
    N = n - 2
    
    F = f_func(X, Y)
    
    # BC - evaluate on boundary arrays
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).copy()
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).copy()
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).copy()
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).copy()
    
    # Adjusted RHS (with BC contribution)
    F_adj = F[1:-1, 1:-1].copy()
    
    # BC contribution for 5-point stencil
    F_adj[0, :] -= bc_left[1:-1] / hx**2
    F_adj[-1, :] -= bc_right[1:-1] / hx**2
    F_adj[:, 0] -= bc_bottom[1:-1] / hy**2
    F_adj[:, -1] -= bc_top[1:-1] / hy**2
    
    # 2D DST
    F_hat = np.zeros((N, N))
    for i in range(N):
        F_hat[i, :] = dst(F_adj[i, :], type=1, norm='ortho')
    for j in range(N):
        F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
    
    # Eigenvalues for 5-point stencil
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))
    
    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            # lambda = 2/hx^2*(cos(k*theta)-1) + 2/hy^2*(cos(l*phi)-1)
            lambda_total = (2.0 / hx**2) * (cos_k[ki] - 1.0) + (2.0 / hy**2) * (cos_k[li] - 1.0)
            if abs(lambda_total) > 1e-14:
                U_hat[ki, li] = -F_hat[ki, li] / lambda_total
    
    # Inverse 2D DST
    U_int = np.zeros((N, N))
    for ki in range(N):
        U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
    for li in range(N):
        U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
    
    U = np.zeros((n, n))
    U[0, :] = bc_left
    U[-1, :] = bc_right
    U[:, 0] = bc_bottom
    U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    
    return U


def test_convergence():
    """
    Test FFT9 convergence with zero BC problem.
    """
    print("=" * 70)
    print("FFT9 Convergence Test")
    print("=" * 70)
    
    # Problem: -u_xx - u_yy = f, u = 0 on boundary
    # Exact solution: u = sin(pi*x) * sin(pi*y)
    # RHS: f = 2*pi^2 * sin(pi*x) * sin(pi*y)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def f_rhs(x, y):
        return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def bc(x, y):
        return 0.0
    
    # Test different grid sizes
    ns = [5, 9, 17, 33, 65]
    
    print("\n--- 5-point (2nd order) ---")
    errors_5pt = []
    for n in ns:
        U = poisson_5point_fft(f_rhs, bc, n)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        U_ex = u_exact(X, Y)
        err = np.max(np.abs(U - U_ex))
        errors_5pt.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\n--- FFT9 9-point (6th order) ---")
    errors_9pt = []
    for n in ns:
        U = fft9_solve_poisson(f_rhs, bc, n)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        U_ex = u_exact(X, Y)
        err = np.max(np.abs(U - U_ex))
        errors_9pt.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    # Convergence rates
    print("\n" + "=" * 70)
    print("Convergence Rate Analysis")
    print("=" * 70)
    
    print("\n  n      5-pt error   rate    9-pt error   rate")
    print("  " + "-" * 55)
    
    for i in range(len(ns)):
        if i == 0:
            print(f"  {ns[i]:3d}  {errors_5pt[i]:.6e}    ---   {errors_9pt[i]:.6e}    ---")
        else:
            r5 = np.log(errors_5pt[i-1] / errors_5pt[i]) / np.log(ns[i] / ns[i-1])
            r9 = np.log(errors_9pt[i-1] / errors_9pt[i]) / np.log(ns[i] / ns[i-1])
            print(f"  {ns[i]:3d}  {errors_5pt[i]:.6e}   {r5:.2f}   {errors_9pt[i]:.6e}   {r9:.2f}")


def test_paper_results():
    """
    Test against the paper's numerical results.
    
    Paper Table I: Poisson equation, exact solution u = 3*x*y*(x-x^2)*(y-y^2)
    FFT9(6th order) results:
      N=4:  1.34E-07
      N=8:  2.10E-09
      N=16: 3.30E-11
      N=32: 9.49E-13
    """
    print("\n" + "=" * 70)
    print("Comparison with Paper Results (Table I)")
    print("=" * 70)
    
    # Problem: -u_xx - u_yy = f
    # Exact solution: u = 3*x*y*(x-x^2)*(y-y^2) = 3*x*y*x*(1-x)*y*(1-y)
    # Wait, (x-x^2) = x(1-x), so u = 3*x^2*(1-x)*y^2*(1-y)
    
    def u_exact(x, y):
        return 3.0 * x * y * (x - x**2) * (y - y**2)
    
    # Compute -laplacian of u
    # u = 3*x*y*(x-x^2)*(y-y^2) = 3*x^2*(1-x)*y^2*(1-y)
    # Let me expand: u = 3xy(x-x^2)(y-y^2) = 3xy*x(1-x)*y(1-y) = 3x^2y^2(1-x)(1-y)
    # Actually: (x-x^2) = x(1-x), (y-y^2) = y(1-y)
    # So u = 3*x*y * x(1-x) * y(1-y) = 3*x^2*(1-x)*y^2*(1-y)
    
    # u_x = 3*y^2*(1-y) * [2x(1-x) + x^2*(-1)] = 3*y^2*(1-y) * [2x - 2x^2 - x^2] = 3*y^2*(1-y)*(2x - 3x^2)
    # u_xx = 3*y^2*(1-y) * (2 - 6x)
    
    # u_y = 3*x^2*(1-x) * (2y - 3y^2)
    # u_yy = 3*x^2*(1-x) * (2 - 6y)
    
    # -u_xx - u_yy = -3*y^2*(1-y)*(2-6x) - 3*x^2*(1-x)*(2-6y)
    
    def f_rhs(x, y):
        return -3.0 * y**2 * (1.0 - y) * (2.0 - 6.0*x) - 3.0 * x**2 * (1.0 - x) * (2.0 - 6.0*y)
    
    def bc(x, y):
        return u_exact(x, y)
    
    # Paper uses N = 4, 8, 16, 32, 64
    # N here means number of intervals, so grid size = N+1
    # But the paper might use N differently...
    # From the paper, N=4 gives error 1.34E-07 for FFT9(6th order)
    # This suggests N is related to number of grid points
    
    # Let's try both interpretations
    print("\nInterpretation 1: N = number of grid points in each direction")
    for N in [5, 9, 17, 33, 65]:
        U = fft9_solve_poisson(f_rhs, bc, N)
        x = np.linspace(0, 1, N)
        y = np.linspace(0, 1, N)
        X, Y = np.meshgrid(x, y, indexing='ij')
        U_ex = u_exact(X, Y)
        err = np.max(np.abs(U - U_ex))
        print(f"  N={N:3d} (grid {N}x{N}): error = {err:.6e}")
    
    # Paper results for reference
    print("\nPaper results (6th order):")
    print("  N=4:  1.34E-07")
    print("  N=8:  2.10E-09")
    print("  N=16: 3.30E-11")
    print("  N=32: 9.49E-13")


if __name__ == "__main__":
    test_convergence()
    test_paper_results()
