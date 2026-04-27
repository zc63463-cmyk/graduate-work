#!/usr/bin/env python3
"""
FFT9 Algorithm - CORRECT Implementation

The 4th-order 9-point compact scheme for Laplacian(u) = g is:
  L_h * u = R_h * g

where:
  L_h = (1/(6h^2)) * [1   4   1;
                       4  -20   4;
                       1   4   1]

  R_h = [0   -1/12   0;
         -1/12  4/3  -1/12;
          0   -1/12   0]

KEY FIX: R_h has NEGATIVE off-center coefficients (-1/12), not positive!
This is because R_h = I + (h^2/12)*Laplacian, and the Laplacian stencil
has negative off-center coefficients.

For -Laplacian(u) = f:
  Laplacian(u) = -f, so g = -f
  L_h * u = R_h * (-f)

Eigenvalues of L_h (for DST diagonalization):
  lambda_L(k,l) = (1/(6h^2)) * [-20 + 8*cos(k*theta) + 8*cos(l*phi) + 4*cos(k*theta)*cos(l*phi)]

Eigenvalues of R_h:
  lambda_R(k,l) = 4/3 - (1/6)*(cos(k*theta) + cos(l*phi))

Solve in frequency domain:
  u_hat(k,l) = lambda_R(k,l) / lambda_L(k,l) * g_hat(k,l)
"""

import numpy as np
from scipy.fft import dst, idst
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time


def fft9_solve(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    FFT9 solver for Poisson equation -Laplacian(u) = f
    with Dirichlet boundary condition u = g on boundary.
    
    Uses 4th-order 9-point compact finite difference scheme.
    """
    # Grid
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    
    h = hx
    assert abs(hx - hy) < 1e-12, "FFT9 requires hx = hy"
    
    N = n - 2  # interior points per direction
    
    # Evaluate f at all grid points
    F = f_func(X, Y)
    
    # For -Laplacian(u) = f: Laplacian(u) = -f, so g = -f
    G = -F  # g = Laplacian(u) = -f
    
    # Apply R_h to g at interior points
    # R_h = [0 -1/12 0; -1/12 4/3 -1/12; 0 -1/12 0]
    Rg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ii = i + 1
            jj = j + 1
            Rg[i, j] = (4.0/3.0) * G[ii, jj] \
                     - (1.0/12.0) * (G[ii-1, jj] + G[ii+1, jj]
                                    + G[ii, jj-1] + G[ii, jj+1])
    
    # BC correction
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)
    
    # L_h at boundary-adjacent points involves BC values
    # Need to subtract their contribution from Rg
    bc_corr = np.zeros((N, N))
    
    # Left boundary: u(0,y) contributes to L_h u(1,y)
    for j in range(N):
        jj = j + 1
        val = 4.0 * bc_left[jj]
        if jj - 1 >= 0:
            val += 1.0 * bc_left[jj-1]
        if jj + 1 < n:
            val += 1.0 * bc_left[jj+1]
        bc_corr[0, j] -= (1.0 / (6.0 * h**2)) * val
    
    # Right boundary
    for j in range(N):
        jj = j + 1
        val = 4.0 * bc_right[jj]
        if jj - 1 >= 0:
            val += 1.0 * bc_right[jj-1]
        if jj + 1 < n:
            val += 1.0 * bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0 / (6.0 * h**2)) * val
    
    # Bottom boundary
    for i in range(N):
        ii = i + 1
        val = 4.0 * bc_bottom[ii]
        if ii - 1 >= 0:
            val += 1.0 * bc_bottom[ii-1]
        if ii + 1 < n:
            val += 1.0 * bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0 / (6.0 * h**2)) * val
    
    # Top boundary
    for i in range(N):
        ii = i + 1
        val = 4.0 * bc_top[ii]
        if ii - 1 >= 0:
            val += 1.0 * bc_top[ii-1]
        if ii + 1 < n:
            val += 1.0 * bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0 / (6.0 * h**2)) * val
    
    Rg += bc_corr
    
    # 2D DST
    Rg_hat = np.zeros((N, N))
    for i in range(N):
        Rg_hat[i, :] = dst(Rg[i, :], type=1, norm='ortho')
    for j in range(N):
        Rg_hat[:, j] = dst(Rg_hat[:, j], type=1, norm='ortho')
    
    # Solve in frequency domain
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))
    
    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            lambda_L = (1.0/(6.0*h**2)) * (
                -20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li] + 4.0*cos_k[ki]*cos_k[li]
            )
            if abs(lambda_L) > 1e-14:
                U_hat[ki, li] = Rg_hat[ki, li] / lambda_L
    
    # Inverse 2D DST
    U_int = np.zeros((N, N))
    for ki in range(N):
        U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
    for li in range(N):
        U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
    
    # Assemble full solution
    U = np.zeros((n, n))
    U[0, :] = bc_left
    U[-1, :] = bc_right
    U[:, 0] = bc_bottom
    U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    
    return U


def poisson_5point(n, f_func, bc_func, sx=1.0, sy=1.0):
    """Standard 5-point FD + FFT solver for -Laplacian(u) = f (2nd order)"""
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    N = n - 2
    
    F = f_func(X, Y)
    
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)
    
    F_adj = F[1:-1, 1:-1].copy()
    F_adj[0, :] -= bc_left[1:-1] / hx**2
    F_adj[-1, :] -= bc_right[1:-1] / hx**2
    F_adj[:, 0] -= bc_bottom[1:-1] / hy**2
    F_adj[:, -1] -= bc_top[1:-1] / hy**2
    
    F_hat = np.zeros((N, N))
    for i in range(N):
        F_hat[i, :] = dst(F_adj[i, :], type=1, norm='ortho')
    for j in range(N):
        F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
    
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))
    
    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            lam = (2.0/hx**2)*(cos_k[ki]-1.0) + (2.0/hy**2)*(cos_k[li]-1.0)
            if abs(lam) > 1e-14:
                U_hat[ki, li] = -F_hat[ki, li] / lam
    
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
    """Test FFT9 convergence with zero BC: u = sin(pi*x)*sin(pi*y)"""
    print("=" * 70)
    print("FFT9 Convergence Test (CORRECT R_h with -1/12 coefficients)")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def f_rhs(x, y):
        return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def bc(x, y):
        return 0.0
    
    ns = [5, 9, 17, 33, 65]
    
    print("\n--- 5-point (2nd order) ---")
    errors_5 = []
    for n in ns:
        U = poisson_5point(n, f_rhs, bc)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y)))
        errors_5.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\n--- FFT9 9-point (4th order, correct R_h) ---")
    errors_9 = []
    for n in ns:
        U = fft9_solve(n, f_rhs, bc)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y)))
        errors_9.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    # Convergence rates
    print("\n" + "=" * 70)
    print("Convergence Rates")
    print("=" * 70)
    print(f"  {'n':>3s}  {'5pt error':>12s}  {'rate':>5s}  {'9pt error':>12s}  {'rate':>5s}")
    print("  " + "-" * 50)
    for i in range(len(ns)):
        if i == 0:
            print(f"  {ns[i]:3d}  {errors_5[i]:12.6e}   ---  {errors_9[i]:12.6e}   ---")
        else:
            r5 = np.log(errors_5[i-1]/errors_5[i]) / np.log(ns[i]/ns[i-1])
            r9 = np.log(errors_9[i-1]/errors_9[i]) / np.log(ns[i]/ns[i-1])
            print(f"  {ns[i]:3d}  {errors_5[i]:12.6e}  {r5:5.2f}  {errors_9[i]:12.6e}  {r9:5.2f}")


def test_paper_problem():
    """Test against paper Table I: u = 3*x*y*(x-x^2)*(y-y^2)"""
    print("\n" + "=" * 70)
    print("Paper Table I Comparison")
    print("Exact: u = 3*x*y*(x-x^2)*(y-y^2)")
    print("=" * 70)
    
    def u_exact(x, y):
        return 3.0 * x * y * (x - x**2) * (y - y**2)
    
    def f_rhs(x, y):
        return -3.0 * (2.0 - 6.0*x) * y**2 * (1.0 - y) - 3.0 * x**2 * (1.0 - x) * (2.0 - 6.0*y)
    
    def bc(x, y):
        return u_exact(x, y)
    
    ns = [5, 9, 17, 33, 65]
    
    print("\n--- FFT9 9-point (4th order, correct R_h) ---")
    errors = []
    for n in ns:
        U = fft9_solve(n, f_rhs, bc)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y)))
        errors.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\nPaper results (6th order):")
    print("  N=4:  1.34E-07")
    print("  N=8:  2.10E-09")
    print("  N=16: 3.30E-11")
    print("  N=32: 9.49E-13")
    
    print("\nConvergence rates:")
    for i in range(1, len(ns)):
        rate = np.log(errors[i-1]/errors[i]) / np.log(ns[i]/ns[i-1])
        print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_convergence()
    test_paper_problem()
