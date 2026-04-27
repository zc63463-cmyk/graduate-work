#!/usr/bin/env python3
"""
FFT9 Algorithm - TRULY CORRECT Implementation

The 4th-order 9-point compact scheme for Laplacian(u) = g is:
  L_h * u = R_h * g

L_h = (1/(6h^2)) * [1   4   1;
                       4  -20   4;
                       1   4   1]

R_h = I + (h^2/12) * Lap_h
    = I + (1/12) * [0  1  0; 1 -4 1; 0  1  0]
    = [0   1/12   0;
       1/12 2/3  1/12;
       0   1/12   0]

KEY: R_h has POSITIVE off-center coefficients (+1/12) and center = 2/3!
Sum of coefficients = 4*(1/12) + 2/3 = 1/3 + 2/3 = 1 ✓

For -Laplacian(u) = f:
  g = Laplacian(u) = -f
  L_h * u = R_h * (-f)
"""

import numpy as np
from scipy.fft import dst, idst
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time


def fft9_solve(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    FFT9 solver for -Laplacian(u) = f with Dirichlet BC.
    Uses 4th-order 9-point compact finite difference scheme.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    h = hx
    N = n - 2
    
    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f
    
    # Apply R_h to g: R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0]
    Rg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ii, jj = i + 1, j + 1
            Rg[i, j] = (2.0/3.0) * G[ii, jj] \
                     + (1.0/12.0) * (G[ii-1, jj] + G[ii+1, jj]
                                    + G[ii, jj-1] + G[ii, jj+1])
    
    # BC correction for L_h
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)
    
    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        # Left boundary
        val = 4.0 * bc_left[jj]
        if jj-1 >= 0: val += 1.0 * bc_left[jj-1]
        if jj+1 < n: val += 1.0 * bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val
        # Right boundary
        val = 4.0 * bc_right[jj]
        if jj-1 >= 0: val += 1.0 * bc_right[jj-1]
        if jj+1 < n: val += 1.0 * bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val
    for i in range(N):
        ii = i + 1
        # Bottom boundary
        val = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val += 1.0 * bc_bottom[ii-1]
        if ii+1 < n: val += 1.0 * bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val
        # Top boundary
        val = 4.0 * bc_top[ii]
        if ii-1 >= 0: val += 1.0 * bc_top[ii-1]
        if ii+1 < n: val += 1.0 * bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val
    
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
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li] + 4.0*cos_k[ki]*cos_k[li])
            if abs(lam_L) > 1e-14:
                U_hat[ki, li] = Rg_hat[ki, li] / lam_L
    
    # Inverse 2D DST
    U_int = np.zeros((N, N))
    for ki in range(N):
        U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
    for li in range(N):
        U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
    
    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    return U


def poisson_5point(n, f_func, bc_func, sx=1.0, sy=1.0):
    """5-point FD + FFT for -Laplacian(u) = f (2nd order)"""
    x = np.linspace(0, sx, n); y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx/(n-1); hy = sy/(n-1); N = n-2
    F = f_func(X, Y)
    bc_l = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_r = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_b = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_t = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)
    F_adj = F[1:-1,1:-1].copy()
    F_adj[0,:] -= bc_l[1:-1]/hx**2; F_adj[-1,:] -= bc_r[1:-1]/hx**2
    F_adj[:,0] -= bc_b[1:-1]/hy**2; F_adj[:,-1] -= bc_t[1:-1]/hy**2
    F_hat = np.zeros((N,N))
    for i in range(N): F_hat[i,:] = dst(F_adj[i,:], type=1, norm='ortho')
    for j in range(N): F_hat[:,j] = dst(F_hat[:,j], type=1, norm='ortho')
    k_v = np.arange(1,N+1); ck = np.cos(k_v*np.pi/(N+1))
    U_hat = np.zeros((N,N))
    for ki in range(N):
        for li in range(N):
            lam = (2.0/hx**2)*(ck[ki]-1.0) + (2.0/hy**2)*(ck[li]-1.0)
            if abs(lam)>1e-14: U_hat[ki,li] = -F_hat[ki,li]/lam
    U_int = np.zeros((N,N))
    for ki in range(N): U_int[ki,:] = idst(U_hat[ki,:], type=1, norm='ortho')
    for li in range(N): U_int[:,li] = idst(U_int[:,li], type=1, norm='ortho')
    U = np.zeros((n,n))
    U[0,:]=bc_l; U[-1,:]=bc_r; U[:,0]=bc_b; U[:,-1]=bc_t
    U[1:-1,1:-1] = U_int
    return U


def test_convergence():
    """Test FFT9 convergence with zero BC"""
    print("=" * 70)
    print("FFT9 Convergence Test (R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0])")
    print("=" * 70)
    
    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc(x, y): return 0.0
    
    ns = [5, 9, 17, 33, 65, 129]
    
    print("\n--- 5-point (2nd order) ---")
    e5 = []
    for n in ns:
        U = poisson_5point(n, f_rhs, bc)
        x = np.linspace(0,1,n); X,Y = np.meshgrid(x,x,indexing='ij')
        err = np.max(np.abs(U - u_exact(X,Y))); e5.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\n--- FFT9 9-point (4th order) ---")
    e9 = []
    for n in ns:
        U = fft9_solve(n, f_rhs, bc)
        x = np.linspace(0,1,n); X,Y = np.meshgrid(x,x,indexing='ij')
        err = np.max(np.abs(U - u_exact(X,Y))); e9.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\n--- Convergence Rates ---")
    print(f"  {'n':>3s}  {'5pt':>12s}  {'rate':>5s}  {'9pt':>12s}  {'rate':>5s}")
    for i in range(len(ns)):
        if i == 0:
            print(f"  {ns[i]:3d}  {e5[i]:12.6e}   ---  {e9[i]:12.6e}   ---")
        else:
            r5 = np.log(e5[i-1]/e5[i])/np.log(ns[i]/ns[i-1])
            r9 = np.log(e9[i-1]/e9[i])/np.log(ns[i]/ns[i-1])
            print(f"  {ns[i]:3d}  {e5[i]:12.6e}  {r5:5.2f}  {e9[i]:12.6e}  {r9:5.2f}")


def test_paper_problem():
    """Test against paper Table I: u = 3*x*y*(x-x^2)*(y-y^2)"""
    print("\n" + "=" * 70)
    print("Paper Table I Comparison")
    print("=" * 70)
    
    def u_exact(x, y): return 3.0*x*y*(x-x**2)*(y-y**2)
    def f_rhs(x, y):
        return -3.0*(2.0-6.0*x)*y**2*(1.0-y) - 3.0*x**2*(1.0-x)*(2.0-6.0*y)
    def bc(x, y): return u_exact(x, y)
    
    ns = [5, 9, 17, 33, 65]
    
    print("\n--- FFT9 9-point (4th order) ---")
    e9 = []
    for n in ns:
        U = fft9_solve(n, f_rhs, bc)
        x = np.linspace(0,1,n); X,Y = np.meshgrid(x,x,indexing='ij')
        err = np.max(np.abs(U - u_exact(X,Y))); e9.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    print("\nPaper results (6th order):")
    print("  N=4: 1.34E-07, N=8: 2.10E-09, N=16: 3.30E-11, N=32: 9.49E-13")
    
    print("\nConvergence rates:")
    for i in range(1, len(ns)):
        rate = np.log(e9[i-1]/e9[i])/np.log(ns[i]/ns[i-1])
        print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_convergence()
    test_paper_problem()
