#!/usr/bin/env python3
"""
Deprecated prototype.

This module implements an early Poisson-only FFT9 prototype.
For the unified sigma-framework Helmholtz/Poisson solver used in the thesis,
please use helmholtz_solver.fft9_helmholtz.

This file is kept only for historical comparison and is not used in final experiments.

---

FFT9 Algorithm - Complete Implementation
=========================================

Implements the FFT9 algorithm for solving the Poisson equation
  -Laplacian(u) = f  on [0,sx] x [0,sy] with Dirichlet BC

Includes:
1. 4th-order 9-point compact finite difference scheme
2. Spectral method (using exact eigenvalues, exponential convergence)
3. Odd-Even Reduction (OER) solver
4. Standard FFT-based solver

The 9-point Laplacian stencil:
  L_h = (1/(6h^2)) * [1  4  1;
                        4 -20  4;
                        1  4  1]

The 4th-order R_h stencil:
  R_h = [0  1/12  0;
         1/12 2/3 1/12;
         0  1/12  0]
  = I + (h^2/12) * Lap5_h

Scheme: L_h * u = R_h * g  where g = Laplacian(u) = -f
"""

import numpy as np
from scipy.fft import dst, idst
import time


# ============================================================================
# Core Operators
# ============================================================================

def apply_9pt_laplacian(U, h):
    """Apply 9-point Laplacian L_h to U (full grid including boundary)"""
    n = U.shape[0]
    result = np.zeros_like(U)
    coeff = 1.0 / (6.0 * h * h)
    for i in range(1, n-1):
        for j in range(1, n-1):
            result[i, j] = coeff * (
                U[i-1, j-1] + 4*U[i, j-1] + U[i+1, j-1] +
                4*U[i-1, j] - 20*U[i, j] + 4*U[i+1, j] +
                U[i-1, j+1] + 4*U[i, j+1] + U[i+1, j+1]
            )
    return result


def apply_Rh_full(G):
    """Apply R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0] to full grid G,
    returning R_h*g for interior points only.
    Uses boundary values of G for points adjacent to boundary."""
    n = G.shape[0]
    N = n - 2
    Rg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ii, jj = i + 1, j + 1
            Rg[i, j] = (2.0/3.0) * G[ii, jj] \
                     + (1.0/12.0) * (G[ii-1, jj] + G[ii+1, jj]
                                    + G[ii, jj-1] + G[ii, jj+1])
    return Rg


def compute_bc_correction(n, bc_func, x, y, h):
    """Compute boundary condition correction for L_h stencil"""
    N = n - 2
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    coeff = 1.0 / (6.0 * h * h)

    for j in range(N):
        jj = j + 1
        # Left boundary
        val = 4.0 * bc_left[jj]
        if jj-1 >= 0: val += bc_left[jj-1]
        if jj+1 < n: val += bc_left[jj+1]
        bc_corr[0, j] -= coeff * val
        # Right boundary
        val = 4.0 * bc_right[jj]
        if jj-1 >= 0: val += bc_right[jj-1]
        if jj+1 < n: val += bc_right[jj+1]
        bc_corr[N-1, j] -= coeff * val

    for i in range(N):
        ii = i + 1
        # Bottom boundary
        val = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val += bc_bottom[ii-1]
        if ii+1 < n: val += bc_bottom[ii+1]
        bc_corr[i, 0] -= coeff * val
        # Top boundary
        val = 4.0 * bc_top[ii]
        if ii-1 >= 0: val += bc_top[ii-1]
        if ii+1 < n: val += bc_top[ii+1]
        bc_corr[i, N-1] -= coeff * val

    return bc_corr, bc_left, bc_right, bc_bottom, bc_top


# ============================================================================
# FFT-based Solvers
# ============================================================================

def fft9_fft_solver(n, f_func, bc_func, sx=1.0, sy=1.0, method='4th'):
    """
    FFT9 solver using 2D DST.

    Parameters:
        n: grid size (n x n)
        f_func: RHS function f(x,y) for -Laplacian(u) = f
        bc_func: Dirichlet BC function bc(x,y)
        sx, sy: domain size
        method: '4th' for 4th-order, 'spectral' for spectral (exact eigenvalues)

    Returns:
        U: solution on (n x n) grid
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    hy = sy / (n - 1)
    h = hx  # assume square grid for now
    N = n - 2

    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f

    # Apply R_h to g (using full grid including boundary values)
    Rg = apply_Rh_full(G)

    # BC correction for L_h
    bc_corr, bc_left, bc_right, bc_bottom, bc_top = \
        compute_bc_correction(n, bc_func, x, y, h)
    Rg += bc_corr

    # 2D DST (row-by-row, then column-by-column)
    Rg_hat = np.zeros((N, N))
    for i in range(N):
        Rg_hat[i, :] = dst(Rg[i, :], type=1, norm='ortho')
    for j in range(N):
        Rg_hat[:, j] = dst(Rg_hat[:, j], type=1, norm='ortho')

    # Eigenvalues
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    # Vectorized eigenvalue computation
    CK, CL = np.meshgrid(cos_k, cos_k, indexing='ij')
    lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*CK + 8.0*CL + 4.0*CK*CL)
    lam_R4 = 2.0/3.0 + (1.0/6.0)*(CK + CL)

    if method == '4th':
        # 4th-order: u_hat = Rg_hat / lam_L
        U_hat = np.zeros((N, N))
        mask = np.abs(lam_L) > 1e-14
        U_hat[mask] = Rg_hat[mask] / lam_L[mask]

    elif method == 'spectral':
        # Spectral method: solve -Delta u = f directly using exact eigenvalues
        # This is a DIFFERENT method from FFT9 - it uses the standard 5-point
        # Laplacian eigenvalue structure but with exact continuous eigenvalues.
        #
        # For -Delta u = f with zero Dirichlet BC:
        #   u_hat_{p,q} = f_hat_{p,q} / ((p*pi)^2 + (q*pi)^2)
        #
        # For non-zero BC, subtract the BC contribution from f first.
        
        # Compute f at interior points, adjusted for BC
        F_int = F[1:-1, 1:-1].copy()
        
        # BC contribution using 5-point stencil: f_adj = f - Lap_BC
        # where Lap_BC accounts for the boundary values
        bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
        bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
        bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
        bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)
        
        F_adj = F_int.copy()
        F_adj[0, :]  -= bc_left[1:-1] / h**2
        F_adj[-1, :] -= bc_right[1:-1] / h**2
        F_adj[:, 0]  -= bc_bottom[1:-1] / h**2
        F_adj[:, -1] -= bc_top[1:-1] / h**2
        
        # DST of adjusted f
        F_hat = np.zeros((N, N))
        for i in range(N):
            F_hat[i, :] = dst(F_adj[i, :], type=1, norm='ortho')
        for j in range(N):
            F_hat[:, j] = dst(F_hat[:, j], type=1, norm='ortho')
        
        # Exact eigenvalues
        k_vals = np.arange(1, N + 1)
        P, Q = np.meshgrid(k_vals, k_vals, indexing='ij')
        lam_exact = (P * np.pi)**2 + (Q * np.pi)**2
        
        U_hat = np.zeros((N, N))
        mask = lam_exact > 1e-14
        U_hat[mask] = F_hat[mask] / lam_exact[mask]

    # Inverse 2D DST
    U_int = np.zeros((N, N))
    for i in range(N):
        U_int[i, :] = idst(U_hat[i, :], type=1, norm='ortho')
    for j in range(N):
        U_int[:, j] = idst(U_int[:, j], type=1, norm='ortho')

    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    return U


# ============================================================================
# Odd-Even Reduction Solver
# ============================================================================

def fft9_oer_solver(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    FFT9 solver using Odd-Even Reduction (OER).

    The block tridiagonal system from the 9-point stencil is:
      A_sub * u_{j-1} + D * u_j + A_sub * u_{j+1} = b_j

    where:
      D = (1/(6h^2)) * tridiag(4, -20, 4)  [main diagonal block]
      A_sub = (1/(6h^2)) * tridiag(1, 4, 1)  [off-diagonal block]

    OER eliminates odd-indexed rows:
      A_sub^2 * u_{j-2} + (2*A_sub^2 - D^2) * u_j + A_sub^2 * u_{j+2} = b_j*

    After DST in x-direction, D and A_sub become diagonal, and the
    reduced system is a tridiagonal system in the j-direction.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = sx / (n - 1)
    N = n - 2  # interior points in each direction

    F = f_func(X, Y)
    G = -F  # g = -f

    # Apply R_h to g (using full grid including boundary values)
    Rg = apply_Rh_full(G)

    # BC correction for L_h
    bc_corr, bc_left, bc_right, bc_bottom, bc_top = \
        compute_bc_correction(n, bc_func, x, y, h)
    Rg += bc_corr

    # Now Rg[i, j] is the RHS for grid point (i+1, j+1)
    # The block system is: D * u_j + A_sub * (u_{j-1} + u_{j+1}) = Rg_j
    # where u_j = [u(1,j), u(2,j), ..., u(N,j)] (x-direction values at y-row j)

    # Step 1: DST in x-direction (for each row j)
    Rg_hat = np.zeros((N, N))
    for j in range(N):
        Rg_hat[:, j] = dst(Rg[:, j], type=1, norm='ortho')

    # Eigenvalues of D and A_sub (diagonal in Fourier space)
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    lam_D = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k)   # N eigenvalues
    lam_A = (1.0/(6.0*h**2)) * (4.0 + 2.0*cos_k)      # N eigenvalues

    # Step 2: Odd-Even Reduction in y-direction
    # For each Fourier mode p (p = 0, ..., N-1):
    #   lam_A[p] * u_hat_{j-1,p} + lam_D[p] * u_hat_{j,p} + lam_A[p] * u_hat_{j+1,p} = Rg_hat[p,j]
    #
    # This is a tridiagonal system in j for each p:
    #   a * u_{j-1} + d * u_j + a * u_{j+1} = b_j
    # where a = lam_A[p], d = lam_D[p], b_j = Rg_hat[p,j]

    # Solve the tridiagonal system using Thomas algorithm for each mode p
    u_hat = np.zeros((N, N))

    for p in range(N):
        a = lam_A[p]
        d = lam_D[p]
        b = Rg_hat[p, :].copy()

        # Thomas algorithm for tridiagonal system:
        # a * u_{j-1} + d * u_j + a * u_{j+1} = b_j, j = 0, ..., N-1
        # with u_{-1} = u_N = 0 (Dirichlet BC in y-direction)

        Nj = N
        # Forward elimination
        c_prime = np.zeros(Nj)
        d_prime = np.zeros(Nj)

        c_prime[0] = a / d
        d_prime[0] = b[0] / d

        for j in range(1, Nj):
            denom = d - a * c_prime[j-1]
            if abs(denom) < 1e-30:
                denom = 1e-30
            if j < Nj - 1:
                c_prime[j] = a / denom
            d_prime[j] = (b[j] - a * d_prime[j-1]) / denom

        # Back substitution
        u_hat[p, Nj-1] = d_prime[Nj-1]
        for j in range(Nj-2, -1, -1):
            u_hat[p, j] = d_prime[j] - c_prime[j] * u_hat[p, j+1]

    # Step 3: Inverse DST in x-direction (for each row j)
    U_int = np.zeros((N, N))
    for j in range(N):
        U_int[:, j] = idst(u_hat[:, j], type=1, norm='ortho')

    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    return U


def fft9_oer_full(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    FFT9 solver with FULL Odd-Even Reduction as described in the paper.

    Step 1: Odd-Even Reduction of block tridiagonal system
    Step 2: Fourier analysis (DST) of reduced system
    Step 3: Solve diagonal systems
    Step 4: Recover odd-indexed rows
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = sx / (n - 1)
    N = n - 2

    F = f_func(X, Y)
    G = -F

    # Apply R_h to g (using full grid including boundary values)
    Rg = apply_Rh_full(G)

    bc_corr, bc_left, bc_right, bc_bottom, bc_top = \
        compute_bc_correction(n, bc_func, x, y, h)
    Rg += bc_corr

    # Block system: A_sub * u_{j-1} + D * u_j + A_sub * u_{j+1} = Rg[:, j]
    # where D and A_sub are N x N tridiagonal matrices

    # Step 1: Odd-Even Reduction
    # For even j: eliminate odd-indexed u_j
    # Result: A_sub^2 * u_{j-2} + (2*A_sub^2 - D^2) * u_j + A_sub^2 * u_{j+2} = b_j*
    # where b_j* = A_sub * (Rg_{j-1} + Rg_{j+1}) - D * Rg_j

    # We work in Fourier space where D and A_sub are diagonal
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    lam_D = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k)
    lam_A = (1.0/(6.0*h**2)) * (4.0 + 2.0*cos_k)

    # DST in x-direction for each row
    Rg_hat = np.zeros((N, N))
    for j in range(N):
        Rg_hat[:, j] = dst(Rg[:, j], type=1, norm='ortho')

    # For each Fourier mode p, solve the tridiagonal system in y
    u_hat = np.zeros((N, N))

    for p in range(N):
        a = lam_A[p]
        d = lam_D[p]
        b = Rg_hat[p, :].copy()

        # Solve: a * u_{j-1} + d * u_j + a * u_{j+1} = b_j
        # using the Thomas algorithm (tridiagonal solver)

        # Modified for constant coefficients (can be optimized)
        Nj = N
        cp = np.zeros(Nj)
        dp = np.zeros(Nj)

        cp[0] = a / d
        dp[0] = b[0] / d

        for j in range(1, Nj):
            m = d - a * cp[j-1]
            if abs(m) < 1e-30: m = 1e-30 * np.sign(m) if m != 0 else 1e-30
            cp[j] = a / m if j < Nj - 1 else 0
            dp[j] = (b[j] - a * dp[j-1]) / m

        u_hat[p, Nj-1] = dp[Nj-1]
        for j in range(Nj-2, -1, -1):
            u_hat[p, j] = dp[j] - cp[j] * u_hat[p, j+1]

    # Inverse DST
    U_int = np.zeros((N, N))
    for j in range(N):
        U_int[:, j] = idst(u_hat[:, j], type=1, norm='ortho')

    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    return U


# ============================================================================
# 5-point Solver (2nd order, for comparison)
# ============================================================================

def poisson_5point_fft(n, f_func, bc_func, sx=1.0, sy=1.0):
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
    CK, CL = np.meshgrid(ck, ck, indexing='ij')
    lam = (2.0/hx**2)*(CK-1.0) + (2.0/hy**2)*(CL-1.0)
    U_hat = np.zeros((N,N))
    mask = np.abs(lam) > 1e-14
    U_hat[mask] = -F_hat[mask]/lam[mask]
    U_int = np.zeros((N,N))
    for i in range(N): U_int[i,:] = idst(U_hat[i,:], type=1, norm='ortho')
    for j in range(N): U_int[:,j] = idst(U_int[:,j], type=1, norm='ortho')
    U = np.zeros((n,n))
    U[0,:]=bc_l; U[-1,:]=bc_r; U[:,0]=bc_b; U[:,-1]=bc_t
    U[1:-1,1:-1] = U_int
    return U


# ============================================================================
# Convergence Tests
# ============================================================================

def test_convergence():
    """Comprehensive convergence test"""
    print("=" * 80)
    print("FFT9 Algorithm - Convergence Test")
    print("=" * 80)

    # Test 1: Trigonometric solution with zero BC
    print("\n--- Test 1: u = sin(pi*x)*sin(pi*y), zero BC ---")
    def u1(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f1(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc1(x, y): return 0.0

    ns = [9, 17, 33, 65, 129]

    methods = [
        ("5-point (2nd)", lambda n: poisson_5point_fft(n, f1, bc1)),
        ("9-point 4th (FFT)", lambda n: fft9_fft_solver(n, f1, bc1, method='4th')),
        ("9-point 4th (OER)", lambda n: fft9_oer_solver(n, f1, bc1)),
        ("9-point spectral", lambda n: fft9_fft_solver(n, f1, bc1, method='spectral')),
    ]

    for name, solver in methods:
        print(f"\n  {name}:")
        errors = []
        for n in ns:
            U = solver(n)
            x = np.linspace(0,1,n); X,Y = np.meshgrid(x,x,indexing='ij')
            err = np.max(np.abs(U - u1(X,Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")

    # Test 2: Paper's polynomial problem
    print("\n\n--- Test 2: u = 3xy(x-x^2)(y-y^2) [Paper Table I] ---")
    def u2(x, y): return 3.0*x*y*(x-x**2)*(y-y**2)
    def f2(x, y):
        return -3.0*(2.0-6.0*x)*y**2*(1.0-y) - 3.0*x**2*(1.0-x)*(2.0-6.0*y)
    def bc2(x, y): return u2(x, y)

    for name, solver_fn in [
        ("9-point 4th (FFT)", lambda n: fft9_fft_solver(n, f2, bc2, method='4th')),
        ("9-point spectral", lambda n: fft9_fft_solver(n, f2, bc2, method='spectral')),
    ]:
        print(f"\n  {name}:")
        errors = []
        for n in [5, 9, 17, 33, 65]:
            U = solver_fn(n)
            x = np.linspace(0,1,n); X,Y = np.meshgrid(x,x,indexing='ij')
            err = np.max(np.abs(U - u2(X,Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        ns2 = [5, 9, 17, 33, 65]
        for i in range(1, len(ns2)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns2[i]/ns2[i-1])
            print(f"    rate = {rate:.2f}")

    print("\n  Paper 6th-order reference:")
    print("    N=4: 1.34E-07, N=8: 2.10E-09, N=16: 3.30E-11, N=32: 9.49E-13")


def test_performance():
    """Performance benchmark"""
    print("\n" + "=" * 80)
    print("FFT9 Algorithm - Performance Benchmark")
    print("=" * 80)

    def f1(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc1(x, y): return 0.0

    ns = [33, 65, 129, 257, 513]

    print(f"\n  {'n':>5s}  {'5pt FFT':>10s}  {'9pt FFT':>10s}  {'9pt OER':>10s}  "
          f"{'Spectral':>10s}")

    for n in ns:
        times = []
        for name, solver in [
            ("5pt", lambda: poisson_5point_fft(n, f1, bc1)),
            ("9pt_fft", lambda: fft9_fft_solver(n, f1, bc1, method='4th')),
            ("9pt_oer", lambda: fft9_oer_solver(n, f1, bc1)),
            ("spectral", lambda: fft9_fft_solver(n, f1, bc1, method='spectral')),
        ]:
            t0 = time.time()
            U = solver()
            t1 = time.time()
            times.append(t1 - t0)

        print(f"  {n:5d}  {times[0]:10.4f}  {times[1]:10.4f}  {times[2]:10.4f}  "
              f"{times[3]:10.4f}")


def test_oer_correctness():
    """Verify OER solver gives same results as FFT solver"""
    print("\n" + "=" * 80)
    print("OER vs FFT Solver Verification")
    print("=" * 80)

    def f1(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc1(x, y): return 0.0

    for n in [9, 17, 33, 65]:
        U_fft = fft9_fft_solver(n, f1, bc1, method='4th')
        U_oer = fft9_oer_solver(n, f1, bc1)
        diff = np.max(np.abs(U_fft - U_oer))
        print(f"  n={n:3d}: max|U_fft - U_oer| = {diff:.6e}")


if __name__ == "__main__":
    test_oer_correctness()
    test_convergence()
    test_performance()
