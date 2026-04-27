#!/usr/bin/env python3
"""
FFT9 Algorithm - 6th Order Scheme Implementation

The 4th-order compact 9-point scheme for -Laplacian(u) = f is:
  L_h * u = R_h^(4) * g    where g = -f

For 6th order, we add a correction term using the discrete biharmonic of g:
  L_h * u = R_h^(4) * g + alpha * h^4 * Delta_h^2 * g

In Fourier space, this becomes:
  lambda_L * u_hat = (lambda_R + alpha * h^4 * lambda_5^2) * g_hat

where lambda_5 is the eigenvalue of the 5-point Laplacian.

The 5x5 physical stencil for Delta_h^2 is:
  (1/h^4) * [0  0  1  0  0;
              0  2 -8  2  0;
              1 -8 20 -8  1;
              0  2 -8  2  0;
              0  0  1  0  0]
"""

import numpy as np
from scipy.fft import dst, idst
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time


def apply_Rh4(G, N):
    """Apply 4th-order R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0] to interior G"""
    # G is (N, N) interior values only
    Rg = (2.0/3.0) * G.copy()
    # Edge neighbors
    Rg[1:, :]  += (1.0/12.0) * G[:-1, :]
    Rg[:-1, :] += (1.0/12.0) * G[1:, :]
    Rg[:, 1:]  += (1.0/12.0) * G[:, :-1]
    Rg[:, :-1] += (1.0/12.0) * G[:, 1:]
    return Rg


def apply_lap5(F_full, h):
    """Apply 5-point Laplacian to F (full grid including boundary)"""
    n = F_full.shape[0]
    N = n - 2
    lap = np.zeros((N, N))
    lap = (F_full[2:, 1:-1] + F_full[:-2, 1:-1] +
           F_full[1:-1, 2:] + F_full[1:-1, :-2] -
           4.0 * F_full[1:-1, 1:-1]) / h**2
    return lap


def apply_lap5_sq_interior(G, h):
    """Apply 5-point Laplacian squared to interior G (with zero BC)"""
    N = G.shape[0]
    # First application: pad with zeros
    G_pad = np.zeros((N+2, N+2))
    G_pad[1:-1, 1:-1] = G
    lap1 = apply_lap5(G_pad, h)
    # Second application: pad with zeros
    lap1_pad = np.zeros((N+2, N+2))
    lap1_pad[1:-1, 1:-1] = lap1
    lap2 = apply_lap5(lap1_pad, h)
    return lap2


def fft9_4th_order(n, f_func, bc_func, sx=1.0, sy=1.0):
    """4th-order FFT9 solver for -Laplacian(u) = f with Dirichlet BC"""
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    h = hx
    N = n - 2

    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f

    # Apply R_h^(4) to g (interior only)
    G_int = G[1:-1, 1:-1]
    Rg = apply_Rh4(G_int, N)

    # BC correction for L_h
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        val_l = 4.0 * bc_left[jj]
        if jj-1 >= 0: val_l += bc_left[jj-1]
        if jj+1 < n: val_l += bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val_l
        val_r = 4.0 * bc_right[jj]
        if jj-1 >= 0: val_r += bc_right[jj-1]
        if jj+1 < n: val_r += bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val_r
    for i in range(N):
        ii = i + 1
        val_b = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val_b += bc_bottom[ii-1]
        if ii+1 < n: val_b += bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val_b
        val_t = 4.0 * bc_top[ii]
        if ii-1 >= 0: val_t += bc_top[ii-1]
        if ii+1 < n: val_t += bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val_t

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
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li]
                                          + 4.0*cos_k[ki]*cos_k[li])
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


def fft9_6th_order_physical(n, f_func, bc_func, sx=1.0, sy=1.0, alpha=1.0/360.0):
    """
    6th-order FFT9 solver using Delta_h^2 correction in physical space.

    Scheme: L_h * u = R_h^(4) * g + alpha * h^4 * Delta_h^2 * g
    where g = -f for -Laplacian(u) = f.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    h = hx
    N = n - 2

    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f

    # Apply R_h^(4) to g (interior only)
    G_int = G[1:-1, 1:-1]
    Rg = apply_Rh4(G_int, N)

    # Apply Delta_h^2 correction: compute Delta_h(g) then Delta_h(Delta_h(g))
    # g on full grid (interior only, zero BC on boundary of g)
    Delta_h_g = apply_lap5_sq_interior(G_int, h)

    # Add correction
    Rg += alpha * h**4 * Delta_h_g

    # BC correction for L_h
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        val_l = 4.0 * bc_left[jj]
        if jj-1 >= 0: val_l += bc_left[jj-1]
        if jj+1 < n: val_l += bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val_l
        val_r = 4.0 * bc_right[jj]
        if jj-1 >= 0: val_r += bc_right[jj-1]
        if jj+1 < n: val_r += bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val_r
    for i in range(N):
        ii = i + 1
        val_b = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val_b += bc_bottom[ii-1]
        if ii+1 < n: val_b += bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val_b
        val_t = 4.0 * bc_top[ii]
        if ii-1 >= 0: val_t += bc_top[ii-1]
        if ii+1 < n: val_t += bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val_t

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
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li]
                                          + 4.0*cos_k[ki]*cos_k[li])
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


def fft9_6th_order_fourier(n, f_func, bc_func, sx=1.0, sy=1.0, alpha=1.0/360.0):
    """
    6th-order FFT9 solver using correction in Fourier space.

    In Fourier space, the 6th-order scheme is:
      lambda_L * u_hat = (lambda_R + alpha * h^4 * lambda_5^2) * g_hat

    where lambda_5 is the eigenvalue of the 5-point Laplacian.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    h = hx
    N = n - 2

    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f

    # Apply R_h^(4) to g (interior only)
    G_int = G[1:-1, 1:-1]
    Rg = apply_Rh4(G_int, N)

    # BC correction for L_h
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        val_l = 4.0 * bc_left[jj]
        if jj-1 >= 0: val_l += bc_left[jj-1]
        if jj+1 < n: val_l += bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val_l
        val_r = 4.0 * bc_right[jj]
        if jj-1 >= 0: val_r += bc_right[jj-1]
        if jj+1 < n: val_r += bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val_r
    for i in range(N):
        ii = i + 1
        val_b = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val_b += bc_bottom[ii-1]
        if ii+1 < n: val_b += bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val_b
        val_t = 4.0 * bc_top[ii]
        if ii-1 >= 0: val_t += bc_top[ii-1]
        if ii+1 < n: val_t += bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val_t

    Rg += bc_corr

    # 2D DST
    Rg_hat = np.zeros((N, N))
    for i in range(N):
        Rg_hat[i, :] = dst(Rg[i, :], type=1, norm='ortho')
    for j in range(N):
        Rg_hat[:, j] = dst(Rg_hat[:, j], type=1, norm='ortho')

    # Also transform g for the correction term
    G_hat = np.zeros((N, N))
    for i in range(N):
        G_hat[i, :] = dst(G_int[i, :], type=1, norm='ortho')
    for j in range(N):
        G_hat[:, j] = dst(G_hat[:, j], type=1, norm='ortho')

    # Solve in frequency domain with 6th-order correction
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li]
                                          + 4.0*cos_k[ki]*cos_k[li])
            # 4th-order R_h eigenvalue
            lam_R = 2.0/3.0 + (1.0/6.0)*(cos_k[ki] + cos_k[li])
            # 5-point Laplacian eigenvalue
            lam_5 = (2.0/h**2)*(cos_k[ki] + cos_k[li] - 2.0)
            # 6th-order correction
            lam_R6 = lam_R + alpha * h**4 * lam_5**2

            if abs(lam_L) > 1e-14:
                U_hat[ki, li] = (Rg_hat[ki, li] + alpha * h**4 * lam_5**2 * G_hat[ki, li]) / lam_L
                # Note: Rg_hat already contains the R_h^(4) * g, so we just add the correction

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


def build_Lh_matrix(N, h):
    """Build L_h = (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1]"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            rows.append(idx); cols.append(idx); vals.append(-20.0/(6.0*h**2))
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(4.0/(6.0*h**2))
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(4.0/(6.0*h**2))
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(4.0/(6.0*h**2))
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(4.0/(6.0*h**2))
            if i>0 and j>0: rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j>0: rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i>0 and j<N-1: rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j<N-1: rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def matrix_6th_order(n, f_func, bc_func, alpha=1.0/360.0, sx=1.0, sy=1.0):
    """
    6th-order solver using sparse matrix (for verification).
    Scheme: L_h * u = R_h^(4) * g + alpha * h^4 * Delta_h^2 * g
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    h = hx
    N = n - 2

    F = f_func(X, Y)
    G = -F

    # Build L_h matrix
    Lh = build_Lh_matrix(N, h)

    # R_h^(4) * g on interior
    G_int = G[1:-1, 1:-1]
    rhs = apply_Rh4(G_int, N)

    # Delta_h^2 * g correction
    delta2_g = apply_lap5_sq_interior(G_int, h)
    rhs += alpha * h**4 * delta2_g

    # BC correction
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        val_l = 4.0 * bc_left[jj]
        if jj-1 >= 0: val_l += bc_left[jj-1]
        if jj+1 < n: val_l += bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val_l
        val_r = 4.0 * bc_right[jj]
        if jj-1 >= 0: val_r += bc_right[jj-1]
        if jj+1 < n: val_r += bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val_r
    for i in range(N):
        ii = i + 1
        val_b = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val_b += bc_bottom[ii-1]
        if ii+1 < n: val_b += bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val_b
        val_t = 4.0 * bc_top[ii]
        if ii-1 >= 0: val_t += bc_top[ii-1]
        if ii+1 < n: val_t += bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val_t

    rhs += bc_corr

    u_int = spsolve(Lh, rhs.flatten()).reshape((N, N))

    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = u_int
    return U


def test_6th_order():
    """Test 6th-order convergence with different correction approaches"""
    print("=" * 80)
    print("FFT9 6th-Order Scheme Test")
    print("=" * 80)

    # Test 1: Smooth trigonometric solution
    print("\n--- Test 1: u = sin(pi*x)*sin(pi*y), zero BC ---")
    def u_exact1(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs1(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc1(x, y): return 0.0

    ns = [9, 17, 33, 65, 129]

    # Test different alpha values for the correction term
    alphas = [1.0/360.0, 1.0/180.0, 1.0/240.0, 1.0/720.0, 1.0/120.0]

    for alpha in alphas:
        print(f"\n  alpha = {alpha:.6f} (1/{1.0/alpha:.0f})")
        errors = []
        for n in ns:
            U = fft9_6th_order_physical(n, f_rhs1, bc1, alpha=alpha)
            x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
            err = np.max(np.abs(U - u_exact1(X, Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")

    # Test 2: Paper's problem (polynomial solution)
    print("\n\n--- Test 2: u = 3xy(x-x^2)(y-y^2) [Paper Table I] ---")
    def u_exact2(x, y): return 3.0*x*y*(x-x**2)*(y-y**2)
    def f_rhs2(x, y):
        return -3.0*(2.0-6.0*x)*y**2*(1.0-y) - 3.0*x**2*(1.0-x)*(2.0-6.0*y)
    def bc2(x, y): return u_exact2(x, y)

    for alpha in alphas:
        print(f"\n  alpha = {alpha:.6f} (1/{1.0/alpha:.0f})")
        errors = []
        for n in [5, 9, 17, 33, 65]:
            U = fft9_6th_order_physical(n, f_rhs2, bc2, alpha=alpha)
            x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
            err = np.max(np.abs(U - u_exact2(X, Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        for i in range(1, len(errors)):
            rate = np.log(errors[i-1]/errors[i])/np.log([5,9,17,33,65][i]/[5,9,17,33,65][i-1])
            print(f"    rate = {rate:.2f}")

    # Reference: paper's 6th-order results
    print("\n  Paper 6th-order reference:")
    print("    N=4: 1.34E-07, N=8: 2.10E-09, N=16: 3.30E-11, N=32: 9.49E-13")


def test_4th_vs_6th():
    """Compare 4th and 6th order schemes"""
    print("\n" + "=" * 80)
    print("4th vs 6th Order Comparison")
    print("=" * 80)

    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc(x, y): return 0.0

    ns = [9, 17, 33, 65, 129]

    print("\n--- 4th order ---")
    e4 = []
    for n in ns:
        U = fft9_4th_order(n, f_rhs, bc)
        x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y))); e4.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")

    print("\n--- 6th order (alpha=1/360) ---")
    e6 = []
    for n in ns:
        U = fft9_6th_order_physical(n, f_rhs, bc, alpha=1.0/360.0)
        x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y))); e6.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")

    print("\n--- Convergence rates ---")
    print(f"  {'n':>3s}  {'4th':>12s}  {'rate':>5s}  {'6th':>12s}  {'rate':>5s}")
    for i in range(len(ns)):
        if i == 0:
            print(f"  {ns[i]:3d}  {e4[i]:12.6e}   ---  {e6[i]:12.6e}   ---")
        else:
            r4 = np.log(e4[i-1]/e4[i])/np.log(ns[i]/ns[i-1])
            r6 = np.log(e6[i-1]/e6[i])/np.log(ns[i]/ns[i-1])
            print(f"  {ns[i]:3d}  {e4[i]:12.6e}  {r4:5.2f}  {e6[i]:12.6e}  {r6:5.2f}")


def test_fourier_correction():
    """Test 6th order using Fourier-space correction"""
    print("\n" + "=" * 80)
    print("6th Order via Fourier-Space Correction")
    print("=" * 80)

    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc(x, y): return 0.0

    ns = [9, 17, 33, 65, 129]

    # Try different corrections in Fourier space
    print("\n--- Correction: alpha * h^4 * lambda_5^2 ---")
    for alpha in [1.0/360.0, 1.0/180.0, 1.0/720.0]:
        print(f"\n  alpha = {alpha:.6f}")
        errors = []
        for n in ns:
            U = fft9_6th_order_fourier(n, f_rhs, bc, alpha=alpha)
            x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
            err = np.max(np.abs(U - u_exact(X, Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")


def test_matrix_verification():
    """Verify 6th order scheme using sparse matrix solver"""
    print("\n" + "=" * 80)
    print("6th Order - Matrix Verification")
    print("=" * 80)

    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc(x, y): return 0.0

    ns = [9, 17, 33, 65]

    for alpha in [1.0/360.0, 1.0/180.0]:
        print(f"\n  alpha = {alpha:.6f} (Matrix solver)")
        errors = []
        for n in ns:
            U = matrix_6th_order(n, f_rhs, bc, alpha=alpha)
            x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
            err = np.max(np.abs(U - u_exact(X, Y)))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_4th_vs_6th()
    test_6th_order()
    test_matrix_verification()
