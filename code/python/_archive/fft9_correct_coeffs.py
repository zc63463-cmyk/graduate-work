#!/usr/bin/env python3
"""
Derive and verify the correct 9-point stencil coefficients.

The key insight: The standard 4th-order 9-point scheme for Laplacian(u) = f is:

L_h * u = R_h * f

where:
  L_h = (1/6) * [1  4  1; 4  -20  4; 1  4  1]   (divided by h^2)
  R_h = [0  1/12  0; 1/12  5/6  1/12; 0  1/12  0]   (NOT divided by h^2)

This is the 4th-order version. The 6th-order version needs additional terms.

The paper's R_h template (1/360)*[0 48 0; 48 0 48; 0 48 0] 
= [0  2/15  0; 2/15  0  2/15; 0  2/15  0]
does NOT have a center term, which seems wrong for a 6th-order scheme.

Let me verify by Taylor expansion.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.fft import dst, idst


def test_correct_9point():
    """
    Test the standard 4th-order 9-point compact scheme.
    
    For Laplacian(u) = f, the 4th-order scheme is:
    (1/(6h^2)) * [1  4  1; 4  -20  4; 1  4  1] u 
    = [0  1/12  0; 1/12  5/6  1/12; 0  1/12  0] f
    """
    print("=" * 70)
    print("Standard 4th-order 9-point compact scheme")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # For -Laplacian(u) = f: f = 2*pi^2 * sin(pi*x)*sin(pi*y)
    # Laplacian(u) = -f
    
    for n in [5, 9, 17, 33, 65]:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # Exact solution at interior points
        U_ex = u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        
        # RHS: -f (for Laplacian(u) = -f)
        neg_f = -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # Apply R_h to (-f) at interior points
        # R_h = [0  1/12  0; 1/12  5/6  1/12; 0  1/12  0]
        Rf = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                Rf[i, j] = (5.0/6.0) * neg_f[ii, jj] \
                         + (1.0/12.0) * (neg_f[ii-1, jj] + neg_f[ii+1, jj] 
                                       + neg_f[ii, jj-1] + neg_f[ii, jj+1])
        
        # Build L_h matrix
        A = build_Lh_matrix(N, h)
        
        # Solve
        u_vec = spsolve(A, Rf.flatten())
        U_num = u_vec.reshape((N, N))
        
        err = np.max(np.abs(U_num - U_ex))
        print(f"  n={n:3d}: error = {err:.6e}")


def test_6th_order_9point():
    """
    Test the 6th-order 9-point compact scheme.
    
    For Laplacian(u) = f with h = hx = hy, the 6th-order scheme is:
    
    (1/(6h^2)) * [1  4  1; 4  -20  4; 1  4  1] u 
    = [0  1/12  0; 1/12  5/6  1/12; 0  1/12  0] f
    + (h^2/12) * [0  0  0; 0  1  0; 0  0  0] * Laplacian(f)
    
    This adds a correction term involving Laplacian(f) to achieve 6th order.
    
    Alternatively, we can write:
    (1/(6h^2)) * [1  4  1; 4  -20  4; 1  4  1] u 
    = R_4th f + (h^2/12) * f_xx + (h^2/12) * f_yy
    
    where R_4th is the 4th-order R_h template.
    
    But computing f_xx and f_yy requires knowing f analytically.
    We can approximate them using finite differences of f.
    """
    print("\n" + "=" * 70)
    print("6th-order 9-point compact scheme (with f-Laplacian correction)")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    for n in [5, 9, 17, 33, 65]:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        U_ex = u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        
        # -f for Laplacian(u) = -f
        neg_f = -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # Apply 4th-order R_h to (-f)
        Rf_4th = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                Rf_4th[i, j] = (5.0/6.0) * neg_f[ii, jj] \
                             + (1.0/12.0) * (neg_f[ii-1, jj] + neg_f[ii+1, jj] 
                                           + neg_f[ii, jj-1] + neg_f[ii, jj+1])
        
        # Compute Laplacian(f) correction term
        # Laplacian(f) = f_xx + f_yy
        # For f = 2*pi^2 * sin(pi*x)*sin(pi*y):
        # f_xx = -2*pi^4 * sin(pi*x)*sin(pi*y) (but this is -f * pi^2)
        # f_yy = -2*pi^4 * sin(pi*x)*sin(pi*y)
        # Laplacian(f) = -4*pi^4 * sin(pi*x)*sin(pi*y) = -2*pi^2 * f
        
        # But we need to approximate Laplacian(f) numerically
        # Use 5-point stencil for Laplacian(f)
        lap_f = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                lap_f[i, j] = (1.0/h**2) * (
                    neg_f[ii-1, jj] + neg_f[ii+1, jj] + neg_f[ii, jj-1] + neg_f[ii, jj+1]
                    - 4.0 * neg_f[ii, jj]
                )
        
        # 6th-order RHS = 4th-order RHS + (h^2/12) * Laplacian(f)
        Rf_6th = Rf_4th + (h**2 / 12.0) * lap_f
        
        # Build L_h matrix and solve
        A = build_Lh_matrix(N, h)
        u_vec = spsolve(A, Rf_6th.flatten())
        U_num = u_vec.reshape((N, N))
        
        err = np.max(np.abs(U_num - U_ex))
        print(f"  n={n:3d}: error = {err:.6e}")


def build_Lh_matrix(N, h):
    """Build the L_h matrix for 9-point stencil"""
    total = N * N
    rows = []
    cols = []
    vals = []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            
            rows.append(idx); cols.append(idx); vals.append(-20.0 / (6.0 * h**2))
            
            if i > 0:
                rows.append(idx); cols.append((i-1)*N+j); vals.append(4.0 / (6.0 * h**2))
            if i < N-1:
                rows.append(idx); cols.append((i+1)*N+j); vals.append(4.0 / (6.0 * h**2))
            if j > 0:
                rows.append(idx); cols.append(i*N+(j-1)); vals.append(4.0 / (6.0 * h**2))
            if j < N-1:
                rows.append(idx); cols.append(i*N+(j+1)); vals.append(4.0 / (6.0 * h**2))
            
            if i > 0 and j > 0:
                rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(1.0 / (6.0 * h**2))
            if i < N-1 and j > 0:
                rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(1.0 / (6.0 * h**2))
            if i > 0 and j < N-1:
                rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(1.0 / (6.0 * h**2))
            if i < N-1 and j < N-1:
                rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(1.0 / (6.0 * h**2))
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def test_fft9_correct():
    """
    Test FFT9 with the CORRECT 4th-order R_h template (including center term).
    """
    print("\n" + "=" * 70)
    print("FFT9 with correct 4th-order R_h (including center 5/6 term)")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    ns = [5, 9, 17, 33, 65]
    errors = []
    
    for n in ns:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # -f (for Laplacian(u) = -f, i.e., -Laplacian(u) = f)
        neg_f = -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # Apply R_h to (-f) - CORRECT VERSION with center term
        # R_h = [0  1/12  0; 1/12  5/6  1/12; 0  1/12  0]
        RF = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                RF[i, j] = (5.0/6.0) * neg_f[ii, jj] \
                         + (1.0/12.0) * (neg_f[ii-1, jj] + neg_f[ii+1, jj] 
                                       + neg_f[ii, jj-1] + neg_f[ii, jj+1])
        
        # 2D DST
        RF_hat = np.zeros((N, N))
        for i in range(N):
            RF_hat[i, :] = dst(RF[i, :], type=1, norm='ortho')
        for j in range(N):
            RF_hat[:, j] = dst(RF_hat[:, j], type=1, norm='ortho')
        
        # Eigenvalues of L_h
        k_vals = np.arange(1, N + 1)
        cos_k = np.cos(k_vals * np.pi / (N + 1))
        
        # Eigenvalues of R_h (with center term)
        # R_h f_{i,j} = (5/6)*f_{i,j} + (1/12)*(f_{i-1,j}+f_{i+1,j}+f_{i,j-1}+f_{i,j+1})
        # For eigenfunction sin(k*i*theta)*sin(l*j*phi):
        # lambda_R(k,l) = 5/6 + (1/12)*(2*cos(k*theta) + 2*cos(l*phi))
        #               = 5/6 + (1/6)*(cos(k*theta) + cos(l*phi))
        
        U_hat = np.zeros((N, N))
        for ki in range(N):
            for li in range(N):
                lambda_A = (1.0/(6.0*h**2)) * (
                    -20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li] + 4.0*cos_k[ki]*cos_k[li]
                )
                # lambda_R = 5.0/6.0 + (1.0/6.0)*(cos_k[ki] + cos_k[li])
                # But we already applied R_h in physical space, so we just divide by lambda_A
                if abs(lambda_A) > 1e-14:
                    U_hat[ki, li] = RF_hat[ki, li] / lambda_A
        
        # Inverse 2D DST
        U_int = np.zeros((N, N))
        for ki in range(N):
            U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
        for li in range(N):
            U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
        
        U = np.zeros((n, n))
        U[1:-1, 1:-1] = U_int
        
        # Error
        U_ex = u_exact(X, Y)
        err = np.max(np.abs(U - U_ex))
        errors.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")
    
    # Convergence rates
    print("\nConvergence rates:")
    for i in range(1, len(ns)):
        rate = np.log(errors[i-1]/errors[i]) / np.log(ns[i]/ns[i-1])
        print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_correct_9point()
    test_6th_order_9point()
    test_fft9_correct()
