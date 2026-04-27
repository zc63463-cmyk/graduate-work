#!/usr/bin/env python3
"""
Verify the 9-point stencil by direct matrix construction and comparison.
This will tell us whether the stencil itself is correct.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.fft import dst, idst

def build_9point_matrix(N, h):
    """
    Build the 9-point finite difference matrix for Laplacian on NxN interior points.
    
    L_h = (1/(6h^2)) * [1  4  1; 4 -20 4; 1  4  1]
    
    This approximates the Laplacian operator.
    """
    total = N * N  # total unknowns
    
    # Build sparse matrix
    rows = []
    cols = []
    vals = []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j  # row index
            
            # Center point: -20/(6h^2)
            rows.append(idx)
            cols.append(idx)
            vals.append(-20.0 / (6.0 * h**2))
            
            # Left neighbor: 4/(6h^2)
            if i > 0:
                rows.append(idx)
                cols.append((i-1) * N + j)
                vals.append(4.0 / (6.0 * h**2))
            
            # Right neighbor: 4/(6h^2)
            if i < N-1:
                rows.append(idx)
                cols.append((i+1) * N + j)
                vals.append(4.0 / (6.0 * h**2))
            
            # Bottom neighbor: 4/(6h^2)
            if j > 0:
                rows.append(idx)
                cols.append(i * N + (j-1))
                vals.append(4.0 / (6.0 * h**2))
            
            # Top neighbor: 4/(6h^2)
            if j < N-1:
                rows.append(idx)
                cols.append(i * N + (j+1))
                vals.append(4.0 / (6.0 * h**2))
            
            # Bottom-left diagonal: 1/(6h^2)
            if i > 0 and j > 0:
                rows.append(idx)
                cols.append((i-1) * N + (j-1))
                vals.append(1.0 / (6.0 * h**2))
            
            # Bottom-right diagonal: 1/(6h^2)
            if i < N-1 and j > 0:
                rows.append(idx)
                cols.append((i+1) * N + (j-1))
                vals.append(1.0 / (6.0 * h**2))
            
            # Top-left diagonal: 1/(6h^2)
            if i > 0 and j < N-1:
                rows.append(idx)
                cols.append((i-1) * N + (j+1))
                vals.append(1.0 / (6.0 * h**2))
            
            # Top-right diagonal: 1/(6h^2)
            if i < N-1 and j < N-1:
                rows.append(idx)
                cols.append((i+1) * N + (j+1))
                vals.append(1.0 / (6.0 * h**2))
    
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))
    return A


def build_Rh_matrix(N, h):
    """
    Build the R_h matrix for the RHS.
    
    R_h = (1/360) * [0 48 0; 48 0 48; 0 48 0]
    
    This approximates the identity operator.
    """
    total = N * N
    
    rows = []
    cols = []
    vals = []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            
            # Left neighbor: 48/360
            if i > 0:
                rows.append(idx)
                cols.append((i-1) * N + j)
                vals.append(48.0 / 360.0)
            
            # Right neighbor: 48/360
            if i < N-1:
                rows.append(idx)
                cols.append((i+1) * N + j)
                vals.append(48.0 / 360.0)
            
            # Bottom neighbor: 48/360
            if j > 0:
                rows.append(idx)
                cols.append(i * N + (j-1))
                vals.append(48.0 / 360.0)
            
            # Top neighbor: 48/360
            if j < N-1:
                rows.append(idx)
                cols.append(i * N + (j+1))
                vals.append(48.0 / 360.0)
    
    B = sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))
    return B


def test_9point_stencil():
    """
    Test the 9-point stencil by solving L_h u = R_h f directly.
    """
    print("=" * 70)
    print("Verify 9-point stencil by direct matrix solve")
    print("=" * 70)
    
    # Problem: -Laplacian(u) = f, u = sin(pi*x)*sin(pi*y) on boundary
    # Equation: Laplacian(u) = -f
    # 9-point: L_h u = R_h * (-f)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # f for -Laplacian(u) = f: f = 2*pi^2 * sin(pi*x)*sin(pi*y)
    # -f for Laplacian(u) = -f: -f = -2*pi^2 * sin(pi*x)*sin(pi*y)
    def neg_f(x, y):
        return -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    for n in [5, 9, 17, 33]:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # Build matrices
        A = build_9point_matrix(N, h)
        B = build_Rh_matrix(N, h)
        
        # RHS vector: -f at interior points
        neg_f_vals = neg_f(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        rhs_vec = B @ neg_f_vals.flatten()
        
        # BC correction
        bc_left = u_exact(x[0], y)
        bc_right = u_exact(x[-1], y)
        bc_bottom = u_exact(x, y[0])
        bc_top = u_exact(x, y[-1])
        
        # For zero BC, no correction needed
        # But let's check: u_exact(0,y) = sin(0)*sin(pi*y) = 0
        # u_exact(1,y) = sin(pi)*sin(pi*y) = 0
        # So BC is zero, no correction needed
        
        # Solve
        u_vec = spsolve(A, rhs_vec)
        U_num = u_vec.reshape((N, N))
        
        # Exact solution at interior points
        U_ex = u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        
        # Error
        err = np.max(np.abs(U_num - U_ex))
        
        # Also compute 5-point error for comparison
        A5 = build_5point_matrix(N, h)
        f_vals = 2.0 * np.pi**2 * np.sin(np.pi * X[1:-1, 1:-1]) * np.sin(np.pi * Y[1:-1, 1:-1])
        u5_vec = spsolve(A5, f_vals.flatten())
        U5_num = u5_vec.reshape((N, N))
        err5 = np.max(np.abs(U5_num - U_ex))
        
        print(f"  n={n:3d}: 5pt error = {err5:.6e}, 9pt error = {err:.6e}")


def build_5point_matrix(N, h):
    """Build standard 5-point Laplacian matrix for -Laplacian(u) = f"""
    total = N * N
    rows = []
    cols = []
    vals = []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            
            rows.append(idx)
            cols.append(idx)
            vals.append(4.0 / h**2)
            
            if i > 0:
                rows.append(idx); cols.append((i-1)*N+j); vals.append(-1.0/h**2)
            if i < N-1:
                rows.append(idx); cols.append((i+1)*N+j); vals.append(-1.0/h**2)
            if j > 0:
                rows.append(idx); cols.append(i*N+(j-1)); vals.append(-1.0/h**2)
            if j < N-1:
                rows.append(idx); cols.append(i*N+(j+1)); vals.append(-1.0/h**2)
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def test_9point_consistency():
    """
    Test that the 9-point stencil is a consistent approximation of Laplacian.
    Apply L_h to a known function and check the result.
    """
    print("\n" + "=" * 70)
    print("9-point stencil consistency check")
    print("=" * 70)
    
    for n in [5, 9, 17, 33, 65]:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # Test function: u = sin(pi*x)*sin(pi*y)
        U = np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # Exact Laplacian: Laplacian(u) = -2*pi^2 * sin(pi*x)*sin(pi*y)
        lap_exact = -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)
        
        # Apply L_h to U at interior points
        lap_num = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                lap_num[i, j] = (1.0 / (6.0 * h**2)) * (
                    U[ii-1, jj-1] + 4*U[ii, jj-1] + U[ii+1, jj-1]
                    + 4*U[ii-1, jj] - 20*U[ii, jj] + 4*U[ii+1, jj]
                    + U[ii-1, jj+1] + 4*U[ii, jj+1] + U[ii+1, jj+1]
                )
        
        # Compare with exact Laplacian at interior points
        lap_ex_int = lap_exact[1:-1, 1:-1]
        err = np.max(np.abs(lap_num - lap_ex_int))
        
        # Also check with R_h applied to Laplacian
        # L_h u should approximate Laplacian(u)
        # But the 6th-order scheme says: L_h u = R_h * Laplacian(u)
        # So we should compare L_h u with R_h * lap_exact
        
        # Apply R_h to lap_exact
        Rlap = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                Rlap[i, j] = (48.0 / 360.0) * (
                    lap_exact[ii-1, jj] + lap_exact[ii+1, jj]
                    + lap_exact[ii, jj-1] + lap_exact[ii, jj+1]
                )
        
        err_Rh = np.max(np.abs(lap_num - Rlap))
        
        print(f"  n={n:3d}: |L_h u - Lap(u)| = {err:.6e}, |L_h u - R_h*Lap(u)| = {err_Rh:.6e}")


def test_fft_vs_matrix():
    """
    Compare FFT9 solver with direct matrix solve.
    """
    print("\n" + "=" * 70)
    print("FFT9 vs Matrix solve comparison")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def f_rhs(x, y):
        return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    def bc(x, y):
        return 0.0
    
    for n in [5, 9]:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # Matrix solve
        A = build_9point_matrix(N, h)
        B = build_Rh_matrix(N, h)
        
        # For -Laplacian(u) = f, we have Laplacian(u) = -f
        # So L_h u = R_h * (-f)
        neg_f_vals = -f_rhs(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        rhs_vec = B @ neg_f_vals.flatten()
        u_mat = spsolve(A, rhs_vec).reshape((N, N))
        
        # FFT solve
        U_fft = np.zeros((n, n))
        F = f_rhs(X, Y)
        
        # Apply -R_h to f (for -Laplacian(u) = f)
        RF = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                ii = i + 1
                jj = j + 1
                RF[i, j] = -(48.0/360.0) * (F[ii-1,jj] + F[ii+1,jj] + F[ii,jj-1] + F[ii,jj+1])
        
        # 2D DST
        RF_hat = np.zeros((N, N))
        for i in range(N):
            RF_hat[i, :] = dst(RF[i, :], type=1, norm='ortho')
        for j in range(N):
            RF_hat[:, j] = dst(RF_hat[:, j], type=1, norm='ortho')
        
        # Solve in frequency domain
        k_vals = np.arange(1, N+1)
        cos_k = np.cos(k_vals * np.pi / (N+1))
        
        U_hat = np.zeros((N, N))
        for ki in range(N):
            for li in range(N):
                lam = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li] + 4.0*cos_k[ki]*cos_k[li])
                if abs(lam) > 1e-14:
                    U_hat[ki, li] = RF_hat[ki, li] / lam
        
        # Inverse DST
        U_int = np.zeros((N, N))
        for ki in range(N):
            U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
        for li in range(N):
            U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')
        
        U_fft[1:-1, 1:-1] = U_int
        
        # Compare
        diff = np.max(np.abs(u_mat - U_int))
        err_mat = np.max(np.abs(u_mat - u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])))
        err_fft = np.max(np.abs(U_int - u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])))
        
        print(f"  n={n}: |matrix - FFT| = {diff:.6e}")
        print(f"         matrix error = {err_mat:.6e}")
        print(f"         FFT error    = {err_fft:.6e}")


if __name__ == "__main__":
    test_9point_stencil()
    test_9point_consistency()
    test_fft_vs_matrix()
