#!/usr/bin/env python3
"""
FINAL verification: Build 9-point scheme with sparse matrix and test.
Key question: Does L_h * u = R_h * g give 4th-order convergence?
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def build_Lh_9point(N, h):
    """Build L_h = (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] as sparse matrix"""
    total = N * N
    rows, cols, vals = [], [], []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # Center: -20
            rows.append(idx); cols.append(idx); vals.append(-20.0 / (6.0 * h**2))
            # Left/Right: 4
            if i > 0:
                rows.append(idx); cols.append((i-1)*N+j); vals.append(4.0 / (6.0 * h**2))
            if i < N-1:
                rows.append(idx); cols.append((i+1)*N+j); vals.append(4.0 / (6.0 * h**2))
            # Top/Bottom: 4
            if j > 0:
                rows.append(idx); cols.append(i*N+(j-1)); vals.append(4.0 / (6.0 * h**2))
            if j < N-1:
                rows.append(idx); cols.append(i*N+(j+1)); vals.append(4.0 / (6.0 * h**2))
            # Diagonals: 1
            if i > 0 and j > 0:
                rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(1.0 / (6.0 * h**2))
            if i < N-1 and j > 0:
                rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(1.0 / (6.0 * h**2))
            if i > 0 and j < N-1:
                rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(1.0 / (6.0 * h**2))
            if i < N-1 and j < N-1:
                rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(1.0 / (6.0 * h**2))
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def build_Rh_5point(N, h):
    """
    Build R_h = I + (h^2/12) * Laplacian
    = I + (1/12) * [0 -1 0; -1 4 -1; 0 -1 0] (5-point)
    
    R_h template: [0 -1/12 0; -1/12 4/3 -1/12; 0 -1/12 0]
    This is the CORRECT R_h for 4th-order compact scheme.
    """
    total = N * N
    rows, cols, vals = [], [], []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # Center: 4/3
            rows.append(idx); cols.append(idx); vals.append(4.0/3.0)
            # Left/Right: -1/12
            if i > 0:
                rows.append(idx); cols.append((i-1)*N+j); vals.append(-1.0/12.0)
            if i < N-1:
                rows.append(idx); cols.append((i+1)*N+j); vals.append(-1.0/12.0)
            # Top/Bottom: -1/12
            if j > 0:
                rows.append(idx); cols.append(i*N+(j-1)); vals.append(-1.0/12.0)
            if j < N-1:
                rows.append(idx); cols.append(i*N+(j+1)); vals.append(-1.0/12.0)
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def build_Rh_9point(N, h):
    """
    Build R_h with 9-point template (including diagonal terms).
    
    For 4th-order: R_h = I + (h^2/12)*Laplacian
    The 9-point version includes diagonal corrections.
    
    R_h = [1/12   -1/6   1/12]
          [ -1/6    4/3   -1/6]
          [1/12   -1/6   1/12]
    
    Wait, this doesn't look right either. Let me think more carefully.
    
    Actually, the standard 4th-order 9-point scheme from the paper uses:
    R_h = (1/360) * [0 48 0; 48 0 48; 0 48 0]
    But this has no center term and sum = 192/360 = 8/15 ≠ 1
    
    Let me try another approach. The 4th-order compact scheme is:
    L_h * u = R_h * g
    where L_h approximates Laplacian with O(h^2) error
    and R_h approximates Identity with O(h^2) error
    but combined they give O(h^4).
    
    The key: T_L = (h^2/12)*Lap^2(u) + O(h^4)
             T_R should equal (h^2/12)*Lap(g) + O(h^4)
    since g = Lap(u), we have Lap(g) = Lap^2(u)
    So T_L = T_R to O(h^2), and the difference is O(h^4). ✓
    
    R_h = I + (h^2/12)*Lap_h  (where Lap_h is the 5-point Laplacian)
    
    In stencil form:
    R_h = I + (1/12)*[0 -1 0; -1 4 -1; 0 -1 0]
        = [0 -1/12 0; -1/12 4/3 -1/12; 0 -1/12 0]  (5 non-zero entries)
    
    This is correct but only has 5 points, not 9.
    """
    # Just use the 5-point version
    return build_Rh_5point(N, h)


def build_A5(N, h):
    """Standard 5-point Laplacian for -Laplacian(u) = f"""
    total = N * N
    rows, cols, vals = [], [], []
    
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            rows.append(idx); cols.append(idx); vals.append(4.0/h**2)
            if i > 0:
                rows.append(idx); cols.append((i-1)*N+j); vals.append(-1.0/h**2)
            if i < N-1:
                rows.append(idx); cols.append((i+1)*N+j); vals.append(-1.0/h**2)
            if j > 0:
                rows.append(idx); cols.append(i*N+(j-1)); vals.append(-1.0/h**2)
            if j < N-1:
                rows.append(idx); cols.append(i*N+(j+1)); vals.append(-1.0/h**2)
    
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def test_all():
    """Test all methods"""
    print("=" * 70)
    print("Complete verification of 9-point vs 5-point methods")
    print("=" * 70)
    
    def u_exact(x, y):
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # For -Laplacian(u) = f:
    def f_rhs(x, y):
        return 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    # g = Laplacian(u) = -f
    def g_func(x, y):
        return -2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    
    print("\nMethod 1: 5-point (2nd order) - A5*u = f")
    print("Method 2: 9-point (4th order compact) - L_h*u = R_h*g")
    print("Method 3: 9-point L_h alone (2nd order) - L_h*u = g (no R_h)")
    print()
    
    ns = [5, 9, 17, 33, 65, 129]
    errors = {1: [], 2: [], 3: []}
    
    for n in ns:
        N = n - 2
        h = 1.0 / (n - 1)
        
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        U_ex = u_exact(X[1:-1, 1:-1], Y[1:-1, 1:-1])
        
        # Method 1: 5-point
        A5 = build_A5(N, h)
        f_vec = f_rhs(X[1:-1, 1:-1], Y[1:-1, 1:-1]).flatten()
        u5 = spsolve(A5, f_vec).reshape((N, N))
        e5 = np.max(np.abs(u5 - U_ex))
        errors[1].append(e5)
        
        # Method 2: 9-point compact (L_h * u = R_h * g)
        Lh = build_Lh_9point(N, h)
        Rh = build_Rh_5point(N, h)
        g_vec = g_func(X[1:-1, 1:-1], Y[1:-1, 1:-1]).flatten()
        rhs2 = Rh @ g_vec
        u9 = spsolve(Lh, rhs2).reshape((N, N))
        e9 = np.max(np.abs(u9 - U_ex))
        errors[2].append(e9)
        
        # Method 3: 9-point L_h only (L_h * u = g, without R_h)
        u9b = spsolve(Lh, g_vec).reshape((N, N))
        e9b = np.max(np.abs(u9b - U_ex))
        errors[3].append(e9b)
        
        print(f"n={n:3d}: 5pt={e5:.6e}, 9pt_compact={e9:.6e}, 9pt_only={e9b:.6e}")
    
    # Convergence rates
    print("\nConvergence rates:")
    print(f"  {'n':>3s}  {'5pt':>8s}  {'9compact':>8s}  {'9only':>8s}")
    for i in range(1, len(ns)):
        r5 = np.log(errors[1][i-1]/errors[1][i]) / np.log(ns[i]/ns[i-1])
        r9c = np.log(errors[2][i-1]/errors[2][i]) / np.log(ns[i]/ns[i-1])
        r9o = np.log(errors[3][i-1]/errors[3][i]) / np.log(ns[i]/ns[i-1])
        print(f"  {ns[i]:3d}  {r5:8.2f}  {r9c:8.2f}  {r9o:8.2f}")


if __name__ == "__main__":
    test_all()
