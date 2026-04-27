#!/usr/bin/env python3
"""
Test the paper's exact 6th-order formula with matrix method.

Paper's formula:
  (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] * u 
  = (1/360) * [0 48 0; 48 0 48; 0 48 0] * f

where f = Laplacian(u) (for Lap(u) = f)
or f = -Laplacian(u) (for -Lap(u) = f, depending on convention)
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def build_Lh(N, h):
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


def build_Rh_paper(N):
    """
    Build paper's R_h = (1/360)*[0 48 0; 48 0 48; 0 48 0]
    This is a 5-point template with NO center term.
    """
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # NO center term!
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(48.0/360.0)
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(48.0/360.0)
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(48.0/360.0)
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(48.0/360.0)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def build_Rh_4th(N):
    """Build 4th-order R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0]"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            rows.append(idx); cols.append(idx); vals.append(2.0/3.0)
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(1.0/12.0)
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(1.0/12.0)
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(1.0/12.0)
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(1.0/12.0)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def build_Rh_4th_neg(N):
    """Build R_h = [0 -1/12 0; -1/12 4/3 -1/12; 0 -1/12 0]"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            rows.append(idx); cols.append(idx); vals.append(4.0/3.0)
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(-1.0/12.0)
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(-1.0/12.0)
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(-1.0/12.0)
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(-1.0/12.0)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def test_all_schemes():
    """Test all R_h variants"""
    print("=" * 80)
    print("Test all R_h variants for -Laplacian(u) = f")
    print("u = sin(pi*x)*sin(pi*y)")
    print("=" * 80)
    
    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    
    ns = [5, 9, 17, 33, 65, 129]
    
    # g = Lap(u) = -f
    # Equation: L_h * u = R_h * g = R_h * (-f)
    
    schemes = [
        ("R_h paper [0 2/15 0; 2/15 0 2/15; 0 2/15 0]", "paper"),
        ("R_h 4th [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0]", "4th"),
        ("R_h neg [0 -1/12 0; -1/12 4/3 -1/12; 0 -1/12 0]", "neg"),
    ]
    
    all_errors = {}
    
    for name, scheme in schemes:
        print(f"\n--- {name} ---")
        errors = []
        for n in ns:
            N = n - 2; h = 1.0/(n-1)
            x = np.linspace(0,1,n); y = np.linspace(0,1,n)
            X, Y = np.meshgrid(x, y, indexing='ij')
            U_ex = u_exact(X[1:-1,1:-1], Y[1:-1,1:-1])
            
            Lh = build_Lh(N, h)
            if scheme == "paper":
                Rh = build_Rh_paper(N)
            elif scheme == "4th":
                Rh = build_Rh_4th(N)
            else:
                Rh = build_Rh_4th_neg(N)
            
            # g = -f at interior points
            g_int = -f_rhs(X[1:-1,1:-1], Y[1:-1,1:-1]).flatten()
            rhs = Rh @ g_int
            
            u_num = spsolve(Lh, rhs).reshape((N, N))
            err = np.max(np.abs(u_num - U_ex))
            errors.append(err)
            print(f"  n={n:3d}: error = {err:.6e}")
        
        all_errors[scheme] = errors
    
    # Convergence rates
    print("\n" + "=" * 80)
    print("Convergence Rates")
    print("=" * 80)
    for name, scheme in schemes:
        errors = all_errors[scheme]
        print(f"\n{name}:")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_all_schemes()
