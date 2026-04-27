#!/usr/bin/env python3
"""Validate GMRES implementation against scipy."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.sparse.linalg import gmres as scipy_gmres
from gmres_solver import gmres, build_helmholtz_matrix

# Build a test problem that requires multiple iterations
n = 33
k2 = 100.0
A, grid_info = build_helmholtz_matrix(n, k2, 'dirichlet')
N = A.shape[0]

# Use a RHS with many Fourier modes
x = np.linspace(0, 1, n)
X, Y = np.meshgrid(x[1:-1], x[1:-1], indexing='ij')
def f_rhs(X, Y):
    # Excites many eigenvectors
    return np.exp(-10 * ((X - 0.5)**2 + (Y - 0.5)**2))
b = f_rhs(X, Y).reshape(-1)

print(f"Matrix size: {N}")
print(f"k2 = {k2}")
print()

# Test 1: No restart (restart >= N)
print("=" * 60)
print("Test 1: GMRES without restart (restart=N)")
print("=" * 60)
x_ours, info = gmres(A, b, tol=1e-10, restart=N, max_iter=N, return_history=True)
x_scipy, flag = scipy_gmres(A, b, rtol=1e-10, restart=N, maxiter=N)
print(f"  Our GMRES:  iters={info['iterations']}, success={info['success']}, final_res={info['residuals'][-1]:.2e}")
print(f"  scipy GMRES: converged={flag==0}, final_res={np.linalg.norm(b - A @ x_scipy):.2e}")
print(f"  Solution diff: {np.linalg.norm(x_ours - x_scipy):.2e}")
print()

# Test 2: With restart (restart=30)
print("=" * 60)
print("Test 2: GMRES with restart m=30")
print("=" * 60)
x_ours2, info2 = gmres(A, b, tol=1e-10, restart=30, max_iter=5000, return_history=True)
x_scipy2, flag2 = scipy_gmres(A, b, rtol=1e-10, restart=30, maxiter=5000)
print(f"  Our GMRES(30):  iters={info2['iterations']}, success={info2['success']}, final_res={info2['residuals'][-1]:.2e}")
print(f"  scipy GMRES(30): converged={flag2==0}, final_res={np.linalg.norm(b - x_scipy2):.2e}")
print(f"  Solution diff: {np.linalg.norm(x_ours2 - x_scipy2):.2e}")
print()

# Test 3: Small problem, no restart, compare exactly
print("=" * 60)
print("Test 3: Small problem, exact comparison")
print("=" * 60)
n_small = 9
k2_small = 10.0
A_small, _ = build_helmholtz_matrix(n_small, k2_small, 'dirichlet')
N_small = A_small.shape[0]
x_s = np.linspace(0, 1, n_small)
X_s, Y_s = np.meshgrid(x_s[1:-1], x_s[1:-1], indexing='ij')
b_small = f_rhs(X_s, Y_s).reshape(-1)

x_ours3, info3 = gmres(A_small, b_small, tol=1e-12, restart=N_small, max_iter=N_small, return_history=True)
x_scipy3, flag3 = scipy_gmres(A_small, b_small, rtol=1e-12, restart=N_small, maxiter=N_small)
res_ours = np.linalg.norm(b_small - A_small @ x_ours3)
res_scipy = np.linalg.norm(b_small - A_small @ x_scipy3)
print(f"  Our GMRES:  iters={info3['iterations']}, res={res_ours:.2e}")
print(f"  scipy GMRES: converged={flag3==0}, res={res_scipy:.2e}")
print(f"  Solution diff: {np.linalg.norm(x_ours3 - x_scipy3):.2e}")

# Detailed residual history
print(f"\n  Our residual history ({len(info3['residuals'])} entries):")
for i, r in enumerate(info3['residuals']):
    print(f"    iter {i}: residual = {r:.6e}")

# Check restart bug: compare restart vs no-restart on same problem
print()
print("=" * 60)
print("Test 4: Restart bug check - m=5 vs m=N")
print("=" * 60)
x_nr, info_nr = gmres(A_small, b_small, tol=1e-10, restart=N_small, max_iter=2*N_small, return_history=True)
x_r5, info_r5 = gmres(A_small, b_small, tol=1e-10, restart=5, max_iter=2*N_small, return_history=True)
res_nr = np.linalg.norm(b_small - A_small @ x_nr)
res_r5 = np.linalg.norm(b_small - A_small @ x_r5)
print(f"  No restart: iters={info_nr['iterations']}, res={res_nr:.2e}")
print(f"  restart=5:  iters={info_r5['iterations']}, res={res_r5:.2e}, success={info_r5['success']}")
print(f"  restart=5 residual history: {info_r5['residuals'][:10]}")
