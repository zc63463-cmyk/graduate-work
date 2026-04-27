#!/usr/bin/env python3
"""Diagnose GMRES convergence issue with different test problems."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from gmres_solver import gmres_helmholtz, gmres, test_problem_dirichlet, build_helmholtz_matrix

k2 = 10.0
n = 33

# DIAGNOSIS 1: Simple eigenfunction test
print("=" * 60)
print("DIAGNOSIS 1: Simple test u = sin(pi*x)*sin(pi*y)")
print("=" * 60)
u_exact, f_rhs, bc = test_problem_dirichlet(k2)
U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-10, restart=30, return_history=True)
print(f"  Iterations: {info['iterations']}")
print(f"  Residuals: {info['residuals'][:5]}")
print(f"  WHY: sin(pi*x)*sin(pi*y) is the (1,1) EIGENFUNCTION of the discrete Laplacian!")
print(f"       b = f is exactly in the direction of eigenvector v_{1,1}")
print(f"       GMRES finds solution in 1 step because K_1 = span{{b}} already contains the answer")
print()

# DIAGNOSIS 2: Higher-frequency mode
print("=" * 60)
print("DIAGNOSIS 2: Single higher mode u = sin(2*pi*x)*sin(3*pi*y)")
print("=" * 60)
def f_high(x, y):
    return ((2*np.pi)**2 + (3*np.pi)**2 + k2) * np.sin(2*np.pi*x) * np.sin(3*np.pi*y)
U2, info2 = gmres_helmholtz(n, f_high, bc, k2=k2, tol=1e-10, restart=30, return_history=True)
print(f"  Iterations: {info2['iterations']}")
print(f"  Number of residuals: {len(info2['residuals'])}")
print(f"  STILL trivial: This is also a single eigenfunction (2,3) of the discrete Laplacian!")
print()

# DIAGNOSIS 3: Multi-mode problem
print("=" * 60)
print("DIAGNOSIS 3: Multi-mode problem (5 modes)")
print("=" * 60)
def f_multimode(x, y):
    f = np.zeros_like(x)
    modes = [(1,1), (1,2), (2,3), (3,1), (4,5)]
    for (p, q) in modes:
        eigval = (p*np.pi)**2 + (q*np.pi)**2
        f += (eigval + k2) * np.sin(p*np.pi*x) * np.sin(q*np.pi*y)
    return f
U3, info3 = gmres_helmholtz(n, f_multimode, bc, k2=k2, tol=1e-10, restart=30, return_history=True)
print(f"  Iterations: {info3['iterations']}")
print(f"  Number of residuals: {len(info3['residuals'])}")
print(f"  First 5 residuals: {info3['residuals'][:5]}")
print(f"  Last 5 residuals: {info3['residuals'][-5:]}")
print(f"  Better but still fast: GMRES only needs ~5 iterations for 5 eigenmodes")
print()

# DIAGNOSIS 4: Random RHS (true test of GMRES)
print("=" * 60)
print("DIAGNOSIS 4: Random RHS (excites ALL eigenvectors)")
print("=" * 60)
A, grid_info = build_helmholtz_matrix(n, k2, 'dirichlet')
np.random.seed(42)
b_random = np.random.randn(A.shape[0])
x_sol, info4 = gmres(A, b_random, tol=1e-10, restart=30, max_iter=5000, return_history=True)
print(f"  Matrix size: {A.shape[0]}")
print(f"  Iterations: {info4['iterations']}")
print(f"  Success: {info4['success']}")
print(f"  Number of residuals: {len(info4['residuals'])}")
if len(info4['residuals']) > 10:
    print(f"  First 5: {info4['residuals'][:5]}")
    print(f"  Last 5: {info4['residuals'][-5:]}")
print()

# DIAGNOSIS 5: Smooth but non-eigenfunction problem
print("=" * 60)
print("DIAGNOSIS 5: Gaussian RHS f = exp(-50*((x-0.5)^2+(y-0.5)^2))")
print("=" * 60)
def f_gaussian(x, y):
    return np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))
for n_test in [33, 65]:
    U5, info5 = gmres_helmholtz(n_test, f_gaussian, bc, k2=k2, tol=1e-10, 
                                 restart=30, max_iter=5000, return_history=True)
    print(f"  n={n_test}: iterations={info5['iterations']}, success={info5['success']}")
    if len(info5['residuals']) > 10:
        print(f"    First 5: {info5['residuals'][:5]}")
        print(f"    Last 5: {info5['residuals'][-5:]}")
print()

# DIAGNOSIS 6: Larger k2 (near resonance)
print("=" * 60)
print("DIAGNOSIS 6: Different k2 values with Gaussian RHS, n=65")
print("=" * 60)
for k2_test in [0, 10, 50, 100, 500, 1000]:
    U6, info6 = gmres_helmholtz(65, f_gaussian, bc, k2=k2_test, tol=1e-10, 
                                 restart=30, max_iter=10000, return_history=True)
    print(f"  k2={k2_test:6.1f}: iters={info6['iterations']:5d}, success={info6['success']}")
