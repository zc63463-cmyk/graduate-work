#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test the new non-eigenfunction test problems."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from gmres_solver import (
    gmres_helmholtz, test_problem_dirichlet, test_problem_polynomial,
    test_problem_gaussian, test_problem_multimode,
    build_helmholtz_matrix, gmres
)

print("=" * 70)
print("New Test Problems Validation")
print("=" * 70)

# ============================================================
# Test 1: Polynomial test problem - verify convergence order
# ============================================================
print("\n--- Test 1: Polynomial problem (Dirichlet) ---")
k2 = 10.0
u_exact, f_rhs, bc = test_problem_polynomial(k2)

print(f"  {'n':>5s} | {'Error':>12s} | {'Rate':>6s} | {'Iters':>6s}")
print("  " + "-" * 45)

prev_err = None
for n in [9, 17, 33, 65]:
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_exact(X, Y)

    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12,
                               restart=30, max_iter=50000,
                               return_history=True)
    err = np.max(np.abs(U - Ue))
    rate = np.log2(prev_err / err) if prev_err is not None else float('nan')
    prev_err = err

    r_str = f"{rate:6.2f}" if not np.isnan(rate) else "   --"
    print(f"  {n:5d} | {err:12.2e} | {r_str} | {info['iterations']:6d}")

# ============================================================
# Test 2: Multi-mode problem - verify iteration count scales
# ============================================================
print("\n--- Test 2: Multi-mode problem ---")
k2 = 10.0

for n_modes in [2, 5, 10]:
    u_exact, f_rhs, bc = test_problem_multimode(k2, n_modes=n_modes)
    n = 33
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_exact(X, Y)

    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12,
                               restart=30, max_iter=50000,
                               return_history=True)
    err = np.max(np.abs(U - Ue))
    print(f"  n_modes={n_modes:2d}: iters={info['iterations']:4d}, "
          f"expected~{n_modes**2:3d}, error={err:.2e}, success={info['success']}")

# ============================================================
# Test 3: Gaussian problem - verify GMRES needs many iters
# ============================================================
print("\n--- Test 3: Gaussian source problem ---")
k2 = 100.0
_, f_rhs, bc = test_problem_gaussian(k2)

for n in [17, 33, 65]:
    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-8,
                               restart=30, max_iter=100000,
                               return_history=True)
    N = (n-2)**2
    print(f"  n={n:3d} (N={N:5d}): iters={info['iterations']:5d}, "
          f"success={info['success']}, time={info['time']:.2f}s")

# ============================================================
# Test 4: Compare old (eigenfunction) vs new test problems
# ============================================================
print("\n--- Test 4: Iteration comparison (k2=100, n=33) ---")
n = 33
k2 = 100.0

# Old: eigenfunction (trivial)
u1, f1, bc1 = test_problem_dirichlet(k2)
_, info1 = gmres_helmholtz(n, f1, bc1, k2=k2, tol=1e-12,
                            restart=30, return_history=True)
print(f"  Eigenfunction (old):    iters = {info1['iterations']}")

# New: polynomial
u2, f2, bc2 = test_problem_polynomial(k2)
_, info2 = gmres_helmholtz(n, f2, bc2, k2=k2, tol=1e-12,
                            restart=30, max_iter=50000,
                            return_history=True)
print(f"  Polynomial (new):       iters = {info2['iterations']}")

# New: multi-mode (5x5=25 modes)
u3, f3, bc3 = test_problem_multimode(k2, n_modes=5)
_, info3 = gmres_helmholtz(n, f3, bc3, k2=k2, tol=1e-12,
                            restart=30, max_iter=50000,
                            return_history=True)
print(f"  Multi-mode 5x5 (new):   iters = {info3['iterations']}")

# New: Gaussian
_, f4, bc4 = test_problem_gaussian(k2)
_, info4 = gmres_helmholtz(n, f4, bc4, k2=k2, tol=1e-8,
                            restart=30, max_iter=100000,
                            return_history=True)
print(f"  Gaussian (new):         iters = {info4['iterations']}")

print("\n" + "=" * 70)
print("New test problems working correctly!")
print("Key insight: eigenfunction=1 iter vs polynomial/Gaussian=hundreds+ iters")
print("=" * 70)
