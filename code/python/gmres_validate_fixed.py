#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate the fixed GMRES implementation against scipy."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from scipy.sparse.linalg import gmres as scipy_gmres
from gmres_solver import gmres, build_helmholtz_matrix

print("=" * 70)
print("GMRES Fix Validation: Our implementation vs scipy")
print("=" * 70)

# Test 1: Small dense system (no restart needed)
print("\n--- Test 1: Small random dense system (no restart) ---")
np.random.seed(42)
N = 20
A_dense = np.random.randn(N, N) + 5 * np.eye(N)
b = np.random.randn(N)

x_ours, info = gmres(A_dense, b, tol=1e-12, restart=N, max_iter=N, return_history=True)
x_scipy, flag = scipy_gmres(A_dense, b, rtol=1e-12, restart=N, maxiter=N)

res_ours = np.linalg.norm(b - A_dense @ x_ours)
res_scipy = np.linalg.norm(b - A_dense @ x_scipy)
diff = np.linalg.norm(x_ours - x_scipy)

print(f"  Our GMRES:  success={info['success']}, iters={info['iterations']}, actual_res={res_ours:.2e}")
print(f"  scipy GMRES: converged={flag==0}, actual_res={res_scipy:.2e}")
print(f"  Solution diff: {diff:.2e}")
print(f"  PASS: {diff < 1e-8}")

# Test 2: Small dense system WITH restart
print("\n--- Test 2: Small random dense system (restart m=5) ---")
np.random.seed(123)
N = 30
A_dense2 = np.random.randn(N, N) + 8 * np.eye(N)
b2 = np.random.randn(N)

x_ours2, info2 = gmres(A_dense2, b2, tol=1e-12, restart=5, max_iter=2000, return_history=True)
x_scipy2, flag2 = scipy_gmres(A_dense2, b2, rtol=1e-12, restart=5, maxiter=2000)

res_ours2 = np.linalg.norm(b2 - A_dense2 @ x_ours2)
res_scipy2 = np.linalg.norm(b2 - A_dense2 @ x_scipy2)
diff2 = np.linalg.norm(x_ours2 - x_scipy2)

print(f"  Our GMRES(5):  iters={info2['iterations']}, actual_res={res_ours2:.2e}")
print(f"  scipy GMRES(5): converged={flag2==0}, actual_res={res_scipy2:.2e}")
print(f"  Solution diff: {diff2:.2e}")
print(f"  PASS: {diff2 < 1e-6}")

# Test 3: Helmholtz sparse system (Dirichlet, k2=100)
print("\n--- Test 3: Helmholtz sparse system (n=33, k2=100, restart m=30) ---")
n = 33
k2 = 100.0
A3, grid_info = build_helmholtz_matrix(n, k2, 'dirichlet')
b3 = np.random.randn(A3.shape[0])

x_ours3, info3 = gmres(A3, b3, tol=1e-10, restart=30, max_iter=5000, return_history=True)
x_scipy3, flag3 = scipy_gmres(A3, b3, rtol=1e-10, restart=30, maxiter=5000)

res_ours3 = np.linalg.norm(b3 - A3 @ x_ours3) / np.linalg.norm(b3)
res_scipy3 = np.linalg.norm(b3 - A3 @ x_scipy3) / np.linalg.norm(b3)
diff3 = np.linalg.norm(x_ours3 - x_scipy3) / np.linalg.norm(x_scipy3)

print(f"  Our GMRES(30):  iters={info3['iterations']}, rel_res={res_ours3:.2e}")
print(f"  scipy GMRES(30): converged={flag3==0}, rel_res={res_scipy3:.2e}")
print(f"  Relative solution diff: {diff3:.2e}")
print(f"  PASS: {diff3 < 1e-4}")

# Test 4: Very small system where we can verify exactly
print("\n--- Test 4: Exact verification (4x4 system, no restart) ---")
A4 = np.array([
    [4, -1, 0, 0],
    [-1, 4, -1, 0],
    [0, -1, 4, -1],
    [0, 0, -1, 4]
], dtype=float)
b4 = np.array([1, 2, 3, 4], dtype=float)

x_exact = np.linalg.solve(A4, b4)
x_ours4, info4 = gmres(A4, b4, tol=1e-14, restart=4, max_iter=4, return_history=True)

err4 = np.linalg.norm(x_ours4 - x_exact)
tracked_res = info4['residuals'][-1]
actual_res = np.linalg.norm(b4 - A4 @ x_ours4)

print(f"  Exact solution: {x_exact}")
print(f"  Our solution:   {x_ours4}")
print(f"  Tracked residual: {tracked_res:.2e}")
print(f"  Actual residual:  {actual_res:.2e}")
print(f"  Error vs exact:   {err4:.2e}")
print(f"  PASS: {err4 < 1e-10}")

# Test 5: Helmholtz with GAUSSIAN RHS (realistic test)
print("\n--- Test 5: Gaussian RHS (realistic problem, restart=30) ---")
n = 33
k2 = 100.0
A5, _ = build_helmholtz_matrix(n, k2, 'dirichlet')
h = 1.0 / (n - 1)
x_grid = np.linspace(0, 1, n)
X, Y = np.meshgrid(x_grid[1:-1], x_grid[1:-1], indexing='ij')

def gaussian_rhs(x, y):
    return 100.0 * np.exp(-((x - 0.3)**2 + (y - 0.7)**2) / 0.01)

f_vals = gaussian_rhs(X, Y)
b5 = f_vals.reshape(-1)

x_ours5, info5 = gmres(A5, b5, tol=1e-8, restart=30, max_iter=10000, return_history=True)
x_scipy5, flag5 = scipy_gmres(A5, b5, rtol=1e-8, restart=30, maxiter=10000)

actual_res5 = np.linalg.norm(b5 - A5 @ x_ours5) / np.linalg.norm(b5)
scipy_res5 = np.linalg.norm(b5 - A5 @ x_scipy5) / np.linalg.norm(b5)
diff5 = np.linalg.norm(x_ours5 - x_scipy5) / max(np.linalg.norm(x_scipy5), 1e-15)

print(f"  Our GMRES:  iters={info5['iterations']}, success={info5['success']}, rel_res={actual_res5:.2e}")
print(f"  scipy GMRES: converged={flag5==0}, rel_res={scipy_res5:.2e}")
print(f"  Relative diff: {diff5:.2e}")
print(f"  PASS: {diff5 < 0.01}")

# Summary
print("\n" + "=" * 70)
all_pass = (diff < 1e-8) and (diff2 < 1e-6) and (diff3 < 1e-4) and (err4 < 1e-10) and (diff5 < 0.01)
print(f"OVERALL: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")
print("=" * 70)
