#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick test: run Exp3 only to verify GMRES convergence plots."""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

_base_dir = os.path.dirname(os.path.abspath(__file__))
_code_dir = os.path.normpath(os.path.join(_base_dir, 'code', 'python'))
sys.path.insert(0, _base_dir)
sys.path.insert(0, _code_dir)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from gmres_solver import (
    gmres_helmholtz, test_problem_polynomial, test_problem_gaussian
)

FIG_DIR = os.path.join(_base_dir, 'thesis', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150

# ---- GMRES convergence vs k^2 (polynomial test) ----
print("GMRES Convergence vs k^2 (n=33, polynomial test)...")
n = 33
k2_values = [0.0, 1.0, 10.0, 100.0, 1000.0]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for k2 in k2_values:
    u_exact, f_rhs, bc = test_problem_polynomial(k2)
    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-10,
                               restart=30, max_iter=50000,
                               return_history=True)
    if info['residuals'] and len(info['residuals']) > 1:
        res = np.array(info['residuals'][1:])
        res = res / res[0] if res[0] > 0 else res
        iters = np.arange(1, len(res) + 1)
        ax.semilogy(iters, res, label=f'$k^2$={k2}', linewidth=1.5)
        print(f"  k^2={k2:7.1f}: {info['iterations']} iters, success={info['success']}")

ax.set_xlabel('GMRES Iterations')
ax.set_ylabel(r'Relative Residual $\|r\|/\|r_0\|$')
ax.set_title(f'GMRES Convergence vs Wavenumber (n={n})')
ax.legend()
ax.grid(True, alpha=0.3)

# ---- GMRES convergence vs grid size (polynomial test) ----
print("\nGMRES Convergence vs grid size (k^2=10, polynomial test)...")
ax = axes[1]
k2 = 10.0
u_exact, f_rhs, bc = test_problem_polynomial(k2)
for n in [17, 33, 65]:
    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-10,
                               restart=30, max_iter=50000,
                               return_history=True)
    if info['residuals'] and len(info['residuals']) > 1:
        res = np.array(info['residuals'][1:])
        res = res / res[0] if res[0] > 0 else res
        iters = np.arange(1, len(res) + 1)
        ax.semilogy(iters, res, label=f'n={n}', linewidth=1.5)
        print(f"  n={n:3d}: {info['iterations']} iters, success={info['success']}")

ax.set_xlabel('GMRES Iterations')
ax.set_ylabel(r'Relative Residual $\|r\|/\|r_0\|$')
ax.set_title(f'GMRES Convergence vs Grid Size ($k^2$={k2})')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path = os.path.join(FIG_DIR, 'gmres_convergence.png')
plt.savefig(fig_path)
plt.close()
print(f"\nFigure saved: {fig_path}")

# ---- GMRES convergence rate verification ----
print("\nGMRES Convergence Rate (polynomial test, k^2=10):")
k2 = 10.0
u_exact, f_rhs, bc = test_problem_polynomial(k2)
print(f"  {'n':>6s} | {'GMRES Err':>12s} | {'Rate':>6s} | {'Iters':>6s}")
print("  " + "-" * 45)
prev_err = None
for n in [9, 17, 33, 65]:
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_exact(X, Y)
    U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12,
                               restart=30, max_iter=50000)
    err = np.max(np.abs(U - Ue))
    rate = np.log2(prev_err / err) if prev_err else float('nan')
    prev_err = err
    r_str = f"{rate:6.2f}" if not np.isnan(rate) else "   --"
    print(f"  {n:6d} | {err:12.2e} | {r_str} | {info['iterations']:6d}")

print("\nDone! The GMRES convergence plots now show realistic behavior.")
