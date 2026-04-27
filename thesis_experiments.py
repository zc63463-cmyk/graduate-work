#!/usr/bin/env python3
"""
Thesis Experiments: Unified Script for Generating All Figures and Tables
========================================================================

Generates all numerical experiment results and plots for the thesis:
  1. Convergence rate verification (2nd-order 5pt, 4th-order FFT9)
  2. Direct solver accuracy comparison (FA/CR/FACR/FFT9)
  3. GMRES convergence behavior analysis
  4. Direct vs iterative solver efficiency comparison
  5. Helmholtz equation with different wavenumbers k
  6. Neumann boundary condition results

All figures saved to thesis/figures/
All tables printed to console (can be copied into LaTeX)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import time
import io

# Fix Windows console encoding
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Add parent directory to path
_base_dir = os.path.dirname(os.path.abspath(__file__))
_code_dir = os.path.normpath(os.path.join(_base_dir, 'code', 'python'))
sys.path.insert(0, _base_dir)
sys.path.insert(0, _code_dir)

from helmholtz_solver import (
    fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz,
    point5_helmholtz, solve_helmholtz,
    _helmholtz_test_problem_dirichlet, _helmholtz_test_problem_neumann
)
from gmres_solver import (
    gmres_helmholtz, test_problem_dirichlet, test_problem_neumann,
    test_problem_polynomial, test_problem_gaussian, test_problem_multimode
)

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'thesis', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Font: try SimHei first for Chinese, fallback to DejaVu Sans
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150
plt.rcParams['savefig.bbox'] = 'tight'

# ============================================================================
# Experiment 1: Convergence Rate Verification
# ============================================================================

def exp1_convergence_rates():
    """Verify 2nd-order (5pt) and 4th-order (FFT9) convergence rates."""
    print("=" * 70)
    print("Exp1: Convergence Rate Verification (Dirichlet BC)")
    print("=" * 70)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
    ns = [9, 17, 33, 65, 129]

    methods = {
        'FA (2nd)': 'fa',
        'CR (2nd)': 'cr',
        'FACR (2nd)': 'facr',
        'FFT9 (4th)': 'fft9',
    }

    results = {m: {'n': [], 'err': [], 'rate': []} for m in methods}

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        for name, method in methods.items():
            U = solve_helmholtz(n, f_rhs, bc, k2=k2, method=method)
            err = np.max(np.abs(U - Ue))
            results[name]['n'].append(n)
            results[name]['err'].append(err)

    # Compute convergence rates
    for name in methods:
        for i in range(1, len(ns)):
            rate = np.log2(results[name]['err'][i-1] / results[name]['err'][i])
            results[name]['rate'].append(rate)

    # Print table
    header = f"{'n':>6s}"
    for name in methods:
        header += f" | {name+' Err':>12s} {name+' Rate':>8s}"
    print(header)
    print("-" * len(header))

    for idx, n in enumerate(ns):
        row = f"{n:6d}"
        for name in methods:
            e = results[name]['err'][idx]
            r = results[name]['rate'][idx-1] if idx > 0 else float('nan')
            r_str = f"{r:6.2f}" if not np.isnan(r) else "   --"
            row += f" | {e:12.2e} {r_str:>8s}"
        print(row)

    # Plot convergence rates
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = {'FA (2nd)': 'o', 'CR (2nd)': 's', 'FACR (2nd)': '^', 'FFT9 (4th)': 'D'}
    for name, mk in markers.items():
        hs = [1.0 / (n - 1) for n in results[name]['n']]
        ax.loglog(hs, results[name]['err'], marker=mk, label=name, linewidth=2)

    # Reference lines
    h_ref = np.array(hs)
    ax.loglog(h_ref, h_ref**2 * results['FA (2nd)']['err'][0] / hs[0]**2,
              'k--', alpha=0.3, label=r'$O(h^2)$')
    ax.loglog(h_ref, h_ref**4 * results['FFT9 (4th)']['err'][0] / hs[0]**4,
              'k:', alpha=0.3, label=r'$O(h^4)$')

    ax.set_xlabel('Grid spacing h')
    ax.set_ylabel(r'Max error $\|u - u_h\|_\infty$')
    ax.set_title(r'Convergence Rates ($k^2=10$, Dirichlet BC)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'convergence_rates.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")

    return results


# ============================================================================
# Experiment 2: Direct Solver Accuracy Comparison
# ============================================================================

def exp2_accuracy_comparison():
    """Compare accuracy of all direct solvers at various grid sizes."""
    print("\n" + "=" * 70)
    print("Exp2: Direct Solver Accuracy Comparison (k^2=10, Dirichlet)")
    print("=" * 70)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
    ns = [9, 17, 33, 65, 129]

    methods = ['fa', 'cr', 'facr', 'fft9']
    method_names = {'fa': 'FA', 'cr': 'CR', 'facr': 'FACR', 'fft9': 'FFT9'}

    print(f"\n  {'n':>6s}", end="")
    for m in methods:
        print(f" | {method_names[m]:>12s}", end="")
    print()
    print("  " + "-" * (8 + 15 * len(methods)))

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)
        print(f"  {n:6d}", end="")
        for m in methods:
            U = solve_helmholtz(n, f_rhs, bc, k2=k2, method=m)
            err = np.max(np.abs(U - Ue))
            print(f" | {err:12.2e}", end="")
        print()


# ============================================================================
# Experiment 3: GMRES Convergence Behavior
# ============================================================================

def exp3_gmres_convergence():
    """Analyze GMRES convergence for different k and grid sizes.

    Uses polynomial test problem (NOT eigenfunction) so that GMRES
    exhibits realistic iterative convergence behavior.
    """
    print("\n" + "=" * 70)
    print("Exp3: GMRES Convergence Behavior Analysis")
    print("  (Using polynomial test problem: u = x^2(1-x)^2 * y^2(1-y)^2)")
    print("=" * 70)

    # 3a: Different k values, fixed grid
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
            res = np.array(info['residuals'][1:])  # Skip initial
            res = res / res[0] if res[0] > 0 else res
            iters = np.arange(1, len(res) + 1)
            ax.semilogy(iters, res, label=f'$k^2$={k2}', linewidth=1.5)

    ax.set_xlabel('GMRES Iterations')
    ax.set_ylabel(r'Relative Residual $\|r\|/\|r_0\|$')
    ax.set_title(f'GMRES Convergence vs Wavenumber (n={n})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3b: Different grid sizes, fixed k
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

    ax.set_xlabel('GMRES Iterations')
    ax.set_ylabel(r'Relative Residual $\|r\|/\|r_0\|$')
    ax.set_title(f'GMRES Convergence vs Grid Size ($k^2$={k2})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'gmres_convergence.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")

    # 3c: GMRES convergence rate verification
    print("\n  GMRES Convergence Rate (Dirichlet, k^2=10, polynomial test):")
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


# ============================================================================
# Experiment 4: Direct vs Iterative Efficiency Comparison
# ============================================================================

def exp4_efficiency_comparison():
    """Compare computation time and accuracy of direct vs iterative solvers."""
    print("\n" + "=" * 70)
    print("Exp4: Direct vs Iterative Efficiency Comparison")
    print("=" * 70)

    k2 = 10.0
    # Direct solvers use eigenfunction test (gives exact discrete solution)
    u_exact_direct, f_rhs_direct, bc_direct = _helmholtz_test_problem_dirichlet(k2)
    # GMRES uses polynomial test (realistic iteration behavior)
    u_exact_gmres, f_rhs_gmres, bc_gmres = test_problem_polynomial(k2)
    ns = [33, 65, 129, 257]

    solvers_direct = [
        ("FA", lambda n: fa_helmholtz(n, f_rhs_direct, bc_direct, k2=k2)),
        ("CR", lambda n: cr_helmholtz(n, f_rhs_direct, bc_direct, k2=k2)),
        ("FACR", lambda n: facr_helmholtz(n, f_rhs_direct, bc_direct, k2=k2)),
        ("FFT9", lambda n: fft9_helmholtz(n, f_rhs_direct, bc_direct, k2=k2)),
    ]
    solvers_iter = [
        ("GMRES(30)", lambda n: gmres_helmholtz(n, f_rhs_gmres, bc_gmres, k2=k2,
                                                  tol=1e-10, restart=30, max_iter=50000)),
    ]

    all_solvers = solvers_direct + solvers_iter

    # Time comparison
    print(f"\n  k^2 = {k2}")
    header = f"  {'n':>6s}"
    for name, _ in all_solvers:
        header += f" | {name+'(ms)':>12s}"
    print(header)
    print("  " + "-" * (8 + 15 * len(all_solvers)))

    times = {name: [] for name, _ in all_solvers}
    errors = {name: [] for name, _ in all_solvers}

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue_direct = u_exact_direct(X, Y)
        Ue_gmres = u_exact_gmres(X, Y)

        row = f"  {n:6d}"
        for name, solver_fn in all_solvers:
            # Warm up
            try:
                _ = solver_fn(n)
            except:
                pass
            # Time it
            t0 = time.time()
            try:
                if name.startswith("GMRES"):
                    U, info = solver_fn(n)
                    Ue = Ue_gmres
                else:
                    U = solver_fn(n)
                    Ue = Ue_direct
                dt = (time.time() - t0) * 1000
                err = np.max(np.abs(U - Ue))
            except Exception as e:
                dt = float('nan')
                err = float('nan')
                print(f"    Error ({name}, n={n}): {e}")

            times[name].append(dt)
            errors[name].append(err)
            row += f" | {dt:12.2f}"
        print(row)

    # Accuracy comparison table
    print(f"\n  Accuracy (Max Error):")
    header = f"  {'n':>6s}"
    for name, _ in all_solvers:
        header += f" | {name:>12s}"
    print(header)
    print("  " + "-" * (8 + 15 * len(all_solvers)))

    for idx, n in enumerate(ns):
        row = f"  {n:6d}"
        for name, _ in all_solvers:
            e = errors[name][idx]
            row += f" | {e:12.2e}"
        print(row)

    # Plot time comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for name, _ in all_solvers:
        ax1.semilogy(ns, times[name], marker='o', label=name, linewidth=2)
    ax1.set_xlabel('Grid size n')
    ax1.set_ylabel('Computation Time (ms)')
    ax1.set_title(f'Time Comparison ($k^2$={k2}, Dirichlet)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    for name, _ in all_solvers:
        hs = [1.0 / (n - 1) for n in ns]
        ax2.loglog(hs, errors[name], marker='o', label=name, linewidth=2)
    ax2.set_xlabel('Grid spacing h')
    ax2.set_ylabel('Max Error')
    ax2.set_title(f'Accuracy Comparison ($k^2$={k2}, Dirichlet)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'efficiency_comparison.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Experiment 5: Different Wavenumbers
# ============================================================================

def exp5_wavenumber_effects():
    """Test solver performance with different Helmholtz wavenumbers k."""
    print("\n" + "=" * 70)
    print("Exp5: Wavenumber Effects on Solver Performance")
    print("=" * 70)

    ns = [33, 65, 129]
    k2_values = [0.0, 1.0, 10.0, 100.0, 1000.0]

    # FA solver accuracy for different k
    print("\n  FA Direct Solver:")
    print(f"  {'k2':>8s}", end="")
    for n in ns:
        print(f" | n={n} Err{' ':>6s}", end="")
    print()
    print("  " + "-" * (10 + 16 * len(ns)))

    for k2 in k2_values:
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
        row = f"  {k2:8.1f}"
        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)
            U = fa_helmholtz(n, f_rhs, bc, k2=k2)
            err = np.max(np.abs(U - Ue))
            row += f" | {err:12.2e}"
        print(row)

    # GMRES iterations for different k (using polynomial test problem)
    print("\n  GMRES Iteration Count (polynomial test problem):")
    print(f"  {'k2':>8s}", end="")
    for n in ns:
        print(f" | n={n} Iters{' ':>2s}", end="")
    print()
    print("  " + "-" * (10 + 14 * len(ns)))

    for k2 in k2_values:
        u_exact, f_rhs, bc = test_problem_polynomial(k2)
        row = f"  {k2:8.1f}"
        for n in ns:
            U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-10,
                                       restart=30, max_iter=50000)
            row += f" | {info['iterations']:12d}"
        print(row)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    n = 65
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

    ax.set_xlabel('GMRES Iterations')
    ax.set_ylabel('Relative Residual')
    ax.set_title(f'GMRES Convergence vs Wavenumber (n={n})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'wavenumber_effects.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Experiment 6: Neumann Boundary Condition
# ============================================================================

def exp6_neumann_bc():
    """Test solvers with Neumann boundary conditions."""
    print("\n" + "=" * 70)
    print("Exp6: Neumann Boundary Condition Test")
    print("=" * 70)

    k2 = 1.0
    ns = [9, 17, 33, 65, 129]

    # FA/CR/FACR for Neumann
    print(f"\n  FFT Direct Methods (Neumann BC, k^2={k2}):")
    print(f"  {'n':>6s} | {'FA Err':>12s} | {'Rate':>6s} | {'CR Err':>12s} | {'Rate':>6s} | {'FACR Err':>12s} | {'Rate':>6s}")
    print("  " + "-" * 85)

    u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)
    prev = {}
    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        U_fa = fa_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
        U_cr = cr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
        U_fc = facr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')

        e_fa = np.max(np.abs(U_fa - Ue))
        e_cr = np.max(np.abs(U_cr - Ue))
        e_fc = np.max(np.abs(U_fc - Ue))

        rates = {}
        for name, err in [('FA', e_fa), ('CR', e_cr), ('FACR', e_fc)]:
            if name in prev:
                rates[name] = np.log2(prev[name] / err)
            else:
                rates[name] = float('nan')
            prev[name] = err

        r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
        r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
        r_fc = f"{rates['FACR']:6.2f}" if not np.isnan(rates['FACR']) else "   --"

        print(f"  {n:6d} | {e_fa:12.2e} | {r_fa} | {e_cr:12.2e} | {r_cr} | {e_fc:12.2e} | {r_fc}")

    # GMRES for Neumann (using polynomial test, NOT eigenfunction)
    print(f"\n  GMRES Iterative (Neumann BC, k^2={k2}, polynomial test):")
    u_exact_poly, f_rhs_poly, bc_poly = test_problem_polynomial(k2)
    print(f"  {'n':>6s} | {'Error':>12s} | {'Rate':>6s} | {'Iters':>6s}")
    print("  " + "-" * 45)
    prev_err = None
    for n in [9, 17, 33, 65]:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact_poly(X, Y)
        U, info = gmres_helmholtz(n, f_rhs_poly, bc_poly, k2=k2,
                                   bc_type='neumann', tol=1e-10,
                                   restart=30, max_iter=50000)
        err = np.max(np.abs(U - Ue))
        rate = np.log2(prev_err / err) if prev_err else float('nan')
        prev_err = err
        r_str = f"{rate:6.2f}" if not np.isnan(rate) else "   --"
        print(f"  {n:6d} | {err:12.2e} | {r_str} | {info['iterations']:6d}")

    # Plot Neumann convergence
    fig, ax = plt.subplots(figsize=(8, 6))
    u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)
    ns_plot = [9, 17, 33, 65, 129]

    err_fa, err_cr, err_fc = [], [], []
    for n in ns_plot:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)
        err_fa.append(np.max(np.abs(fa_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann') - Ue)))
        err_cr.append(np.max(np.abs(cr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann') - Ue)))
        err_fc.append(np.max(np.abs(facr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann') - Ue)))

    hs = [1.0 / (n - 1) for n in ns_plot]
    ax.loglog(hs, err_fa, 'o-', label='FA', linewidth=2)
    ax.loglog(hs, err_cr, 's-', label='CR', linewidth=2)
    ax.loglog(hs, err_fc, '^-', label='FACR', linewidth=2)

    # Reference line O(h^2)
    h_ref = np.array(hs)
    ax.loglog(h_ref, h_ref**2 * err_fa[0] / hs[0]**2, 'k--', alpha=0.3, label=r'$O(h^2)$')

    ax.set_xlabel('Grid spacing h')
    ax.set_ylabel(r'Max error $\|u - u_h\|_\infty$')
    ax.set_title(f'Neumann BC Convergence ($k^2$={k2})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'neumann_convergence.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Run All Experiments
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  Thesis Numerical Experiments - Full Results")
    print("  FFT-based Fast Solvers for Poisson/Helmholtz Equations")
    print("  vs GMRES Iterative Method Comparison")
    print("=" * 70 + "\n")

    results = {}

    results['exp1'] = exp1_convergence_rates()
    exp2_accuracy_comparison()
    exp3_gmres_convergence()
    exp4_efficiency_comparison()
    exp5_wavenumber_effects()
    exp6_neumann_bc()

    print("\n" + "=" * 70)
    print("All experiments complete! Figures saved to thesis/figures/")
    print("=" * 70)
