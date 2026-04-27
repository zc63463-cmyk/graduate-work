#!/usr/bin/env python3
"""
Advanced Thesis Experiments - Enhanced Visualizations for Master's Thesis
========================================================================
Adds: error heatmaps, solution surface plots, condition number analysis,
      FACR(l) parameter study, GMRES restart sensitivity, complex test problems,
      complexity scaling verification.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import sys, os, time, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

_base_dir = os.path.dirname(os.path.abspath(__file__))
_code_dir = os.path.normpath(os.path.join(_base_dir, 'code', 'python'))
sys.path.insert(0, _base_dir)
sys.path.insert(0, _code_dir)

from helmholtz_solver import (
    fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz,
    solve_helmholtz, _helmholtz_test_problem_dirichlet, _helmholtz_test_problem_neumann
)
from gmres_solver import gmres_helmholtz, test_problem_dirichlet, test_problem_neumann

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'thesis', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150
plt.rcParams['savefig.bbox'] = 'tight'

# ============================================================================
# Exp7: Error Heatmaps & Solution Surface
# ============================================================================
def exp7_error_visualization():
    """Generate error heatmaps and solution surface plots."""
    print("=" * 70)
    print("Exp7: Error Heatmaps & Solution Surface Visualization")
    print("=" * 70)

    k2 = 10.0
    n = 65
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_exact(X, Y)

    # Solve with FA (2nd order) and FFT9 (4th order)
    U_fa = solve_helmholtz(n, f_rhs, bc, k2=k2, method='fa')
    U_fft9 = solve_helmholtz(n, f_rhs, bc, k2=k2, method='fft9')

    err_fa = np.abs(U_fa - Ue)
    err_fft9 = np.abs(U_fft9 - Ue)

    # --- Figure 1: Solution + Error heatmaps (2x2) ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    c = ax.pcolormesh(X, Y, Ue, cmap='RdBu_r', shading='auto')
    fig.colorbar(c, ax=ax)
    ax.set_title(r'Exact Solution $u = \sin(\pi x)\sin(\pi y)$')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_aspect('equal')

    ax = axes[0, 1]
    c = ax.pcolormesh(X, Y, U_fa, cmap='RdBu_r', shading='auto')
    fig.colorbar(c, ax=ax)
    ax.set_title(r'FA Solution (2nd order, $n=65$)')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_aspect('equal')

    ax = axes[1, 0]
    c = ax.pcolormesh(X, Y, np.log10(err_fa + 1e-20), cmap='hot_r', shading='auto')
    fig.colorbar(c, ax=ax, label=r'$\log_{10}$ |error|')
    ax.set_title(r'FA Error ($\max$ = {:.2e})'.format(err_fa.max()))
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_aspect('equal')

    ax = axes[1, 1]
    c = ax.pcolormesh(X, Y, np.log10(err_fft9 + 1e-20), cmap='hot_r', shading='auto')
    fig.colorbar(c, ax=ax, label=r'$\log_{10}$ |error|')
    ax.set_title(r'FFT9 Error ($\max$ = {:.2e})'.format(err_fft9.max()))
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_aspect('equal')

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'error_heatmaps.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")

    # --- Figure 2: 3D solution surface ---
    fig = plt.figure(figsize=(14, 5))
    ax = fig.add_subplot(131, projection='3d')
    ax.plot_surface(X, Y, Ue, cmap='RdBu_r', alpha=0.8)
    ax.set_title('Exact Solution')
    ax.set_xlabel('x'); ax.set_ylabel('y')

    ax = fig.add_subplot(132, projection='3d')
    ax.plot_surface(X, Y, U_fa, cmap='RdBu_r', alpha=0.8)
    ax.set_title('FA Solution')
    ax.set_xlabel('x'); ax.set_ylabel('y')

    ax = fig.add_subplot(133, projection='3d')
    ax.plot_surface(X, Y, err_fa * 1e4, cmap='hot', alpha=0.8)
    ax.set_title('FA Error (x1e4)')
    ax.set_xlabel('x'); ax.set_ylabel('y')

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'solution_surface.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")

    # --- Figure 3: Error heatmaps at different grid sizes ---
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    for idx, n in enumerate([17, 33, 65, 129]):
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)
        U = solve_helmholtz(n, f_rhs, bc, k2=k2, method='fa')
        err = np.abs(U - Ue)
        ax = axes[idx]
        c = ax.pcolormesh(X, Y, np.log10(err + 1e-20), cmap='hot_r',
                          shading='auto', vmin=-6, vmax=-2)
        ax.set_title(f'n={n}, max={err.max():.1e}')
        ax.set_aspect('equal')
        ax.set_xlabel('x'); ax.set_ylabel('y')
    fig.colorbar(c, ax=axes[-1], label=r'$\log_{10}$ |error|')
    plt.suptitle(r'FA Error Distribution vs Grid Size ($k^2=10$)', y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'error_grid_refinement.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")


# ============================================================================
# Exp8: Condition Number Analysis
# ============================================================================
def exp8_condition_number():
    """Analyze condition numbers for different k and n."""
    print("\n" + "=" * 70)
    print("Exp8: Condition Number Analysis")
    print("=" * 70)

    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh

    ns = [9, 17, 33, 65]
    k2_values = [0, 1, 10, 100, 1000]

    print(f"\n  {'n':>6s}", end="")
    for k2 in k2_values:
        print(f" | k2={k2:>5.0f}", end="")
    print()
    print("  " + "-" * (8 + 12 * len(k2_values)))

    for n in ns:
        N = n - 2
        h = 1.0 / (n - 1)
        row = f"  {n:6d}"
        for k2 in k2_values:
            # Build tridiagonal T
            diag_T = np.full(N, 4.0/h**2 + k2)
            off_T = np.full(N-1, -1.0/h**2)
            T = np.diag(diag_T) + np.diag(off_T, 1) + np.diag(off_T, -1)

            # Build full 2D matrix (small enough for dense)
            I_N = np.eye(N)
            A = np.kron(I_N, T) + np.kron(T, I_N) - k2 * np.kron(I_N, I_N) + k2 * np.eye(N**2)

            # Condition number (use eigenvalues for symmetric positive definite)
            eigs = np.linalg.eigvalsh(A)
            cond = eigs[-1] / max(eigs[0], 1e-15)
            row += f" | {cond:10.1f}"
        print(row)

    # Plot condition number vs h
    fig, ax = plt.subplots(figsize=(8, 6))
    ns_fine = [9, 17, 33, 65, 129]
    for k2 in [0, 10, 100]:
        conds = []
        for n in ns_fine:
            N = n - 2
            h = 1.0 / (n - 1)
            # Analytical: lambda_max / lambda_min
            lam_min = (2/h**2) * (2 - 2*np.cos(np.pi/(N+1))) + k2
            lam_max = (2/h**2) * (2 + 2) + k2  # = 8/h^2 + k^2
            conds.append(lam_max / max(lam_min, 1e-15))
        hs = [1.0/(n-1) for n in ns_fine]
        ax.loglog(hs, conds, 'o-', label=f'$k^2={k2}$', linewidth=2)

    # Reference O(h^{-2})
    h_ref = np.array(hs)
    ax.loglog(h_ref, h_ref**(-2) * conds[0] * hs[0]**2, 'k--', alpha=0.3, label=r'$O(h^{-2})$')

    ax.set_xlabel('Grid spacing h')
    ax.set_ylabel(r'Condition number $\kappa(A)$')
    ax.set_title('Condition Number vs Grid Spacing')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'condition_number.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Exp9: FACR(l) Parameter Study
# ============================================================================
def exp9_facr_parameter():
    """Study FACR(l) performance for different l values."""
    print("\n" + "=" * 70)
    print("Exp9: FACR(l) Parameter Study")
    print("=" * 70)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
    ns = [33, 65, 129, 257, 513]
    l_values = [0, 1, 2, 3, 4, 5, 6]

    print(f"\n  Computation time (ms):")
    header = f"  {'n':>6s}"
    for l in l_values:
        header += f" | l={l:>2d}"
    print(header)
    print("  " + "-" * (8 + 9 * len(l_values)))

    times = {l: [] for l in l_values}

    for n in ns:
        row = f"  {n:6d}"
        for l in l_values:
            # Warm up
            try:
                _ = facr_helmholtz(n, f_rhs, bc, k2=k2, facr_l=l)
            except:
                pass
            t0 = time.time()
            try:
                U = facr_helmholtz(n, f_rhs, bc, k2=k2, facr_l=l)
                dt = (time.time() - t0) * 1000
            except Exception as e:
                dt = float('nan')
            times[l].append(dt)
            row += f" | {dt:6.1f}"
        print(row)

    # Find optimal l for each n
    print(f"\n  Optimal l for each n:")
    for idx, n in enumerate(ns):
        best_l = min(l_values, key=lambda l: times[l][idx] if not np.isnan(times[l][idx]) else 1e9)
        print(f"    n={n}: l_opt={best_l} (time={times[best_l][idx]:.1f} ms)")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    for l in l_values:
        valid = [(n, t) for n, t in zip(ns, times[l]) if not np.isnan(t)]
        if valid:
            ns_v, ts_v = zip(*valid)
            ax.semilogy(ns_v, ts_v, 'o-', label=f'FACR(l={l})', linewidth=1.5)

    ax.set_xlabel('Grid size n')
    ax.set_ylabel('Computation Time (ms)')
    ax.set_title(r'FACR(l) Performance vs CR Steps ($k^2=10$)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'facr_parameter.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Exp10: GMRES Restart Parameter Sensitivity
# ============================================================================
def exp10_gmres_restart():
    """Study GMRES convergence with different restart parameters."""
    print("\n" + "=" * 70)
    print("Exp10: GMRES Restart Parameter Sensitivity")
    print("=" * 70)

    k2 = 10.0
    n = 65
    restart_values = [10, 20, 30, 50, 100]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 10a: Convergence curves for different restart
    ax = axes[0]
    u_exact, f_rhs, bc = test_problem_dirichlet(k2)
    for restart in restart_values:
        U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-10,
                                   restart=restart, return_history=True)
        if info['residuals'] and len(info['residuals']) > 1:
            res = np.array(info['residuals'][1:])
            res = res / res[0] if res[0] > 0 else res
            iters = np.arange(1, len(res) + 1)
            ax.semilogy(iters, res, label=f'm={restart}', linewidth=1.5)

    ax.set_xlabel('GMRES Iterations')
    ax.set_ylabel(r'Relative Residual $\|r\|/\|r_0\|$')
    ax.set_title(f'GMRES(m) Restart Sensitivity (n={n}, $k^2$={k2})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 10b: Iteration count vs restart for different n
    ax = axes[1]
    for n in [33, 65, 129]:
        u_exact, f_rhs, bc = test_problem_dirichlet(k2)
        iters_list = []
        for restart in restart_values:
            U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-8,
                                       restart=restart)
            iters_list.append(info['iterations'])
        ax.plot(restart_values, iters_list, 'o-', label=f'n={n}', linewidth=2)

    ax.set_xlabel('Restart Parameter m')
    ax.set_ylabel('Total Iterations')
    ax.set_title(f'GMRES(m) Iterations vs Restart ($k^2$={k2})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'gmres_restart.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")


# ============================================================================
# Exp11: Complexity Scaling Verification
# ============================================================================
def exp11_complexity_scaling():
    """Verify theoretical O(N^2 log N) scaling empirically."""
    print("\n" + "=" * 70)
    print("Exp11: Complexity Scaling Verification")
    print("=" * 70)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
    ns = [33, 65, 129, 257, 513, 1025]
    n_runs = 3  # average over runs

    methods = {
        'FA': lambda n: fa_helmholtz(n, f_rhs, bc, k2=k2),
        'CR': lambda n: cr_helmholtz(n, f_rhs, bc, k2=k2),
        'FACR': lambda n: facr_helmholtz(n, f_rhs, bc, k2=k2),
    }

    times = {m: [] for m in methods}

    for n in ns:
        for name, solver_fn in methods.items():
            # Warm up
            try: _ = solver_fn(n)
            except: pass
            dt_list = []
            for _ in range(n_runs):
                t0 = time.time()
                try:
                    _ = solver_fn(n)
                    dt = (time.time() - t0) * 1000
                except:
                    dt = float('nan')
                dt_list.append(dt)
            times[name].append(np.nanmean(dt_list))

    # Print table
    print(f"\n  {'n':>6s}", end="")
    for name in methods:
        print(f" | {name+'(ms)':>10s}", end="")
    print()
    print("  " + "-" * (8 + 13 * len(methods)))
    for idx, n in enumerate(ns):
        row = f"  {n:6d}"
        for name in methods:
            row += f" | {times[name][idx]:10.1f}"
        print(row)

    # Plot with theoretical reference
    fig, ax = plt.subplots(figsize=(8, 6))
    N_vals = [(n-2)**2 for n in ns]

    for name in methods:
        ax.loglog(N_vals, times[name], 'o-', label=name, linewidth=2)

    # Reference O(N^2 log N)
    N_ref = np.array(N_vals, dtype=float)
    ref_time = times['FA'][0] * (N_ref * np.log(N_ref)) / (N_vals[0]**2 * np.log(N_vals[0]**2))
    ax.loglog(N_ref, ref_time, 'k--', alpha=0.3, label=r'$O(N^2 \log N)$')

    ax.set_xlabel(r'DOFs $N^2$')
    ax.set_ylabel('Computation Time (ms)')
    ax.set_title(r'Complexity Scaling Verification ($k^2=10$)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig_path = os.path.join(FIG_DIR, 'complexity_scaling.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Exp12: Complex Test Problem (non-separable)
# ============================================================================
def exp12_complex_test():
    """Test with a more complex exact solution (multiple frequency components)."""
    print("\n" + "=" * 70)
    print("Exp12: Complex Test Problem (Multiple Frequencies)")
    print("=" * 70)

    k2 = 10.0

    # Complex exact solution with multiple modes
    def u_complex(X, Y):
        return (np.sin(np.pi*X) * np.sin(np.pi*Y)
                + 0.5 * np.sin(2*np.pi*X) * np.sin(3*np.pi*Y)
                + 0.3 * np.sin(5*np.pi*X) * np.sin(4*np.pi*Y))

    def f_complex(X, Y):
        u = u_complex(X, Y)
        laplacian_u = (-(np.pi**2 + (2*np.pi)**2) * np.sin(np.pi*X) * np.sin(np.pi*Y)
                       - 0.5 * ((2*np.pi)**2 + (3*np.pi)**2) * np.sin(2*np.pi*X) * np.sin(3*np.pi*Y)
                       - 0.3 * ((5*np.pi)**2 + (4*np.pi)**2) * np.sin(5*np.pi*X) * np.sin(4*np.pi*Y))
        return -laplacian_u + k2 * u

    def bc_complex(*args):
        return 0.0  # homogeneous Dirichlet

    ns = [9, 17, 33, 65, 129]

    print(f"\n  {'n':>6s} | {'FA Err':>12s} | {'Rate':>6s} | {'FFT9 Err':>12s} | {'Rate':>6s}")
    print("  " + "-" * 60)

    prev_fa = prev_9 = None
    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_complex(X, Y)

        U_fa = fa_helmholtz(n, f_complex, bc_complex, k2=k2)
        U_fft9 = fft9_helmholtz(n, f_complex, bc_complex, k2=k2)

        e_fa = np.max(np.abs(U_fa - Ue))
        e_9 = np.max(np.abs(U_fft9 - Ue))

        r_fa = f"{np.log2(prev_fa/e_fa):6.2f}" if prev_fa else "   --"
        r_9 = f"{np.log2(prev_9/e_9):6.2f}" if prev_9 else "   --"
        prev_fa, prev_9 = e_fa, e_9

        print(f"  {n:6d} | {e_fa:12.2e} | {r_fa} | {e_9:12.2e} | {r_9}")

    # Error heatmap for complex problem
    n = 65
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_complex(X, Y)
    U_fa = fa_helmholtz(n, f_complex, bc_complex, k2=k2)
    U_fft9 = fft9_helmholtz(n, f_complex, bc_complex, k2=k2)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    ax = axes[0]
    c = ax.pcolormesh(X, Y, Ue, cmap='RdBu_r', shading='auto')
    fig.colorbar(c, ax=ax)
    ax.set_title('Exact (complex)')
    ax.set_aspect('equal')

    ax = axes[1]
    c = ax.pcolormesh(X, Y, np.log10(np.abs(U_fa - Ue) + 1e-20), cmap='hot_r', shading='auto')
    fig.colorbar(c, ax=ax, label=r'$\log_{10}$ |err|')
    ax.set_title(f'FA Error (max={np.max(np.abs(U_fa-Ue)):.1e})')
    ax.set_aspect('equal')

    ax = axes[2]
    c = ax.pcolormesh(X, Y, np.log10(np.abs(U_fft9 - Ue) + 1e-20), cmap='hot_r', shading='auto')
    fig.colorbar(c, ax=ax, label=r'$\log_{10}$ |err|')
    ax.set_title(f'FFT9 Error (max={np.max(np.abs(U_fft9-Ue)):.1e})')
    ax.set_aspect('equal')

    plt.suptitle(r'Complex Test Problem ($u = \sin\pi x\sin\pi y + 0.5\sin 2\pi x\sin 3\pi y + 0.3\sin 5\pi x\sin 4\pi y$)', y=1.02)
    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'complex_test_heatmap.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"\n  Figure saved: {fig_path}")


# ============================================================================
# Exp13: Wavenumber k Effect - Comprehensive (with condition number heatmap)
# ============================================================================
def exp13_wavenumber_comprehensive():
    """Comprehensive wavenumber study with accuracy heatmap."""
    print("\n" + "=" * 70)
    print("Exp13: Comprehensive Wavenumber Study")
    print("=" * 70)

    ns = [17, 33, 65, 129]
    k2_values = [0, 1, 5, 10, 20, 50, 100, 200, 500, 1000]

    # FA accuracy matrix
    err_matrix = np.zeros((len(k2_values), len(ns)))
    for i, k2 in enumerate(k2_values):
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
        for j, n in enumerate(ns):
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)
            U = fa_helmholtz(n, f_rhs, bc, k2=k2)
            err_matrix[i, j] = np.max(np.abs(U - Ue))

    # Heatmap of log10(error)
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(ns, k2_values, np.log10(err_matrix + 1e-20),
                       cmap='YlOrRd', shading='auto')
    fig.colorbar(c, ax=ax, label=r'$\log_{10}$ (max error)')
    ax.set_xlabel('Grid size n')
    ax.set_ylabel(r'$k^2$')
    ax.set_title(r'FA Solver Accuracy vs Grid Size and Wavenumber')
    ax.set_yscale('log')

    fig_path = os.path.join(FIG_DIR, 'wavenumber_heatmap.png')
    plt.savefig(fig_path)
    plt.close()
    print(f"  Figure saved: {fig_path}")


# ============================================================================
# Run All Advanced Experiments
# ============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("  Advanced Thesis Experiments - Master's Level")
    print("=" * 70 + "\n")

    exp7_error_visualization()
    exp8_condition_number()
    exp9_facr_parameter()
    exp10_gmres_restart()
    exp11_complexity_scaling()
    exp12_complex_test()
    exp13_wavenumber_comprehensive()

    print("\n" + "=" * 70)
    print("All advanced experiments complete!")
    print(f"Figures saved to {FIG_DIR}")
    print("=" * 70)
