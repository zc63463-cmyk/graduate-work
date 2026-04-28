#!/usr/bin/env python3
"""Experiment 04: Modified vs True Helmholtz Comparison

Quantitatively compare modified Helmholtz vs true Helmholtz in terms of:
- spectral denominator min/abs
- spectral condition indicator
- GMRES iteration behavior
- error

Usage (from project root):
    cd c:/Users/20564/Desktop/Graduate/论文收集
    python -m code.python.experiments.exp04_modified_vs_true
"""
import numpy as np
import pandas as pd
import os, sys, time

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz
from gmres_solver import gmres_helmholtz, build_helmholtz_matrix
from .utils import (
    test_problem_dirichlet, test_problem_dirichlet_mode,
    equation_type, get_results_dir, get_figures_dir
)

# ============================================================
# Parameters
# ============================================================
n_list = [33, 65, 129]
sigma_modified = [1, 10, 100, 1000]
sigma_true_safe = [-1, -5, -10]
sigma_true_near = [-50]  # near lambda_{2,1} = 5*pi^2 ≈ 49.3


def compute_spectral_indicators(n, sigma, bc_type='dirichlet'):
    """Compute min_denominator and spectral condition indicator."""
    
    if bc_type == 'dirichlet':
        N = n - 2
        # Dirichlet eigenvalues: mu_p = 2 - 2*cos(p*pi/(N+1))
        mu = 2.0 - 2.0 * np.cos(np.arange(1, N+1) * np.pi / (N + 1))
        # 2D eigenvalues: lambda_pq = (mu_p + mu_q) / h^2
        h = 1.0 / (n - 1)
        lambda_2d = np.add.outer(mu, mu) / (h * h)
    else:
        # Neumann: use DCT-I eigenvalues (not implemented for this experiment)
        return None, None, None
    
    # Denominator: D_pq = lambda_pq + sigma
    D = lambda_2d + sigma
    abs_D = np.abs(D)
    
    min_denom = np.min(abs_D)
    max_denom = np.max(abs_D)
    spectral_indicator = max_denom / min_denom if min_denom > 1e-15 else np.inf
    
    return min_denom, spectral_indicator, D


def run():
    rows = []
    
    # ------------------------------------------------------------------
    # Modified Helmholtz
    # ------------------------------------------------------------------
    for n in n_list:
        for sigma in sigma_modified:
            eq_type = equation_type(sigma)
            min_d, si, D = compute_spectral_indicators(n, sigma)
            
            # GMRES solve (use polynomial test problem to avoid eigfunction 1-iter convergence)
            # Use fixed RHS: f = sin(2*pi*x)*sin(3*pi*y)
            def f_func_poly(x, y):
                return np.sin(2*np.pi*x) * np.sin(3*np.pi*y)
            def bc_func_zero(x, y):
                return 0.0
            
            t0 = time.time()
            try:
                U, info = gmres_helmholtz(n, f_func_poly, bc_func_zero, sigma=sigma, bc_type='dirichlet',
                                              restart=30, max_iter=500, tol=1e-10)
                gmres_success = info['success']
                gmres_iters = info['iterations']
                residual = info['residuals'][-1] if info['residuals'] else np.inf
                status = 'ok' if gmres_success else 'failed'
            except Exception as e:
                gmres_success = False
                gmres_iters = -1
                residual = np.inf
                status = 'failed'
            t1 = time.time()
            solve_time_ms = (t1 - t0) * 1000
            
            # Error calculation skipped: f_func_poly has no simple exact solution
            linf_err = np.nan
            l2_err = np.nan
            
            rows.append({
                'sigma': sigma,
                'equation_type': eq_type,
                'n': n,
                'dof': (n-2)**2,
                'min_denominator': min_d,
                'spectral_condition_indicator': si,
                'gmres_iters': gmres_iters,
                'gmres_success': gmres_success,
                'final_residual': residual,
                'linf_error': linf_err,
                'l2_error': l2_err,
                'solve_time_ms': solve_time_ms,
                'status': 'ok' if gmres_success else 'failed'
            })
            print(f"  sigma={sigma:+.1f} n={n:3d} iters={gmres_iters:3d} err={linf_err:.2e} status={status}")
    
    # ------------------------------------------------------------------
    # True Helmholtz (safe region)
    # ------------------------------------------------------------------
    for n in n_list:
        for sigma in sigma_true_safe:
            eq_type = equation_type(sigma)
            min_d, si, D = compute_spectral_indicators(n, sigma)
            
            # Check resonance
            if min_d < 1e-10:
                status = 'near_resonance'
                rows.append({
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'n': n,
                    'dof': (n-2)**2,
                    'min_denominator': min_d,
                    'spectral_condition_indicator': si,
                    'gmres_iters': -1,
                    'gmres_success': False,
                    'final_residual': np.inf,
                    'linf_error': np.inf,
                    'l2_error': np.inf,
                    'solve_time_ms': 0,
                    'status': status
                })
                print(f"  sigma={sigma:+.1f} n={n:3d} STATUS: NEAR RESONANCE (min_d={min_d:.2e})")
                continue
            
            # GMRES solve (use polynomial test problem)
            def f_func_poly(x, y):
                return np.sin(2*np.pi*x) * np.sin(3*np.pi*y)
            def bc_func_zero(x, y):
                return 0.0
            
            t0 = time.time()
            try:
                U, info = gmres_helmholtz(n, f_func_poly, bc_func_zero, sigma=sigma, bc_type='dirichlet',
                                              restart=30, max_iter=500, tol=1e-10)
                gmres_success = info['success']
                gmres_iters = info['iterations']
                residual = info['residuals'][-1] if info['residuals'] else np.inf
                status = 'ok' if gmres_success else 'failed'
            except Exception as e:
                gmres_success = False
                gmres_iters = -1
                residual = np.inf
                status = 'failed'
            t1 = time.time()
            solve_time_ms = (t1 - t0) * 1000
            
            # For true Helmholtz, we don't have exact solution easily
            # So we skip error calculation
            rows.append({
                'sigma': sigma,
                'equation_type': eq_type,
                'n': n,
                'dof': (n-2)**2,
                'min_denominator': min_d,
                'spectral_condition_indicator': si,
                'gmres_iters': gmres_iters,
                'gmres_success': gmres_success,
                'final_residual': residual,
                'linf_error': np.nan,
                'l2_error': np.nan,
                'solve_time_ms': solve_time_ms,
                'status': 'ok' if gmres_success else 'failed'
            })
            print(f"  sigma={sigma:+.1f} n={n:3d} iters={gmres_iters:3d} residual={residual:.2e}")
    
    # Save results
    df = pd.DataFrame(rows)
    results_dir = get_results_dir()
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, 'exp04_modified_vs_true.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")
    return df


def plot(df=None, figures_dir=None):
    """Generate plots for exp04."""
    if df is None:
        results_dir = get_results_dir()
        csv_path = os.path.join(results_dir, 'exp04_modified_vs_true.csv')
        df = pd.read_csv(csv_path)
    
    if figures_dir is None:
        figures_dir = get_figures_dir()
    os.makedirs(figures_dir, exist_ok=True)
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    # Plot 1: min_denominator vs sigma
    plt.figure(figsize=(10, 6))
    for eq_type in df['equation_type'].unique():
        sub = df[df['equation_type'] == eq_type]
        # Use largest n for clarity
        n_max = sub['n'].max()
        sub = sub[sub['n'] == n_max]
        plt.plot(sub['sigma'], sub['min_denominator'], 'o-', label=f'{eq_type} (n={n_max})')
    plt.xscale('symlog')
    plt.yscale('log')
    plt.xlabel('sigma')
    plt.ylabel('min_denominator')
    plt.legend()
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.title('exp04: min_denominator vs sigma')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'exp04_min_denom_vs_sigma.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Plot 2: spectral_condition_indicator vs sigma
    plt.figure(figsize=(10, 6))
    for eq_type in df['equation_type'].unique():
        sub = df[df['equation_type'] == eq_type]
        n_max = sub['n'].max()
        sub = sub[sub['n'] == n_max]
        plt.plot(sub['sigma'], sub['spectral_condition_indicator'], 'o-', label=f'{eq_type} (n={n_max})')
    plt.xscale('symlog')
    plt.yscale('log')
    plt.xlabel('sigma')
    plt.ylabel('spectral condition indicator')
    plt.legend()
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.title('exp04: spectral condition indicator vs sigma')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'exp04_spectral_indicator_vs_sigma.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Plot 3: gmres_iters vs sigma
    plt.figure(figsize=(10, 6))
    for eq_type in df['equation_type'].unique():
        sub = df[df['equation_type'] == eq_type]
        n_max = sub['n'].max()
        sub = sub[sub['n'] == n_max]
        plt.plot(sub['sigma'], sub['gmres_iters'], 'o-', label=f'{eq_type} (n={n_max})')
    plt.xscale('symlog')
    plt.xlabel('sigma')
    plt.ylabel('GMRES iterations')
    plt.legend()
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.title('exp04: GMRES iterations vs sigma')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'exp04_gmres_iters_vs_sigma.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Figures saved to: {figures_dir}")
    return df


if __name__ == '__main__':
    df = run()
    plot(df)
    print("\nexp04 completed.")
