#!/usr/bin/env python3
"""Experiment 04: Modified vs True Helmholtz Comparison

Quantitatively compare modified Helmholtz vs true Helmholtz in terms of:
- spectral denominator min/abs
- spectral condition indicator
- GMRES iteration behavior (with non-eigenfunction Gaussian RHS)

Note: Previous versions used f=sin(2πx)sin(3πy) which is a Laplacian eigenfunction,
causing GMRES to converge in 1 iteration regardless of sigma. This version uses
a Gaussian RHS (non-eigenfunction) so GMRES iterations are representative.

Usage (from project root):
    cd c:/Users/20564/Desktop/Graduate/论文收集
    python -m code.python.experiments.exp04_modified_vs_true
"""
import numpy as np
import pandas as pd
import os, sys, time

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz
from gmres_solver import gmres_helmholtz
from .utils import (
    test_problem_dirichlet, test_problem_dirichlet_mode,
    test_problem_gaussian_rhs,
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


def _gmres_solve(n, f_rhs, bc_func, sigma, bc_type='dirichlet'):
    """Run GMRES and return summary dict."""
    t0 = time.time()
    try:
        U, info = gmres_helmholtz(n, f_rhs, bc_func, sigma=sigma, bc_type=bc_type,
                                      restart=30, max_iter=500, tol=1e-10)
        gmres_success = info['success']
        gmres_iters = info['iterations']
        residual = info['residuals'][-1] if info['residuals'] else np.inf
        # Compute relative residual using RHS norm from the sparse matrix build
        from gmres_solver import build_helmholtz_matrix, _build_rhs
        A_mat, grid_info = build_helmholtz_matrix(n, sigma=sigma, bc_type=bc_type)
        bc_x, bc_y = grid_info['bc_x'], grid_info['bc_y']
        b_vec = _build_rhs(n, f_rhs, bc_func, sigma, bc_x, bc_y)
        norm_b = np.linalg.norm(b_vec)
        rel_residual = residual / norm_b if norm_b > 0 else np.inf
        status = 'ok' if gmres_success else 'failed'
    except Exception as e:
        gmres_success = False
        gmres_iters = -1
        residual = np.inf
        rel_residual = np.inf
        status = f'error: {e}'
    t1 = time.time()
    solve_time_ms = (t1 - t0) * 1000
    
    return {
        'gmres_success': gmres_success,
        'gmres_iters': gmres_iters,
        'final_abs_residual': residual,
        'final_relative_residual': rel_residual,
        'solve_time_ms': solve_time_ms,
        'status': status,
    }


def run():
    rows = []
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    
    # Use Gaussian RHS (non-eigenfunction) for GMRES iterations
    _, f_gaussian, bc_zero = test_problem_gaussian_rhs(0.0)  # sigma unused for RHS shape
    
    def bc_func_zero(x, y):
        return 0.0
    
    # ------------------------------------------------------------------
    # Modified Helmholtz
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Modified Helmholtz (sigma > 0)")
    print("=" * 70)
    for n in n_list:
        for sigma in sigma_modified:
            eq_type = equation_type(sigma)
            min_d, si, D = compute_spectral_indicators(n, sigma)
            
            # GMRES with Gaussian RHS (non-eigenfunction)
            gmres_result = _gmres_solve(n, f_gaussian, bc_func_zero, sigma, bc_type='dirichlet')
            
            rows.append({
                'sigma': sigma,
                'equation_type': eq_type,
                'rhs_type': 'gaussian',
                'is_eigenmode_rhs': False,
                'n': n,
                'dof': (n-2)**2,
                'min_denominator': min_d,
                'spectral_condition_indicator': si,
                'gmres_iters': gmres_result['gmres_iters'],
                'gmres_success': gmres_result['gmres_success'],
                'final_abs_residual': gmres_result['final_abs_residual'],
                'final_relative_residual': gmres_result['final_relative_residual'],
                'solve_time_ms': gmres_result['solve_time_ms'],
                'status': gmres_result['status'],
            })
            print(f"  sigma={sigma:+.1f} n={n:3d} iters={gmres_result['gmres_iters']:3d} "
                  f"abs_res={gmres_result['final_abs_residual']:.2e} status={gmres_result['status']}")
    
    # ------------------------------------------------------------------
    # True Helmholtz (safe region)
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("True Helmholtz (safe region, sigma < 0)")
    print("=" * 70)
    for n in n_list:
        for sigma in sigma_true_safe:
            eq_type = equation_type(sigma)
            min_d, si, D = compute_spectral_indicators(n, sigma)
            
            # Check resonance
            if min_d < 1e-10:
                rows.append({
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'rhs_type': 'gaussian',
                    'is_eigenmode_rhs': False,
                    'n': n,
                    'dof': (n-2)**2,
                    'min_denominator': min_d,
                    'spectral_condition_indicator': si,
                    'gmres_iters': -1,
                    'gmres_success': False,
                    'final_abs_residual': np.inf,
                    'final_relative_residual': np.inf,
                    'solve_time_ms': 0,
                    'status': 'near_resonance',
                })
                print(f"  sigma={sigma:+.1f} n={n:3d} STATUS: NEAR RESONANCE (min_d={min_d:.2e})")
                continue
            
            gmres_result = _gmres_solve(n, f_gaussian, bc_func_zero, sigma, bc_type='dirichlet')
            
            rows.append({
                'sigma': sigma,
                'equation_type': eq_type,
                'rhs_type': 'gaussian',
                'is_eigenmode_rhs': False,
                'n': n,
                'dof': (n-2)**2,
                'min_denominator': min_d,
                'spectral_condition_indicator': si,
                'gmres_iters': gmres_result['gmres_iters'],
                'gmres_success': gmres_result['gmres_success'],
                'final_abs_residual': gmres_result['final_abs_residual'],
                'final_relative_residual': gmres_result['final_relative_residual'],
                'solve_time_ms': gmres_result['solve_time_ms'],
                'status': gmres_result['status'],
            })
            print(f"  sigma={sigma:+.1f} n={n:3d} iters={gmres_result['gmres_iters']:3d} "
                  f"abs_res={gmres_result['final_abs_residual']:.2e}")
    
    # ------------------------------------------------------------------
    # True Helmholtz (near-resonance)
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("True Helmholtz (near-resonance)")
    print("=" * 70)
    for n in n_list:
        for sigma in sigma_true_near:
            eq_type = equation_type(sigma)
            min_d, si, D = compute_spectral_indicators(n, sigma)
            
            # Find the resonant mode (closest eigenvalue)
            h = 1.0 / (n - 1)
            N = n - 2
            mu = 2.0 - 2.0 * np.cos(np.arange(1, N+1) * np.pi / (N + 1))
            lambda_2d = np.add.outer(mu, mu) / (h * h)
            D_all = lambda_2d + sigma
            abs_D = np.abs(D_all)
            idx_min = np.unravel_index(np.argmin(abs_D), abs_D.shape)
            resonant_p = idx_min[0] + 1  # 1-indexed
            resonant_q = idx_min[1] + 1
            
            if min_d < 1e-12:
                rows.append({
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'rhs_type': 'gaussian',
                    'is_eigenmode_rhs': False,
                    'n': n,
                    'dof': (n-2)**2,
                    'min_denominator': min_d,
                    'spectral_condition_indicator': si,
                    'resonant_p': resonant_p,
                    'resonant_q': resonant_q,
                    'gmres_iters': -1,
                    'gmres_success': False,
                    'final_abs_residual': np.inf,
                    'final_relative_residual': np.inf,
                    'solve_time_ms': 0,
                    'status': 'resonance',
                })
                print(f"  sigma={sigma:+.1f} n={n:3d} STATUS: RESONANCE "
                      f"(min_d={min_d:.2e}, mode=({resonant_p},{resonant_q}))")
                continue
            
            gmres_result = _gmres_solve(n, f_gaussian, bc_func_zero, sigma, bc_type='dirichlet')
            
            rows.append({
                'sigma': sigma,
                'equation_type': eq_type,
                'rhs_type': 'gaussian',
                'is_eigenmode_rhs': False,
                'n': n,
                'dof': (n-2)**2,
                'min_denominator': min_d,
                'spectral_condition_indicator': si,
                'resonant_p': resonant_p,
                'resonant_q': resonant_q,
                'gmres_iters': gmres_result['gmres_iters'],
                'gmres_success': gmres_result['gmres_success'],
                'final_abs_residual': gmres_result['final_abs_residual'],
                'final_relative_residual': gmres_result['final_relative_residual'],
                'solve_time_ms': gmres_result['solve_time_ms'],
                'status': 'near_resonance' if min_d < 1.0 else gmres_result['status'],
            })
            print(f"  sigma={sigma:+.1f} n={n:3d} iters={gmres_result['gmres_iters']:3d} "
                  f"min_d={min_d:.2e} mode=({resonant_p},{resonant_q}) "
                  f"abs_res={gmres_result['final_abs_residual']:.2e}")
    
    # Save results
    df = pd.DataFrame(rows)
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, 'exp04_modified_vs_true.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")
    
    # Generate plots
    plot(df, figures_dir)
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
    plt.title('exp04: min_denominator vs sigma (Gaussian RHS)')
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
    plt.title('exp04: spectral condition indicator vs sigma (Gaussian RHS)')
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'exp04_spectral_indicator_vs_sigma.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Plot 3: GMRES iterations vs sigma (non-eigenfunction RHS)
    # Only plot cases where gmres_iters > 0 (not resonance/excluded)
    df_plot = df[df['gmres_iters'] > 0].copy()
    if len(df_plot) > 0:
        plt.figure(figsize=(10, 6))
        for eq_type in df_plot['equation_type'].unique():
            sub = df_plot[df_plot['equation_type'] == eq_type]
            n_max = sub['n'].max()
            sub = sub[sub['n'] == n_max].sort_values('sigma')
            if len(sub) > 0:
                plt.plot(sub['sigma'], sub['gmres_iters'], 'o-', label=f'{eq_type} (n={n_max})')
        plt.xscale('symlog')
        plt.xlabel('sigma')
        plt.ylabel('GMRES iterations')
        plt.legend()
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.title('exp04: GMRES iterations vs sigma (Gaussian RHS, non-eigenfunction)')
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'exp04_gmres_iters_vs_sigma.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Figures saved to: {figures_dir}")
    return df


if __name__ == '__main__':
    df = run()
    print("\nexp04 completed.")
