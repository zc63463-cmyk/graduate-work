#!/usr/bin/env python3
"""Experiment 05: True Helmholtz Near-Resonance Scan

Systematically scan sigma near discrete Laplacian eigenvalues to demonstrate:
- min_denominator decreases as sigma approaches -lambda_{p,q}^h
- solution norm amplifies near resonance
- GMRES iterations increase or convergence slows near resonance

Usage (from project root):
    cd c:/Users/20564/Desktop/Graduate/论文收集
    python -m code.python.experiments.exp05_true_helmholtz_resonance
"""
import numpy as np
import pandas as pd
import os, sys, time, warnings

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from helmholtz_solver import fa_helmholtz
from gmres_solver import gmres_helmholtz, build_helmholtz_matrix
from .utils import (
    test_problem_gaussian_rhs, equation_type, get_results_dir, get_figures_dir
)


def compute_discrete_eigenvalue(n, p, q, sx=1.0, sy=1.0):
    """Compute discrete 5-point Laplacian eigenvalue lambda_{p,q}^h.

    lambda_{p,q}^h = (4/h^2) * [sin^2(p*pi/(2*N+2)) + sin^2(q*pi/(2*N+2))]
    where N = n-2 interior points.
    """
    N = n - 2
    h = sx / (n - 1)
    mu_p = (4.0 / (h * h)) * np.sin(p * np.pi / (2 * (N + 1)))**2
    mu_q = (4.0 / (h * h)) * np.sin(q * np.pi / (2 * (N + 1)))**2
    return mu_p + mu_q


def run():
    """Run experiment 05: near-resonance scan."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    rows = []
    
    n = 65
    N = n - 2
    sx, sy = 1.0, 1.0
    h = sx / (n - 1)
    
    # Target mode for resonance
    p_target, q_target = 2, 1
    lambda_target = compute_discrete_eigenvalue(n, p_target, q_target, sx, sy)
    print(f"n={n}, target mode ({p_target},{q_target}), "
          f"lambda_{{p,q}}^h = {lambda_target:.6f}")
    
    # Delta scan: sigma = -(lambda_target + delta)
    # Negative delta: sigma > -lambda (slightly away from resonance)
    # Positive delta: sigma < -lambda (past resonance)
    deltas = [-1e-1, -5e-2, -1e-2, -5e-3, -1e-3, -5e-4, -1e-4,
              1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
    
    # Gaussian RHS (non-eigenfunction)
    _, f_gaussian, bc_zero = test_problem_gaussian_rhs(0.0)
    def bc_func_zero(x, y):
        return 0.0
    
    print("\n" + "=" * 80)
    print(f"Near-resonance scan: sigma = -(lambda_{{{p_target},{q_target}}}^h + delta)")
    print(f"  lambda_{{{p_target},{q_target}}}^h = {lambda_target:.6f}")
    print("=" * 80)
    
    for delta in deltas:
        sigma = -(lambda_target + delta)
        eq_type = equation_type(sigma)
        
        # Compute min_denominator
        mu = 2.0 - 2.0 * np.cos(np.arange(1, N+1) * np.pi / (N + 1))
        lambda_2d = np.add.outer(mu, mu) / (h * h)
        D = lambda_2d + sigma
        abs_D = np.abs(D)
        min_denom = np.min(abs_D)
        idx_min = np.unravel_index(np.argmin(abs_D), abs_D.shape)
        resonant_p = idx_min[0] + 1
        resonant_q = idx_min[1] + 1
        spectral_cond = np.max(abs_D) / min_denom if min_denom > 1e-15 else np.inf
        
        # FA solve for solution norm
        t0 = time.time()
        try:
            U_fa = fa_helmholtz(n, f_gaussian, bc_func_zero, sigma=sigma,
                                bc_type='dirichlet', sx=sx, sy=sy)
            sol_l2 = np.sqrt(np.mean(U_fa**2))
            sol_linf = np.max(np.abs(U_fa))
            fa_success = True
        except Exception as e:
            sol_l2 = np.nan
            sol_linf = np.nan
            fa_success = False
        t1 = time.time()
        fa_time_ms = (t1 - t0) * 1000
        
        # GMRES solve
        t0 = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                U_gm, info = gmres_helmholtz(n, f_gaussian, bc_func_zero, sigma=sigma,
                                              bc_type='dirichlet', restart=30,
                                              max_iter=1000, tol=1e-10)
            gmres_success = info['success']
            gmres_iters = info['iterations']
            abs_res = info['residuals'][-1] if info['residuals'] else np.inf
        except Exception as e:
            gmres_success = False
            gmres_iters = -1
            abs_res = np.inf
        t1 = time.time()
        gmres_time_ms = (t1 - t0) * 1000
        
        row = {
            'n': n, 'n_interior': N, 'h': h,
            'p_target': p_target, 'q_target': q_target,
            'lambda_h': lambda_target,
            'delta': delta,
            'sigma': sigma,
            'equation_type': eq_type,
            'min_denominator': min_denom,
            'spectral_condition': spectral_cond,
            'resonant_p': resonant_p,
            'resonant_q': resonant_q,
            'solution_l2': sol_l2,
            'solution_linf': sol_linf,
            'gmres_iters': gmres_iters,
            'gmres_success': gmres_success,
            'gmres_abs_residual': abs_res,
            'fa_time_ms': fa_time_ms,
            'gmres_time_ms': gmres_time_ms,
        }
        rows.append(row)
        
        print(f"  delta={delta:+10.1e} sigma={sigma:+10.4f} min_d={min_denom:.2e} "
              f"||u||_inf={sol_linf:.2e} gmres_it={gmres_iters:4d} "
              f"res={abs_res:.2e}")
    
    # Also add a few safe sigma values as baseline comparison
    baseline_sigmas = [-1, -5, -10, -20]
    print("\nBaseline (safe) sigma values:")
    for sigma in baseline_sigmas:
        eq_type = equation_type(sigma)
        D = lambda_2d + sigma
        abs_D = np.abs(D)
        min_denom = np.min(abs_D)
        spectral_cond = np.max(abs_D) / min_denom if min_denom > 1e-15 else np.inf
        
        try:
            U_fa = fa_helmholtz(n, f_gaussian, bc_func_zero, sigma=sigma,
                                bc_type='dirichlet', sx=sx, sy=sy)
            sol_l2 = np.sqrt(np.mean(U_fa**2))
            sol_linf = np.max(np.abs(U_fa))
        except:
            sol_l2 = np.nan
            sol_linf = np.nan
        
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                U_gm, info = gmres_helmholtz(n, f_gaussian, bc_func_zero, sigma=sigma,
                                              bc_type='dirichlet', restart=30,
                                              max_iter=1000, tol=1e-10)
            gmres_iters = info['iterations']
            abs_res = info['residuals'][-1] if info['residuals'] else np.inf
        except:
            gmres_iters = -1
            abs_res = np.inf
        
        rows.append({
            'n': n, 'n_interior': N, 'h': h,
            'p_target': p_target, 'q_target': q_target,
            'lambda_h': lambda_target,
            'delta': np.nan,  # not a near-resonance point
            'sigma': sigma,
            'equation_type': eq_type,
            'min_denominator': min_denom,
            'spectral_condition': spectral_cond,
            'resonant_p': np.nan, 'resonant_q': np.nan,
            'solution_l2': sol_l2,
            'solution_linf': sol_linf,
            'gmres_iters': gmres_iters,
            'gmres_success': True,
            'gmres_abs_residual': abs_res,
            'fa_time_ms': np.nan, 'gmres_time_ms': np.nan,
        })
        print(f"  sigma={sigma:+10.1f} min_d={min_denom:.2e} "
              f"||u||_inf={sol_linf:.2e} gmres_it={gmres_iters:4d}")
    
    # Save
    df = pd.DataFrame(rows)
    os.makedirs(results_dir, exist_ok=True)
    csv_path = os.path.join(results_dir, 'exp05_resonance.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")
    
    plot(df, figures_dir)
    return df


def plot(df=None, figures_dir=None):
    """Generate near-resonance plots."""
    if df is None:
        results_dir = get_results_dir()
        csv_path = os.path.join(results_dir, 'exp05_resonance.csv')
        df = pd.read_csv(csv_path)
    if figures_dir is None:
        figures_dir = get_figures_dir()
    os.makedirs(figures_dir, exist_ok=True)
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    # Only plot near-resonance points (where delta is not NaN)
    df_nr = df.dropna(subset=['delta']).sort_values('delta')
    
    if len(df_nr) == 0:
        print("No near-resonance data to plot.")
        return df
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot 1: min_denominator vs |delta|
    ax = axes[0]
    ax.semilogy(np.abs(df_nr['delta']), df_nr['min_denominator'], 'bo-', markersize=6)
    ax.set_xlabel('$|\\delta|$')
    ax.set_ylabel('$d_{\\min}$')
    ax.set_title('Min denominator vs $|\\delta|$')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: solution norm vs |delta|
    ax = axes[1]
    ax.semilogy(np.abs(df_nr['delta']), df_nr['solution_linf'], 'ro-', markersize=6)
    ax.set_xlabel('$|\\delta|$')
    ax.set_ylabel('$\\|u\\|_\\infty$')
    ax.set_title('Solution norm vs $|\\delta|$')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: GMRES iterations vs |delta|
    ax = axes[2]
    df_gmres = df_nr[df_nr['gmres_iters'] > 0]
    if len(df_gmres) > 0:
        ax.plot(np.abs(df_gmres['delta']), df_gmres['gmres_iters'], 'go-', markersize=6)
    ax.set_xlabel('$|\\delta|$')
    ax.set_ylabel('GMRES iterations')
    ax.set_title('GMRES iterations vs $|\\delta|$')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Near-Resonance Scan: $\\sigma = -(\\lambda_{{2,1}}^h + \\delta)$, n={df_nr["n"].iloc[0]}',
                 fontsize=12)
    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'exp05_resonance.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {fig_path}")
    plt.close()
    
    return df


if __name__ == '__main__':
    df = run()
    print("\nexp05 completed.")
