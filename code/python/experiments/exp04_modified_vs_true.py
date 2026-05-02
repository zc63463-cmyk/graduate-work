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
import os, sys, time, argparse, shutil

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz
from gmres_solver import gmres_helmholtz
from gmres_solver import build_helmholtz_matrix, _build_rhs, gmres
try:
    from .utils import (
        test_problem_dirichlet, test_problem_dirichlet_mode,
        test_problem_gaussian_rhs,
        equation_type, get_results_dir, get_figures_dir
    )
except ImportError:  # Support direct execution with PYTHONPATH=code/python.
    from experiments.utils import (
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

PATCH4_CONDITION_CSV = "exp04_condition_check.csv"
PATCH4_HISTORY_CSV = "exp04_gmres_history.csv"
PATCH4_CONDITION_FIG = "exp04_condition_check"
PATCH4_HISTORY_FIG = "exp04_gmres_history"
PATCH4_HISTORY_SIGMAS = [10, 1000, -10, -50]
PATCH4_HISTORY_N = 129
PATCH4_RESTART = 30
PATCH4_TOL = 1e-10
PATCH4_MAX_ITER = 500


def compute_spectral_indicators(n, sigma, bc_type='dirichlet'):
    """Compute min_denominator and spectral condition indicator.

    For the Dirichlet five-point discretization, DST-I orthogonally
    diagonalizes A = Q.T @ diag(d_pq) @ Q. Hence, for nonsingular A,
    cond_2(A) = max(abs(d_pq)) / min(abs(d_pq)).

    This equivalence relies on orthogonal diagonalization / symmetric normal
    structure. Do not treat this denominator indicator as cond_2(A) for the
    original nonsymmetric Neumann ghost-point matrix, mixed-boundary systems,
    non-orthogonal diagonalizations, or general non-normal matrices.
    """
    
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


def dirichlet_denominator_grid(n, sigma):
    """Return Dirichlet five-point eigenvalues and denominators."""
    N = n - 2
    h = 1.0 / (n - 1)
    modes = np.arange(1, N + 1)
    mu = 2.0 - 2.0 * np.cos(modes * np.pi / (N + 1))
    lambdas = np.add.outer(mu, mu) / h**2
    denominators = lambdas + sigma
    return modes, lambdas, denominators, h


def denominator_summary(n, sigma):
    """Compute spectral denominator summary for the Dirichlet five-point system."""
    modes, lambdas, denominators, h = dirichlet_denominator_grid(n, sigma)
    abs_d = np.abs(denominators)
    idx = np.unravel_index(np.argmin(abs_d), abs_d.shape)
    d_min_abs = float(abs_d[idx])
    d_max_abs = float(np.max(abs_d))
    indicator = d_max_abs / d_min_abs if d_min_abs > 0 else np.inf
    return {
        "N": n - 2,
        "h": h,
        "lambda_min": float(np.min(lambdas)),
        "lambda_max": float(np.max(lambdas)),
        "d_min_abs": d_min_abs,
        "d_max_abs": d_max_abs,
        "nearest_p": int(idx[0] + 1),
        "nearest_q": int(idx[1] + 1),
        "nearest_lambda": float(lambdas[idx]),
        "spectral_indicator": indicator,
    }


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
        from matplotlib.lines import Line2D

        line_handles = []
        capped_handle_added = False
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for idx, eq_type in enumerate(df_plot['equation_type'].unique()):
            sub = df_plot[df_plot['equation_type'] == eq_type]
            n_max = sub['n'].max()
            sub = sub[sub['n'] == n_max].sort_values('sigma')
            if len(sub) > 0:
                color = colors[idx % len(colors)]
                success = sub['gmres_success'].astype(str).str.lower().eq('true')
                capped = (~success) | sub['status'].isin(['failed', 'near_resonance'])

                plt.plot(sub['sigma'], sub['gmres_iters'], '-', color=color, lw=1.8)
                plt.scatter(sub.loc[~capped, 'sigma'], sub.loc[~capped, 'gmres_iters'],
                            marker='o', s=46, color=color, edgecolor='white', zorder=3)
                if capped.any():
                    label = 'capped / not converged' if not capped_handle_added else None
                    plt.scatter(sub.loc[capped, 'sigma'], sub.loc[capped, 'gmres_iters'],
                                marker='x', s=80, color='red', linewidths=2.0,
                                zorder=4, label=label)
                    capped_handle_added = True

                line_handles.append(Line2D([0], [0], color=color, lw=1.8,
                                           label=f'{eq_type} (n={n_max})'))
        cap_line = plt.axhline(501, color='0.35', lw=1.0, ls='--',
                               label='501 = capped under max_iter=500')
        plt.xscale('symlog')
        plt.xlabel('sigma')
        plt.ylabel('GMRES iterations')
        status_handles = [
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='0.3',
                   markeredgecolor='white', markersize=7, label='converged'),
            Line2D([0], [0], marker='x', linestyle='None', color='red',
                   markersize=8, markeredgewidth=2.0, label='capped / not converged'),
            cap_line,
        ]
        plt.legend(handles=line_handles + status_handles, fontsize=9)
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.title('exp04: GMRES iterations vs sigma (Gaussian RHS, non-eigenfunction)')
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'exp04_gmres_iters_vs_sigma.png'), dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Figures saved to: {figures_dir}")
    return df


def _copy_to_thesis(figures_dir, filenames):
    thesis_figures = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "thesis", "figures")
    )
    os.makedirs(thesis_figures, exist_ok=True)
    copied = []
    for name in filenames:
        src = os.path.join(figures_dir, name)
        dst = os.path.join(thesis_figures, name)
        shutil.copy2(src, dst)
        copied.append(dst)
    return copied


def run_condition_check():
    """Dense cond_2 check for the Dirichlet five-point spectral indicator.

    This check verifies a special orthogonally diagonalized symmetric system;
    it is not a claim about Neumann ghost-point, mixed-boundary, or non-normal
    matrices.
    """
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    rows = []
    n_list_dense = [17, 33]
    sigma_list = [-50, -10, -5, -1, 1, 10, 100, 1000]
    for n in n_list_dense:
        for sigma in sigma_list:
            A, _ = build_helmholtz_matrix(n, sigma=sigma, bc_type="dirichlet")
            dense_cond2 = float(np.linalg.cond(A.toarray(), 2))
            summary = denominator_summary(n, sigma)
            spectral_indicator = summary["spectral_indicator"]
            rel_diff = abs(dense_cond2 - spectral_indicator) / dense_cond2
            rows.append({
                "n": n,
                "N": summary["N"],
                "h": summary["h"],
                "sigma": sigma,
                "equation_type": equation_type(sigma),
                "lambda_min": summary["lambda_min"],
                "lambda_max": summary["lambda_max"],
                "d_min_abs": summary["d_min_abs"],
                "d_max_abs": summary["d_max_abs"],
                "nearest_p": summary["nearest_p"],
                "nearest_q": summary["nearest_q"],
                "nearest_lambda": summary["nearest_lambda"],
                "spectral_indicator": spectral_indicator,
                "dense_cond2": dense_cond2,
                "relative_difference": rel_diff,
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, PATCH4_CONDITION_CSV)
    df.to_csv(csv_path, index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    markers = {17: "o", 33: "s"}
    for n in n_list_dense:
        sub = df[df["n"] == n]
        ax.scatter(
            sub["spectral_indicator"],
            sub["dense_cond2"],
            marker=markers[n],
            s=48,
            alpha=0.85,
            label=f"n={n}",
        )
    lo = min(df["spectral_indicator"].min(), df["dense_cond2"].min())
    hi = max(df["spectral_indicator"].max(), df["dense_cond2"].max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="y=x")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("spectral denominator indicator")
    ax.set_ylabel(r"dense $\mathrm{cond}_2(A)$")
    ax.set_title("exp04 condition check: DST denominator vs dense cond2")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.legend(fontsize=9)
    fig.tight_layout()

    png = os.path.join(figures_dir, f"{PATCH4_CONDITION_FIG}.png")
    pdf = os.path.join(figures_dir, f"{PATCH4_CONDITION_FIG}.pdf")
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _copy_to_thesis(figures_dir, [f"{PATCH4_CONDITION_FIG}.png", f"{PATCH4_CONDITION_FIG}.pdf"])
    print(f"Condition check CSV saved to: {csv_path}")
    print(f"Condition check figures saved to: {png}, {pdf}")
    return df


def _filter_gmres_residuals(residuals, restart, reported_iterations, max_iter):
    """Remove restart-boundary recomputed residuals from gmres() history."""
    filtered = []
    for idx, value in enumerate(residuals):
        if idx > 0 and idx % (restart + 1) == 0:
            continue
        filtered.append(float(value))

    max_reported = min(reported_iterations, max_iter)
    if len(filtered) > max_reported + 1:
        filtered = filtered[:max_reported + 1]

    rows = [(i, value) for i, value in enumerate(filtered)]
    if reported_iterations > max_iter and rows[-1][0] != reported_iterations:
        rows.append((reported_iterations, rows[-1][1]))
    return rows


def run_gmres_history():
    """Generate per-iteration GMRES residual histories for representative sigma."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    _, f_gaussian, bc_zero = test_problem_gaussian_rhs(0.0)
    rows = []
    for sigma in PATCH4_HISTORY_SIGMAS:
        A, grid_info = build_helmholtz_matrix(
            PATCH4_HISTORY_N,
            sigma=sigma,
            bc_type="dirichlet",
        )
        b = _build_rhs(
            PATCH4_HISTORY_N,
            f_gaussian,
            bc_zero,
            sigma,
            grid_info["bc_x"],
            grid_info["bc_y"],
        )
        _, info = gmres(
            A,
            b,
            tol=PATCH4_TOL,
            max_iter=PATCH4_MAX_ITER,
            restart=PATCH4_RESTART,
            return_history=True,
        )
        summary = denominator_summary(PATCH4_HISTORY_N, sigma)
        success = bool(info["success"])
        if success:
            status = "ok"
        elif sigma < 0 and summary["d_min_abs"] < 1.0:
            status = "near_resonance"
        else:
            status = "failed"

        norm_b = np.linalg.norm(b)
        residual_rows = _filter_gmres_residuals(
            info["residuals"],
            PATCH4_RESTART,
            int(info["iterations"]),
            PATCH4_MAX_ITER,
        )
        final_abs = float(info["final_abs_residual"])
        final_rel = float(info["final_relative_residual"])
        for iteration, abs_residual in residual_rows:
            rel_residual = abs_residual / norm_b if norm_b > 0 else np.inf
            rows.append({
                "sigma": sigma,
                "equation_type": equation_type(sigma),
                "n": PATCH4_HISTORY_N,
                "N": PATCH4_HISTORY_N - 2,
                "iteration": iteration,
                "restart": PATCH4_RESTART,
                "restart_cycle": iteration // PATCH4_RESTART,
                "abs_residual": abs_residual,
                "rel_residual": rel_residual,
                "gmres_success": success,
                "status": status,
                "tol": PATCH4_TOL,
                "tol_type": "absolute_residual",
                "max_iter": PATCH4_MAX_ITER,
                "final_abs_residual": final_abs,
                "final_rel_residual": final_rel,
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, PATCH4_HISTORY_CSV)
    df.to_csv(csv_path, index=False)

    plot_gmres_history(df, figures_dir)
    print(f"GMRES history CSV saved to: {csv_path}")
    return df


def plot_gmres_history(df, figures_dir):
    """Plot absolute and relative residual history for exp04."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    colors = {
        10: "#1f77b4",
        1000: "#2ca02c",
        -10: "#ff7f0e",
        -50: "#d62728",
    }

    for sigma in PATCH4_HISTORY_SIGMAS:
        sub = df[df["sigma"] == sigma].sort_values("iteration")
        if sub.empty:
            continue
        success = bool(sub["gmres_success"].iloc[-1])
        status = str(sub["status"].iloc[-1])
        label = f"sigma={sigma:+g}"
        if not success:
            label += f" ({status})"
        color = colors[sigma]
        axes[0].plot(sub["iteration"], sub["abs_residual"], color=color, lw=1.8, label=label)
        axes[1].plot(sub["iteration"], sub["rel_residual"], color=color, lw=1.8, label=label)
        if not success:
            for ax, col in [(axes[0], "abs_residual"), (axes[1], "rel_residual")]:
                ax.scatter(
                    sub["iteration"].iloc[-1],
                    sub[col].iloc[-1],
                    marker="x",
                    s=70,
                    color="red",
                    linewidths=2.0,
                    zorder=4,
                )

    axes[0].axhline(PATCH4_TOL, color="0.2", lw=1.1, ls="--", label="absolute tolerance")
    axes[0].set_title("(a) Absolute residual history")
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel("absolute residual")
    axes[1].set_title("(b) Relative residual history")
    axes[1].set_xlabel("iteration")
    axes[1].set_ylabel("relative residual")
    for ax in axes:
        ax.set_yscale("log")
        ax.grid(True, which="both", ls="--", alpha=0.35)
    axes[1].legend(fontsize=8, loc="best")
    axes[0].legend(fontsize=8, loc="best")

    capped_handle = Line2D([0], [0], marker="x", linestyle="None", color="red",
                           markersize=8, markeredgewidth=2.0,
                           label="capped / not converged marker")
    handles, labels = axes[1].get_legend_handles_labels()
    if "capped / not converged marker" not in labels:
        handles.append(capped_handle)
        labels.append("capped / not converged marker")
        axes[1].legend(handles, labels, fontsize=8, loc="best")

    png = os.path.join(figures_dir, f"{PATCH4_HISTORY_FIG}.png")
    pdf = os.path.join(figures_dir, f"{PATCH4_HISTORY_FIG}.pdf")
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _copy_to_thesis(figures_dir, [f"{PATCH4_HISTORY_FIG}.png", f"{PATCH4_HISTORY_FIG}.pdf"])
    print(f"GMRES history figures saved to: {png}, {pdf}")


def run_patch4_outputs():
    """Generate Patch4 supplemental condition and residual-history outputs."""
    condition_df = run_condition_check()
    history_df = run_gmres_history()
    return condition_df, history_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run exp04 or Patch4 supplemental outputs.")
    parser.add_argument(
        "--regenerate-main",
        action="store_true",
        help="Also regenerate historical exp04_modified_vs_true.csv and legacy figures.",
    )
    args = parser.parse_args()
    if args.regenerate_main:
        df = run()
        print("\nexp04 main experiment completed.")
    run_patch4_outputs()
    print("\nexp04 Patch4 outputs completed.")
