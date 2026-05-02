#!/usr/bin/env python3
"""Experiment 05: True Helmholtz Near-Resonance Scan

Systematically scan sigma near discrete Laplacian eigenvalues to demonstrate:
- min_denominator decreases as sigma approaches -lambda_{p,q}^h
- solution norm amplifies near resonance
- GMRES iterations increase or convergence slows near resonance

Usage (from project root):
    cd code/python
    python -m experiments.exp05_true_helmholtz_resonance
"""
import numpy as np
import pandas as pd
import os, sys, time, warnings, shutil
from pathlib import Path
from scipy.fft import dstn, idstn

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from helmholtz_solver import fa_helmholtz
from gmres_solver import gmres_helmholtz, build_helmholtz_matrix, _build_rhs, gmres
try:
    from .utils import (
        test_problem_gaussian_rhs, equation_type, get_results_dir, get_figures_dir
    )
except ImportError:  # Support direct execution with PYTHONPATH=code/python.
    from experiments.utils import (
        test_problem_gaussian_rhs, equation_type, get_results_dir, get_figures_dir
    )


PATCH5_MULTIMODE_CSV = "exp05_multimode_resonance.csv"
PATCH5_GMRES_HISTORY_CSV = "exp05_resonance_gmres_history.csv"
PATCH5_SUMMARY_FIG = "exp05_multimode_resonance_summary"
PATCH5_PROJECTION_FIG = "exp05_dominant_mode_projection"
PATCH5_HISTORY_FIG = "exp05_resonance_gmres_history"
PATCH5_N = 65
PATCH5_DELTAS = [1e-1, 1e-2, 1e-3, 1e-4]
PATCH5_HISTORY_DELTAS = [1e-1, 1e-3, 1e-4]
PATCH5_TARGET_SETS = [
    ("mode_11", [(1, 1)]),
    ("mode_23_32", [(2, 3), (3, 2)]),
    ("mode_33", [(3, 3)]),
]
PATCH5_RESTART = 30
PATCH5_TOL = 1e-10
PATCH5_MAX_ITER = 1000


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


def compute_dirichlet_spectrum(n, sx=1.0):
    """Return the full Dirichlet 5-point Laplacian spectrum on an n x n grid."""
    N = n - 2
    h = sx / (n - 1)
    modes = np.arange(1, N + 1)
    mu = (4.0 / (h * h)) * np.sin(modes * np.pi / (2 * (N + 1)))**2
    return np.add.outer(mu, mu), h


def _copy_to_thesis(fig_path):
    repo_root = Path(__file__).resolve().parents[3]
    thesis_figures = repo_root / "thesis" / "figures"
    thesis_figures.mkdir(parents=True, exist_ok=True)
    target = thesis_figures / Path(fig_path).name
    shutil.copy2(fig_path, target)
    return target


def _target_modal_amplitude(U, p, q):
    """Return the orthonormal DST-I coefficient for target Dirichlet mode."""
    from scipy.fft import dst

    interior = U[1:-1, 1:-1]
    coeffs = dst(dst(interior, type=1, norm='ortho', axis=0),
                 type=1, norm='ortho', axis=1)
    return abs(coeffs[p - 1, q - 1])


def _target_modes_string(target_modes):
    return ";".join(f"({p},{q})" for p, q in target_modes)


def _normalized_offcenter_gaussian(n):
    """Return normalized off-center Gaussian RHS and its interior array."""
    x = np.linspace(0.0, 1.0, n)
    y = np.linspace(0.0, 1.0, n)
    X, Y = np.meshgrid(x[1:-1], y[1:-1], indexing="ij")
    raw = np.exp(-40.0 * ((X - 0.37) ** 2 + (Y - 0.41) ** 2))
    norm = np.linalg.norm(raw)
    if norm == 0:
        raise ValueError("Gaussian RHS has zero norm")
    scale = 1.0 / norm

    def f_rhs(xv, yv):
        return scale * np.exp(-40.0 * ((xv - 0.37) ** 2 + (yv - 0.41) ** 2))

    def bc_zero(xv, yv):
        return 0.0

    return f_rhs, bc_zero, raw * scale


def _modal_solution_for_target(n, target_modes, delta):
    """Compute DST modal solution and dominant projection for a target set."""
    lambda_2d, h = compute_dirichlet_spectrum(n)
    f_rhs, bc_zero, f_interior = _normalized_offcenter_gaussian(n)
    f_hat = dstn(f_interior, type=1, norm="ortho")

    target_lambdas = np.array([lambda_2d[p - 1, q - 1] for p, q in target_modes], dtype=float)
    lambda_target = float(target_lambdas[0])
    target_lambda_spread = float(np.max(target_lambdas) - np.min(target_lambdas))
    if target_lambda_spread > 1e-9:
        raise ValueError(f"Target set {target_modes} is not degenerate enough")

    kappa2 = lambda_target + delta
    sigma = -kappa2
    denominators = lambda_2d - kappa2
    abs_d = np.abs(denominators)
    idx_min = np.unravel_index(np.argmin(abs_d), abs_d.shape)

    u_hat_full = f_hat / denominators
    u_hat_dom = np.zeros_like(u_hat_full)
    fhat_target_sq = 0.0
    uhat_target_sq = 0.0
    for p, q in target_modes:
        i, j = p - 1, q - 1
        u_hat_dom[i, j] = u_hat_full[i, j]
        fhat_target_sq += float(np.abs(f_hat[i, j]) ** 2)
        uhat_target_sq += float(np.abs(u_hat_full[i, j]) ** 2)

    u_full = idstn(u_hat_full, type=1, norm="ortho")
    u_dom = idstn(u_hat_dom, type=1, norm="ortho")
    diff = u_full - u_dom

    full_norm = float(np.linalg.norm(u_full))
    dom_norm = float(np.linalg.norm(u_dom))
    target_uhat_norm = float(np.sqrt(uhat_target_sq))
    predicted_target_uhat_norm = float(np.sqrt(fhat_target_sq) / abs(delta))
    modal_rel_error = (
        abs(target_uhat_norm - predicted_target_uhat_norm) / predicted_target_uhat_norm
        if predicted_target_uhat_norm > 0 else np.nan
    )
    total_energy = float(np.sum(np.abs(u_hat_full) ** 2))
    target_energy = float(uhat_target_sq)
    target_energy_fraction = target_energy / total_energy if total_energy > 0 else np.nan
    if full_norm > 0 and dom_norm > 0:
        shape_corr = float(abs(np.vdot(u_full, u_dom)) / (full_norm * dom_norm))
    else:
        shape_corr = np.nan
    dominant_rel_error = float(np.linalg.norm(diff) / full_norm) if full_norm > 0 else np.nan

    return {
        "n": n,
        "N": n - 2,
        "h": h,
        "f_rhs": f_rhs,
        "bc_zero": bc_zero,
        "lambda_2d": lambda_2d,
        "lambda_target": lambda_target,
        "target_lambda_spread": target_lambda_spread,
        "kappa2": kappa2,
        "sigma": sigma,
        "denominators": denominators,
        "d_min_abs": float(abs_d[idx_min]),
        "nearest_p": int(idx_min[0] + 1),
        "nearest_q": int(idx_min[1] + 1),
        "num_near_zero": int(np.sum(abs_d <= max(10.0 * abs(delta), 1e-12))),
        "u_full": u_full,
        "u_dom": u_dom,
        "diff": diff,
        "u_linf": float(np.max(np.abs(u_full))),
        "u_l2": float(np.sqrt(np.mean(u_full**2))),
        "target_uhat_norm": target_uhat_norm,
        "predicted_target_uhat_norm": predicted_target_uhat_norm,
        "modal_amplitude_relative_error": modal_rel_error,
        "target_energy_fraction": target_energy_fraction,
        "shape_correlation": shape_corr,
        "dominant_relative_error": dominant_rel_error,
    }


def _gmres_for_sigma(n, f_rhs, bc_zero, sigma, d_min_abs):
    """Run unpreconditioned GMRES(30) with exp05 settings."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            _, info = gmres_helmholtz(
                n, f_rhs, bc_zero, sigma=sigma, bc_type="dirichlet",
                restart=PATCH5_RESTART, max_iter=PATCH5_MAX_ITER, tol=PATCH5_TOL,
            )
        success = bool(info["success"])
        status = "ok" if success else ("near_resonance" if d_min_abs < 1.0 else "failed")
        return {
            "gmres_iterations": int(info["iterations"]),
            "gmres_success": success,
            "gmres_status": status,
            "gmres_final_abs_residual": float(info["final_abs_residual"]),
            "gmres_final_rel_residual": float(info["final_relative_residual"]),
        }
    except Exception as exc:
        return {
            "gmres_iterations": -1,
            "gmres_success": False,
            "gmres_status": f"error: {exc}",
            "gmres_final_abs_residual": np.inf,
            "gmres_final_rel_residual": np.inf,
        }


def _full_grid_from_interior(interior):
    n = interior.shape[0] + 2
    full = np.zeros((n, n))
    full[1:-1, 1:-1] = interior
    return full


def run():
    """Run experiment 05: near-resonance scan."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    rows = []
    
    n = 65
    N = n - 2
    sx, sy = 1.0, 1.0
    h = sx / (n - 1)
    
    # Representative mode from the degenerate modal pair (2,3)/(3,2).
    p_target, q_target = 2, 3
    lambda_target = compute_discrete_eigenvalue(n, p_target, q_target, sx, sy)
    print(f"n={n}, modal pair ({p_target},{q_target})/({q_target},{p_target}), "
          f"lambda_{{p,q}}^h = {lambda_target:.6f}")
    
    # Delta scan: kappa^2 = lambda_target + delta, sigma = -kappa^2.
    deltas = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
    
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
            target_amp = _target_modal_amplitude(U_fa, p_target, q_target)
            fa_success = True
        except Exception as e:
            sol_l2 = np.nan
            sol_linf = np.nan
            target_amp = np.nan
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
            'target_mode_abs': target_amp,
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
            gmres_success = info.get('success', False)
        except:
            gmres_iters = -1
            abs_res = np.inf
            gmres_success = False
        
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
            'gmres_success': gmres_success,
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
    plot_summary(df, figures_dir, f_gaussian, bc_func_zero)
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
    
    p_target = int(df_nr["p_target"].iloc[0])
    q_target = int(df_nr["q_target"].iloc[0])
    plt.suptitle(f'Near-Resonance Scan: $\\sigma = -(\\lambda_{{{p_target},{q_target}}}^h + \\delta)$, n={df_nr["n"].iloc[0]}',
                 fontsize=12)
    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'exp05_resonance.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {fig_path}")
    thesis_path = _copy_to_thesis(fig_path)
    print(f"Figure copied to {thesis_path}")
    plt.close()
    
    return df


def plot_summary(df=None, figures_dir=None, f_func=None, bc_func=None):
    """Generate a 2x2 near-resonance summary figure for Chapter 6."""
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

    df_nr = df.dropna(subset=['delta']).copy()
    df_nr = df_nr[df_nr['delta'] > 0].sort_values('delta', ascending=False)
    if len(df_nr) == 0:
        print("No positive near-resonance data to plot.")
        return df

    n = int(df_nr['n'].iloc[0])
    p_target = int(df_nr['p_target'].iloc[0])
    q_target = int(df_nr['q_target'].iloc[0])
    deltas = np.asarray(df_nr['delta'], dtype=float)
    solution_l2 = np.asarray(df_nr['solution_l2'], dtype=float)
    target_amp = np.asarray(df_nr['target_mode_abs'], dtype=float)
    gmres_success = df_nr['gmres_success'].map(
        lambda value: str(value).lower() == 'true'
    ).to_numpy()

    if f_func is None or bc_func is None:
        _, f_func, bc_func = test_problem_gaussian_rhs(0.0)

    delta_wide = 1e-1
    delta_near = 1e-4
    lambda_target = float(df_nr['lambda_h'].iloc[0])
    lambda_2d, _ = compute_dirichlet_spectrum(n)
    kappa2_wide = lambda_target + delta_wide
    kappa2_near = lambda_target + delta_near
    min_denom_wide = float(np.min(np.abs(lambda_2d - kappa2_wide)))
    min_denom_near = float(np.min(np.abs(lambda_2d - kappa2_near)))
    sigma_wide = -kappa2_wide
    sigma_near = -kappa2_near
    U_wide = fa_helmholtz(n, f_func, bc_func, sigma=sigma_wide,
                          bc_type='dirichlet', sx=1.0, sy=1.0)
    U_near = fa_helmholtz(n, f_func, bc_func, sigma=sigma_near,
                          bc_type='dirichlet', sx=1.0, sy=1.0)

    max_wide = np.max(np.abs(U_wide))
    max_near = np.max(np.abs(U_near))
    extent = [0.0, 1.0, 0.0, 1.0]
    gmres_abs_residual = np.asarray(df_nr['gmres_abs_residual'], dtype=float)
    gmres_tol = 1e-10

    print(f"Detuning comparison field: delta={delta_wide:.1e}, "
          f"kappa2={kappa2_wide:.6f}, min_denom={min_denom_wide:.6e}")
    print(f"Detuning comparison field: delta={delta_near:.1e}, "
          f"kappa2={kappa2_near:.6f}, "
          f"min_denom={min_denom_near:.6e}")

    plt.rcParams.update({
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'legend.fontsize': 8,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
    })
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.0), constrained_layout=True)

    ax = axes[0, 0]
    ax.loglog(deltas, solution_l2, 'o-', label=r'$\|u\|_2$ (RMS)')
    ax.loglog(deltas, target_amp, 's-', label=fr'$|\hat u_{{{p_target},{q_target}}}|$')
    finite = np.isfinite(target_amp) & (target_amp > 0)
    if np.any(finite):
        c_ref = target_amp[finite][0] * deltas[finite][0]
        ax.loglog(deltas, c_ref / deltas, 'k--', linewidth=1.0,
                  label=r'reference $C/|\delta|$')
    ax.set_xlabel(r'$|\delta|$')
    ax.set_ylabel('amplification')
    ax.set_title('(a) resonance amplification')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8)
    ax.invert_xaxis()

    ax = axes[0, 1]
    if np.any(gmres_success):
        ax.loglog(deltas[gmres_success], gmres_abs_residual[gmres_success],
                  'o-', label='converged')
    if np.any(~gmres_success):
        ax.loglog(deltas[~gmres_success], gmres_abs_residual[~gmres_success],
                  'x', color='tab:red', markersize=8,
                  label='capped / not converged')
    ax.axhline(gmres_tol, color='0.35', linestyle='--', linewidth=1.0,
               label='absolute tolerance')
    ax.set_xlabel(r'$|\delta|$')
    ax.set_ylabel('final absolute residual')
    ax.set_title('(b) unpreconditioned GMRES residual')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8)
    ax.invert_xaxis()

    ax = axes[1, 0]
    im0 = ax.imshow(U_wide.T, origin='lower', extent=extent, aspect='equal',
                    cmap='RdBu_r')
    ax.set_title(f'(c) raw field, delta={delta_wide:.0e}\n'
                 f'max |u|={max_wide:.2e}; min |d_pq|={min_denom_wide:.2e}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    cbar0 = fig.colorbar(im0, ax=ax, shrink=0.88, pad=0.02)
    cbar0.set_label('u')

    ax = axes[1, 1]
    im1 = ax.imshow(U_near.T, origin='lower', extent=extent, aspect='equal',
                    cmap='RdBu_r')
    ax.set_title(f'(d) raw field, delta={delta_near:.0e}\n'
                 f'max |u|={max_near:.2e}; min |d_pq|={min_denom_near:.2e}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    cbar1 = fig.colorbar(im1, ax=ax, shrink=0.88, pad=0.02)
    cbar1.set_label('u')
    fig.suptitle(
        fr'True Helmholtz near-resonance summary: modal pair ({p_target},{q_target})/({q_target},{p_target}), n={n}',
        fontsize=11,
    )

    fig_path = os.path.join(figures_dir, 'exp05_near_resonance_summary.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"Summary figure saved to {fig_path}")
    thesis_path = _copy_to_thesis(fig_path)
    print(f"Summary figure copied to {thesis_path}")
    plt.close(fig)

    return df


def run_multimode_resonance():
    """Run Patch5 multimode near-resonance modal projection scan."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    rows = []
    for mode_label, target_modes in PATCH5_TARGET_SETS:
        for delta in PATCH5_DELTAS:
            modal = _modal_solution_for_target(PATCH5_N, target_modes, delta)
            gmres_result = _gmres_for_sigma(
                PATCH5_N,
                modal["f_rhs"],
                modal["bc_zero"],
                modal["sigma"],
                modal["d_min_abs"],
            )
            rows.append({
                "mode_label": mode_label,
                "target_modes": _target_modes_string(target_modes),
                "n": modal["n"],
                "N": modal["N"],
                "h": modal["h"],
                "delta": delta,
                "kappa2": modal["kappa2"],
                "sigma": modal["sigma"],
                "lambda_target": modal["lambda_target"],
                "target_lambda_spread": modal["target_lambda_spread"],
                "d_min_abs": modal["d_min_abs"],
                "nearest_p": modal["nearest_p"],
                "nearest_q": modal["nearest_q"],
                "num_near_zero": modal["num_near_zero"],
                "u_linf": modal["u_linf"],
                "u_l2": modal["u_l2"],
                "target_uhat_norm": modal["target_uhat_norm"],
                "predicted_target_uhat_norm": modal["predicted_target_uhat_norm"],
                "modal_amplitude_relative_error": modal["modal_amplitude_relative_error"],
                "target_energy_fraction": modal["target_energy_fraction"],
                "shape_correlation": modal["shape_correlation"],
                "dominant_relative_error": modal["dominant_relative_error"],
                **gmres_result,
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, PATCH5_MULTIMODE_CSV)
    df.to_csv(csv_path, index=False)
    print(f"Patch5 multimode CSV saved to {csv_path}")
    plot_multimode_summary(df, figures_dir)
    plot_dominant_mode_projection(figures_dir)
    return df


def plot_multimode_summary(df, figures_dir):
    """Plot Patch5 2x2 multimode near-resonance summary."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)
    mode_labels = [label for label, _ in PATCH5_TARGET_SETS]
    display_labels = {
        "mode_11": r"$(1,1)$",
        "mode_23_32": r"$(2,3)/(3,2)$",
        "mode_33": r"$(3,3)$",
    }
    color_map = {
        "mode_11": "tab:blue",
        "mode_23_32": "tab:green",
        "mode_33": "tab:purple",
    }

    ax = axes[0, 0]
    pair_amp = _target_pair_modal_amplitudes(PATCH5_N, PATCH5_DELTAS)
    ax.loglog(
        pair_amp["delta"], pair_amp["amp_23"], "o-",
        color="tab:green", linewidth=2.0, markersize=6,
        label=r"$|\hat u_{2,3}|$",
    )
    ax.loglog(
        pair_amp["delta"], pair_amp["amp_32"], "s-",
        color="teal", linewidth=2.0, markersize=6,
        label=r"$|\hat u_{3,2}|$",
    )
    c_ref = float(pair_amp["amp_pair_norm"].iloc[0] * pair_amp["delta"].iloc[0])
    ax.loglog(pair_amp["delta"], c_ref / pair_amp["delta"], "k--", lw=1.0,
              label=r"reference $C/|\delta|$")
    for col, color, yoff in [("amp_23", "tab:green", 12), ("amp_32", "teal", -18)]:
        first = pair_amp.iloc[0]
        last = pair_amp.iloc[-1]
        ratio = float(last[col] / first[col])
        ax.annotate(
            f"{ratio:.0f}x",
            xy=(last["delta"], last[col]),
            xytext=(12, yoff),
            textcoords="offset points",
            fontsize=8,
            color=color,
            arrowprops=dict(arrowstyle="->", color=color, lw=0.8),
        )
    ax.set_xlabel(r"$|\delta|$ (approaching resonance $\rightarrow$)")
    ax.set_ylabel("modal coefficient magnitude")
    ax.set_title(r"(a) target pair amplification")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.invert_xaxis()
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    energy_min = []
    for mode_label in mode_labels:
        sub = df[df["mode_label"] == mode_label].sort_values("delta", ascending=False)
        energy_min.append(float(sub["target_energy_fraction"].min()))
        marker = "D" if mode_label == "mode_23_32" else "o"
        ax.semilogx(
            sub["delta"], sub["target_energy_fraction"], marker=marker, linestyle="-",
            linewidth=2.2 if mode_label == "mode_23_32" else 1.5,
            markersize=6 if mode_label == "mode_23_32" else 5,
            color=color_map[mode_label],
            label=display_labels[mode_label],
        )
    ax.set_xlabel(r"$|\delta|$ (approaching resonance $\rightarrow$)")
    ax.set_ylabel("target subspace energy fraction")
    ax.set_title("(b) target energy fraction (near 1)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.invert_xaxis()
    ax.set_ylim(max(0.99, min(energy_min) - 0.002), 1.0005)
    ax.text(
        0.04, 0.10, "all values near 1",
        transform=ax.transAxes, fontsize=8,
        bbox=dict(facecolor="white", edgecolor="0.6", alpha=0.85),
    )
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    corr_min = []
    for mode_label in mode_labels:
        sub = df[df["mode_label"] == mode_label].sort_values("delta", ascending=False)
        corr_min.append(float(sub["shape_correlation"].min()))
        marker = "D" if mode_label == "mode_23_32" else "o"
        ax.semilogx(
            sub["delta"], sub["shape_correlation"], marker=marker, linestyle="-",
            linewidth=2.2 if mode_label == "mode_23_32" else 1.5,
            markersize=6 if mode_label == "mode_23_32" else 5,
            color=color_map[mode_label],
            label=display_labels[mode_label],
        )
    ax.set_xlabel(r"$|\delta|$ (approaching resonance $\rightarrow$)")
    ax.set_ylabel("shape correlation")
    ax.set_title("(c) dominant projection correlation (near 1)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.invert_xaxis()
    ax.set_ylim(max(0.99, min(corr_min) - 0.002), 1.0005)
    ax.text(
        0.04, 0.10, "full solution already\nmatches target subspace",
        transform=ax.transAxes, fontsize=8,
        bbox=dict(facecolor="white", edgecolor="0.6", alpha=0.85),
    )
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    for mode_label in mode_labels:
        sub = df[df["mode_label"] == mode_label].sort_values("delta", ascending=False)
        success = sub["gmres_success"].astype(str).str.lower().eq("true")
        ax.semilogx(
            sub["delta"], sub["gmres_iterations"],
            "-", color=color_map[mode_label], lw=1.5, alpha=0.8,
            label=display_labels[mode_label],
        )
        ax.scatter(
            sub.loc[success, "delta"],
            sub.loc[success, "gmres_iterations"],
            marker="o", color=color_map[mode_label], edgecolor="black",
            linewidth=0.4, s=42, zorder=3,
        )
        ax.scatter(
            sub.loc[~success, "delta"],
            sub.loc[~success, "gmres_iterations"],
            marker="x", color="tab:red", s=70, linewidths=1.8,
            zorder=4,
        )
    ax.set_xlabel(r"$|\delta|$ (approaching resonance $\rightarrow$)")
    ax.set_ylabel("GMRES iterations")
    ax.set_title("(d) GMRES capped iterations (all hit 1001)")
    ax.grid(True, which="both", ls="--", alpha=0.35)
    ax.invert_xaxis()
    ax.set_ylim(990, 1012)
    ax.text(
        0.04, 0.10, "detailed stagnation\nshown in residual history",
        transform=ax.transAxes, fontsize=8,
        bbox=dict(facecolor="white", edgecolor="0.6", alpha=0.85),
    )
    ax.scatter([], [], marker="x", color="tab:red", s=70, linewidths=1.8,
               label="capped / not converged")
    ax.legend(fontsize=8)

    fig.suptitle("True Helmholtz modal amplification and GMRES capped behavior, Gaussian RHS", fontsize=12)
    png = os.path.join(figures_dir, f"{PATCH5_SUMMARY_FIG}.png")
    pdf = os.path.join(figures_dir, f"{PATCH5_SUMMARY_FIG}.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _copy_to_thesis(png)
    _copy_to_thesis(pdf)
    print(f"Patch5 multimode summary saved to {png}, {pdf}")


def plot_dominant_mode_projection(figures_dir):
    """Plot full/dominant/difference fields for (2,3)/(3,2), delta=1e-4."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modal = _modal_solution_for_target(PATCH5_N, [(2, 3), (3, 2)], 1e-4)
    U_full = _full_grid_from_interior(modal["u_full"])
    U_dom = _full_grid_from_interior(modal["u_dom"])
    U_diff = _full_grid_from_interior(modal["diff"])
    vmax = max(np.max(np.abs(U_full)), np.max(np.abs(U_dom)))
    diff_vmax = np.max(np.abs(U_diff))
    extent = [0.0, 1.0, 0.0, 1.0]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), constrained_layout=True)
    panels = [
        (U_full, "full solution", vmax),
        (U_dom, "dominant projection", vmax),
        (U_diff, "full - dominant", diff_vmax),
    ]
    for ax, (field, title, local_vmax) in zip(axes, panels):
        im = ax.imshow(field.T, origin="lower", extent=extent, cmap="RdBu_r",
                       vmin=-local_vmax, vmax=local_vmax, aspect="equal")
        ax.set_title(f"{title}\nmax |u|={np.max(np.abs(field)):.2e}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, shrink=0.88, pad=0.02)
    fig.suptitle("Near-resonance dominant modal projection: target (2,3)/(3,2), delta=1e-4",
                 fontsize=11)
    png = os.path.join(figures_dir, f"{PATCH5_PROJECTION_FIG}.png")
    pdf = os.path.join(figures_dir, f"{PATCH5_PROJECTION_FIG}.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _copy_to_thesis(png)
    _copy_to_thesis(pdf)
    print(f"Patch5 dominant projection saved to {png}, {pdf}")


def _filter_gmres_residuals(residuals, restart, reported_iterations, max_iter):
    """Remove restart-boundary recomputed residuals while preserving capped count."""
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


def run_resonance_gmres_history():
    """Generate Patch5 residual histories for the (2,3)/(3,2) resonance branch."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    f_rhs, bc_zero, _ = _normalized_offcenter_gaussian(PATCH5_N)
    rows = []
    target_modes = [(2, 3), (3, 2)]
    lambda_2d, _ = compute_dirichlet_spectrum(PATCH5_N)
    lambda_target = float(lambda_2d[1, 2])
    spread = abs(float(lambda_2d[1, 2] - lambda_2d[2, 1]))
    if spread > 1e-9:
        raise ValueError("(2,3)/(3,2) target pair is unexpectedly non-degenerate")

    for delta in PATCH5_HISTORY_DELTAS:
        kappa2 = lambda_target + delta
        sigma = -kappa2
        A, grid_info = build_helmholtz_matrix(PATCH5_N, sigma=sigma, bc_type="dirichlet")
        b = _build_rhs(PATCH5_N, f_rhs, bc_zero, sigma, grid_info["bc_x"], grid_info["bc_y"])
        _, info = gmres(
            A, b, tol=PATCH5_TOL, max_iter=PATCH5_MAX_ITER,
            restart=PATCH5_RESTART, return_history=True,
        )
        success = bool(info["success"])
        status = "ok" if success else "near_resonance"
        norm_b = np.linalg.norm(b)
        for iteration, abs_residual in _filter_gmres_residuals(
            info["residuals"], PATCH5_RESTART, int(info["iterations"]), PATCH5_MAX_ITER
        ):
            rows.append({
                "target_modes": _target_modes_string(target_modes),
                "delta": delta,
                "n": PATCH5_N,
                "iteration": iteration,
                "restart": PATCH5_RESTART,
                "restart_cycle": iteration // PATCH5_RESTART,
                "abs_residual": abs_residual,
                "rel_residual": abs_residual / norm_b if norm_b > 0 else np.inf,
                "gmres_success": success,
                "status": status,
                "tol": PATCH5_TOL,
                "tol_type": "absolute_residual",
                "max_iter": PATCH5_MAX_ITER,
                "final_abs_residual": float(info["final_abs_residual"]),
                "final_rel_residual": float(info["final_relative_residual"]),
            })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, PATCH5_GMRES_HISTORY_CSV)
    df.to_csv(csv_path, index=False)
    print(f"Patch5 GMRES history CSV saved to {csv_path}")
    plot_resonance_gmres_history(df, figures_dir)
    return df


def _target_pair_energy_fractions(n, deltas):
    """Return individual modal energy fractions for the (2,3)/(3,2) pair."""
    lambda_2d, _ = compute_dirichlet_spectrum(n)
    _, _, f_interior = _normalized_offcenter_gaussian(n)
    f_hat = dstn(f_interior, type=1, norm="ortho")
    lambda_target = float(lambda_2d[1, 2])

    rows = []
    for delta in deltas:
        denominators = lambda_2d - (lambda_target + delta)
        u_hat = f_hat / denominators
        total_energy = float(np.sum(np.abs(u_hat) ** 2))
        e23 = float(np.abs(u_hat[1, 2]) ** 2 / total_energy)
        e32 = float(np.abs(u_hat[2, 1]) ** 2 / total_energy)
        rows.append({
            "delta": float(delta),
            "mode_23": e23,
            "mode_32": e32,
            "combined": e23 + e32,
        })
    return pd.DataFrame(rows).sort_values("delta", ascending=False)


def _target_pair_modal_amplitudes(n, deltas):
    """Return |u_hat| values for the two modes in the (2,3)/(3,2) pair."""
    lambda_2d, _ = compute_dirichlet_spectrum(n)
    _, _, f_interior = _normalized_offcenter_gaussian(n)
    f_hat = dstn(f_interior, type=1, norm="ortho")
    lambda_target = float(lambda_2d[1, 2])

    rows = []
    base_amp_23 = None
    base_amp_32 = None
    for delta in sorted(deltas, reverse=True):
        denominators = lambda_2d - (lambda_target + delta)
        u_hat = f_hat / denominators
        amp_23 = float(abs(u_hat[1, 2]))
        amp_32 = float(abs(u_hat[2, 1]))
        if base_amp_23 is None:
            base_amp_23 = amp_23
            base_amp_32 = amp_32
        rows.append({
            "delta": float(delta),
            "amp_23": amp_23,
            "amp_32": amp_32,
            "amp_pair_norm": float(np.sqrt(amp_23**2 + amp_32**2)),
            "ratio_23": amp_23 / base_amp_23,
            "ratio_32": amp_32 / base_amp_32,
        })
    return pd.DataFrame(rows)


def plot_resonance_gmres_history(df, figures_dir):
    """Plot Patch5 near-resonance GMRES residual history."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for delta in PATCH5_HISTORY_DELTAS:
        sub = df[df["delta"] == delta].sort_values("iteration")
        success = bool(sub["gmres_success"].iloc[-1])
        label = f"delta={delta:.0e}" + ("" if success else " (capped)")
        axes[0].plot(sub["iteration"], sub["abs_residual"], lw=1.8, label=label)
        axes[1].plot(sub["iteration"], sub["rel_residual"], lw=1.8, label=label)
        if not success:
            axes[0].scatter(sub["iteration"].iloc[-1], sub["abs_residual"].iloc[-1],
                            marker="x", color="tab:red", s=70, linewidths=2.0, zorder=4)
            axes[1].scatter(sub["iteration"].iloc[-1], sub["rel_residual"].iloc[-1],
                            marker="x", color="tab:red", s=70, linewidths=2.0, zorder=4)
    axes[0].axhline(PATCH5_TOL, color="0.25", ls="--", lw=1.1, label="absolute tolerance")
    axes[0].set_title("(a) Absolute residual")
    axes[0].set_ylabel("absolute residual")
    axes[1].set_title("(b) Relative residual")
    axes[1].set_ylabel("relative residual")
    for ax in axes:
        ax.set_xlabel("iteration")
        ax.set_yscale("log")
        ax.grid(True, which="both", ls="--", alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle("Near-resonance GMRES(30) residual history: target (2,3)/(3,2)",
                 fontsize=11)
    png = os.path.join(figures_dir, f"{PATCH5_HISTORY_FIG}.png")
    pdf = os.path.join(figures_dir, f"{PATCH5_HISTORY_FIG}.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    _copy_to_thesis(png)
    _copy_to_thesis(pdf)
    print(f"Patch5 GMRES history figure saved to {png}, {pdf}")


def run_patch5_outputs():
    multimode_df = run_multimode_resonance()
    history_df = run_resonance_gmres_history()
    return multimode_df, history_df


if __name__ == '__main__':
    df = run()
    run_patch5_outputs()
    print("\nexp05 completed.")
