"""
Experiment 03: Neumann and mixed boundary verification (V2).

This script keeps the original Experiment 4 scope in Chapter 6:
five-point FA/CR/FACR-like solvers under Neumann and mixed boundary
conditions.  It adds L2, boundary residual, zero-mode, and compatibility
diagnostics, plus clearer summary and field/error visualizations.
"""

import os
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz
from .utils import (
    compute_convergence_rate,
    equation_type,
    get_figures_dir,
    get_results_dir,
    test_problem_mixed_dn,
    test_problem_mixed_nd,
    test_problem_neumann,
    trapezoid_weights,
    weighted_mean_full_grid,
)


SOLVERS = {
    'FA': fa_helmholtz,
    'CR': cr_helmholtz,
    'FACR': facr_helmholtz,
}

METHOD_COLORS = {
    'FA': 'C0',
    'CR': 'C1',
    'FACR': 'C2',
}

NS = [9, 17, 33, 65, 129]
SX = 1.0
SY = 1.0

CASES = [
    {
        'sub_exp': 'a_pure_neumann_sigma1',
        'display': 'Pure Neumann, sigma=1',
        'sigma': 1.0,
        'bc_type': 'neumann',
        'problem': 'neumann',
        'zero_mode': False,
        'plot_in_summary': True,
    },
    {
        'sub_exp': 'a_pure_neumann_sigma10',
        'display': 'Pure Neumann, sigma=10',
        'sigma': 10.0,
        'bc_type': 'neumann',
        'problem': 'neumann',
        'zero_mode': False,
        'plot_in_summary': True,
    },
    {
        'sub_exp': 'b_mixed_ND_sigma10',
        'display': 'Mixed (N,D), sigma=10',
        'sigma': 10.0,
        'bc_type': ('N', 'D'),
        'problem': 'mixed_nd',
        'zero_mode': False,
        'plot_in_summary': True,
    },
    {
        'sub_exp': 'c_poisson_neumann_sigma0',
        'display': 'Neumann Poisson, sigma=0',
        'sigma': 0.0,
        'bc_type': 'neumann',
        'problem': 'neumann',
        'zero_mode': True,
        'plot_in_summary': True,
    },
    {
        'sub_exp': 'd_mixed_DN_sigma10',
        'display': 'Mixed (D,N), sigma=10 smoke check',
        'sigma': 10.0,
        'bc_type': ('D', 'N'),
        'problem': 'mixed_dn',
        'zero_mode': False,
        'plot_in_summary': False,
    },
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _copy_to_thesis(fig_path: str) -> None:
    target_dir = _repo_root() / 'thesis' / 'figures'
    target_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(fig_path, target_dir / Path(fig_path).name)


def _grid(n: int):
    x = np.linspace(0.0, SX, n)
    y = np.linspace(0.0, SY, n)
    return np.meshgrid(x, y, indexing='ij')


def _problem_data(case, n: int):
    X, Y = _grid(n)
    sigma = case['sigma']
    if case['problem'] == 'neumann':
        u_exact, f_fun, bc = test_problem_neumann(sigma, SX, SY)
    elif case['problem'] == 'mixed_nd':
        u_exact, f_fun, bc = test_problem_mixed_nd(sigma, SX, SY)
    elif case['problem'] == 'mixed_dn':
        u_exact, f_fun, bc = test_problem_mixed_dn(sigma, SX, SY)
    else:
        raise ValueError(f"Unknown problem type: {case['problem']}")
    U_exact = u_exact(X, Y)
    f_rhs = f_fun(X, Y)
    return X, Y, U_exact, f_rhs, f_fun, bc


def _rms_error(err: np.ndarray) -> float:
    return float(np.sqrt(np.mean(err**2)))


def _boundary_flux_linf(U: np.ndarray, h: float, bc_type) -> float:
    """Second-order one-sided derivative residual on homogeneous Neumann edges."""
    residuals = []

    if bc_type == 'neumann' or (isinstance(bc_type, tuple) and bc_type[0] == 'N'):
        residuals.append((-3.0 * U[0, :] + 4.0 * U[1, :] - U[2, :]) / (2.0 * h))
        residuals.append((3.0 * U[-1, :] - 4.0 * U[-2, :] + U[-3, :]) / (2.0 * h))

    if bc_type == 'neumann' or (isinstance(bc_type, tuple) and bc_type[1] == 'N'):
        residuals.append((-3.0 * U[:, 0] + 4.0 * U[:, 1] - U[:, 2]) / (2.0 * h))
        residuals.append((3.0 * U[:, -1] - 4.0 * U[:, -2] + U[:, -3]) / (2.0 * h))

    if not residuals:
        return np.nan
    return float(max(np.max(np.abs(r)) for r in residuals))


def _boundary_dirichlet_linf(U: np.ndarray, U_exact: np.ndarray, bc_type) -> float:
    """Maximum boundary value residual on Dirichlet edges."""
    residuals = []

    if isinstance(bc_type, tuple) and bc_type[0] == 'D':
        residuals.append(U[0, :] - U_exact[0, :])
        residuals.append(U[-1, :] - U_exact[-1, :])

    if isinstance(bc_type, tuple) and bc_type[1] == 'D':
        residuals.append(U[:, 0] - U_exact[:, 0])
        residuals.append(U[:, -1] - U_exact[:, -1])

    if not residuals:
        return np.nan
    return float(max(np.max(np.abs(r)) for r in residuals))


def _trapezoid_integral(F: np.ndarray) -> float:
    wx = trapezoid_weights(F.shape[0])
    wy = trapezoid_weights(F.shape[1])
    h = SX / (F.shape[0] - 1)
    return float(h * h * np.sum(wx[:, None] * wy[None, :] * F))


def _solve_case(case, n: int, method: str, solver):
    X, Y, U_exact, f_rhs, f_fun, bc = _problem_data(case, n)
    h = SX / (n - 1)
    U_num = solver(n, f_fun, bc, sigma=case['sigma'], bc_type=case['bc_type'], sx=SX, sy=SY)

    if U_num.shape != U_exact.shape:
        raise ValueError(
            f"{method} returned shape {U_num.shape}, expected {U_exact.shape} "
            f"for {case['sub_exp']} at n={n}"
        )

    raw_error = U_num - U_exact
    mean_offset = 0.0
    if case['zero_mode']:
        mean_offset = weighted_mean_full_grid(raw_error, SX, SY)
        error = raw_error - mean_offset
    else:
        error = raw_error

    linf_error = float(np.max(np.abs(error)))
    l2_error = _rms_error(error)

    pure_neumann = case['bc_type'] == 'neumann'
    weighted_mean_solution = (
        weighted_mean_full_grid(U_num, SX, SY)
        if pure_neumann else np.nan
    )
    compatibility_integral = (
        _trapezoid_integral(f_rhs)
        if pure_neumann else np.nan
    )

    return {
        'method': method,
        'sigma': case['sigma'],
        'equation_type': equation_type(case['sigma']),
        'bc_type': str(case['bc_type']),
        'sub_exp': case['sub_exp'],
        'display': case['display'],
        'n': n,
        'h': h,
        'linf_error': linf_error,
        'l2_error': l2_error,
        'observed_order_linf': np.nan,
        'observed_order_l2': np.nan,
        'boundary_flux_linf': _boundary_flux_linf(U_num, h, case['bc_type']),
        'boundary_dirichlet_linf': _boundary_dirichlet_linf(
            U_num, U_exact, case['bc_type']
        ),
        'mean_offset': float(mean_offset),
        'weighted_mean_solution': float(weighted_mean_solution),
        'compatibility_integral': float(compatibility_integral),
        # Backward-compatible aliases for older tables/scripts.
        'error': linf_error,
        'rate': np.nan,
    }


def run_experiment() -> pd.DataFrame:
    print("Running Experiment 03 V2: Neumann and mixed boundary verification")
    rows = []

    for case in CASES:
        print(f"\n{case['display']}")
        for method, solver in SOLVERS.items():
            for n in NS:
                row = _solve_case(case, n, method, solver)
                rows.append(row)
                print(
                    f"  {method:4s} n={n:3d}: "
                    f"Linf={row['linf_error']:.3e}, "
                    f"L2={row['l2_error']:.3e}, "
                    f"flux={row['boundary_flux_linf']:.3e}"
                )

    df = pd.DataFrame(rows)

    for (sub_exp, method), group in df.groupby(['sub_exp', 'method']):
        idx = group.sort_values('n').index
        linf_rates = compute_convergence_rate(df.loc[idx, 'linf_error'].to_numpy())
        l2_rates = compute_convergence_rate(df.loc[idx, 'l2_error'].to_numpy())
        df.loc[idx, 'observed_order_linf'] = linf_rates
        df.loc[idx, 'observed_order_l2'] = l2_rates
        df.loc[idx, 'rate'] = linf_rates

    out_csv = os.path.join(get_results_dir(), 'exp03_neumann_mixed.csv')
    df.to_csv(out_csv, index=False)
    print(f"\nSaved results to {out_csv}")
    return df


def _positive(values, floor=1e-18):
    arr = np.asarray(values, dtype=float)
    return np.where(np.isfinite(arr), np.maximum(np.abs(arr), floor), np.nan)


def _add_h2_reference(ax, h_values, y_values):
    h = np.asarray(h_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    mask = np.isfinite(h) & np.isfinite(y) & (h > 0.0) & (y > 0.0)
    if not np.any(mask):
        return
    h0 = np.max(h[mask])
    y0 = y[mask][np.argmax(h[mask])]
    h_ref = np.array([np.min(h[mask]), np.max(h[mask])])
    ax.loglog(h_ref, y0 * (h_ref / h0) ** 2, 'k:', linewidth=1.2, label='$O(h^2)$')


def _plot_case_convergence(ax, df: pd.DataFrame, sub_exp: str, title: str):
    for method, color in METHOD_COLORS.items():
        data = df[(df['sub_exp'] == sub_exp) & (df['method'] == method)].sort_values('h')
        ax.loglog(
            data['h'], data['linf_error'],
            marker='o', linestyle='-', color=color, linewidth=1.2,
            label=f'{method} $L_\\infty$'
        )
        ax.loglog(
            data['h'], data['l2_error'],
            marker='s', linestyle='--', color=color, linewidth=1.2,
            label=f'{method} $L_2$'
        )

    ref = df[(df['sub_exp'] == sub_exp) & (df['method'] == 'FA')].sort_values('h')
    _add_h2_reference(ax, ref['h'], ref['linf_error'])
    ax.set_title(title)
    ax.set_xlabel('$h$')
    ax.set_ylabel('error')
    ax.grid(True, which='both', alpha=0.3)


def plot_legacy_linf(df: pd.DataFrame):
    """Keep the original figure path for backward compatibility."""
    cases = [case for case in CASES if case['plot_in_summary']]
    fig, axes = plt.subplots(1, len(cases), figsize=(4.0 * len(cases), 3.4))
    if len(cases) == 1:
        axes = [axes]

    for ax, case in zip(axes, cases):
        for method, color in METHOD_COLORS.items():
            data = df[
                (df['sub_exp'] == case['sub_exp'])
                & (df['method'] == method)
            ].sort_values('h')
            ax.loglog(data['h'], data['linf_error'], 'o-', color=color, label=method)

        ref = df[
            (df['sub_exp'] == case['sub_exp'])
            & (df['method'] == 'FA')
        ].sort_values('h')
        _add_h2_reference(ax, ref['h'], ref['linf_error'])
        ax.set_title(case['display'])
        ax.set_xlabel('$h$')
        ax.set_ylabel('$L_\\infty$ error')
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(fontsize=7)

    fig.suptitle('Experiment 03: Neumann/Mixed boundary convergence')
    fig.tight_layout()
    fig_path = os.path.join(get_figures_dir(), 'exp03_neumann_mixed.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    _copy_to_thesis(fig_path)
    print(f"Saved legacy figure to {fig_path}")


def plot_summary(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)

    core_cases = [
        ('a_pure_neumann_sigma1', 'Pure N $\\sigma=1$'),
        ('a_pure_neumann_sigma10', 'Pure N $\\sigma=10$'),
        ('b_mixed_ND_sigma10', 'Mixed (N,D)'),
        ('c_poisson_neumann_sigma0', 'N Poisson'),
    ]

    ax = axes[0]
    for sub_exp, label in core_cases:
        data = df[
            (df['sub_exp'] == sub_exp)
            & (df['method'] == 'FA')
        ].sort_values('h')
        ax.loglog(
            data['h'], data['linf_error'],
            marker='o', linewidth=1.4, label=label
        )
    ref = df[
        (df['sub_exp'] == 'a_pure_neumann_sigma1')
        & (df['method'] == 'FA')
    ].sort_values('h')
    _add_h2_reference(ax, ref['h'], ref['linf_error'])
    ax.set_title('(a) Linf convergence, FA solver')
    ax.set_xlabel('$h$')
    ax.set_ylabel('$L_\\infty$ error')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[1]
    flux_ref = df[
        (df['sub_exp'] == 'a_pure_neumann_sigma1')
        & (df['method'] == 'FA')
    ].sort_values('h')
    ax.loglog(
        flux_ref['h'], _positive(flux_ref['boundary_flux_linf']),
        marker='o', linewidth=1.3, label='Neumann flux (all cases coincide)'
    )
    mixed_nd = df[
        (df['sub_exp'] == 'b_mixed_ND_sigma10')
        & (df['method'] == 'FA')
    ].sort_values('h')
    ax.loglog(
        mixed_nd['h'], _positive(mixed_nd['boundary_dirichlet_linf']),
        marker='x', linestyle='--', linewidth=1.3,
        label='Mixed (N,D) Dirichlet'
    )
    ax.set_title('(b) Boundary constraint residuals')
    ax.set_xlabel('$h$')
    ax.set_ylabel('boundary residual')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[2]
    for method, color in METHOD_COLORS.items():
        data = df[
            (df['sub_exp'] == 'c_poisson_neumann_sigma0')
            & (df['method'] == method)
        ].sort_values('h')
        ax.loglog(
            data['h'], _positive(data['weighted_mean_solution']),
            marker='o', linestyle='-', color=color, linewidth=1.2,
            label=f'{method} $|\\bar u_w|$'
        )
        ax.loglog(
            data['h'], _positive(data['mean_offset']),
            marker='s', linestyle='--', color=color, linewidth=1.2,
            label=f'{method} removed offset'
        )
    ax.set_title('(c) Neumann Poisson zero-mode check')
    ax.set_xlabel('$h$')
    ax.set_ylabel('absolute value')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=6, ncol=2)

    fig.suptitle('Neumann/Mixed boundary verification summary', fontsize=15)
    fig_path = os.path.join(get_figures_dir(), 'exp03_neumann_mixed_summary.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    _copy_to_thesis(fig_path)
    print(f"Saved summary figure to {fig_path}")


def _solve_for_field(case, n: int, method='FA'):
    solver = SOLVERS[method]
    X, Y, U_exact, f_rhs, f_fun, bc = _problem_data(case, n)
    U_num = solver(n, f_fun, bc, sigma=case['sigma'], bc_type=case['bc_type'], sx=SX, sy=SY)
    if case['zero_mode']:
        U_num = U_num - weighted_mean_full_grid(U_num - U_exact, SX, SY)
    err = U_num - U_exact
    return X, Y, U_exact, U_num, err


def plot_fields(n: int = 65):
    field_cases = [
        next(case for case in CASES if case['sub_exp'] == 'a_pure_neumann_sigma1'),
        next(case for case in CASES if case['sub_exp'] == 'b_mixed_ND_sigma10'),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.2), constrained_layout=True)

    for row, case in enumerate(field_cases):
        _, _, U_exact, U_num, err = _solve_for_field(case, n, method='FA')
        vmin = min(float(U_exact.min()), float(U_num.min()))
        vmax = max(float(U_exact.max()), float(U_num.max()))
        log_err = np.log10(np.abs(err) + 1e-16)
        max_err = float(np.max(np.abs(err)))

        titles = [
            f"{case['display']} exact",
            f"{case['display']} FA numerical",
            f"{case['display']} log10 error\nmax={max_err:.2e}",
        ]
        arrays = [U_exact, U_num, log_err]
        cmaps = ['viridis', 'viridis', 'magma']

        for col in range(3):
            ax = axes[row, col]
            if col < 2:
                im = ax.imshow(
                    arrays[col].T, origin='lower', extent=[0, 1, 0, 1],
                    cmap=cmaps[col], vmin=vmin, vmax=vmax, aspect='equal'
                )
            else:
                im = ax.imshow(
                    arrays[col].T, origin='lower', extent=[0, 1, 0, 1],
                    cmap=cmaps[col], aspect='equal'
                )
            ax.set_title(titles[col], fontsize=9)
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(f'Neumann/Mixed representative fields and errors, n={n}', fontsize=14)
    fig_path = os.path.join(get_figures_dir(), 'exp03_neumann_mixed_fields.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    _copy_to_thesis(fig_path)
    print(f"Saved field figure to {fig_path}")


def main():
    df = run_experiment()
    plot_legacy_linf(df)
    plot_summary(df)
    plot_fields(n=65)
    return df


if __name__ == '__main__':
    main()
