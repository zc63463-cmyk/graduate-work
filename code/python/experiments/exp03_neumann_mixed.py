#!/usr/bin/env python3
"""Experiment 03: Neumann & Mixed Boundary Conditions

Sub-experiments:
  (a) Pure Neumann modified H. (σ=+1, +10)
  (b) Mixed BC (N,D): u = cos(πx)sin(πy), σ=+10
  (c) Pure Neumann Poisson (σ=0): compatibility condition + zero mean (trapezoid-weighted)
  (d) Mixed BC (D,N): u = sin(πx)cos(πy), σ=+10  [smoke test]
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz
from .utils import (
    get_results_dir, get_figures_dir,
    test_problem_neumann, test_problem_mixed_nd, test_problem_mixed_dn,
    compute_convergence_rate, equation_type, weighted_mean_full_grid
)


def run():
    """Run experiment 03: Neumann & mixed BC verification."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    rows = []

    ns = [9, 17, 33, 65, 129]
    sx, sy = 1.0, 1.0

    # ---- Sub-experiment (a): Pure Neumann modified H. ----
    for sigma in [1.0, 10.0]:
        eq_type = equation_type(sigma)
        u_exact, f_rhs, bc = test_problem_neumann(sigma, sx, sy)

        print(f"\n{'='*70}")
        print(f"Sub-exp (a): Pure Neumann  σ={sigma:+.1f}  ({eq_type})")
        print(f"{'='*70}")

        for method_name, solver_fn in [('FA', fa_helmholtz), ('CR', cr_helmholtz),
                                         ('FACR', facr_helmholtz)]:
            errors = []
            for n in ns:
                x = np.linspace(0, sx, n)
                X, Y = np.meshgrid(x, x, indexing='ij')
                Ue = u_exact(X, Y)

                try:
                    U = solver_fn(n, f_rhs, bc, sigma=sigma, bc_type='neumann', sx=sx, sy=sy)
                    err = np.max(np.abs(U - Ue))
                except Exception as e:
                    print(f"  {method_name:>5s} n={n:3d}: ERROR {e}")
                    err = np.nan

                errors.append(err)
                rows.append({
                    'method': method_name,
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'bc_type': 'neumann',
                    'sub_exp': 'a_pure_neumann',
                    'n': n,
                    'h': sx / (n - 1),
                    'error': err,
                })

            rates = compute_convergence_rate(errors)
            for i in range(len(ns)):
                rows[-len(ns) + i]['rate'] = rates[i]

            print(f"\n  {method_name}:")
            print(f"  {'n':>5s} | {'error':>12s} | {'rate':>6s}")
            print(f"  {'-'*30}")
            for i, n in enumerate(ns):
                r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
                e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
                print(f"  {n:5d} | {e_str} | {r_str}")

    # ---- Sub-experiment (b): Mixed BC (N,D) ----
    sigma = 10.0
    eq_type = equation_type(sigma)
    u_exact, f_rhs, bc = test_problem_mixed_nd(sigma, sx, sy)

    print(f"\n{'='*70}")
    print(f"Sub-exp (b): Mixed (N,D)  σ={sigma:+.1f}  ({eq_type})")
    print(f"{'='*70}")

    for method_name, solver_fn in [('FA', fa_helmholtz), ('CR', cr_helmholtz),
                                     ('FACR', facr_helmholtz)]:
        errors = []
        for n in ns:
            x = np.linspace(0, sx, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U = solver_fn(n, f_rhs, bc, sigma=sigma, bc_type=('N', 'D'), sx=sx, sy=sy)
                err = np.max(np.abs(U - Ue))
            except Exception as e:
                print(f"  {method_name:>5s} n={n:3d}: ERROR {e}")
                err = np.nan

            errors.append(err)
            rows.append({
                'method': method_name,
                'sigma': sigma,
                'equation_type': eq_type,
                'bc_type': 'N,D',
                'sub_exp': 'b_mixed_ND',
                'n': n,
                'h': sx / (n - 1),
                'error': err,
            })

        rates = compute_convergence_rate(errors)
        for i in range(len(ns)):
            rows[-len(ns) + i]['rate'] = rates[i]

        print(f"\n  {method_name}:")
        print(f"  {'n':>5s} | {'error':>12s} | {'rate':>6s}")
        print(f"  {'-'*30}")
        for i, n in enumerate(ns):
            r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
            e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
            print(f"  {n:5d} | {e_str} | {r_str}")

    # ---- Sub-experiment (c): Pure Neumann Poisson ----
    sigma = 0.0
    eq_type = 'Poisson'
    u_exact, f_rhs, bc = test_problem_neumann(sigma, sx, sy)

    print(f"\n{'='*70}")
    print(f"Sub-exp (c): Pure Neumann Poisson  σ={sigma:.1f}  ({eq_type})")
    print(f"  (compatibility condition + trapezoid-weighted zero-mean check)")
    print(f"{'='*70}")

    for method_name, solver_fn in [('FA', fa_helmholtz), ('CR', cr_helmholtz),
                                     ('FACR', facr_helmholtz)]:
        errors = []
        for n in ns:
            x = np.linspace(0, sx, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U = solver_fn(n, f_rhs, bc, sigma=sigma, bc_type='neumann', sx=sx, sy=sy)
                # For Neumann Poisson, solution has arbitrary constant offset.
                # Use trapezoid-weighted mean (consistent with DCT-I orthogonality).
                mean_diff = weighted_mean_full_grid(U - Ue, sx, sy)
                err = np.max(np.abs(U - Ue - mean_diff))
                zero_mean = abs(mean_diff)
            except Exception as e:
                print(f"  {method_name:>5s} n={n:3d}: ERROR {e}")
                err = np.nan
                zero_mean = np.nan

            errors.append(err)
            rows.append({
                'method': method_name,
                'sigma': sigma,
                'equation_type': eq_type,
                'bc_type': 'neumann',
                'sub_exp': 'c_poisson_neumann',
                'n': n,
                'h': sx / (n - 1),
                'error': err,
                'mean_offset': zero_mean,
            })

        rates = compute_convergence_rate(errors)
        for i in range(len(ns)):
            rows[-len(ns) + i]['rate'] = rates[i]

        print(f"\n  {method_name}:")
        print(f"  {'n':>5s} | {'error':>12s} | {'rate':>6s} | {'mean_off':>12s}")
        print(f"  {'-'*50}")
        for i, n in enumerate(ns):
            r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
            e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
            m_str = f"{rows[-len(ns)+i].get('mean_offset', np.nan):12.2e}"
            print(f"  {n:5d} | {e_str} | {r_str} | {m_str}")

    # ---- Sub-experiment (d): Mixed BC (D,N) smoke test ----
    sigma = 10.0
    eq_type = equation_type(sigma)
    u_exact, f_rhs, bc = test_problem_mixed_dn(sigma, sx, sy)

    print(f"\n{'='*70}")
    print(f"Sub-exp (d): Mixed (D,N)  σ={sigma:+.1f}  ({eq_type})")
    print(f"  (smoke test: x=Dirichlet, y=Neumann)")
    print(f"{'='*70}")

    for method_name, solver_fn in [('FA', fa_helmholtz), ('CR', cr_helmholtz),
                                     ('FACR', facr_helmholtz)]:
        errors = []
        for n in ns:
            x = np.linspace(0, sx, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U = solver_fn(n, f_rhs, bc, sigma=sigma, bc_type=('D', 'N'), sx=sx, sy=sy)
                err = np.max(np.abs(U - Ue))
            except Exception as e:
                print(f"  {method_name:>5s} n={n:3d}: ERROR {e}")
                err = np.nan

            errors.append(err)
            rows.append({
                'method': method_name,
                'sigma': sigma,
                'equation_type': eq_type,
                'bc_type': 'D,N',
                'sub_exp': 'd_mixed_DN',
                'n': n,
                'h': sx / (n - 1),
                'error': err,
            })

        rates = compute_convergence_rate(errors)
        for i in range(len(ns)):
            rows[-len(ns) + i]['rate'] = rates[i]

        print(f"\n  {method_name}:")
        print(f"  {'n':>5s} | {'error':>12s} | {'rate':>6s}")
        print(f"  {'-'*30}")
        for i, n in enumerate(ns):
            r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
            e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
            print(f"  {n:5d} | {e_str} | {r_str}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, 'exp03_neumann_mixed.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")

    plot(df, figures_dir)
    return df


def plot(df=None, figures_dir=None):
    """Generate convergence plots for Neumann & mixed BC."""
    if df is None:
        results_dir = get_results_dir()
        df = pd.read_csv(os.path.join(results_dir, 'exp03_neumann_mixed.csv'))
    if figures_dir is None:
        figures_dir = get_figures_dir()

    sub_exps = sorted(df['sub_exp'].unique())
    markers = {'FA': 'o', 'CR': 's', 'FACR': '^'}
    colors  = {'FA': 'C0', 'CR': 'C1', 'FACR': 'C2'}

    fig, axes = plt.subplots(1, len(sub_exps), figsize=(5 * len(sub_exps), 5), squeeze=False)

    for idx, sub_exp in enumerate(sub_exps):
        ax = axes[0, idx]
        sub = df[df['sub_exp'] == sub_exp].dropna(subset=['error'])
        if len(sub) == 0:
            continue

        sigma = sub['sigma'].iloc[0]
        bc_type = sub['bc_type'].iloc[0]
        eq_type = sub['equation_type'].iloc[0]

        for method in sub['method'].unique():
            md = sub[sub['method'] == method].sort_values('h')
            if len(md) > 1:
                ax.loglog(md['h'], md['error'], marker=markers.get(method, 'o'),
                          label=method, color=colors.get(method, 'k'), linewidth=1.5)

        h_vals = sub['h'].sort_values().unique()
        if len(h_vals) >= 2:
            h_ref = np.array([h_vals[0], h_vals[-1]])
            ax.loglog(h_ref, 0.5 * h_ref**2, 'k--', alpha=0.3, label='$O(h^2)$')

        ax.set_xlabel('$h$')
        ax.set_ylabel('$\\|u - u_h\\|_\\infty$')
        ax.set_title(f'{sub_exp}: $\\sigma={sigma}$, {bc_type}')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'exp03_neumann_mixed.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {fig_path}")
    plt.close()


if __name__ == '__main__':
    run()
