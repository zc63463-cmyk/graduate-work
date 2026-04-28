#!/usr/bin/env python3
"""Experiment 01: Convergence Order Verification (Upgraded: 3-track)

Verify theoretical convergence rates for Poisson (σ=0), modified H. (σ=+10),
and true H. (σ=-5) — the signature three-track comparison.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
from helmholtz_solver import fa_helmholtz, cr_helmholtz, facr_helmholtz, fft9_helmholtz
from .utils import (
    get_results_dir, get_figures_dir,
    test_problem_dirichlet, test_problem_dirichlet_mode,
    compute_convergence_rate, equation_type
)


def run():
    """Run experiment 01: convergence order verification."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    rows = []

    # Three tracks: Poisson, modified H., true H.
    tracks = [
        (0.0,  'dirichlet', 'homogeneous', test_problem_dirichlet),
        (10.0, 'dirichlet', 'homogeneous', test_problem_dirichlet),
        (-5.0, 'dirichlet', 'nonhom', lambda sigma, sx=1.0, sy=1.0:
            test_problem_dirichlet_mode(sigma, m=2, n=3, sx=sx, sy=sy)),
        # Note: true H. uses sin(2πx)sin(3πy) with λ_{2,3}=13π²≈128.3, σ=-5 is safe
    ]

    methods = [
        ('FA',   'fa',   lambda n, f, bc, **kw: fa_helmholtz(n, f, bc, **kw)),
        ('CR',   'cr',   lambda n, f, bc, **kw: cr_helmholtz(n, f, bc, **kw)),
        ('FACR', 'facr', lambda n, f, bc, **kw: facr_helmholtz(n, f, bc, **kw)),
        ('FFT9', 'fft9', lambda n, f, bc, **kw: fft9_helmholtz(n, f, bc, **kw)),
    ]

    ns = [9, 17, 33, 65, 129]
    sx, sy = 1.0, 1.0

    for sigma, bc_type, bc_sub, problem_fn in tracks:
        eq_type = equation_type(sigma)
        u_exact, f_rhs, bc = problem_fn(sigma, sx, sy)

        print(f"\n{'='*70}")
        print(f"σ = {sigma:+.1f}  ({eq_type})  bc_subtype={bc_sub}")
        print(f"{'='*70}")

        for method_name, method_key, solver_fn in methods:
            errors = []
            for n in ns:
                x = np.linspace(0, sx, n)
                X, Y = np.meshgrid(x, x, indexing='ij')
                Ue = u_exact(X, Y)

                try:
                    U = solver_fn(n, f_rhs, bc, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)
                    err = np.max(np.abs(U - Ue))
                except Exception as e:
                    print(f"  {method_name:>5s} n={n:3d}: ERROR {e}")
                    err = np.nan

                errors.append(err)
                rows.append({
                    'method': method_name,
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'bc_subtype': bc_sub,
                    'n': n,
                    'h': sx / (n - 1),
                    'error': err,
                })

            rates = compute_convergence_rate(errors)
            for i, n in enumerate(ns):
                if i < len(rows):
                    rows[-len(ns) + i]['rate'] = rates[i]

            # Print table
            print(f"\n  {method_name}:")
            print(f"  {'n':>5s} | {'h':>10s} | {'error':>12s} | {'rate':>6s}")
            print(f"  {'-'*45}")
            for i, n in enumerate(ns):
                h_val = sx / (n - 1)
                r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
                e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
                print(f"  {n:5d} | {h_val:10.2e} | {e_str} | {r_str}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, 'exp01_convergence.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")

    # --- Plot ---
    plot(df, figures_dir)
    return df


def plot(df=None, figures_dir=None):
    """Generate convergence rate plots from CSV data."""
    if df is None:
        results_dir = get_results_dir()
        df = pd.read_csv(os.path.join(results_dir, 'exp01_convergence.csv'))
    if figures_dir is None:
        figures_dir = get_figures_dir()

    sigmas = sorted(df['sigma'].unique())
    n_sigmas = len(sigmas)

    fig, axes = plt.subplots(1, n_sigmas, figsize=(5 * n_sigmas, 5), squeeze=False)
    markers = {'FA': 'o', 'CR': 's', 'FACR': '^', 'FFT9': 'D'}
    colors  = {'FA': 'C0', 'CR': 'C1', 'FACR': 'C2', 'FFT9': 'C3'}

    for idx, sigma in enumerate(sigmas):
        ax = axes[0, idx]
        sub = df[df['sigma'] == sigma].dropna(subset=['error'])
        eq_type = sub['equation_type'].iloc[0] if len(sub) > 0 else ''

        for method in sub['method'].unique():
            md = sub[sub['method'] == method].sort_values('h')
            if len(md) > 1:
                ax.loglog(md['h'], md['error'], marker=markers.get(method, 'o'),
                          label=method, color=colors.get(method, 'k'), linewidth=1.5)

        # Reference lines
        h_vals = sub['h'].sort_values().unique()
        if len(h_vals) >= 2:
            h_ref = np.array([h_vals[0], h_vals[-1]])
            e_ref = 0.5 * h_ref**2
            ax.loglog(h_ref, e_ref, 'k--', alpha=0.3, label='$O(h^2)$')
            e_ref4 = 0.05 * h_ref**4
            ax.loglog(h_ref, e_ref4, 'k:', alpha=0.3, label='$O(h^4)$')

        ax.set_xlabel('$h$')
        ax.set_ylabel('$\\|u - u_h\\|_\\infty$')
        ax.set_title(f'$\\sigma = {sigma:+.0f}$ ({eq_type})')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'exp01_convergence.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {fig_path}")
    plt.close()


if __name__ == '__main__':
    run()
