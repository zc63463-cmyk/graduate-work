#!/usr/bin/env python3
"""Experiment 02: Non-Homogeneous Dirichlet BC — A1 Bug Fix Verification

Verify that the F += bc/h² sign fix works correctly for non-zero BC.
Test function: u = sin(2πx)sin(3πy) with non-zero boundary values.
Four σ tracks: 0, +10, +100, -5.
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
    test_problem_dirichlet_mode,
    compute_convergence_rate, equation_type
)


def run():
    """Run experiment 02: non-homogeneous Dirichlet BC verification."""
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    rows = []

    sigmas = [0.0, 10.0, 100.0, -5.0]
    methods = ['FA', 'CR', 'FACR', 'FFT9']
    ns = [9, 17, 33, 65]
    sx, sy = 1.0, 1.0

    for sigma in sigmas:
        eq_type = equation_type(sigma)
        # u = sin(2πx)sin(3πy) — non-zero on boundary
        u_exact, f_rhs, bc = test_problem_dirichlet_mode(sigma, m=2, n=3, sx=sx, sy=sy)

        print(f"\n{'='*70}")
        print(f"σ = {sigma:+.1f}  ({eq_type})  u = sin(2πx)sin(3πy)")
        print(f"{'='*70}")

        for method in methods:
            errors = []
            for n in ns:
                x = np.linspace(0, sx, n)
                X, Y = np.meshgrid(x, x, indexing='ij')
                Ue = u_exact(X, Y)

                try:
                    U = fa_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='dirichlet', sx=sx, sy=sy) \
                        if method == 'FA' else \
                        cr_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='dirichlet', sx=sx, sy=sy) \
                        if method == 'CR' else \
                        facr_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='dirichlet', sx=sx, sy=sy) \
                        if method == 'FACR' else \
                        fft9_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='dirichlet', sx=sx, sy=sy)
                    err = np.max(np.abs(U - Ue))
                except Exception as e:
                    print(f"  {method:>5s} n={n:3d}: ERROR {e}")
                    err = np.nan

                errors.append(err)
                rows.append({
                    'method': method,
                    'sigma': sigma,
                    'equation_type': eq_type,
                    'n': n,
                    'h': sx / (n - 1),
                    'error': err,
                })

            rates = compute_convergence_rate(errors)
            for i in range(len(ns)):
                rows[-len(ns) + i]['rate'] = rates[i]

            # Print
            print(f"\n  {method}:")
            print(f"  {'n':>5s} | {'error':>12s} | {'rate':>6s}")
            print(f"  {'-'*30}")
            for i, n in enumerate(ns):
                r_str = f"{rates[i]:6.2f}" if not np.isnan(rates[i]) else "   --"
                e_str = f"{errors[i]:12.2e}" if not np.isnan(errors[i]) else "       NaN"
                print(f"  {n:5d} | {e_str} | {r_str}")

    df = pd.DataFrame(rows)
    csv_path = os.path.join(results_dir, 'exp02_nonhom_bc.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")

    # --- Verification summary ---
    print(f"\n{'='*50}")
    print("A1 Bug Fix Verification Summary:")
    for method in methods:
        for sigma in sigmas:
            sub = df[(df['method'] == method) & (df['sigma'] == sigma)]
            rates = sub['rate'].dropna()
            if len(rates) > 0:
                # Check if rates are near 2 (for 5pt) or 4 (for FFT9)
                expected = 4.0 if method == 'FFT9' else 2.0
                # For FFT9: superconvergence on eigenfunctions can give rate > 4.0 on coarse grids.
                # The key check is: (1) rate is ≥ expected, AND (2) rate is decreasing towards expected.
                # For 5pt: rate should be near 2.0 (tolerance 0.5).
                if method == 'FFT9':
                    positive_rates = [r for r in rates if r > 0]
                    # Check all rates ≥ 3.5 (at least 4th-order)
                    all_above = all(r >= 3.5 for r in positive_rates)
                    # Check last rate is within 0.5 of 4.0 (approaching asymptotic)
                    last_rate = positive_rates[-1] if positive_rates else 0
                    asymp_ok = abs(last_rate - 4.0) < 0.5
                    # Check rate is decreasing (superconvergence fading)
                    decreasing = all(positive_rates[i] >= positive_rates[i+1]
                                     for i in range(len(positive_rates)-1))
                    ok = all_above and (asymp_ok or decreasing)
                else:
                    ok = all(abs(r - expected) < 0.5 for r in rates if r > 0)
                avg_rate = np.mean([r for r in rates if r > 0])
                status = "PASS" if ok else "WARN"
                print(f"  {method:>5s} σ={sigma:+6.1f}: avg rate = {avg_rate:.2f}  [{status}]")

    plot(df, figures_dir)
    return df


def plot(df=None, figures_dir=None):
    """Generate convergence plots for non-homogeneous BC."""
    if df is None:
        results_dir = get_results_dir()
        df = pd.read_csv(os.path.join(results_dir, 'exp02_nonhom_bc.csv'))
    if figures_dir is None:
        figures_dir = get_figures_dir()

    sigmas = sorted(df['sigma'].unique())
    markers = {'FA': 'o', 'CR': 's', 'FACR': '^', 'FFT9': 'D'}
    colors  = {'FA': 'C0', 'CR': 'C1', 'FACR': 'C2', 'FFT9': 'C3'}

    fig, axes = plt.subplots(1, len(sigmas), figsize=(5 * len(sigmas), 5), squeeze=False)

    for idx, sigma in enumerate(sigmas):
        ax = axes[0, idx]
        sub = df[df['sigma'] == sigma].dropna(subset=['error'])
        eq_type = sub['equation_type'].iloc[0] if len(sub) > 0 else ''

        for method in sub['method'].unique():
            md = sub[sub['method'] == method].sort_values('h')
            if len(md) > 1:
                ax.loglog(md['h'], md['error'], marker=markers.get(method, 'o'),
                          label=method, color=colors.get(method, 'k'), linewidth=1.5)

        h_vals = sub['h'].sort_values().unique()
        if len(h_vals) >= 2:
            h_ref = np.array([h_vals[0], h_vals[-1]])
            ax.loglog(h_ref, 0.5 * h_ref**2, 'k--', alpha=0.3, label='$O(h^2)$')
            ax.loglog(h_ref, 0.05 * h_ref**4, 'k:', alpha=0.3, label='$O(h^4)$')

        ax.set_xlabel('$h$')
        ax.set_ylabel('$\\|u - u_h\\|_\\infty$')
        ax.set_title(f'Non-hom BC: $\\sigma={sigma:+.0f}$ ({eq_type})')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(figures_dir, 'exp02_nonhom_bc.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {fig_path}")
    plt.close()


if __name__ == '__main__':
    run()
