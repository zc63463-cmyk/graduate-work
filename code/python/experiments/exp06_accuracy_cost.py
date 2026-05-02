#!/usr/bin/env python3
"""Experiment 06 supplement: accuracy-cost comparison.

This script compares the current implementations on one smooth homogeneous
Dirichlet manufactured solution. It is an implementation-level baseline:
FA/CR/FACR-like and GMRES solve the five-point second-order system, while FFT9
solves the Dirichlet compact fourth-order system.
"""

from __future__ import annotations

import os
import shutil
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

from gmres_solver import (  # noqa: E402
    _assemble_solution,
    _build_rhs,
    build_helmholtz_matrix,
    gmres,
)
from helmholtz_solver import (  # noqa: E402
    cr_helmholtz,
    fa_helmholtz,
    facr_helmholtz,
    fft9_helmholtz,
)

try:
    from .utils import get_figures_dir, get_results_dir  # noqa: E402
except ImportError:  # Support direct execution as a file with PYTHONPATH=code/python.
    from experiments.utils import get_figures_dir, get_results_dir  # noqa: E402


SX = 1.0
SY = 1.0
SIGMA_LIST = [0.0, 10.0]
N_LIST_DIRECT = [33, 65, 129, 257, 513]
N_LIST_GMRES = [33, 65, 129]
WARMUP = 2
REPEAT = 7

GMRES_RESTART = 30
GMRES_TOL_TYPE = "absolute_residual"
GMRES_TOL = 1e-10
GMRES_MAXITER = 2000

CSV_NAME = "exp06_accuracy_cost.csv"
ERROR_TIME_FIG = "exp06_accuracy_cost_error_time"
TIME_SCALING_FIG = "exp06_time_scaling"


def u_exact(x, y):
    """u=x^2(1-x)^2 y^2(1-y)^2."""
    px = x**2 * (1.0 - x) ** 2
    py = y**2 * (1.0 - y) ** 2
    return px * py


def f_rhs_factory(sigma):
    """RHS for (-Delta + sigma)u=f with the polynomial exact solution."""

    def f_rhs(x, y):
        px = x**2 * (1.0 - x) ** 2
        py = y**2 * (1.0 - y) ** 2
        d2px = 2.0 - 12.0 * x + 12.0 * x**2
        d2py = 2.0 - 12.0 * y + 12.0 * y**2
        return -(d2px * py + px * d2py) + sigma * px * py

    return f_rhs


def bc_zero(x, y):
    return 0.0


def full_grid(n):
    x = np.linspace(0.0, SX, n)
    y = np.linspace(0.0, SY, n)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return X, Y


def error_metrics(U, n):
    X, Y = full_grid(n)
    err = U - u_exact(X, Y)
    return float(np.max(np.abs(err))), float(np.sqrt(np.mean(err**2)))


def summarize_times(times):
    arr = np.asarray(times, dtype=float)
    q25, q75 = np.percentile(arr, [25, 75])
    return float(np.median(arr)), float(q75 - q25)


def time_solver_call(solver_fn):
    times = []
    last = None
    for _ in range(WARMUP):
        last = solver_fn()
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        last = solver_fn()
        times.append(time.perf_counter() - t0)
    return last, *summarize_times(times)


def direct_solver_registry():
    return {
        "FA": ("five_point_second_order", fa_helmholtz),
        "CR": ("five_point_second_order", cr_helmholtz),
        "FACR": ("five_point_second_order", facr_helmholtz),
        "FFT9": ("fft9_compact_fourth_order_dirichlet", fft9_helmholtz),
    }


def run_direct_method(method, discretization, solver, n, sigma):
    f_rhs = f_rhs_factory(sigma)

    def solve_once():
        return solver(
            n,
            f_rhs,
            bc_zero,
            sigma=sigma,
            bc_type="dirichlet",
            sx=SX,
            sy=SY,
        )

    U, t_med, t_iqr = time_solver_call(solve_once)
    linf, l2 = error_metrics(U, n)
    return make_row(
        method=method,
        discretization=discretization,
        n=n,
        sigma=sigma,
        time_s_median=t_med,
        time_s_iqr=t_iqr,
        linf_error=linf,
        l2_error=l2,
        timing_scope="solver_call_including_internal_rhs_and_transform_setup",
    )


def run_gmres_method(n, sigma):
    f_rhs = f_rhs_factory(sigma)
    A, grid_info = build_helmholtz_matrix(
        n, sigma=sigma, bc_type="dirichlet", sx=SX, sy=SY
    )
    b = _build_rhs(n, f_rhs, bc_zero, sigma, grid_info["bc_x"], grid_info["bc_y"], SX, SY)
    norm_b = float(np.linalg.norm(b))
    last_info = None
    last_u = None
    times = []

    def solve_once():
        return gmres(
            A,
            b,
            tol=GMRES_TOL,
            max_iter=GMRES_MAXITER,
            restart=GMRES_RESTART,
            return_history=True,
        )

    for _ in range(WARMUP):
        last_u, last_info = solve_once()
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        last_u, last_info = solve_once()
        times.append(time.perf_counter() - t0)

    t_med, t_iqr = summarize_times(times)
    U = _assemble_solution(last_u, grid_info, bc_zero, SX, SY)
    linf, l2 = error_metrics(U, n)
    final_abs = float(last_info["final_abs_residual"])
    final_rel = final_abs / norm_b if norm_b > 0 else np.inf

    return make_row(
        method="GMRES30",
        discretization="five_point_second_order",
        n=n,
        sigma=sigma,
        time_s_median=t_med,
        time_s_iqr=t_iqr,
        linf_error=linf,
        l2_error=l2,
        gmres_restart=GMRES_RESTART,
        gmres_tol_type=GMRES_TOL_TYPE,
        gmres_tol=GMRES_TOL,
        gmres_iterations=min(int(last_info["iterations"]), GMRES_MAXITER),
        gmres_converged=bool(last_info["success"]),
        gmres_final_absres=final_abs,
        gmres_final_relres=float(final_rel),
        timing_scope="gmres_core_solve_only_matrix_and_rhs_setup_excluded",
    )


def make_row(
    method,
    discretization,
    n,
    sigma,
    time_s_median,
    time_s_iqr,
    linf_error,
    l2_error,
    timing_scope,
    gmres_restart=np.nan,
    gmres_tol_type="",
    gmres_tol=np.nan,
    gmres_iterations=np.nan,
    gmres_converged="",
    gmres_final_absres=np.nan,
    gmres_final_relres=np.nan,
):
    N = n - 2
    h = SX / (n - 1)
    return {
        "case": f"polynomial_dirichlet_sigma_{sigma:g}",
        "sigma": sigma,
        "n": n,
        "N": N,
        "h": h,
        "dof": N * N,
        "method": method,
        "discretization": discretization,
        "time_s_median": time_s_median,
        "time_s_iqr": time_s_iqr,
        "linf_error": linf_error,
        "l2_error": l2_error,
        "gmres_restart": gmres_restart,
        "gmres_tol_type": gmres_tol_type,
        "gmres_tol": gmres_tol,
        "gmres_iterations": gmres_iterations,
        "gmres_converged": gmres_converged,
        "gmres_final_absres": gmres_final_absres,
        "gmres_final_relres": gmres_final_relres,
        "timing_scope": timing_scope,
        "warmup": WARMUP,
        "repeat": REPEAT,
    }


def copy_to_thesis_figures(*paths):
    root = Path(__file__).resolve().parents[3]
    thesis_dir = root / "thesis" / "figures"
    thesis_dir.mkdir(parents=True, exist_ok=True)
    copied = []
    for path in paths:
        src = Path(path)
        dst = thesis_dir / src.name
        shutil.copy2(src, dst)
        copied.append(dst)
    return copied


def positive(series):
    return np.maximum(np.asarray(series, dtype=float), 1e-300)


def plot_error_time(df, figures_dir):
    from matplotlib.lines import Line2D
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(1, len(SIGMA_LIST), figsize=(13.6, 4.9), sharey=True)
    if len(SIGMA_LIST) == 1:
        axes = [axes]

    markers = {"FA": "o", "CR": "s", "FACR": "^", "FFT9": "D", "GMRES30": "X"}
    colors = {"FA": "C0", "CR": "C1", "FACR": "C2", "FFT9": "C3", "GMRES30": "C4"}
    labels = {
        "FA": "FA (5-point FFT)",
        "CR": "CR (5-point CR)",
        "FACR": "FACR-like (5-point)",
        "FFT9": "FFT9 (4th-order)",
        "GMRES30": "GMRES30 (unpreconditioned)",
    }
    method_order = ["FA", "CR", "FACR", "FFT9", "GMRES30"]
    dof_norm = LogNorm(vmin=float(df["dof"].min()), vmax=float(df["dof"].max()))
    cmap = plt.get_cmap("viridis")
    scatter_for_cbar = None

    for ax, sigma in zip(axes, SIGMA_LIST):
        sub = df[df["sigma"] == sigma]
        for method in method_order:
            md = sub[sub["method"] == method].sort_values("time_s_median")
            if md.empty:
                continue
            ax.loglog(
                positive(md["time_s_median"]),
                positive(md["linf_error"]),
                color=colors.get(method, None),
                linewidth=1.5,
                alpha=0.65,
            )
            scatter_for_cbar = ax.scatter(
                positive(md["time_s_median"]),
                positive(md["linf_error"]),
                c=md["dof"],
                cmap=cmap,
                norm=dof_norm,
                marker=markers.get(method, "o"),
                s=58,
                edgecolors=colors.get(method, "black"),
                linewidths=1.2,
                zorder=4,
            )
            if method == "GMRES30":
                failed = md[md["gmres_converged"].astype(str).str.lower() == "false"]
                if not failed.empty:
                    ax.scatter(
                        positive(failed["time_s_median"]),
                        positive(failed["linf_error"]),
                        marker="x",
                        s=90,
                        color="black",
                        linewidths=2.0,
                        label="GMRES30 not converged",
                        zorder=5,
                    )

        ax.set_title(f"$\\sigma={sigma:g}$")
        ax.set_xlabel("median solve time (s)")
        ax.grid(True, which="both", ls="--", alpha=0.35)
    axes[0].set_ylabel("$L_\\infty$ error")
    handles = [
        Line2D(
            [0],
            [0],
            marker=markers[method],
            linestyle="-",
            color=colors[method],
            markerfacecolor="white",
            markeredgecolor=colors[method],
            markeredgewidth=1.2,
            linewidth=1.5,
            markersize=7,
            label=labels[method],
        )
        for method in method_order
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="x",
            linestyle="None",
            color="black",
            markeredgewidth=2.0,
            markersize=8,
            label="GMRES30 not converged",
        )
    )
    axes[-1].legend(handles=handles, fontsize=7.5, loc="best")
    if scatter_for_cbar is not None:
        fig.subplots_adjust(left=0.08, right=0.86, bottom=0.16, top=0.82, wspace=0.12)
        cax = fig.add_axes([0.89, 0.18, 0.018, 0.62])
        cbar = fig.colorbar(scatter_for_cbar, cax=cax)
        cbar.set_label("unknowns $N^2$ (log scale)")
    fig.suptitle(
        "Log-log accuracy-cost comparison on polynomial Dirichlet manufactured solution"
    )

    png = Path(figures_dir) / f"{ERROR_TIME_FIG}.png"
    pdf = Path(figures_dir) / f"{ERROR_TIME_FIG}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_time_scaling(df, figures_dir):
    fig, axes = plt.subplots(1, len(SIGMA_LIST), figsize=(12, 4.8), sharey=True)
    if len(SIGMA_LIST) == 1:
        axes = [axes]

    markers = {"FA": "o", "CR": "s", "FACR": "^", "FFT9": "D", "GMRES30": "X"}
    colors = {"FA": "C0", "CR": "C1", "FACR": "C2", "FFT9": "C3", "GMRES30": "C4"}
    labels = {
        "FA": "FA (5-point FFT)",
        "CR": "CR (5-point CR)",
        "FACR": "FACR-like (5-point)",
        "FFT9": "FFT9 (4th-order)",
        "GMRES30": "GMRES30 (unpreconditioned)",
    }

    for ax, sigma in zip(axes, SIGMA_LIST):
        sub = df[df["sigma"] == sigma]
        for method in ["FA", "CR", "FACR", "FFT9", "GMRES30"]:
            md = sub[sub["method"] == method].sort_values("dof")
            if md.empty:
                continue
            ax.loglog(
                md["dof"],
                positive(md["time_s_median"]),
                marker=markers.get(method, "o"),
                color=colors.get(method, None),
                linewidth=1.5,
                label=labels.get(method, method),
            )

        fa = sub[sub["method"] == "FA"].sort_values("dof")
        if len(fa) >= 2:
            x = np.asarray(fa["dof"], dtype=float)
            n_int = np.sqrt(x)
            ref = x * np.log2(np.maximum(n_int, 2.0))
            scale = float(fa["time_s_median"].iloc[0]) / ref[0]
            ax.loglog(x, scale * ref, "k--", alpha=0.55, label="$O(N^2\\log N)$ ref.")

        ax.set_title(f"$\\sigma={sigma:g}$")
        ax.set_xlabel("unknowns $N^2$")
        ax.grid(True, which="both", ls="--", alpha=0.35)
    axes[0].set_ylabel("median solve-call time (s)")
    axes[-1].legend(fontsize=8)
    fig.suptitle("Time scaling of current implementations")
    fig.tight_layout()

    png = Path(figures_dir) / f"{TIME_SCALING_FIG}.png"
    pdf = Path(figures_dir) / f"{TIME_SCALING_FIG}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot(df=None, figures_dir=None):
    if df is None:
        df = pd.read_csv(Path(get_results_dir()) / CSV_NAME)
    if figures_dir is None:
        figures_dir = get_figures_dir()
    Path(figures_dir).mkdir(parents=True, exist_ok=True)

    paths = []
    paths.extend(plot_error_time(df, figures_dir))
    paths.extend(plot_time_scaling(df, figures_dir))
    copy_to_thesis_figures(*paths)
    return paths


def run():
    rows = []
    solvers = direct_solver_registry()

    for sigma in SIGMA_LIST:
        print(f"\n=== sigma={sigma:g} ===")
        for method, (disc, solver) in solvers.items():
            for n in N_LIST_DIRECT:
                print(f"  {method:7s} n={n:4d}", flush=True)
                rows.append(run_direct_method(method, disc, solver, n, sigma))

        for n in N_LIST_GMRES:
            print(f"  {'GMRES30':7s} n={n:4d}", flush=True)
            rows.append(run_gmres_method(n, sigma))

    df = pd.DataFrame(rows)
    results_dir = Path(get_results_dir())
    csv_path = results_dir / CSV_NAME
    df.to_csv(csv_path, index=False)
    print(f"\nCSV saved to {csv_path}")

    paths = plot(df, get_figures_dir())
    for path in paths:
        print(f"Figure saved to {path}")
    print("Figures copied to thesis/figures")
    return df


if __name__ == "__main__":
    run()
