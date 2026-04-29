#!/usr/bin/env python3
"""Complex manufactured solution visualization supplement.

This script is a visualization supplement for Chapter 6, not a new core
experiment section. It first verifies grid consistency, boundary values, and
convergence behavior, then generates trustworthy field and convergence figures.
"""

import csv
import os
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
from helmholtz_solver import fa_helmholtz, fft9_helmholtz

from .utils import compute_convergence_rate, get_figures_dir, get_results_dir


SIGMA = 10.0
SX = 1.0
SY = 1.0
PLOT_N = 129
CONVERGENCE_NS = [33, 65, 129, 257]
SANITY_N = 65

CSV_NAME = "exp06_complex_manufactured_convergence.csv"
FIELD_FIG_NAME = "exp06_complex_manufactured_fields.png"
CONVERGENCE_FIG_NAME = "exp06_complex_manufactured_convergence.png"


def u_exact(x, y):
    """Analytical multimode manufactured solution."""
    return (
        np.sin(np.pi * x) * np.sin(np.pi * y)
        + 0.5 * np.sin(2.0 * np.pi * x) * np.sin(3.0 * np.pi * y)
        + 0.3 * np.sin(5.0 * np.pi * x) * np.sin(4.0 * np.pi * y)
    )


def f_rhs(x, y):
    """RHS for (-Delta + sigma)u=f, with per-mode eigenvalues."""
    return (
        (2.0 * np.pi**2 + SIGMA) * np.sin(np.pi * x) * np.sin(np.pi * y)
        + 0.5
        * (13.0 * np.pi**2 + SIGMA)
        * np.sin(2.0 * np.pi * x)
        * np.sin(3.0 * np.pi * y)
        + 0.3
        * (41.0 * np.pi**2 + SIGMA)
        * np.sin(5.0 * np.pi * x)
        * np.sin(4.0 * np.pi * y)
    )


def bc_zero(x, y):
    """Homogeneous Dirichlet boundary value."""
    return 0.0


def full_grid(n):
    x = np.linspace(0.0, SX, n)
    y = np.linspace(0.0, SY, n)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return x, y, X, Y


def discrete_l2(error, h):
    return h * np.sqrt(np.sum(error**2))


def boundary_max_abs(U):
    return max(
        float(np.max(np.abs(U[0, :]))),
        float(np.max(np.abs(U[-1, :]))),
        float(np.max(np.abs(U[:, 0]))),
        float(np.max(np.abs(U[:, -1]))),
    )


def solve_on_grid(n):
    x, y, X, Y = full_grid(n)
    U_exact = u_exact(X, Y)
    U_fa = fa_helmholtz(
        n, f_rhs, bc_zero, sigma=SIGMA, bc_type="dirichlet", sx=SX, sy=SY
    )
    U_fft9 = fft9_helmholtz(
        n, f_rhs, bc_zero, sigma=SIGMA, bc_type="dirichlet", sx=SX, sy=SY
    )
    validate_grid_consistency(n, U_exact, U_fa, U_fft9)
    validate_plot_grid_orientation(x, y, U_exact)
    return x, y, U_exact, U_fa, U_fft9


def validate_grid_consistency(n, U_exact, U_fa, U_fft9):
    expected_shape = (n, n)
    arrays = {
        "exact": U_exact,
        "FA": U_fa,
        "FFT9": U_fft9,
    }
    for name, U in arrays.items():
        if U.shape != expected_shape:
            raise RuntimeError(f"{name} shape {U.shape} != expected {expected_shape}")

    interior_shape = U_exact[1:-1, 1:-1].shape
    if interior_shape != (n - 2, n - 2):
        raise RuntimeError(
            f"Interior shape {interior_shape} != expected {(n - 2, n - 2)}"
        )

    exact_boundary = boundary_max_abs(U_exact)
    fa_boundary = boundary_max_abs(U_fa)
    fft9_boundary = boundary_max_abs(U_fft9)
    if exact_boundary > 1e-12:
        raise RuntimeError(f"Exact boundary is not zero: {exact_boundary:.3e}")
    if fa_boundary > 1e-12:
        raise RuntimeError(f"FA boundary is not zero: {fa_boundary:.3e}")
    if fft9_boundary > 1e-12:
        raise RuntimeError(f"FFT9 boundary is not zero: {fft9_boundary:.3e}")


def validate_plot_grid_orientation(x, y, U_exact):
    """Check that U[i,j] represents u(x_i,y_j), then imshow uses U.T."""
    n = U_exact.shape[0]
    i = min(3, n - 2)
    j = min(5, n - 2)
    expected = u_exact(x[i], y[j])
    if not np.isclose(U_exact[i, j], expected, rtol=0.0, atol=1e-14):
        raise RuntimeError("Grid indexing mismatch: U[i,j] is not u(x_i,y_j).")
    if U_exact.T.shape != (len(y), len(x)):
        raise RuntimeError("Plot transpose check failed for imshow orientation.")


def compute_error_metrics(U, U_exact, h):
    error = U - U_exact
    return {
        "linf": float(np.max(np.abs(error))),
        "l2": float(discrete_l2(error, h)),
    }


def run_sanity_check():
    x, y, U_exact, U_fa, U_fft9 = solve_on_grid(SANITY_N)
    h = SX / (SANITY_N - 1)
    fa = compute_error_metrics(U_fa, U_exact, h)
    fft9 = compute_error_metrics(U_fft9, U_exact, h)

    if fa["linf"] >= 0.1:
        raise RuntimeError(
            f"Sanity check failed: FA max error is suspiciously large "
            f"({fa['linf']:.3e})."
        )
    if fft9["linf"] >= 0.1:
        raise RuntimeError(
            f"Sanity check failed: FFT9 max error is suspiciously large "
            f"({fft9['linf']:.3e})."
        )
    if fft9["linf"] >= fa["linf"]:
        raise RuntimeError(
            f"Sanity check failed: FFT9 max error ({fft9['linf']:.3e}) "
            f"is not smaller than FA ({fa['linf']:.3e})."
        )

    return {
        "n": SANITY_N,
        "h": h,
        "fa_linf": fa["linf"],
        "fa_l2": fa["l2"],
        "fft9_linf": fft9["linf"],
        "fft9_l2": fft9["l2"],
    }


def run_convergence_study():
    rows = []
    fa_linf_errors = []
    fa_l2_errors = []
    fft9_linf_errors = []
    fft9_l2_errors = []

    for n in CONVERGENCE_NS:
        h = SX / (n - 1)
        _, _, U_exact, U_fa, U_fft9 = solve_on_grid(n)
        fa = compute_error_metrics(U_fa, U_exact, h)
        fft9 = compute_error_metrics(U_fft9, U_exact, h)

        fa_linf_errors.append(fa["linf"])
        fa_l2_errors.append(fa["l2"])
        fft9_linf_errors.append(fft9["linf"])
        fft9_l2_errors.append(fft9["l2"])
        rows.append(
            {
                "n": n,
                "h": h,
                "fa_linf": fa["linf"],
                "fa_l2": fa["l2"],
                "fft9_linf": fft9["linf"],
                "fft9_l2": fft9["l2"],
            }
        )

    rate_columns = {
        "fa_order_linf": compute_convergence_rate(fa_linf_errors),
        "fa_order_l2": compute_convergence_rate(fa_l2_errors),
        "fft9_order_linf": compute_convergence_rate(fft9_linf_errors),
        "fft9_order_l2": compute_convergence_rate(fft9_l2_errors),
    }
    for i, row in enumerate(rows):
        for key, values in rate_columns.items():
            row[key] = values[i]

    fit_orders = {
        "fa_linf": fit_order([r["h"] for r in rows], fa_linf_errors),
        "fa_l2": fit_order([r["h"] for r in rows], fa_l2_errors),
        "fft9_linf": fit_order([r["h"] for r in rows], fft9_linf_errors),
        "fft9_l2": fit_order([r["h"] for r in rows], fft9_l2_errors),
    }
    validate_convergence(fit_orders)
    return rows, fit_orders


def fit_order(h_values, errors):
    h_values = np.asarray(h_values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    mask = (h_values > 0.0) & (errors > 0.0)
    return float(np.polyfit(np.log(h_values[mask]), np.log(errors[mask]), 1)[0])


def validate_convergence(fit_orders):
    checks = [
        ("FA Linf", fit_orders["fa_linf"], 1.8, 2.2),
        ("FA L2", fit_orders["fa_l2"], 1.8, 2.2),
        ("FFT9 Linf", fit_orders["fft9_linf"], 3.7, 4.3),
        ("FFT9 L2", fit_orders["fft9_l2"], 3.7, 4.3),
    ]
    for label, value, lo, hi in checks:
        if not (lo <= value <= hi):
            raise RuntimeError(
                f"Convergence check failed for {label}: order={value:.3f}, "
                f"expected in [{lo:.1f}, {hi:.1f}]."
            )


def save_convergence_csv(rows):
    results_dir = Path(get_results_dir())
    csv_path = results_dir / CSV_NAME
    fieldnames = [
        "n",
        "h",
        "fa_linf",
        "fa_l2",
        "fft9_linf",
        "fft9_l2",
        "fa_order_linf",
        "fa_order_l2",
        "fft9_order_linf",
        "fft9_order_l2",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return csv_path


def copy_to_thesis(fig_path):
    repo_root = Path(__file__).resolve().parents[3]
    thesis_figures = repo_root / "thesis" / "figures"
    thesis_figures.mkdir(parents=True, exist_ok=True)
    target = thesis_figures / Path(fig_path).name
    shutil.copy2(fig_path, target)
    return target


def plot_fields():
    x, y, U_exact, U_fa, U_fft9 = solve_on_grid(PLOT_N)
    h = SX / (PLOT_N - 1)
    fa_error = np.abs(U_fa - U_exact)
    fft9_error = np.abs(U_fft9 - U_exact)
    fa_metrics = compute_error_metrics(U_fa, U_exact, h)
    fft9_metrics = compute_error_metrics(U_fft9, U_exact, h)

    field_vmin = float(min(np.min(U_exact), np.min(U_fa), np.min(U_fft9)))
    field_vmax = float(max(np.max(U_exact), np.max(U_fa), np.max(U_fft9)))
    log_fa_error = np.log10(fa_error + 1e-16)
    log_fft9_error = np.log10(fft9_error + 1e-16)
    err_vmin = float(min(np.min(log_fa_error), np.min(log_fft9_error)))
    err_vmax = float(max(np.max(log_fa_error), np.max(log_fft9_error)))
    extent = [0.0, SX, 0.0, SY]
    j_mid = int(np.argmin(np.abs(y - 0.5)))

    fig, axes = plt.subplots(2, 3, figsize=(14, 8.2))
    field_panels = [
        (axes[0, 0], U_exact, "(a) exact complex field"),
        (axes[0, 1], U_fa, "(b) FA numerical solution"),
        (axes[0, 2], U_fft9, "(c) FFT9 numerical solution"),
    ]
    for ax, data, title in field_panels:
        im = ax.imshow(
            data.T,
            origin="lower",
            extent=extent,
            aspect="equal",
            vmin=field_vmin,
            vmax=field_vmax,
        )
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    error_panels = [
        (
            axes[1, 0],
            log_fa_error,
            f"(d) FA error (max={fa_metrics['linf']:.2e})",
        ),
        (
            axes[1, 1],
            log_fft9_error,
            f"(e) FFT9 error (max={fft9_metrics['linf']:.2e})",
        ),
    ]
    for ax, data, title in error_panels:
        im = ax.imshow(
            data.T,
            origin="lower",
            extent=extent,
            aspect="equal",
            vmin=err_vmin,
            vmax=err_vmax,
        )
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("log10(error + 1e-16)")

    ax = axes[1, 2]
    ax.plot(x, U_exact[:, j_mid], label="exact", linewidth=1.8)
    ax.plot(x, U_fa[:, j_mid], "--", label="FA", linewidth=1.4)
    ax.plot(x, U_fft9[:, j_mid], ":", label="FFT9", linewidth=1.8)
    ax.set_title("(f) line cut at y = 0.5")
    ax.set_xlabel("x")
    ax.set_ylabel("u(x, 0.5)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.suptitle(
        "Complex multimode manufactured solution, sigma=10, n=129", fontsize=12
    )
    fig.tight_layout()

    fig_path = Path(get_figures_dir()) / FIELD_FIG_NAME
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path, copy_to_thesis(fig_path)


def plot_convergence(rows):
    h = np.array([row["h"] for row in rows])
    fa_linf = np.array([row["fa_linf"] for row in rows])
    fa_l2 = np.array([row["fa_l2"] for row in rows])
    fft9_linf = np.array([row["fft9_linf"] for row in rows])
    fft9_l2 = np.array([row["fft9_l2"] for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    plot_error_panel(axes[0], h, fa_linf, fft9_linf, "Linf error")
    plot_error_panel(axes[1], h, fa_l2, fft9_l2, "L2 error")
    fig.suptitle("Complex manufactured solution convergence")
    fig.tight_layout()

    fig_path = Path(get_figures_dir()) / CONVERGENCE_FIG_NAME
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path, copy_to_thesis(fig_path)


def plot_error_panel(ax, h, fa_error, fft9_error, title):
    ax.loglog(h, fa_error, "o-", label="FA")
    ax.loglog(h, fft9_error, "s-", label="FFT9")

    h_ref = np.array([np.min(h), np.max(h)])
    fa_scale = fa_error[-1] / (h[-1] ** 2)
    fft9_scale = fft9_error[-1] / (h[-1] ** 4)
    ax.loglog(h_ref, fa_scale * h_ref**2, "k--", alpha=0.35, label="O(h^2)")
    ax.loglog(h_ref, fft9_scale * h_ref**4, "k:", alpha=0.45, label="O(h^4)")

    ax.set_xlabel("h")
    ax.set_ylabel("error")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()


def run():
    print("Complex manufactured solution supplement")
    print("PDE: (-Delta + sigma)u=f, sigma=10, homogeneous Dirichlet")
    print("Modes: (1,1), (2,3), (5,4)")

    sanity = run_sanity_check()
    print(
        "Sanity n={n}: FA Linf={fa_linf:.3e}, FA L2={fa_l2:.3e}, "
        "FFT9 Linf={fft9_linf:.3e}, FFT9 L2={fft9_l2:.3e}".format(**sanity)
    )

    rows, fit_orders = run_convergence_study()
    csv_path = save_convergence_csv(rows)
    field_fig, field_thesis = plot_fields()
    conv_fig, conv_thesis = plot_convergence(rows)

    print(f"CSV saved to {csv_path}")
    print(f"Field figure saved to {field_fig}")
    print(f"Field figure copied to {field_thesis}")
    print(f"Convergence figure saved to {conv_fig}")
    print(f"Convergence figure copied to {conv_thesis}")
    print(
        "Fitted orders: FA Linf={fa_linf:.3f}, FA L2={fa_l2:.3f}, "
        "FFT9 Linf={fft9_linf:.3f}, FFT9 L2={fft9_l2:.3f}".format(
            **fit_orders
        )
    )

    return {
        "sanity": sanity,
        "rows": rows,
        "fit_orders": fit_orders,
        "csv_path": csv_path,
        "field_fig": field_fig,
        "field_thesis": field_thesis,
        "conv_fig": conv_fig,
        "conv_thesis": conv_thesis,
    }


if __name__ == "__main__":
    run()
