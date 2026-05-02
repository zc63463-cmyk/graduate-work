#!/usr/bin/env python3
"""Generate visualization-only supplement figures for Chapter 6.

This script does not create new experiment CSV files. It only produces two PNG
figures used as visual supplements for existing experiments.
"""

import shutil
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from helmholtz_solver import fft9_helmholtz
from .utils import get_figures_dir, test_problem_nonhom_dirichlet


def _copy_to_thesis(fig_path):
    repo_root = Path(__file__).resolve().parents[3]
    thesis_figures = repo_root / "thesis" / "figures"
    thesis_figures.mkdir(parents=True, exist_ok=True)
    target = thesis_figures / Path(fig_path).name
    shutil.copy2(fig_path, target)
    return target


def generate_temperature_field_comparison():
    sigma = 10.0
    n_grid = 129
    sx = sy = 1.0
    u_exact, f_rhs, bc = test_problem_nonhom_dirichlet(
        sigma, m=2, n=3, sx=sx, sy=sy
    )

    x = np.linspace(0.0, sx, n_grid)
    y = np.linspace(0.0, sy, n_grid)
    X, Y = np.meshgrid(x, y, indexing="ij")
    U_exact = u_exact(X, Y)
    U_num = fft9_helmholtz(
        n_grid, f_rhs, bc, sigma=sigma, bc_type="dirichlet", sx=sx, sy=sy
    )
    abs_err = np.abs(U_num - U_exact)
    log_err = np.log10(abs_err + 1e-16)
    max_err = np.max(abs_err)
    field_vmin = min(np.min(U_exact), np.min(U_num))
    field_vmax = max(np.max(U_exact), np.max(U_num))

    j_mid = int(np.argmin(np.abs(y - 0.5)))

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    extent = [0.0, sx, 0.0, sy]

    im0 = axes[0, 0].imshow(
        U_exact.T,
        origin="lower",
        extent=extent,
        aspect="equal",
        vmin=field_vmin,
        vmax=field_vmax,
    )
    axes[0, 0].set_title("(a) analytical temperature field")
    axes[0, 0].set_xlabel("x")
    axes[0, 0].set_ylabel("y")
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)

    im1 = axes[0, 1].imshow(
        U_num.T,
        origin="lower",
        extent=extent,
        aspect="equal",
        vmin=field_vmin,
        vmax=field_vmax,
    )
    axes[0, 1].set_title("(b) FFT9 numerical solution")
    axes[0, 1].set_xlabel("x")
    axes[0, 1].set_ylabel("y")
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)

    im2 = axes[1, 0].imshow(log_err.T, origin="lower", extent=extent, aspect="equal")
    axes[1, 0].set_title(
        f"(c) log10(|u_h - u_exact| + 1e-16), max={max_err:.2e}"
    )
    axes[1, 0].set_xlabel("x")
    axes[1, 0].set_ylabel("y")
    cbar = fig.colorbar(im2, ax=axes[1, 0], fraction=0.046, pad=0.04)
    cbar.set_label("log10(error + 1e-16)")

    axes[1, 1].plot(x, U_exact[:, j_mid], label="exact", linewidth=1.8)
    axes[1, 1].plot(x, U_num[:, j_mid], "--", label="FFT9", linewidth=1.5)
    axes[1, 1].set_title("(d) line cut at y = 0.5")
    axes[1, 1].set_xlabel("x")
    axes[1, 1].set_ylabel("u(x, 0.5)")
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()

    fig.suptitle(
        "Nonhomogeneous Dirichlet analytical temperature field, sigma=10, n=129",
        fontsize=12,
    )
    fig.tight_layout()

    figures_dir = Path(get_figures_dir())
    fig_path = figures_dir / "exp02_temperature_field_comparison.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    thesis_path = _copy_to_thesis(fig_path)
    return fig_path, thesis_path


def _dirichlet_laplacian_eigenvalues(n_grid, sx=1.0):
    N = n_grid - 2
    h = sx / (n_grid - 1)
    modes = np.arange(1, N + 1)
    mu = (4.0 / (h * h)) * np.sin(modes * np.pi / (2.0 * (N + 1))) ** 2
    return np.add.outer(mu, mu), modes, h


def generate_denominator_heatmap():
    n_grid = 65
    p_target = 2
    q_target = 3
    delta = 1e-3
    lambda_2d, modes, h = _dirichlet_laplacian_eigenvalues(n_grid)
    lambda_target = lambda_2d[p_target - 1, q_target - 1]
    kappa2 = lambda_target + delta
    denominator = lambda_2d - kappa2
    abs_denominator = np.abs(denominator)
    resonance_strength = -np.log10(abs_denominator)

    min_value = np.min(abs_denominator)
    min_locations = np.argwhere(np.isclose(abs_denominator, min_value, rtol=1e-10))
    min_modes = [
        (int(modes[p_idx]), int(modes[q_idx])) for p_idx, q_idx in min_locations
    ]
    min_modes_text = ", ".join(f"({p},{q})" for p, q in min_modes)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    N = n_grid - 2
    im = ax.imshow(
        resonance_strength,
        origin="lower",
        extent=[0.5, N + 0.5, 0.5, N + 0.5],
        aspect="equal",
    )
    ax.set_title(
        "true Helmholtz near-resonance denominator heatmap\n"
        "selected modal pair (2,3)/(3,2), delta=1e-3"
    )
    ax.set_xlabel("q mode")
    ax.set_ylabel("p mode")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("resonance strength = -log10|d_pq|")

    for p_idx, q_idx in min_locations:
        p_mode = int(modes[p_idx])
        q_mode = int(modes[q_idx])
        ax.scatter(
            q_mode,
            p_mode,
            s=90,
            facecolors="none",
            edgecolors="red",
            linewidths=1.6,
        )

    ax.text(
        0.98,
        0.98,
        f"n={n_grid}, h={h:.5f}, min |d_pq|={min_value:.1e}\n"
        f"min modes: {min_modes_text}\n"
        "larger value = stronger near-resonance",
        transform=ax.transAxes,
        fontsize=8,
        color="white",
        ha="right",
        va="top",
        bbox=dict(facecolor="black", alpha=0.45, edgecolor="none"),
    )

    inset_n = 10
    inset = ax.inset_axes([0.57, 0.08, 0.38, 0.38])
    inset.imshow(
        resonance_strength[:inset_n, :inset_n],
        origin="lower",
        extent=[0.5, inset_n + 0.5, 0.5, inset_n + 0.5],
        aspect="equal",
        vmin=np.min(resonance_strength),
        vmax=np.max(resonance_strength),
    )
    inset.set_title("p,q <= 10", fontsize=8)
    inset.set_xlabel("q", fontsize=7)
    inset.set_ylabel("p", fontsize=7)
    inset.tick_params(axis="both", labelsize=7)
    for p_mode, q_mode in min_modes:
        if p_mode <= inset_n and q_mode <= inset_n:
            inset.scatter(
                q_mode,
                p_mode,
                s=55,
                facecolors="none",
                edgecolors="red",
                linewidths=1.2,
            )

    fig.tight_layout()

    figures_dir = Path(get_figures_dir())
    fig_path = figures_dir / "exp05_denominator_heatmap.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    thesis_path = _copy_to_thesis(fig_path)
    return fig_path, thesis_path


def main():
    generated = [
        generate_temperature_field_comparison(),
        generate_denominator_heatmap(),
    ]
    for code_path, thesis_path in generated:
        print(f"Saved {code_path}")
        print(f"Copied to {thesis_path}")


if __name__ == "__main__":
    main()
