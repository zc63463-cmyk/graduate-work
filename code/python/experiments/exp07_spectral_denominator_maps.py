#!/usr/bin/env python3
"""Experiment 07 supplement: spectral denominator risk maps.

Visualizes the five-point Dirichlet frequency denominator

    d_{p,q} = lambda_{p,q} + sigma

through the small-denominator risk

    R_{p,q} = -log10(|d_{p,q}|)

for modified Helmholtz, true Helmholtz away from resonance, and true Helmholtz
near a degenerate modal pair. The output supports the Chapter 4 spectral
classification and small-denominator discussion.

The risk map is not a condition-number plot. Only for the Dirichlet five-point
system, which is orthogonally diagonalized by DST-I, does the denominator ratio
max(abs(d_pq)) / min(abs(d_pq)) coincide with cond_2(A) for nonsingular A.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

try:
    from .utils import get_figures_dir, get_results_dir  # noqa: E402
except ImportError:  # Support direct execution as a file with PYTHONPATH=code/python.
    from experiments.utils import get_figures_dir, get_results_dir  # noqa: E402


N_GRID = 129
CSV_NAME = "exp07_spectral_denominator_summary.csv"
FIG_STEM = "exp07_spectral_denominator_heatmaps"
RISK_EPS = 1e-16
LOW_MODE_LIMIT = 12
SORTED_COUNT = 20


def dirichlet_spectrum(n=N_GRID):
    """Return p/q mesh and five-point Dirichlet eigenvalues."""
    N = n - 2
    h = 1.0 / (n - 1)
    modes = np.arange(1, N + 1)
    mu = 2.0 - 2.0 * np.cos(np.pi * modes / (N + 1))
    P, Q = np.meshgrid(modes, modes, indexing="ij")
    lam = (mu[:, None] + mu[None, :]) / h**2
    return modes, P, Q, lam, h


def value_at_mode(lam, p, q):
    """Access lambda_{p,q} with 1-based p,q."""
    return float(lam[p - 1, q - 1])


def build_cases(lam):
    lambda_12 = value_at_mode(lam, 1, 2)
    lambda_22 = value_at_mode(lam, 2, 2)
    kappa2_safe = 0.5 * (lambda_12 + lambda_22)
    lambda_target = value_at_mode(lam, 2, 3)
    delta = 1e-2
    return [
        {
            "case": "modified_sigma_plus_100",
            "title": "modified",
            "subtitle": "all positive, low small-denominator risk",
            "sigma": 100.0,
            "kappa2": np.nan,
            "target_modes": "",
            "delta": np.nan,
            "markers": [],
        },
        {
            "case": "true_away_from_resonance",
            "title": "true-away",
            "subtitle": "sign-changing but low risk",
            "sigma": -kappa2_safe,
            "kappa2": kappa2_safe,
            "target_modes": "",
            "delta": np.nan,
            "markers": [],
        },
        {
            "case": "true_near_resonance_23_32",
            "title": "true-near",
            "subtitle": "hotspots at (2,3)/(3,2)",
            "sigma": -(lambda_target + delta),
            "kappa2": lambda_target + delta,
            "target_modes": "(2,3);(3,2)",
            "delta": delta,
            "markers": [(2, 3), (3, 2)],
        },
    ]


def summarize_case(case, lam, d, h):
    abs_d = np.abs(d)
    risk = -np.log10(np.maximum(abs_d, RISK_EPS))
    nearest_idx = np.unravel_index(np.argmin(abs_d), abs_d.shape)
    nearest_p = int(nearest_idx[0] + 1)
    nearest_q = int(nearest_idx[1] + 1)
    nearest_lambda = float(lam[nearest_idx])
    d_min_abs = float(np.min(abs_d))
    d_max_abs = float(np.max(abs_d))
    delta = case["delta"]
    near_zero_tol = max(1e-10, 10.0 * abs(delta)) if np.isfinite(delta) else 1e-8

    return {
        "case": case["case"],
        "n": N_GRID,
        "N": N_GRID - 2,
        "h": h,
        "sigma": float(case["sigma"]),
        "kappa2": case["kappa2"],
        "lambda_min": float(np.min(lam)),
        "lambda_max": float(np.max(lam)),
        "d_min": float(np.min(d)),
        "d_max": float(np.max(d)),
        "d_min_abs": d_min_abs,
        "d_max_abs": d_max_abs,
        "spectral_indicator": d_max_abs / d_min_abs,
        "nearest_p": nearest_p,
        "nearest_q": nearest_q,
        "nearest_lambda": nearest_lambda,
        "distance_to_nearest_lambda": abs(float(case["kappa2"]) - nearest_lambda)
        if np.isfinite(case["kappa2"])
        else np.nan,
        "num_positive": int(np.sum(d > 0.0)),
        "num_negative": int(np.sum(d < 0.0)),
        "num_near_zero": int(np.sum(abs_d <= near_zero_tol)),
        "target_modes": case["target_modes"],
        "delta": delta,
        "risk_min": float(np.min(risk)),
        "risk_max": float(np.max(risk)),
        "risk_at_nearest_mode": float(risk[nearest_idx]),
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


def sorted_denominators(d, count=SORTED_COUNT):
    """Return the smallest |d_{p,q}| values and their 1-based mode indices."""
    abs_d = np.abs(d)
    order = np.argsort(abs_d, axis=None)[:count]
    indices = np.array(np.unravel_index(order, abs_d.shape)).T
    modes = [(int(i + 1), int(j + 1)) for i, j in indices]
    values = np.array([float(abs_d[i, j]) for i, j in indices])
    return values, modes


def plot_heatmaps(cases, lam, figures_dir):
    fig = plt.figure(figsize=(15, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 0.9])
    heatmap_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    sort_ax = fig.add_subplot(gs[1, :])

    risk_maps = []
    for case in cases:
        d = lam + case["sigma"]
        risk_maps.append(-np.log10(np.maximum(np.abs(d), RISK_EPS)))

    zoom = LOW_MODE_LIMIT
    vmin = min(float(np.min(risk[:zoom, :zoom])) for risk in risk_maps)
    vmax = max(float(np.max(risk[:zoom, :zoom])) for risk in risk_maps)
    zoom_modes = np.arange(1, zoom + 1)
    Xz, Yz = np.meshgrid(zoom_modes, zoom_modes, indexing="xy")
    heatmap_image = None

    for ax, case, risk in zip(heatmap_axes, cases, risk_maps):
        d = lam + case["sigma"]
        d_zoom = d[:zoom, :zoom]
        heatmap_image = ax.imshow(
            risk[:zoom, :zoom].T,
            origin="lower",
            extent=[0.5, zoom + 0.5, 0.5, zoom + 0.5],
            aspect="auto",
            cmap="inferno",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(case["title"])
        ax.text(
            0.04,
            0.94,
            f"{case['subtitle']}\nmin |d|={np.min(np.abs(d)):.2e}",
            transform=ax.transAxes,
            va="top",
            fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="0.6", alpha=0.82, boxstyle="round,pad=0.25"),
        )
        ax.set_xlabel("p")
        ax.set_ylabel("q")
        ax.set_xlim(0.5, zoom + 0.5)
        ax.set_ylim(0.5, zoom + 0.5)
        ax.set_xticks([1, 3, 6, 9, 12])
        ax.set_yticks([1, 3, 6, 9, 12])
        if float(np.min(d_zoom)) < 0.0 < float(np.max(d_zoom)):
            ax.contour(Xz, Yz, d_zoom.T, levels=[0.0], colors="black", linewidths=1.0)
        for p, q in case["markers"]:
            ax.plot(p, q, marker="o", ms=11, mfc="#2ecc40", mec="black", mew=1.2)
        if case["markers"]:
            ax.text(
                0.05,
                0.78,
                "green markers:\n(2,3), (3,2)",
                transform=ax.transAxes,
                fontsize=9,
                bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
            )
        elif float(np.min(d_zoom)) < 0.0 < float(np.max(d_zoom)):
            ax.text(
                0.05,
                0.08,
                "interpolated sign\nboundary $d=0$",
                transform=ax.transAxes,
                fontsize=9,
                bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
            )

    cbar = fig.colorbar(heatmap_image, ax=heatmap_axes, shrink=0.92, location="right")
    cbar.set_label(r"risk = $-\log_{10}|d_{p,q}|$")

    colors = {
        "modified_sigma_plus_100": "#1f77b4",
        "true_away_from_resonance": "#2ca02c",
        "true_near_resonance_23_32": "#d62728",
    }
    labels = {
        "modified_sigma_plus_100": "modified",
        "true_away_from_resonance": "true-away",
        "true_near_resonance_23_32": "true-near",
    }
    ranks = np.arange(1, SORTED_COUNT + 1)
    near_modes = []
    near_values = []
    for case in cases:
        d = lam + case["sigma"]
        values, modes = sorted_denominators(d)
        sort_ax.plot(
            ranks,
            values,
            marker="o",
            lw=1.8,
            color=colors[case["case"]],
            label=labels[case["case"]],
        )
        if case["case"] == "true_near_resonance_23_32":
            near_modes = modes[:2]
            near_values = values[:2]
            sort_ax.scatter(
                ranks[:2],
                values[:2],
                s=180,
                marker="*",
                color="#2ecc40",
                edgecolor="black",
                linewidth=1.0,
                zorder=5,
                label="near targets (2,3)/(3,2)",
            )

    if len(near_values) >= 2:
        sort_ax.annotate(
            r"(2,3)/(3,2), $|d|\approx 10^{-2}$",
            xy=(1.5, float(np.mean(near_values[:2]))),
            xytext=(3.1, float(np.mean(near_values[:2])) * 2.8),
            arrowprops=dict(arrowstyle="->", lw=1.0),
            fontsize=10,
            bbox=dict(facecolor="white", edgecolor="0.5", alpha=0.88, boxstyle="round,pad=0.25"),
        )

    sort_ax.set_yscale("log")
    sort_ax.axhline(1e-2, color="0.35", ls="--", lw=1.0, alpha=0.85, label=r"$|d|=10^{-2}$")
    sort_ax.set_xlim(0.7, SORTED_COUNT + 0.6)
    sort_ax.set_xticks([1, 5, 10, 15, 20])
    sort_ax.set_xlabel("rank k")
    sort_ax.set_ylabel(r"sorted $|d_{p,q}|$")
    sort_ax.set_title("smallest spectral denominators")
    sort_ax.grid(True, which="both", ls=":", alpha=0.45)
    sort_ax.legend(loc="best", ncol=4)

    fig.suptitle(
        "Five-point Dirichlet small-denominator risk "
        r"$R_{p,q}=-\log_{10}|d_{p,q}|$, n=129"
    )
    png = Path(figures_dir) / f"{FIG_STEM}.png"
    pdf = Path(figures_dir) / f"{FIG_STEM}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    copy_to_thesis_figures(png, pdf)
    return png, pdf


def run():
    _, _, _, lam, h = dirichlet_spectrum(N_GRID)
    cases = build_cases(lam)
    rows = []
    for case in cases:
        d = lam + case["sigma"]
        rows.append(summarize_case(case, lam, d, h))

    df = pd.DataFrame(rows)
    results_dir = Path(get_results_dir())
    csv_path = results_dir / CSV_NAME
    df.to_csv(csv_path, index=False)

    figures_dir = Path(get_figures_dir())
    png, pdf = plot_heatmaps(cases, lam, figures_dir)
    print(f"CSV saved to {csv_path}")
    print(f"Figure saved to {png}")
    print(f"Figure saved to {pdf}")
    print("Figures copied to thesis/figures")
    return df


if __name__ == "__main__":
    run()
