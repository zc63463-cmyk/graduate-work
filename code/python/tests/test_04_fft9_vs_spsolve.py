#!/usr/bin/env python3
"""
Test 04: FFT9 vs Compact Sparse Direct Solver
===============================================

Verifies that the FFT9 (4th-order compact) solver produces the same
discrete solution as solving the compact sparse system directly.

The FFT9 solver internally solves:
    (L_h - sigma * R_h) * u = -(R_h * f) + bc_correction

where L_h is the 9-point Laplacian and R_h is the averaging operator.
This is equivalent to:
    (-L_h + sigma * R_h) * u = R_h * f

We build the sparse matrix for (L_h - sigma * R_h) and the matching RHS,
then compare with the FFT9 frequency-domain solver.

Acceptance: ||u_fft9 - u_spsolve||_inf <= 1e-10
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import spsolve
from helmholtz_solver import (
    fft9_helmholtz, fft9_oer_helmholtz,
    apply_Rh_full, compute_bc_correction_9pt
)


def _build_fft9_compact_system(n, f_func, bc_func, sigma=0.0, sx=1.0, sy=1.0):
    """
    Build the FFT9 compact sparse system matching the solver formulation:

        (L_h - sigma * R_h) * u_int = -(R_h * f) + bc_correction

    This is exactly what the FFT9 solver computes in the frequency domain.

    Returns A (sparse), b (vector).
    """
    h = sx / (n - 1)
    h2 = h * h
    N = n - 2
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')

    # Evaluate f on full grid
    F = np.asarray(f_func(X, Y), dtype=float)

    # Apply R_h to g = -f (matching solver: Rg = R_h * (-f))
    G = -F
    Rg = apply_Rh_full(G)  # shape (N, N)

    # Boundary correction (from the solver's own function)
    bc_corr, _, _, _, _ = compute_bc_correction_9pt(
        n, bc_func, x, y, h, bc_x='D', bc_y='D', sigma=sigma)
    Rg += bc_corr

    # Build sparse matrix for (L_h - sigma * R_h)
    # L_h stencil: (1/(6h²)) * [1 4 1; 4 -20 4; 1 4 1]
    # R_h stencil: [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0]
    # Combined (L_h - sigma * R_h):
    #   Center:      -20/(6h²) - sigma*(2/3)
    #   Direct nbr:   4/(6h²) - sigma*(1/12)
    #   Diagonal nbr: 1/(6h²)  (R_h has no diagonal entries)

    coeff_L = 1.0 / (6.0 * h2)

    center_val = -20.0 * coeff_L - sigma * (2.0 / 3.0)
    neighbor_val = 4.0 * coeff_L - sigma / 12.0
    diag_neighbor_val = 1.0 * coeff_L

    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j

            # Center
            rows.append(idx)
            cols.append(idx)
            vals.append(center_val)

            # Left (i-1, j) — direct neighbor
            if i > 0:
                rows.append(idx)
                cols.append((i - 1) * N + j)
                vals.append(neighbor_val)

            # Right (i+1, j) — direct neighbor
            if i < N - 1:
                rows.append(idx)
                cols.append((i + 1) * N + j)
                vals.append(neighbor_val)

            # Bottom (i, j-1) — direct neighbor
            if j > 0:
                rows.append(idx)
                cols.append(i * N + (j - 1))
                vals.append(neighbor_val)

            # Top (i, j+1) — direct neighbor
            if j < N - 1:
                rows.append(idx)
                cols.append(i * N + (j + 1))
                vals.append(neighbor_val)

            # Diagonal neighbors
            if i > 0 and j > 0:  # bottom-left
                rows.append(idx)
                cols.append((i - 1) * N + (j - 1))
                vals.append(diag_neighbor_val)

            if i > 0 and j < N - 1:  # top-left
                rows.append(idx)
                cols.append((i - 1) * N + (j + 1))
                vals.append(diag_neighbor_val)

            if i < N - 1 and j > 0:  # bottom-right
                rows.append(idx)
                cols.append((i + 1) * N + (j - 1))
                vals.append(diag_neighbor_val)

            if i < N - 1 and j < N - 1:  # top-right
                rows.append(idx)
                cols.append((i + 1) * N + (j + 1))
                vals.append(diag_neighbor_val)

    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N * N, N * N))
    b = Rg.reshape(-1)

    return A, b


class TestFFT9vsSpsolve:
    """Test FFT9 solver matches compact sparse direct solver."""

    @pytest.mark.parametrize("sigma", [0.0, 1.0, 10.0])
    @pytest.mark.parametrize("n", [9, 17])
    def test_fft9_4th_vs_spsolve(self, n, sigma):
        """FFT9 4th-order should match sparse solve of compact system."""
        def u_exact(x, y):
            return np.sin(np.pi * x) * np.sin(np.pi * y)

        def f_rhs(x, y):
            return (2 * np.pi**2 + sigma) * np.sin(np.pi * x) * np.sin(np.pi * y)

        def bc(x, y):
            return 0.0

        # FFT9 solver
        U_fft9 = fft9_helmholtz(n, f_rhs, bc, sigma=sigma, method='4th')

        # Sparse direct solver
        A, b = _build_fft9_compact_system(n, f_rhs, bc, sigma=sigma)
        u_sparse = spsolve(A, b)
        U_sparse = np.zeros((n, n))
        U_sparse[1:-1, 1:-1] = u_sparse.reshape(n - 2, n - 2)

        # Compare interior
        diff = np.max(np.abs(U_fft9[1:-1, 1:-1] - U_sparse[1:-1, 1:-1]))
        assert diff < 1e-10, \
            f"n={n}, sigma={sigma}: ||U_fft9 - U_sparse||_inf = {diff:.2e}"

    def test_fft9_nonhom_bc(self):
        """FFT9 with non-homogeneous BC should match sparse solver."""
        sigma = 5.0
        n = 17

        def u_exact(x, y):
            return np.sin(2 * np.pi * x) * np.sin(3 * np.pi * y) + 1.0

        def f_rhs(x, y):
            lam = (2 * np.pi)**2 + (3 * np.pi)**2
            return (lam + sigma) * np.sin(2 * np.pi * x) * np.sin(3 * np.pi * y) + sigma

        def bc(x, y):
            return 1.0  # u|∂Ω = 1

        U_fft9 = fft9_helmholtz(n, f_rhs, bc, sigma=sigma, method='4th')

        A, b = _build_fft9_compact_system(n, f_rhs, bc, sigma=sigma)
        u_sparse = spsolve(A, b)
        U_sparse = np.zeros((n, n))
        U_sparse[1:-1, 1:-1] = u_sparse.reshape(n - 2, n - 2)

        diff = np.max(np.abs(U_fft9[1:-1, 1:-1] - U_sparse[1:-1, 1:-1]))
        assert diff < 1e-10, \
            f"FFT9 non-hom BC: ||U_fft9 - U_sparse||_inf = {diff:.2e}"

    @pytest.mark.parametrize("n", [9, 17])
    def test_fft9_oer_vs_spsolve(self, n):
        """FFT9 OER solver should match sparse solver."""
        sigma = 0.0

        def u_exact(x, y):
            return np.sin(np.pi * x) * np.sin(np.pi * y)

        def f_rhs(x, y):
            return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)

        def bc(x, y):
            return 0.0

        U_oer = fft9_oer_helmholtz(n, f_rhs, bc, sigma=sigma)

        A, b = _build_fft9_compact_system(n, f_rhs, bc, sigma=sigma)
        u_sparse = spsolve(A, b)
        U_sparse = np.zeros((n, n))
        U_sparse[1:-1, 1:-1] = u_sparse.reshape(n - 2, n - 2)

        diff = np.max(np.abs(U_oer[1:-1, 1:-1] - U_sparse[1:-1, 1:-1]))
        assert diff < 1e-10, \
            f"n={n}: ||U_oer - U_sparse||_inf = {diff:.2e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
