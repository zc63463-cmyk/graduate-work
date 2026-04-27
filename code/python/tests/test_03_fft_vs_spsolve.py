#!/usr/bin/env python3
"""
Test 03: FFT Direct Solver vs scipy.sparse.linalg.spsolve
==========================================================

Verifies that FFT-based solvers solve the SAME discrete system
as a sparse direct solver. This is the fundamental consistency check.

For 5-point Dirichlet:
    Build the sparse matrix A_h and RHS b using the SAME stencil,
    solve with spsolve, compare with FA/CR/FACR.

For 5-point Neumann (modified Helmholtz, k²>0):
    Build ghost-point Neumann matrix, solve with spsolve,
    compare with FA/CR/FACR (after d_scale rescaling).

Acceptance: ||u_fft - u_spsolve||_inf <= 1e-10
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import spsolve
from helmholtz_solver import (
    fa_helmholtz, cr_helmholtz, facr_helmholtz,
    _helmholtz_test_problem_dirichlet,
    _helmholtz_test_problem_neumann,
)
from gmres_solver import build_helmholtz_matrix


def _build_5pt_dirichlet_system(n, f_func, bc_func, sigma=0.0, sx=1.0, sy=1.0):
    """
    Build the 5-point Dirichlet sparse system manually.

    Returns A (sparse), b (vector), and grid info.
    """
    h = sx / (n - 1)
    h2 = h * h
    N = n - 2
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)

    # Build sparse matrix
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # Diagonal: 4/h² + sigma
            rows.append(idx)
            cols.append(idx)
            vals.append(4.0 / h2 + sigma)

            # Left neighbor
            if i > 0:
                rows.append(idx)
                cols.append((i - 1) * N + j)
                vals.append(-1.0 / h2)

            # Right neighbor
            if i < N - 1:
                rows.append(idx)
                cols.append((i + 1) * N + j)
                vals.append(-1.0 / h2)

            # Bottom neighbor
            if j > 0:
                rows.append(idx)
                cols.append(i * N + (j - 1))
                vals.append(-1.0 / h2)

            # Top neighbor
            if j < N - 1:
                rows.append(idx)
                cols.append(i * N + (j + 1))
                vals.append(-1.0 / h2)

    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N * N, N * N))

    # Build RHS
    x_int = x[1:-1]
    y_int = y[1:-1]
    X, Y = np.meshgrid(x_int, y_int, indexing='ij')
    F = np.asarray(f_func(X, Y), dtype=float)

    # Add boundary contributions (F += bc/h²)
    bc_l = bc_func(x[0], y)
    bc_r = bc_func(x[-1], y)
    bc_b = bc_func(x, y[0])
    bc_t = bc_func(x, y[-1])

    bc_l = np.broadcast_to(np.atleast_1d(bc_l), (n,)).astype(float)
    bc_r = np.broadcast_to(np.atleast_1d(bc_r), (n,)).astype(float)
    bc_b = np.broadcast_to(np.atleast_1d(bc_b), (n,)).astype(float)
    bc_t = np.broadcast_to(np.atleast_1d(bc_t), (n,)).astype(float)

    F[0, :] += bc_l[1:-1] / h2
    F[-1, :] += bc_r[1:-1] / h2
    F[:, 0] += bc_b[1:-1] / h2
    F[:, -1] += bc_t[1:-1] / h2

    b = F.reshape(-1)
    return A, b


class TestFFTvsSpsolveDirichlet:
    """Test 5-point FFT solvers match spsolve for Dirichlet BC."""

    @pytest.mark.parametrize("sigma", [0.0, 1.0, 10.0, 100.0])
    @pytest.mark.parametrize("n", [17, 33])
    def test_fa_dirichlet(self, n, sigma):
        """FA solver should match spsolve for Dirichlet BC."""
        k2 = sigma  # for test problem (modified Helmholtz)
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

        # FFT solver
        U_fa = fa_helmholtz(n, f_rhs, bc, sigma=sigma)

        # Sparse direct solver
        A, b = _build_5pt_dirichlet_system(n, f_rhs, bc, sigma=sigma)
        u_sparse = spsolve(A, b)
        U_sparse = np.zeros((n, n))
        U_sparse[1:-1, 1:-1] = u_sparse.reshape(n - 2, n - 2)
        # Fill boundary
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        U_sparse[0, :] = bc(x[0], y) if callable(bc) else 0.0
        U_sparse[-1, :] = bc(x[-1], y) if callable(bc) else 0.0
        U_sparse[:, 0] = bc(x, y[0]) if callable(bc) else 0.0
        U_sparse[:, -1] = bc(x, y[-1]) if callable(bc) else 0.0

        # Compare interior only (boundaries are exact for homogeneous BC)
        diff = np.max(np.abs(U_fa[1:-1, 1:-1] - U_sparse[1:-1, 1:-1]))
        assert diff < 1e-10, \
            f"n={n}, sigma={sigma}: ||U_fa - U_sparse||_inf = {diff:.2e}"

    @pytest.mark.parametrize("method", ['cr', 'facr'])
    def test_cr_facr_dirichlet(self, method):
        """CR and FACR should match FA (and thus spsolve)."""
        sigma = 10.0
        n = 33
        k2 = sigma
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

        U_fa = fa_helmholtz(n, f_rhs, bc, sigma=sigma)

        if method == 'cr':
            U_method = cr_helmholtz(n, f_rhs, bc, sigma=sigma)
        else:
            U_method = facr_helmholtz(n, f_rhs, bc, sigma=sigma)

        diff = np.max(np.abs(U_fa - U_method))
        assert diff < 1e-10, \
            f"method={method}: ||U_fa - U_method||_inf = {diff:.2e}"


class TestFFTvsSpsolveNeumann:
    """Test 5-point FFT solvers match spsolve for Neumann BC (modified Helmholtz)."""

    @pytest.mark.parametrize("sigma", [1.0, 10.0])
    @pytest.mark.parametrize("n", [9, 17])
    def test_fa_neumann_modified(self, n, sigma):
        """
        FA Neumann solver should match spsolve for modified Helmholtz.

        Note: The ghost-point Neumann matrix is asymmetric, but the
        symmetrized version S = D^{-1} G D gives the same solution
        after d_scale rescaling. We compare with the unsymmetric sparse solve.
        """
        k2 = sigma
        u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

        # FFT solver (uses symmetrized DCT-I)
        U_fa = fa_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='neumann')

        # Sparse direct solver (uses unsymmetric ghost-point matrix)
        from gmres_solver import build_helmholtz_matrix, _build_rhs
        A, info = build_helmholtz_matrix(n, sigma=sigma, bc_type='neumann')
        b = _build_rhs(n, f_rhs, bc, k2=sigma, bc_x='N', bc_y='N')
        u_sparse = spsolve(A, b)
        U_sparse = u_sparse.reshape(n, n)

        # Compare (both should solve the same continuous problem)
        diff = np.max(np.abs(U_fa - U_sparse))
        assert diff < 1e-8, \
            f"n={n}, sigma={sigma}: ||U_fa - U_sparse||_inf = {diff:.2e}"


class TestFFTvsSpsolveTrueHelmholtz:
    """Test FFT solver vs spsolve for true Helmholtz (away from resonance)."""

    def test_true_helmholtz_dirichlet(self):
        """True Helmholtz (sigma < 0) with Dirichlet BC away from resonance."""
        sigma = -5.0  # true Helmholtz, well away from eigenvalues
        n = 33
        k2 = -sigma  # for test problem construction

        # Use the standard Dirichlet test problem
        def u_exact(x, y):
            return np.sin(np.pi * x) * np.sin(np.pi * y)

        def f_rhs(x, y):
            return (2 * np.pi**2 + sigma) * np.sin(np.pi * x) * np.sin(np.pi * y)

        def bc(x, y):
            return 0.0

        # FFT solver
        U_fa = fa_helmholtz(n, f_rhs, bc, sigma=sigma)

        # Sparse direct solver
        A, b = _build_5pt_dirichlet_system(n, f_rhs, bc, sigma=sigma)
        u_sparse = spsolve(A, b)
        U_sparse_int = u_sparse.reshape(n - 2, n - 2)

        diff = np.max(np.abs(U_fa[1:-1, 1:-1] - U_sparse_int))
        assert diff < 1e-10, \
            f"sigma={sigma}: ||U_fa - U_sparse||_inf = {diff:.2e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
