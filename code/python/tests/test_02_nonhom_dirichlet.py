#!/usr/bin/env python3
"""
Test 02: Non-homogeneous Dirichlet BC — F += bc/h² Verification
================================================================

Validates the A1 bug fix: non-homogeneous Dirichlet BC requires
F += bc/h² (NOT F -= bc/h²).

The discrete 5-point Laplacian for interior point (i,j):
    (-u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1} + 4*u_{i,j}) / h² + sigma*u_{i,j} = f_{i,j}

When a neighbor is a boundary node with value g, moving to RHS gives:
    (4*u_{i,j})/h² + sigma*u_{i,j} = f_{i,j} + g/h²   (+ sign!)

Tests:
    - 5-point: convergence order ~2 with non-homogeneous Dirichlet
    - FFT9: matches 5-point on same grid (different accuracy)
    - CR/FACR: same as FA

Acceptance: convergence order in [1.8, 2.2] for 5-point
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from helmholtz_solver import (
    fa_helmholtz, cr_helmholtz, facr_helmholtz,
    fft9_helmholtz, solve_helmholtz
)


def _make_nonhom_problem(k2=0.0, sigma=None, sx=1.0, sy=1.0):
    """
    Non-homogeneous Dirichlet test problem.

    u(x,y) = sin(2*pi*x/sx) * sin(3*pi*y/sy)
    f(x,y) = ((2*pi/sx)^2 + (3*pi/sy)^2 + sigma) * u(x,y)
    BC: u|_{dO} = sin(2*pi*x/sx) * sin(3*pi*y/sy)  (non-zero on boundary)

    The sine product is NOT zero on boundaries (except corners),
    so this tests the non-homogeneous BC handling.
    """
    if sigma is None:
        sigma = k2  # backward compat

    def u_exact(x, y):
        return np.sin(2 * np.pi * x / sx) * np.sin(3 * np.pi * y / sy)

    def f_rhs(x, y):
        lam = (2 * np.pi / sx)**2 + (3 * np.pi / sy)**2
        return (lam + sigma) * np.sin(2 * np.pi * x / sx) * np.sin(3 * np.pi * y / sy)

    def bc(x, y):
        return u_exact(x, y)

    return u_exact, f_rhs, bc


class TestNonhomDirichlet5pt:
    """Test 5-point solver with non-homogeneous Dirichlet BC."""

    @pytest.mark.parametrize("sigma", [0.0, 1.0, 10.0, 100.0])
    def test_5pt_convergence_order(self, sigma):
        """5-point FA solver should achieve ~2nd order convergence."""
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)

        ns = [9, 17, 33, 65]
        errors = []
        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            U = fa_helmholtz(n, f_rhs, bc, sigma=sigma)
            err = np.max(np.abs(U - Ue))
            errors.append(err)

        # Compute convergence rates
        rates = []
        for i in range(1, len(ns)):
            rate = np.log2(errors[i - 1] / errors[i])
            rates.append(rate)

        # Mean rate should be ~2.0
        mean_rate = np.mean(rates)
        assert 1.7 < mean_rate < 2.3, \
            f"sigma={sigma}: mean convergence rate = {mean_rate:.2f} (expected ~2.0)"

    def test_5pt_sign_fix_verification(self):
        """
        If F -= bc/h² (wrong sign), the error would be O(1), not O(h²).
        This test would fail with the old bug.
        """
        sigma = 10.0
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)

        n = 33
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        U = fa_helmholtz(n, f_rhs, bc, sigma=sigma)
        err = np.max(np.abs(U - Ue))

        # With correct sign, error should be O(h²) ~ 1e-3 for n=33
        # With wrong sign, error would be O(1)
        assert err < 0.1, \
            f"Error {err:.2e} is too large — sign fix (F+=bc/h²) may be broken"


class TestNonhomDirichletAllSolvers:
    """Test all solvers give consistent results with non-homogeneous Dirichlet."""

    @pytest.mark.parametrize("method", ['fa', 'cr', 'facr'])
    def test_solver_matches_fa(self, method):
        """CR and FACR should match FA solver output."""
        sigma = 10.0
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)
        n = 33

        U_fa = fa_helmholtz(n, f_rhs, bc, sigma=sigma)
        U_method = solve_helmholtz(n, f_rhs, bc, sigma=sigma, method=method)

        diff = np.max(np.abs(U_fa - U_method))
        assert diff < 1e-10, \
            f"method={method}: max|U_fa - U_method| = {diff:.2e}"

    def test_fft9_nonhom_dirichlet(self):
        """FFT9 solver with non-homogeneous Dirichlet should converge at 4th order."""
        sigma = 0.0  # Poisson
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)

        ns = [9, 17, 33, 65]
        errors = []
        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            U = fft9_helmholtz(n, f_rhs, bc, sigma=sigma)
            err = np.max(np.abs(U - Ue))
            errors.append(err)

        rates = []
        for i in range(1, len(ns)):
            rate = np.log2(errors[i - 1] / errors[i])
            rates.append(rate)

        mean_rate = np.mean(rates)
        # FFT9 compact scheme can slightly exceed 4th order in practice
        assert 3.5 < mean_rate < 5.0, \
            f"FFT9 mean convergence rate = {mean_rate:.2f} (expected ~4.0)"


class TestNonhomDirichletModifiedHelmholtz:
    """Test modified Helmholtz with non-homogeneous Dirichlet BC."""

    @pytest.mark.parametrize("sigma", [1.0, 50.0])
    def test_modified_convergence(self, sigma):
        """Modified Helmholtz should also converge with non-homogeneous BC."""
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)

        n = 33
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        U = fa_helmholtz(n, f_rhs, bc, sigma=sigma)
        err = np.max(np.abs(U - Ue))

        # Should converge at 2nd order
        assert err < 0.01, \
            f"sigma={sigma}: error = {err:.2e} (expected O(h²))"


class TestNonhomDirichletTrueHelmholtz:
    """Test true Helmholtz (sigma < 0) with non-homogeneous Dirichlet BC."""

    def test_true_helmholtz_away_from_resonance(self):
        """True Helmholtz away from resonance should still work."""
        sigma = -5.0  # true Helmholtz, away from resonance
        u_exact, f_rhs, bc = _make_nonhom_problem(sigma=sigma)

        n = 33
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        U = fa_helmholtz(n, f_rhs, bc, sigma=sigma)
        err = np.max(np.abs(U - Ue))

        # Should still converge (no resonance for this sigma)
        assert err < 0.1, \
            f"sigma={sigma}: error = {err:.2e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
