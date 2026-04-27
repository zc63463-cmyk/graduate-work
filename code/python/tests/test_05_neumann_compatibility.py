#!/usr/bin/env python3
"""
Test 05: Neumann Compatibility Condition & Weighted Mean
========================================================

For the pure Neumann Poisson equation (sigma=0), the discrete Laplacian
has a zero eigenvalue (constant mode). The system is solvable only if
the compatibility condition is satisfied:

    integral_Omega f dOmega = 0

In the discrete setting with ghost-point DCT-I:
    h² * sum_{i,j} w_x[i] * w_y[j] * f_{i,j} = 0

where w = trapezoidal weights: [1/2, 1, 1, ..., 1, 1/2].

After solving, the solution should have zero weighted mean:
    h² * sum_{i,j} w_x[i] * w_y[j] * u_{i,j} ≈ 0

For modified Helmholtz (sigma > 0), the (0,0) eigenvalue is sigma > 0,
so no compatibility issue arises.

Tests:
    - Compatible RHS: solution has weighted mean ~0
    - Incompatible RHS: warning should be raised
    - Modified Helmholtz Neumann: no compatibility issue
    - u_{0,0} mode = 0 in solution for Poisson

Acceptance: weighted mean < 1e-12 (after projection)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
import warnings
from helmholtz_solver import (
    fa_helmholtz, cr_helmholtz, facr_helmholtz,
    check_neumann_compatibility, _neumann_d_scale,
    _helmholtz_test_problem_neumann,
)


def _trapezoid_weights(n):
    """Trapezoidal quadrature weights: [1/2, 1, ..., 1, 1/2]."""
    w = np.ones(n)
    w[0] = 0.5
    w[-1] = 0.5
    return w


def _weighted_mean(u, h, bc_x='N', bc_y='N'):
    """
    Compute discrete weighted mean of u.

    For Neumann BC: trapezoidal weights.
    For Dirichlet: uniform weights (no boundary).
    """
    n = u.shape[0]
    if bc_x == 'N' and bc_y == 'N':
        wx = _trapezoid_weights(n)
        wy = _trapezoid_weights(n)
        W = wx[:, None] * wy[None, :]
        return h * h * np.sum(W * u) / (1.0)  # just the weighted sum * h²
    return 0.0


class TestNeumannPoissonCompatibility:
    """Test Neumann Poisson compatibility condition."""

    @pytest.mark.parametrize("n", [9, 17, 33])
    def test_compatible_rhs_zero_mean(self, n):
        """Compatible RHS (cos mode) should give solution with ~zero mean."""
        k2 = 0.0  # Poisson
        u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

        U = fa_helmholtz(n, f_rhs, bc, sigma=0.0, bc_type='neumann')

        # Check weighted mean of solution
        h = 1.0 / (n - 1)
        wx = _trapezoid_weights(n)
        wy = _trapezoid_weights(n)
        W = wx[:, None] * wy[None, :]
        weighted_mean_u = h * h * np.sum(W * U)
        # The exact solution cos(pi*x)*cos(pi*y) has nonzero integral,
        # but the discrete Neumann solution is only unique up to a constant.
        # The FFT solver naturally sets the (0,0) mode to zero.
        # So the weighted mean should be approximately zero.
        assert abs(weighted_mean_u) < 1e-10, \
            f"n={n}: weighted mean of Neumann Poisson solution = {weighted_mean_u:.2e}"

    def test_incompatible_rhs_warning(self):
        """Incompatible RHS for Neumann Poisson should raise a warning."""
        n = 17

        # f = 1 (nonzero integral, incompatible with Neumann Poisson)
        def f_rhs(x, y):
            return np.ones_like(x)

        def bc(x, y):
            return 0.0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            U = fa_helmholtz(n, f_rhs, bc, sigma=0.0, bc_type='neumann')

            # Should have raised a compatibility warning
            compat_warnings = [x for x in w if "compatibility" in str(x.message).lower()
                               or "兼容" in str(x.message)]
            assert len(compat_warnings) > 0, \
                "Expected compatibility warning for incompatible Neumann Poisson RHS"

    @pytest.mark.parametrize("n", [9, 17])
    def test_check_neumann_compatibility_function(self, n):
        """Test the check_neumann_compatibility() utility function."""
        h = 1.0 / (n - 1)

        # Compatible: f = cos(pi*x)*cos(pi*y) with Neumann Poisson
        # This f has integral = 0 (by orthogonality with constant)
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        F_compat = np.cos(np.pi * X) * np.cos(np.pi * Y)

        is_compat, mean_f = check_neumann_compatibility(F_compat, h, 'N', 'N', 0.0)
        assert is_compat, f"Compatible RHS flagged as incompatible: mean_f = {mean_f:.2e}"

        # Incompatible: f = 1
        F_incompat = np.ones((n, n))
        is_compat2, mean_f2 = check_neumann_compatibility(F_incompat, h, 'N', 'N', 0.0)
        assert not is_compat2, f"Incompatible RHS not detected: mean_f = {mean_f2:.2e}"


class TestNeumannModifiedHelmholtz:
    """Test that modified Helmholtz (sigma>0) Neumann has no compatibility issue."""

    @pytest.mark.parametrize("sigma", [1.0, 10.0, 100.0])
    def test_modified_no_compatibility_issue(self, sigma):
        """Modified Helmholtz Neumann should work even with nonzero-mean RHS."""
        n = 17

        # f = 1 (nonzero integral, but sigma > 0 so no compatibility issue)
        def f_rhs(x, y):
            return np.ones_like(x)

        def bc(x, y):
            return 0.0

        # Should NOT raise compatibility warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            U = fa_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='neumann')

            compat_warnings = [x for x in w if "compatibility" in str(x.message).lower()]
            assert len(compat_warnings) == 0, \
                f"sigma={sigma}: unexpected compatibility warning for modified Helmholtz"

        # Solution should be valid (not NaN/Inf)
        assert np.all(np.isfinite(U)), f"sigma={sigma}: solution contains NaN/Inf"


class TestNeumannZeroMode:
    """Test that the (0,0) mode is properly handled."""

    @pytest.mark.parametrize("n", [9, 17, 33])
    def test_poisson_zero_mode_is_zero(self, n):
        """For Neumann Poisson, the (0,0) Fourier coefficient should be ~0."""
        k2 = 0.0
        u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

        U = fa_helmholtz(n, f_rhs, bc, sigma=0.0, bc_type='neumann')

        # The DCT-I coefficient for mode (0,0) corresponds to the
        # weighted mean under the d_scale symmetrization.
        # U = U_tilde * d_sx * d_sy, so U_tilde = U / d_sx / d_sy
        # The (0,0) mode being zero means the solver set it to zero
        # (since the denominator lam_0 + lam_0 + 0 = 0 for Poisson).
        d_sx = _neumann_d_scale(n)
        d_sy = _neumann_d_scale(n)
        U_tilde = U / d_sx[:, None] / d_sy[None, :]

        # Verify via 2D DCT-I with correct axis specification
        from scipy.fft import dct
        u_hat = dct(dct(U_tilde, type=1, norm='ortho', axis=0),
                     type=1, norm='ortho', axis=1)
        u_hat_00 = u_hat[0, 0]

        assert abs(u_hat_00) < 1e-10, \
            f"n={n}: (0,0) mode = {u_hat_00:.2e} (should be ~0 for Neumann Poisson)"


class TestNeumannConvergence:
    """Test Neumann solver convergence rate."""

    @pytest.mark.parametrize("sigma", [1.0, 10.0])
    def test_neumann_convergence_rate(self, sigma):
        """Modified Helmholtz Neumann should converge at ~2nd order."""
        k2 = sigma
        u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

        ns = [9, 17, 33, 65]
        errors = []
        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            U = fa_helmholtz(n, f_rhs, bc, sigma=sigma, bc_type='neumann')
            err = np.max(np.abs(U - Ue))
            errors.append(err)

        rates = []
        for i in range(1, len(ns)):
            rate = np.log2(errors[i - 1] / errors[i])
            rates.append(rate)

        mean_rate = np.mean(rates)
        assert 1.7 < mean_rate < 2.3, \
            f"sigma={sigma}: Neumann convergence rate = {mean_rate:.2f} (expected ~2.0)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
