#!/usr/bin/env python3
"""
Test 06: Modified vs True Helmholtz — Spectral Properties
==========================================================

This test validates the key distinction between modified and true Helmholtz:

Modified Helmholtz (-∇² + κ²)u = f  (σ = +κ² > 0):
    - Frequency domain denominator: λ_{p,q} + κ² > 0 always
    - Matrix is SPD (positive definite)
    - No resonance possible

True Helmholtz (-∇² - κ²)u = f  (σ = -κ² < 0):
    - Frequency domain denominator: λ_{p,q} - κ² can be zero
    - Matrix is indefinite when κ² > λ_min
    - Resonance occurs when κ² ≈ λ_{p,q} for some mode (p,q)

Tests:
    - Modified: all denominators strictly positive
    - True: min |denominator| can approach zero
    - Resonance detection: warning for true H., not for modified
    - Poisson + Neumann zero mode: not flagged as resonance

Acceptance: modified denominators all > 0; true near-resonance detected
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
import warnings
from helmholtz_solver import (
    check_resonance, _eigenvalues_dirichlet, _eigenvalues_neumann,
    fa_helmholtz, _helmholtz_test_problem_dirichlet,
)


class TestModifiedHelmholtzSpectrum:
    """Test spectral properties of modified Helmholtz (sigma > 0)."""

    @pytest.mark.parametrize("sigma", [1.0, 10.0, 100.0])
    def test_dirichlet_denominators_all_positive(self, sigma):
        """For modified Helmholtz, all λ_{p,q} + σ should be strictly positive."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_dirichlet(N, h)

        LXX, LYY = np.meshgrid(lam_x, lam_y, indexing='ij')
        denom = LXX + LYY + sigma

        assert np.all(denom > 0), \
            f"sigma={sigma}: some denominators are non-positive! min = {np.min(denom):.2e}"

    @pytest.mark.parametrize("sigma", [1.0, 10.0])
    def test_neumann_denominators_all_positive(self, sigma):
        """For modified Helmholtz Neumann, (0,0) mode gives κ² > 0."""
        n = 17
        N = n  # Neumann uses full grid
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_neumann(N, h)
        lam_y, _, _ = _eigenvalues_neumann(N, h)

        LXX, LYY = np.meshgrid(lam_x, lam_y, indexing='ij')
        denom = LXX + LYY + sigma

        # (0,0) mode: 0 + 0 + σ = σ > 0
        assert denom[0, 0] == sigma, \
            f"(0,0) denominator = {denom[0,0]:.2e}, expected σ = {sigma}"
        assert np.all(denom > 0), \
            f"sigma={sigma}: some Neumann denominators are non-positive!"

    def test_modified_no_resonance_warning(self):
        """Modified Helmholtz should never trigger resonance warning."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_dirichlet(N, h)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            is_res, min_denom = check_resonance(100.0, lam_x, lam_y, 'D', 'D')

            assert not is_res, "Modified Helmholtz flagged as resonant!"
            resonance_warnings = [x for x in w if "resonance" in str(x.message).lower()]
            assert len(resonance_warnings) == 0, \
                "Unexpected resonance warning for modified Helmholtz"


class TestTrueHelmholtzSpectrum:
    """Test spectral properties of true Helmholtz (sigma < 0)."""

    def test_true_helmholtz_can_be_indefinite(self):
        """For true Helmholtz with large κ², some denominators can be negative."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_neumann(N, h)

        # Use sigma = -200 (true Helmholtz with large κ²)
        sigma = -200.0
        LXX, LYY = np.meshgrid(lam_x, lam_y, indexing='ij')
        denom = LXX + LYY + sigma

        # Some denominators should be negative (matrix is indefinite)
        assert np.any(denom < 0), \
            f"sigma={sigma}: expected some negative denominators for true Helmholtz"

    def test_near_resonance_detection(self):
        """True Helmholtz near resonance should be detected."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_dirichlet(N, h)

        # Choose sigma ≈ -(λ_{1,1} + λ_{1,1}) = -2*λ_{1,1}
        lam_11 = lam_x[0] + lam_y[0]  # smallest 2D eigenvalue
        sigma_near_resonance = -lam_11 + 1e-6  # just above resonance

        is_res, min_denom = check_resonance(sigma_near_resonance, lam_x, lam_y, 'D', 'D')

        # Should detect near-resonance
        assert min_denom < 1.0, \
            f"Near resonance not detected: min_denom = {min_denom:.2e}"

    def test_resonance_warning_for_true_helmholtz(self):
        """True Helmholtz at resonance should trigger a warning."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_dirichlet(N, h)

        # Choose sigma exactly at eigenvalue
        lam_11 = lam_x[0] + lam_y[0]
        sigma_resonant = -lam_11

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            is_res, min_denom = check_resonance(sigma_resonant, lam_x, lam_y, 'D', 'D')

            assert is_res, "True Helmholtz resonance not detected"
            res_warnings = [x for x in w if "resonance" in str(x.message).lower()
                            or "Resonance" in str(x.message)]
            assert len(res_warnings) > 0, \
                "No resonance warning for true Helmholtz at resonance"

    def test_min_denom_approaches_zero(self):
        """As σ → -λ_{p,q}, the minimum denominator should approach zero."""
        n = 33
        N = n - 2
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_dirichlet(N, h)
        lam_y, _, _ = _eigenvalues_dirichlet(N, h)

        lam_11 = lam_x[0] + lam_y[0]

        # Scan sigma approaching resonance
        deltas = [1.0, 0.1, 0.01, 0.001]
        min_denoms = []
        for delta in deltas:
            sigma = -lam_11 + delta
            _, min_d = check_resonance(sigma, lam_x, lam_y, 'D', 'D')
            min_denoms.append(min_d)

        # min_denom should decrease as delta → 0
        for i in range(1, len(min_denoms)):
            assert min_denoms[i] < min_denoms[i - 1], \
                f"min_denom not decreasing: {min_denoms}"


class TestPoissonZeroModeNotResonance:
    """Test that Poisson + Neumann zero mode is NOT flagged as resonance."""

    def test_poisson_neumann_no_resonance_warning(self):
        """Poisson (σ=0) + Neumann should not trigger resonance warning for zero mode."""
        n = 17
        N = n
        h = 1.0 / (n - 1)

        lam_x, _, _ = _eigenvalues_neumann(N, h)
        lam_y, _, _ = _eigenvalues_neumann(N, h)

        # σ=0: (0,0) mode gives denom = 0, but this is expected
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            is_res, min_denom = check_resonance(0.0, lam_x, lam_y, 'N', 'N')

            # min_denom should be 0 (zero mode), but no resonance warning
            # because Poisson + Neumann is a known case
            res_warnings = [x for x in w if "resonance" in str(x.message).lower()]
            # The code should NOT warn for Poisson + Neumann zero mode
            assert len(res_warnings) == 0, \
                "Poisson + Neumann zero mode should not trigger resonance warning"


class TestSigmaParameterFramework:
    """Test the unified σ parameter framework."""

    @pytest.mark.parametrize("sigma_val,sigma_name", [
        (0.0, "Poisson"),
        (10.0, "modified_Helmholtz"),
        (-5.0, "true_Helmholtz"),
    ])
    def test_sigma_consistency(self, sigma_val, sigma_name):
        """Verify that sigma parameter is correctly used in solvers."""
        n = 33

        def u_exact(x, y):
            return np.sin(np.pi * x) * np.sin(np.pi * y)

        def f_rhs(x, y):
            return (2 * np.pi**2 + sigma_val) * np.sin(np.pi * x) * np.sin(np.pi * y)

        def bc(x, y):
            return 0.0

        # Away from resonance
        if sigma_val < 0 and abs(sigma_val) < 2 * np.pi**2:
            # This sigma is below the first eigenvalue, matrix is still SPD
            pass
        elif sigma_val < 0:
            pytest.skip("sigma at resonance, skip")

        U = fa_helmholtz(n, f_rhs, bc, sigma=sigma_val)

        # For Poisson and modified Helmholtz, check convergence
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)
        err = np.max(np.abs(U - Ue))

        # Should be O(h²) accuracy
        if sigma_val >= 0:
            assert err < 0.01, \
                f"{sigma_name} (sigma={sigma_val}): error = {err:.2e}"

    def test_k2_backward_compat(self):
        """k2 parameter should be equivalent to sigma=+k2."""
        n = 33
        k2_val = 10.0

        def u_exact(x, y):
            return np.sin(np.pi * x) * np.sin(np.pi * y)

        def f_rhs(x, y):
            return (2 * np.pi**2 + k2_val) * np.sin(np.pi * x) * np.sin(np.pi * y)

        def bc(x, y):
            return 0.0

        U_k2 = fa_helmholtz(n, f_rhs, bc, k2=k2_val)
        U_sigma = fa_helmholtz(n, f_rhs, bc, sigma=k2_val)

        diff = np.max(np.abs(U_k2 - U_sigma))
        assert diff < 1e-14, \
            f"k2 and sigma not equivalent: diff = {diff:.2e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
