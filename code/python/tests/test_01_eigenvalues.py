#!/usr/bin/env python3
"""
Test 01: Eigenvalue Verification
=================================

Verify that DST-I and DCT-I eigenvalues correctly diagonalize
the discrete Laplacian sub-matrices.

For Dirichlet BC (DST-I):
    B = tridiag(-1, 2, -1)  (1D sub-matrix, N=n-2 interior points)
    Eigenvalues: mu_p = 2 - 2*cos(p*pi/(N+1)), p=1,...,N
    Eigenvectors: v_p(k) = sqrt(2/(N+1)) * sin(p*k*pi/(N+1))
    Test: ||B*v_p - mu_p*v_p||_inf < tol

For Neumann BC (DCT-I, ghost-point):
    Symmetrized matrix S = D^{-1} G D, where G is the ghost-point Neumann matrix
    Eigenvalues: lam_k = (2/h^2)(1 - cos(k*pi/(n-1))), k=0,...,n-1
    Eigenvectors: C_{k,j} (DCT-I columns)
    Test: ||S*e_k - lam_k*e_k||_inf < tol

Acceptance: ||Bv - mu*v||_inf < 1e-12
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest
from scipy.fft import dst, idst, dct, idct


class TestDirichletEigenvalues:
    """Test DST-I eigenvalues against explicit B = tridiag(-1,2,-1)."""

    @pytest.mark.parametrize("N", [5, 10, 20, 50])
    def test_bv_equals_mu_v(self, N):
        """Verify B*v_p = mu_p * v_p for all modes p."""
        # Build B = tridiag(-1, 2, -1) of size N x N
        B = np.zeros((N, N))
        for i in range(N):
            B[i, i] = 2.0
            if i > 0:
                B[i, i - 1] = -1.0
            if i < N - 1:
                B[i, i + 1] = -1.0

        # Eigenvalues
        p = np.arange(1, N + 1)
        mu = 2.0 - 2.0 * np.cos(p * np.pi / (N + 1))

        # Eigenvectors from DST-I (type=1, unnormalized)
        # v_p(k) = sin(p * k * pi / (N+1)), k = 1,...,N
        for pi in range(N):
            v = np.sin((pi + 1) * p[:1] * np.pi / (N + 1))  # wrong
            # Correct: v_p(k) for fixed p, k=1..N
            k = np.arange(1, N + 1)
            v_p = np.sin((pi + 1) * k * np.pi / (N + 1))
            # Normalize
            v_p = v_p / np.linalg.norm(v_p)

            residual = B @ v_p - mu[pi] * v_p
            assert np.max(np.abs(residual)) < 1e-12, \
                f"N={N}, p={pi+1}: ||Bv-mu*v||_inf = {np.max(np.abs(residual)):.2e}"

    @pytest.mark.parametrize("N", [5, 10, 20])
    def test_dst_orthogonality(self, N):
        """Verify DST-I transform matrix is orthogonal: S^T S = I."""
        # Build DST-I matrix
        k = np.arange(1, N + 1)
        p = np.arange(1, N + 1)
        K, P = np.meshgrid(k, p, indexing='ij')
        S = np.sqrt(2.0 / (N + 1)) * np.sin(K * P * np.pi / (N + 1))

        # Check orthogonality
        STS = S.T @ S
        I = np.eye(N)
        assert np.max(np.abs(STS - I)) < 1e-14, \
            f"N={N}: ||S^T S - I||_inf = {np.max(np.abs(STS - I)):.2e}"

    @pytest.mark.parametrize("n", [9, 17, 33])
    def test_eigenvalue_formula_matches_fft(self, n):
        """Verify our eigenvalue formula matches scipy DST-I diagonalization."""
        N = n - 2
        h = 1.0 / (n - 1)

        # Build B matrix
        B = np.zeros((N, N))
        for i in range(N):
            B[i, i] = 2.0
            if i > 0:
                B[i, i - 1] = -1.0
            if i < N - 1:
                B[i, i + 1] = -1.0

        # Our formula
        p = np.arange(1, N + 1)
        mu_formula = 2.0 - 2.0 * np.cos(p * np.pi / (N + 1))

        # Eigenvalues from numpy (sorted)
        mu_numpy = np.sort(np.linalg.eigvalsh(B))

        # Compare
        np.testing.assert_allclose(mu_formula, mu_numpy, atol=1e-10,
                                   err_msg=f"n={n}: eigenvalue formula mismatch")


class TestNeumannEigenvalues:
    """Test DCT-I eigenvalues against the symmetrized Neumann matrix."""

    @pytest.mark.parametrize("n", [9, 17, 33])
    def test_symmetrized_matrix_eigenvalues(self, n):
        """Verify eigenvalues of symmetrized Neumann matrix match formula."""
        h = 1.0 / (n - 1)
        h2 = h * h

        # Build ghost-point Neumann matrix G (n x n, asymmetric)
        # Row i: (u_{i-1} - 2*u_i + u_{i+1}) / h^2  with ghost points
        # At i=0: u_{-1} = u_1 => (u_1 - 2*u_0 + u_1)/h^2 = (-2*u_0 + 2*u_1)/h^2
        # At i=n-1: u_n = u_{n-2} => (-2*u_{n-1} + 2*u_{n-2})/h^2
        G = np.zeros((n, n))
        for i in range(n):
            G[i, i] = 2.0 / h2
            if i == 0:
                G[i, 1] = -2.0 / h2  # ghost: u_{-1} = u_1 => coeff to u_1 is -2/h²
            else:
                G[i, i - 1] = -1.0 / h2
            if i == n - 1:
                G[i, n - 2] = -2.0 / h2  # ghost: u_n = u_{n-2}
            else:
                if i == 0:
                    pass  # already set above (G[0,1] = -2/h²)
                else:
                    G[i, i + 1] = -1.0 / h2

        # Symmetrize: S = D^{-1} G D, D = diag(sqrt(2),1,...,1,sqrt(2))
        d_scale = np.ones(n)
        d_scale[0] = np.sqrt(2.0)
        d_scale[-1] = np.sqrt(2.0)
        D = np.diag(d_scale)
        Dinv = np.diag(1.0 / d_scale)
        S = Dinv @ G @ D

        # Check symmetry
        assert np.max(np.abs(S - S.T)) < 1e-12, \
            f"n={n}: S is not symmetric, ||S-S^T|| = {np.max(np.abs(S - S.T)):.2e}"

        # Our eigenvalue formula
        k = np.arange(0, n)
        lam_formula = (2.0 / h2) * (1.0 - np.cos(k * np.pi / (n - 1)))

        # Eigenvalues from numpy (sorted)
        lam_numpy = np.sort(np.linalg.eigvalsh(S))

        # Compare
        np.testing.assert_allclose(np.sort(lam_formula), lam_numpy, atol=1e-8,
                                   err_msg=f"n={n}: Neumann eigenvalue formula mismatch")

    @pytest.mark.parametrize("n", [9, 17, 33])
    def test_dct1_orthogonality(self, n):
        """Verify DCT-I transform matrix is orthogonal: C^T C = I (with alpha weights)."""
        # Build DCT-I matrix with dual-endpoint weights alpha_k * alpha_j
        alpha = np.ones(n)
        alpha[0] = 1.0 / np.sqrt(2.0)
        alpha[-1] = 1.0 / np.sqrt(2.0)

        C = np.zeros((n, n))
        for k in range(n):
            for j in range(n):
                C[k, j] = alpha[k] * alpha[j] * np.sqrt(2.0 / (n - 1)) * \
                           np.cos(k * j * np.pi / (n - 1))

        # Check orthogonality
        CTC = C.T @ C
        I = np.eye(n)
        assert np.max(np.abs(CTC - I)) < 1e-13, \
            f"n={n}: ||C^T C - I||_inf = {np.max(np.abs(CTC - I)):.2e}"


class TestEigenvalueBounds:
    """Test eigenvalue bounds and condition number formula."""

    @pytest.mark.parametrize("n", [17, 33, 65])
    def test_dirichlet_lambda_min(self, n):
        """Verify lambda_min ~ 2*pi^2 for 2D Dirichlet BC (5-point Laplacian)."""
        N = n - 2
        h = 1.0 / (n - 1)
        p = np.arange(1, N + 1)
        mu = 2.0 - 2.0 * np.cos(p * np.pi / (N + 1))
        # 2D minimum eigenvalue: lambda_{1,1} = (mu_1 + mu_1) / h^2
        lam_min = 2 * mu[0] / h**2

        # Should approach 2*pi^2 as h->0
        expected = 2 * np.pi**2
        rel_err = abs(lam_min - expected) / expected
        assert rel_err < 0.05, \
            f"n={n}: lambda_min={lam_min:.4f}, expected~{expected:.4f}, rel_err={rel_err:.4f}"

    @pytest.mark.parametrize("n", [17, 33, 65])
    def test_dirichlet_condition_number(self, n):
        """Verify kappa(A) ~ 4/(pi^2 * h^2) for 5-point Dirichlet Poisson."""
        N = n - 2
        h = 1.0 / (n - 1)
        p = np.arange(1, N + 1)
        mu = 2.0 - 2.0 * np.cos(p * np.pi / (N + 1))
        # 2D eigenvalues: (mu_p + mu_q) / h^2
        lam_2d_min = 2 * mu[0] / h**2  # p=q=1
        lam_2d_max = 2 * mu[-1] / h**2  # p=q=N

        kappa_actual = lam_2d_max / lam_2d_min
        kappa_formula = 4.0 / (np.pi**2 * h**2)
        ratio = kappa_actual / kappa_formula
        # Should be close to 1 (converges as h->0)
        assert 0.8 < ratio < 1.2, \
            f"n={n}: kappa_actual/kappa_formula = {ratio:.4f}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
