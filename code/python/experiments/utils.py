#!/usr/bin/env python3
"""Shared utilities for thesis experiments."""

import numpy as np
import os, sys

# Ensure parent directory is on path for imports
_CODE_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
if _CODE_DIR not in sys.path:
    sys.path.insert(0, _CODE_DIR)

_RESULTS_DIR = os.path.join(os.path.dirname(__file__), 'results')
_FIGURES_DIR = os.path.join(os.path.dirname(__file__), 'figures')


def get_results_dir():
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    return _RESULTS_DIR

def get_figures_dir():
    os.makedirs(_FIGURES_DIR, exist_ok=True)
    return _FIGURES_DIR


# ============================================================================
# Test Problems
# ============================================================================

def test_problem_dirichlet(sigma, sx=1.0, sy=1.0):
    """Homogeneous Dirichlet: u = sin(πx/sx) sin(πy/sy).
    f = ((π/sx)² + (π/sy)² + σ) u
    """
    def u_exact(x, y):
        return np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)
    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + sigma) * \
               np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)
    def bc(x, y):
        return 0.0
    return u_exact, f_rhs, bc


def test_problem_dirichlet_mode(sigma, m=2, n=3, sx=1.0, sy=1.0):
    """Dirichlet with mode (m,n): u = sin(mπx/sx) sin(nπy/sy).
    f = ((mπ/sx)² + (nπ/sy)² + σ) u
    BC is zero when m,n are integers (homogeneous Dirichlet on [0,sx]×[0,sy]).
    """
    def u_exact(x, y):
        return np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy)
    def f_rhs(x, y):
        return ((m * np.pi / sx)**2 + (n * np.pi / sy)**2 + sigma) * \
               np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy)
    def bc(x, y):
        return u_exact(x, y)
    return u_exact, f_rhs, bc


def test_problem_nonhom_dirichlet(sigma, m=2, n=3, sx=1.0, sy=1.0):
    """Non-homogeneous Dirichlet: u = sin(mπx/sx) sin(nπy/sy) + 1.
    
    The constant shift +1 makes boundary values non-zero: g_D = 1 on ∂Ω.
    
    PDE: (-Δ + σ)u = f, where:
      f = ((mπ/sx)² + (nπ/sy)² + σ) sin(mπx/sx) sin(nπy/sy) + σ
    
    Note: The Laplacian of the constant 1 is zero, so f ≠ ((mπ/sx)² + (nπ/sy)² + σ) * u.
    Only the σ·1 term contributes to f from the constant part.
    
    Good for verifying non-homogeneous Dirichlet boundary correction (A1 bug fix).
    """
    def u_exact(x, y):
        return np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy) + 1.0
    def f_rhs(x, y):
        return ((m * np.pi / sx)**2 + (n * np.pi / sy)**2 + sigma) * \
               np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy) + sigma
    def bc(x, y):
        return 1.0  # u|∂Ω = sin(mπx/sx)sin(nπy/sy) + 1 = 0 + 1 = 1
    return u_exact, f_rhs, bc


def test_problem_neumann(sigma, sx=1.0, sy=1.0):
    """Homogeneous Neumann: u = cos(πx/sx) cos(πy/sy).
    f = ((π/sx)² + (π/sy)² + σ) u
    """
    def u_exact(x, y):
        return np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)
    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + sigma) * \
               np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)
    def bc(x, y):
        return 0.0  # du/dn = 0
    return u_exact, f_rhs, bc


def test_problem_mixed_nd(sigma, sx=1.0, sy=1.0):
    """Mixed BC (x=Neumann, y=Dirichlet): u = cos(πx/sx) sin(πy/sy).
    f = ((π/sx)² + (π/sy)² + σ) u
    """
    def u_exact(x, y):
        return np.cos(np.pi * x / sx) * np.sin(np.pi * y / sy)
    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + sigma) * \
               np.cos(np.pi * x / sx) * np.sin(np.pi * y / sy)
    def bc(x, y):
        return 0.0
    return u_exact, f_rhs, bc


def test_problem_mixed_dn(sigma, sx=1.0, sy=1.0):
    """Mixed BC (x=Dirichlet, y=Neumann): u = sin(πx/sx) cos(πy/sy).
    f = ((π/sx)² + (π/sy)² + σ) u
    """
    def u_exact(x, y):
        return np.sin(np.pi * x / sx) * np.cos(np.pi * y / sy)
    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + sigma) * \
               np.sin(np.pi * x / sx) * np.cos(np.pi * y / sy)
    def bc(x, y):
        return 0.0
    return u_exact, f_rhs, bc


def test_problem_polynomial(sigma, sx=1.0, sy=1.0):
    """Non-eigenfunction Dirichlet: u = x²(sx-x)² y²(sy-y)².
    Good for GMRES convergence testing.
    """
    def u_exact(x, y):
        px = x**2 * (sx - x)**2
        py = y**2 * (sy - y)**2
        return px * py
    def f_rhs(x, y):
        px = x**2 * (sx - x)**2
        py = y**2 * (sy - y)**2
        d2px = 2*(sx - x)**2 - 8*x*(sx - x) + 2*x**2
        d2py = 2*(sy - y)**2 - 8*y*(sy - y) + 2*y**2
        return -(d2px * py + px * d2py) + sigma * px * py
    def bc(x, y):
        return 0.0
    return u_exact, f_rhs, bc


# ============================================================================
# Computation Helpers
# ============================================================================

def compute_convergence_rate(errors):
    """Compute log2 convergence rates from a list of errors."""
    rates = [np.nan]
    for i in range(1, len(errors)):
        if errors[i] > 0 and errors[i-1] > 0:
            rates.append(np.log2(errors[i-1] / errors[i]))
        else:
            rates.append(np.nan)
    return rates


def equation_type(sigma):
    """Return equation type label from sigma."""
    if sigma == 0:
        return 'Poisson'
    elif sigma > 0:
        return 'modified_H'
    else:
        return 'true_H'


def compute_max_error(n, u_solver, u_exact_func, sx=1.0, sy=1.0):
    """Compute max-norm error between solver output and exact solution."""
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    Ue = u_exact_func(X, Y)
    return np.max(np.abs(u_solver - Ue))
