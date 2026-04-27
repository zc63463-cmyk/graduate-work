#!/usr/bin/env python3
"""
Helmholtz Equation Fast Solvers
================================

Unified solver suite for the generalized Helmholtz equation:
    (-nabla^2 + sigma) u = f   on [0, sx] x [0, sy]

where sigma is the Helmholtz parameter:
    sigma > 0 : modified Helmholtz (Yukawa / screened Poisson)
                (-nabla^2 + k^2) u = f,  sigma = +k^2
    sigma < 0 : true Helmholtz (wave / resonance)
                (-nabla^2 - kappa^2) u = f,  sigma = -kappa^2
    sigma = 0 : Poisson equation  -nabla^2 u = f

Supports Dirichlet, Neumann, and mixed boundary conditions.

Methods:
    1. FA  (Fourier Analysis)    -- 2D DST/DCT, O(N^2 log N)
    2. CR  (Cyclic Reduction)    -- 1D transform + Thomas, O(N^2 log N)
    3. FACR(l)                   -- CR + FA hybrid
    4. FFT9 (compact 4th-order)  -- 9-point stencil + DST/DCT
    5. 5-point (2nd-order)       -- baseline for comparison

When sigma=0, all solvers reduce to the Poisson equation solvers.

Grid conventions:
    - Dirichlet BC: unknowns at INTERIOR nodes (n-2 points), DST-I
    - Neumann BC:   unknowns at ALL nodes (n points) using ghost-point
      approach with symmetrization. The ghost-point Neumann Laplacian
      matrix G is asymmetric; we symmetrize via D = diag(sqrt(2),1,...,1,sqrt(2)),
      giving S = D^{-1} G D which IS diagonalized by DCT-I.

      Eigenvalues: lam_k = (2/h^2)(1 - cos(pi*k/(n-1))), k=0,...,n-1
      Symmetrization scaling: d_scale = [sqrt(2), 1, ..., 1, sqrt(2)]
      RHS scaling:  F_tilde = h^2 * F / d_scale_x / d_scale_y
      Solution:     U = U_tilde * d_scale_x * d_scale_y

Key references:
    - Strang (1999), "The Discrete Cosine Transform", SIAM Review
    - Boisvert (1987), Algorithm 651: HFFT, ACM TOMS
    - Swarztrauber (1977), FACR algorithm, SIAM Review
"""

import numpy as np
from scipy.fft import dst, idst, dct, idct
import time
import warnings


# ============================================================================
# PDE Parameters
# ============================================================================

def _resolve_sigma(sigma=None, k2=None):
    """
    Resolve the sigma parameter from either sigma or k2.

    The unified equation is: (-nabla^2 + sigma) u = f

    Parameters
    ----------
    sigma : float, optional
        Direct sigma parameter. If provided, takes precedence.
    k2 : float, optional
        Backward-compatible k^2 parameter (implies sigma = +k^2,
        i.e., modified Helmholtz).

    Returns
    -------
    sigma : float
        The resolved sigma value.

    Raises
    ------
    ValueError
        If neither sigma nor k2 is provided.
    """
    if sigma is not None:
        return sigma
    elif k2 is not None:
        return k2  # modified Helmholtz: sigma = +k^2
    else:
        return 0.0  # default: Poisson


# ============================================================================
# Transform & Eigenvalue Utilities
# ============================================================================

def _get_transform(bc_type, axis='x'):
    """
    Return forward/inverse transform functions for given BC type.

    For Dirichlet: DST type=1 on N=n-2 interior points
    For Neumann:   DCT type=1 on N=n full-grid points (ghost-point approach)
    """
    if bc_type == 'D':
        def fwd(x, **kw): return dst(x, type=1, norm='ortho', **kw)
        def inv(x, **kw): return idst(x, type=1, norm='ortho', **kw)
        return fwd, inv, 'DST-1'
    elif bc_type == 'N':
        def fwd(x, **kw): return dct(x, type=1, norm='ortho', **kw)
        def inv(x, **kw): return idct(x, type=1, norm='ortho', **kw)
        return fwd, inv, 'DCT-1'
    else:
        raise ValueError(f"Unknown BC type '{bc_type}' for {axis}-axis. Use 'D' or 'N'.")


def _eigenvalues_dirichlet(N, h):
    """
    Eigenvalues for Dirichlet BC (DST-1 modes).

    lam_p = (2/h^2)(1 - cos(p*pi/(N+1))),  p = 1, ..., N

    For 5-point stencil diagonal: d_p = 4 - 2cos(p*pi/(N+1))
    """
    p = np.arange(1, N + 1)
    cos_p = np.cos(p * np.pi / (N + 1))
    lam = (2.0 / h**2) * (1.0 - cos_p)
    d_eig = 4.0 - 2.0 * cos_p  # for CR/FACR diagonal (h^2 * lam)
    return lam, d_eig, cos_p


def _eigenvalues_neumann(N, h):
    """
    Eigenvalues for Neumann BC (DCT-1 modes, ghost-point approach).

    N = n (number of grid nodes including boundaries).

    lam_k = (2/h^2)(1 - cos(k*pi/(N-1))),  k = 0, 1, ..., N-1

    Note: k=0 gives lam=0 (constant mode). For Helmholtz with k^2>0,
    the effective eigenvalue is k^2>0, so no singularity.

    For CR/FACR: the per-mode diagonal of the h^2-scaled symmetrized
    Neumann Laplacian is d_eig = 2(1-cos(k*pi/(N-1))).
    The full y-tridiagonal diagonal per mode p is: 2 + d_eig[p] + h^2*k^2
    = 4 - 2*cos(p*pi/(N-1)) + h^2*k^2

    Parameters:
        N: number of grid nodes (= n)
        h: grid spacing
    """
    k = np.arange(0, N)  # k = 0, 1, ..., N-1
    cos_k = np.cos(k * np.pi / (N - 1))
    lam = (2.0 / h**2) * (1.0 - cos_k)
    d_eig = 2.0 - 2.0 * cos_k  # h^2 * lam
    return lam, d_eig, cos_k


def _get_eigenvalues(N, h, bc_type):
    """Get eigenvalues based on boundary condition type."""
    if bc_type == 'D':
        return _eigenvalues_dirichlet(N, h)
    elif bc_type == 'N':
        return _eigenvalues_neumann(N, h)
    else:
        raise ValueError(f"Unknown bc_type '{bc_type}'")


def _get_grid_dims(n, bc_x, bc_y):
    """
    Get transform grid dimensions based on BC types.

    Dirichlet: N = n-2 (interior nodes)
    Neumann:   N = n   (full grid including boundaries, ghost-point approach)

    Returns (Nx, Ny).
    """
    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2
    return Nx, Ny


def _neumann_d_scale(n):
    """
    Return the symmetrization scaling vector for Neumann BC.

    d_scale = [sqrt(2), 1, 1, ..., 1, sqrt(2)]  (length n)

    The ghost-point Neumann matrix G is asymmetric. It can be symmetrized
    via D = diag(d_scale), giving S = D^{-1} G D which is diagonalized
    by DCT-I.
    """
    d = np.ones(n)
    d[0] = np.sqrt(2.0)
    d[-1] = np.sqrt(2.0)
    return d


def _normalize_bc_type(bc_type):
    """
    Normalize bc_type to a tuple (bc_x, bc_y).

    Parameters:
        bc_type: 'dirichlet', 'neumann', or tuple like ('D', 'N', 'D', 'N')

    Returns:
        (bc_x, bc_y) where bc_x is 'D' or 'N', bc_y is 'D' or 'N'
    """
    if isinstance(bc_type, str):
        bc_type_lower = bc_type.lower()
        if bc_type_lower in ('dirichlet', 'd'):
            return ('D', 'D')
        elif bc_type_lower in ('neumann', 'n'):
            return ('N', 'N')
        else:
            raise ValueError(f"Unknown bc_type string '{bc_type}'")
    elif isinstance(bc_type, (tuple, list)):
        if len(bc_type) == 4:
            bc_x_left, bc_x_right, bc_y_bottom, bc_y_top = bc_type
            if bc_x_left != bc_x_right:
                raise NotImplementedError(
                    f"Mixed x-BC ({bc_x_left}/{bc_x_right}) not supported.")
            if bc_y_bottom != bc_y_top:
                raise NotImplementedError(
                    f"Mixed y-BC ({bc_y_bottom}/{bc_y_top}) not supported.")
            return (bc_x_left, bc_y_bottom)
        elif len(bc_type) == 2:
            return tuple(bc_type)
        else:
            raise ValueError(f"bc_type tuple must have 2 or 4 elements, got {len(bc_type)}")
    else:
        raise TypeError(f"bc_type must be str or tuple, got {type(bc_type)}")


def check_resonance(sigma, lam_x, lam_y, bc_x='D', bc_y='D', tol=1e-10):
    """
    Check if sigma causes resonance (denominator near zero).

    For modified Helmholtz (sigma > 0): lam + sigma > 0 always, no resonance.
    For true Helmholtz (sigma < 0): lam + sigma can be zero when
        sigma = -(lam_x + lam_y) for some mode (p,q).
    For Poisson (sigma = 0): Neumann has zero mode (expected, not resonance).
    """
    LXX, LYY = np.meshgrid(lam_x, lam_y, indexing='ij')
    ltot = LXX + LYY + sigma
    min_denom = np.min(np.abs(ltot))
    is_resonant = min_denom < tol

    if is_resonant:
        # Distinguish: Poisson + Neumann zero mode (expected) vs true resonance
        if sigma == 0.0 and (bc_x == 'N' or bc_y == 'N'):
            # This is the expected zero mode for Neumann Poisson, not a resonance
            pass  # Don't warn for expected zero mode
        elif sigma < 0:
            warnings.warn(
                f"Resonance detected (true Helmholtz, sigma={sigma:.6e})! "
                f"Min |denom| = {min_denom:.2e}. Solution may be inaccurate.",
                RuntimeWarning)
        else:
            warnings.warn(
                f"Near-singular system (sigma={sigma:.6e})! "
                f"Min |denom| = {min_denom:.2e}. Check problem setup.",
                RuntimeWarning)
    return is_resonant, min_denom


def check_neumann_compatibility(F_adj, h, bc_x, bc_y, sigma, tol=1e-8):
    """
    Check Neumann Poisson compatibility condition: integral of f must be zero.

    For pure Neumann BC with sigma=0, the discrete Laplacian has a zero
    eigenvalue (constant mode). The system is solvable only if the
    discrete compatibility condition is satisfied: sum(f * d_scale) ≈ 0.

    For sigma != 0, the (0,0) eigenvalue is sigma, so no compatibility
    issue arises (the system is non-singular as long as sigma != 0).

    Parameters
    ----------
    F_adj : ndarray
        RHS array (already scaled by d_scale^{-1} if Neumann)
    h : float
        Grid spacing
    bc_x, bc_y : str
        Boundary condition types ('D' or 'N')
    sigma : float
        Helmholtz parameter
    tol : float
        Tolerance for compatibility check

    Returns
    -------
    is_compatible : bool
        True if compatibility condition is satisfied
    mean_f : float
        Discrete mean of f (weighted by d_scale and h^2)
    """
    if sigma != 0.0:
        return True, 0.0  # Non-zero sigma: no compatibility issue

    has_neumann = (bc_x == 'N' or bc_y == 'N')
    if not has_neumann:
        return True, 0.0  # No Neumann BC: no compatibility issue

    # For Neumann BC, F_adj is already divided by d_scale.
    # The compatibility condition is: h^2 * sum(F_adj * d_scale_x * d_scale_y) ≈ 0
    # Since F_adj = F / d_scale, we need: h^2 * sum(F) ≈ 0
    # But F_adj = F / d_scale, so sum(F) = sum(F_adj * d_scale)
    # The discrete integral is: h^2 * sum(F) = h^2 * sum(F_adj * d_scale)

    if bc_x == 'N' and bc_y == 'N':
        d_sx = _neumann_d_scale(F_adj.shape[0])
        d_sy = _neumann_d_scale(F_adj.shape[1])
        weighted_sum = np.sum(F_adj * d_sx[:, None] * d_sy[None, :])
    elif bc_x == 'N':
        d_sx = _neumann_d_scale(F_adj.shape[0])
        weighted_sum = np.sum(F_adj * d_sx[:, None])
    elif bc_y == 'N':
        d_sy = _neumann_d_scale(F_adj.shape[1])
        weighted_sum = np.sum(F_adj * d_sy[None, :])
    else:
        weighted_sum = np.sum(F_adj)

    discrete_integral = h * h * weighted_sum
    is_compatible = abs(discrete_integral) < tol * max(1.0, np.max(np.abs(F_adj)))

    if not is_compatible:
        warnings.warn(
            f"Neumann Poisson compatibility condition violated! "
            f"Discrete integral of f = {discrete_integral:.6e} (should be ~0). "
            f"Solution will have arbitrary constant offset. "
            f"Consider projecting out the (0,0) mode or using sigma > 0.",
            RuntimeWarning)

    return is_compatible, discrete_integral


# ============================================================================
# Grid Construction & RHS Building
# ============================================================================

def _make_bc_array(bc_func, coord_arr1, coord_arr2, n):
    """Evaluate BC function and ensure it returns a length-n array."""
    val = bc_func(coord_arr1, coord_arr2)
    arr = np.atleast_1d(np.asarray(np.squeeze(val), dtype=float))
    if arr.shape[0] == 1:
        return np.full(n, float(arr[0]))
    return arr


def _build_rhs_and_grid(n, f_func, bc_func, bc_x, bc_y, sx=1.0, sy=1.0):
    """
    Build RHS array on the appropriate grid for given BC types.

    For Dirichlet: f evaluated at interior nodes, with boundary terms subtracted.
    For Neumann:   f evaluated at ALL nodes (including boundaries), scaled by
                   d_scale^{-1} for symmetrization. No h^2 factor included --
                   it is added in CR/FACR solvers. Ghost-point approach makes
                   Neumann BC implicit in the DCT-I basis.

    Returns:
        F_adj:  RHS array of shape (Nx, Ny), already scaled for transforms
        bc_values: dict with boundary values (for Dirichlet) or None (for Neumann)
        grids: (x_node, y_node) grid arrays
        d_info: dict with d_scale_x, d_scale_y (for Neumann rescaling after solve)
    """
    h = sx / (n - 1)
    x_node = np.linspace(0, sx, n)
    y_node = np.linspace(0, sy, n)

    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2

    d_info = {}

    # --- Build RHS for each direction ---
    if bc_x == 'N' and bc_y == 'N':
        # Pure Neumann: evaluate f on full n x n node grid
        X, Y = np.meshgrid(x_node, y_node, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        d_scale_x = _neumann_d_scale(n)
        d_scale_y = _neumann_d_scale(n)
        d_info['d_scale_x'] = d_scale_x
        d_info['d_scale_y'] = d_scale_y

        # Symmetrization scaling: F_tilde = F / d_scale_x / d_scale_y
        # NOTE: h^2 is NOT included here; it's added in CR/FACR Thomas solve step.
        # FA solver uses original-scale eigenvalues (lam + k^2) so no h^2 needed.
        F_adj = F / d_scale_x[:, None] / d_scale_y[None, :]

        bc_values = {'bc_l': None, 'bc_r': None, 'bc_b': None, 'bc_t': None}

    elif bc_x == 'D' and bc_y == 'D':
        # Pure Dirichlet: evaluate f at interior nodes, subtract boundary contributions
        x_int = x_node[1:-1]
        y_int = y_node[1:-1]
        X, Y = np.meshgrid(x_int, y_int, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        bc_values = {}

        # Left/right boundary (x-direction)
        bc_l_full = _make_bc_array(bc_func, x_node[0], y_node, n)
        bc_r_full = _make_bc_array(bc_func, x_node[-1], y_node, n)
        F[0, :] += bc_l_full[1:-1] / h**2
        F[-1, :] += bc_r_full[1:-1] / h**2
        bc_values['bc_l'] = bc_l_full
        bc_values['bc_r'] = bc_r_full

        # Bottom/top boundary (y-direction)
        bc_b_full = _make_bc_array(bc_func, x_node, y_node[0], n)
        bc_t_full = _make_bc_array(bc_func, x_node, y_node[-1], n)
        F[:, 0] += bc_b_full[1:-1] / h**2
        F[:, -1] += bc_t_full[1:-1] / h**2
        bc_values['bc_b'] = bc_b_full
        bc_values['bc_t'] = bc_t_full

        F_adj = F  # Dirichlet: no d_scale needed

    elif bc_x == 'N' and bc_y == 'D':
        # Mixed: x=Neumann (full n nodes), y=Dirichlet (n-2 interior)
        X, Y = np.meshgrid(x_node, y_node[1:-1], indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        d_scale_x = _neumann_d_scale(n)
        d_info['d_scale_x'] = d_scale_x

        # Y-boundary corrections (Dirichlet)
        bc_b_full = _make_bc_array(bc_func, x_node, y_node[0], n)
        bc_t_full = _make_bc_array(bc_func, x_node, y_node[-1], n)
        # Standard Dirichlet corrections (bc_b/bc_t evaluated on n-point grid,
        # F has n rows in x and n-2 cols in y)
        # bc_b_full has length n (matches x-dimension); use all n entries
        F[:, 0] += bc_b_full / h**2
        F[:, -1] += bc_t_full / h**2

        bc_values = {'bc_l': None, 'bc_r': None,
                     'bc_b': bc_b_full, 'bc_t': bc_t_full}

        # Scale by d_scale_x^{-1} for x-Neumann symmetrization (no h^2)
        F_adj = F / d_scale_x[:, None]

    elif bc_x == 'D' and bc_y == 'N':
        # Mixed: x=Dirichlet (n-2 interior), y=Neumann (full n nodes)
        X, Y = np.meshgrid(x_node[1:-1], y_node, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        d_scale_y = _neumann_d_scale(n)
        d_info['d_scale_y'] = d_scale_y

        # X-boundary corrections (Dirichlet)
        bc_l_full = _make_bc_array(bc_func, x_node[0], y_node, n)
        bc_r_full = _make_bc_array(bc_func, x_node[-1], y_node, n)
        # Standard Dirichlet corrections
        F[0, :] += bc_l_full / h**2
        F[-1, :] += bc_r_full / h**2

        bc_values = {'bc_l': bc_l_full, 'bc_r': bc_r_full,
                     'bc_b': None, 'bc_t': None}

        # Scale by d_scale_y^{-1} for y-Neumann symmetrization (no h^2)
        F_adj = F / d_scale_y[None, :]

    else:
        raise ValueError(f"Unexpected BC combination: ({bc_x}, {bc_y})")

    return F_adj, bc_values, (x_node, y_node), d_info


# ============================================================================
# Thomas Algorithm
# ============================================================================

def thomas(a, d, b_vec):
    """
    Thomas algorithm for: a*u_{j-1} + d*u_j + a*u_{j+1} = b_j
    Constant off-diagonal coefficients.
    """
    M = len(b_vec)
    if M == 0: return np.array([])
    if M == 1: return np.array([b_vec[0]/d]) if abs(d)>1e-30 else np.zeros(1)

    cp, dp = np.zeros(M), np.zeros(M)
    cp[0] = a/d if abs(d)>1e-30 else 0.0
    dp[0] = b_vec[0]/d if abs(d)>1e-30 else 0.0

    for j in range(1, M):
        denom = d - a*cp[j-1]
        if abs(denom) < 1e-30: denom = 1e-30*(1 if denom>=0 else -1)
        cp[j] = a/denom if j < M-1 else 0.0
        dp[j] = (b_vec[j] - a*dp[j-1])/denom

    u = np.zeros(M)
    u[M-1] = dp[M-1]
    for j in range(M-2, -1, -1):
        u[j] = dp[j] - cp[j]*u[j+1]
    return u


def thomas_sym(d_vec, e_vec, b_vec):
    """
    Thomas algorithm for symmetric tridiagonal:
      e[j-1]*u_{j-1} + d[j]*u_j + e[j]*u_{j+1} = b[j]
    """
    M = len(b_vec)
    if M == 0: return np.array([])
    if M == 1:
        return np.array([b_vec[0]/d_vec[0]]) if abs(d_vec[0])>1e-30 else np.zeros(1)

    cp = np.zeros(M)
    dp = np.zeros(M)

    cp[0] = e_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0
    dp[0] = b_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0

    for j in range(1, M):
        c_j = e_vec[j-1]
        denom = d_vec[j] - c_j*cp[j-1]
        if abs(denom) < 1e-30: denom = 1e-30*(1 if denom>=0 else -1)
        cp[j] = (e_vec[j]/denom) if j < M-1 else 0.0
        dp[j] = (b_vec[j] - c_j*dp[j-1])/denom

    u = np.zeros(M)
    u[M-1] = dp[M-1]
    for j in range(M-2, -1, -1):
        u[j] = dp[j] - cp[j]*u[j+1]
    return u


# ============================================================================
# 2D Transform Helpers
# ============================================================================

def _transform_2d(F, fwd_x, fwd_y):
    """Apply 2D separable transform: x-direction then y-direction.
    
    fwd_x is applied to columns (x-lines), fwd_y to rows (y-lines).
    For indexing='ij': F[i,:] = y-line at x_i, F[:,j] = x-line at y_j.
    """
    Nx, Ny = F.shape
    Fh = np.zeros((Nx, Ny))
    # First: transform along x (columns) using fwd_x
    for j in range(Ny):
        Fh[:, j] = fwd_x(F[:, j])
    # Then: transform along y (rows) using fwd_y
    for i in range(Nx):
        Fh[i, :] = fwd_y(Fh[i, :])
    return Fh


def _itransform_2d(Fh, inv_x, inv_y):
    """Apply 2D separable inverse transform: y-direction then x-direction."""
    Nx, Ny = Fh.shape
    Ui = np.zeros((Nx, Ny))
    # First: inverse transform along y (rows) using inv_y
    for i in range(Nx):
        Ui[i, :] = inv_y(Fh[i, :])
    # Then: inverse transform along x (columns) using inv_x
    for j in range(Ny):
        Ui[:, j] = inv_x(Ui[:, j])
    return Ui


def _transform_1d_x(F, fwd_x):
    """Apply 1D transform along x-direction (each column)."""
    Nx, Ny = F.shape
    Fh = np.zeros((Nx, Ny))
    for j in range(Ny):
        Fh[:, j] = fwd_x(F[:, j])
    return Fh


def _itransform_1d_x(Fh, inv_x):
    """Apply 1D inverse transform along x-direction."""
    Nx, Ny = Fh.shape
    Ui = np.zeros((Nx, Ny))
    for j in range(Ny):
        Ui[:, j] = inv_x(Fh[:, j])
    return Ui


# ============================================================================
# Neumann off-diagonal vector (for CR/FACR y-tridiagonal)
# ============================================================================

def _neumann_e_vec(N):
    """
    Return the off-diagonal vector for the symmetrized Neumann
    y-tridiagonal system.

    The symmetrized Neumann matrix S_y has off-diagonal:
      e[0] = e[-1] = -sqrt(2)   (boundary rows)
      e[1:-1] = -1              (interior rows)

    Parameters:
        N: size of the system (= n, number of grid nodes)

    Returns:
        e_vec: array of length N-1
    """
    e_vec = np.full(N - 1, -1.0)
    if N > 1:
        e_vec[0] = -np.sqrt(2.0)
        e_vec[-1] = -np.sqrt(2.0)
    return e_vec


# ============================================================================
# Method 1: FA (Fourier Analysis) -- 2D Transform Direct Solve
# ============================================================================

def fa_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0,
                  sigma=None):
    """
    2D-DST/DCT direct Helmholtz solver. O(N^2 log N).

    Solves: (-nabla^2 + sigma) u = f

    Use sigma for the unified framework, or k2 for backward compatibility
    (k2 implies sigma = +k^2, i.e., modified Helmholtz).

    Dirichlet: DST-I on interior nodes (n-2 points)
    Neumann:   DCT-I on full grid (n points) with ghost-point symmetrization
    """
    sigma = _resolve_sigma(sigma, k2)
    bc_x, bc_y = _normalize_bc_type(bc_type)
    Nx, Ny = _get_grid_dims(n, bc_x, bc_y)
    h = sx / (n - 1)

    F_adj, bc_values, grids, d_info = _build_rhs_and_grid(
        n, f_func, bc_func, bc_x, bc_y, sx, sy)

    # Get transforms
    fwd_x, inv_x, _ = _get_transform(bc_x, 'x')
    fwd_y, inv_y, _ = _get_transform(bc_y, 'y')

    # Get eigenvalues
    lam_x, d_x, _ = _get_eigenvalues(Nx, h, bc_x)
    lam_y, d_y, _ = _get_eigenvalues(Ny, h, bc_y)

    # Check resonance
    check_resonance(sigma, lam_x, lam_y, bc_x, bc_y)

    # Check Neumann Poisson compatibility condition
    if sigma == 0.0 and (bc_x == 'N' or bc_y == 'N'):
        check_neumann_compatibility(F_adj, h, bc_x, bc_y, sigma)

    # 2D transform
    Fh = _transform_2d(F_adj, fwd_x, fwd_y)

    # Eigenvalue matrix + sigma
    LXX, LYY = np.meshgrid(lam_x, lam_y, indexing='ij')
    ltot = LXX + LYY + sigma

    # Solve in Fourier space
    Uh = np.zeros((Nx, Ny))
    mask = np.abs(ltot) > 1e-14
    Uh[mask] = Fh[mask] / ltot[mask]

    # Inverse transform
    Ui = _itransform_2d(Uh, inv_x, inv_y)

    # Assemble full n x n grid
    U = _assemble_solution(Ui, n, bc_values, bc_x, bc_y, grids, d_info)
    return U


# ============================================================================
# Method 2: CR (Cyclic Reduction) -- 1D Transform + Thomas
# ============================================================================

def cr_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0,
                  sigma=None):
    """
    CR solver: 1D transform in x + Thomas in y. O(N^2 log N).

    Solves: (-nabla^2 + sigma) u = f

    For each Fourier mode p in x, solve a y-direction tridiagonal system:
      Dirichlet: tridiag(-1, d_p+h^2*sigma, -1)
      Neumann:   tridiag with e_vec=[-sqrt2, -1,..., -1, -sqrt2], d=2+d_eig[p]+h^2*sigma
    """
    sigma = _resolve_sigma(sigma, k2)
    bc_x, bc_y = _normalize_bc_type(bc_type)
    Nx, Ny = _get_grid_dims(n, bc_x, bc_y)
    h = sx / (n - 1)

    F_adj, bc_values, grids, d_info = _build_rhs_and_grid(
        n, f_func, bc_func, bc_x, bc_y, sx, sy)

    # Get transforms for x-direction
    fwd_x, inv_x, _ = _get_transform(bc_x, 'x')

    # Get eigenvalues
    _, d_x, _ = _get_eigenvalues(Nx, h, bc_x)

    # Check resonance
    lam_x, _, _ = _get_eigenvalues(Nx, h, bc_x)
    lam_y, _, _ = _get_eigenvalues(Ny, h, bc_y)
    check_resonance(sigma, lam_x, lam_y, bc_x, bc_y)

    # 1D transform in x
    Fh = _transform_1d_x(F_adj, fwd_x)

    h2 = h * h

    # Per-mode diagonal depends on bc_x (which determines whether d_eig includes +2)
    # The y-tridiagonal diagonal is: 2 + h^2*lam_x[p] + h^2*sigma
    # For Dirichlet x: d_eig = 2 + h^2*lam_x, so diagonal = d_eig + h^2*sigma
    # For Neumann x:   d_eig = h^2*lam_x, so diagonal = 2 + d_eig + h^2*sigma
    if bc_x == 'D':
        dp = d_x + h2 * sigma  # d_eig_x already includes +2
    else:
        dp = 2.0 + d_x + h2 * sigma  # need to add +2 for y-diagonal

    # Tridiagonal structure depends on bc_y
    if bc_y == 'D':
        # Dirichlet: constant off-diag = -1
        Uh = np.zeros((Nx, Ny))
        for p in range(Nx):
            Uh[p, :] = thomas(-1.0, dp[p], h2 * Fh[p, :])
    else:
        # Neumann: e_vec = [-sqrt2, -1, ..., -sqrt2]
        e_vec = _neumann_e_vec(Ny)  # (Ny-1,)
        Uh = np.zeros((Nx, Ny))
        for p in range(Nx):
            d_vec = np.full(Ny, dp[p])
            Uh[p, :] = thomas_sym(d_vec, e_vec, h2 * Fh[p, :])

    # Inverse transform in x
    Ui = _itransform_1d_x(Uh, inv_x)

    # Assemble
    U = _assemble_solution(Ui, n, bc_values, bc_x, bc_y, grids, d_info)
    return U


# ============================================================================
# Method 3: FACR(l) -- Vectorized CR + FA
# ============================================================================

def facr_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0, l=None,
                    sigma=None):
    """
    FACR(l) hybrid algorithm for Helmholtz equation.

    Solves: (-nabla^2 + sigma) u = f

    Algorithm:
      1. 1D transform in x -> N independent tridiagonal systems
      2. l steps of CR on each y-line system (vectorized across modes)
      3. Solve small reduced systems
      4. Back-substitute to recover all y-line values
      5. Inverse transform in x
    """
    sigma = _resolve_sigma(sigma, k2)
    bc_x, bc_y = _normalize_bc_type(bc_type)
    Nx, Ny = _get_grid_dims(n, bc_x, bc_y)
    h = sx / (n - 1)

    if Ny <= 2 or Nx <= 2:
        return fa_helmholtz(n, f_func, bc_func, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)

    F_adj, bc_values, grids, d_info = _build_rhs_and_grid(
        n, f_func, bc_func, bc_x, bc_y, sx, sy)

    h2 = h * h

    if l is None:
        l = max(0, int(np.round(np.log2(max(1, np.log2(Ny))))))
        l = min(l, int(np.log2(Ny)) - 1)
    if l <= 0:
        return fa_helmholtz(n, f_func, bc_func, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)

    # Get transforms for x-direction
    fwd_x, inv_x, _ = _get_transform(bc_x, 'x')

    # Get eigenvalues for resonance check
    lam_x, _, _ = _get_eigenvalues(Nx, h, bc_x)
    lam_y, _, _ = _get_eigenvalues(Ny, h, bc_y)
    check_resonance(sigma, lam_x, lam_y, bc_x, bc_y)

    # Phase 1: 1D transform in x
    Fh = _transform_1d_x(F_adj, fwd_x)

    # Phase 2: Vectorized CR forward reduction
    # Build the y-direction tridiagonal system per mode
    _, d_eig_x, _ = _get_eigenvalues(Nx, h, bc_x)

    # Per-mode diagonal depends on bc_x (same logic as CR solver)
    if bc_x == 'D':
        dp = d_eig_x + h2 * sigma  # d_eig_x already includes +2
    else:
        dp = 2.0 + d_eig_x + h2 * sigma  # need to add +2 for y-diagonal

    # Tridiagonal structure depends on bc_y
    if bc_y == 'D':
        # Dirichlet: constant off-diag -1
        cur_d = np.tile(dp.reshape(-1, 1), (1, Ny))  # (Nx, Ny)
        cur_e = np.full((Nx, Ny - 1), -1.0)
    else:
        # Neumann: variable off-diag [-sqrt2, -1, ..., -sqrt2]
        cur_d = np.tile(dp.reshape(-1, 1), (1, Ny))  # (Nx, Ny)
        e_vec_base = _neumann_e_vec(Ny)  # (Ny-1,)
        cur_e = np.tile(e_vec_base.reshape(1, -1), (Nx, 1))  # (Nx, Ny-1)

    cur_b = h2 * Fh  # (Nx, Ny) -- h^2 scaled for Thomas/CR

    levels = []

    for step in range(l):
        cur_M = cur_d.shape[1]
        if cur_M <= 2:
            break

        levels.append((cur_d.copy(), cur_e.copy(), cur_b.copy()))

        M_new = (cur_M + 1) // 2
        even_idx = np.arange(0, cur_M, 2)

        d_new = cur_d[:, even_idx].copy()
        b_new = cur_b[:, even_idx].copy()

        # Left odd neighbor
        has_left = even_idx > 0
        left_odd = even_idx[has_left] - 1
        d_new[:, has_left] -= cur_e[:, left_odd]**2 / cur_d[:, left_odd]
        b_new[:, has_left] -= (cur_e[:, left_odd] / cur_d[:, left_odd]) * cur_b[:, left_odd]

        # Right odd neighbor
        has_right = even_idx < cur_M - 1
        right_odd = even_idx[has_right] + 1
        e_at_even = cur_e[:, even_idx[has_right]]
        d_new[:, has_right] -= e_at_even**2 / cur_d[:, right_odd]
        b_new[:, has_right] -= (e_at_even / cur_d[:, right_odd]) * cur_b[:, right_odd]

        # Reduced off-diagonal
        if M_new > 1:
            odd_between = np.arange(1, cur_M, 2)[:M_new - 1]
            e_new = -cur_e[:, odd_between - 1] * cur_e[:, odd_between] / cur_d[:, odd_between]
        else:
            e_new = np.empty((Nx, 0))

        cur_d = d_new
        cur_e = e_new
        cur_b = b_new

    # Phase 3: Solve reduced systems
    M_red = cur_d.shape[1]
    u_cur = np.zeros((Nx, M_red))

    for p in range(Nx):
        if M_red == 1:
            u_cur[p, 0] = cur_b[p, 0] / cur_d[p, 0] if abs(cur_d[p, 0]) > 1e-30 else 0.0
        elif cur_e.shape[1] > 0:
            u_cur[p, :] = thomas_sym(cur_d[p], cur_e[p], cur_b[p])
        else:
            u_cur[p, :M_red] = cur_b[p, :M_red] / cur_d[p, :M_red]

    # Phase 4: Vectorized back-substitution
    for rev in reversed(range(len(levels))):
        d_lv, e_lv, b_lv = levels[rev]
        M_lv = d_lv.shape[1]
        M_even = u_cur.shape[1]

        u_full = np.zeros((Nx, M_lv))
        u_full[:, :2 * M_even:2] = u_cur

        # Interior odd rows
        n_interior_odd = min(M_even - 1, M_lv // 2)
        if n_interior_odd > 0:
            odd_int = np.arange(1, 2 * n_interior_odd + 1, 2)
            u_full[:, odd_int] = (b_lv[:, odd_int]
                - e_lv[:, odd_int - 1] * u_full[:, odd_int - 1]
                - e_lv[:, odd_int] * u_full[:, odd_int + 1]) / d_lv[:, odd_int]

        # Last odd row
        if M_lv % 2 == 0:
            j_last = M_lv - 1
            u_full[:, j_last] = (b_lv[:, j_last] - e_lv[:, j_last - 1] * u_full[:, j_last - 1]) / d_lv[:, j_last]

        u_cur = u_full

    Uh = u_cur

    # Phase 5: Inverse transform in x
    Ui = _itransform_1d_x(Uh, inv_x)

    # Assemble
    U = _assemble_solution(Ui, n, bc_values, bc_x, bc_y, grids, d_info)
    return U


# ============================================================================
# Method 4: FFT9 -- 4th-Order Compact Finite Difference
# ============================================================================

def apply_Rh_full(G):
    """
    Apply R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0] to full grid G.
    Returns R_h*g for interior points only.
    """
    n = G.shape[0]
    N = n - 2
    Rg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ii, jj = i + 1, j + 1
            Rg[i, j] = (2.0 / 3.0) * G[ii, jj] \
                     + (1.0 / 12.0) * (G[ii - 1, jj] + G[ii + 1, jj]
                                      + G[ii, jj - 1] + G[ii, jj + 1])
    return Rg


def compute_bc_correction_9pt(n, bc_func, x, y, h, bc_x='D', bc_y='D', sigma=0.0):
    """
    Compute BC correction for FFT9 solver: (-L_h + sigma R_h)u = R_h f.

    For non-homogeneous Dirichlet BC, boundary values must be moved from
    LHS to RHS. Two contributions:
      1. L_h boundary: stencil [1,4,1;4,-20,4;1,4,1]/(6h²) → subtracted (negative)
      2. -sigma R_h boundary: stencil [0,1/12,0;1/12,2/3,1/12;0,1/12,0] → sign depends on sigma

    Corner diagonal neighbors (coeff 1/(6h²) in L_h) are double-counted by
    the x and y boundary loops; a correction is added back at the end.
    R_h has no diagonal stencil entries, so no corner double-counting there.
    """
    N = n - 2
    h2 = h * h
    coeff_L = 1.0 / (6.0 * h2)  # L_h stencil denominator
    coeff_R = 1.0 / 12.0          # R_h off-diagonal stencil coefficient
    bc_corr = np.zeros((N, N))

    bc_l = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float).copy()
    bc_r = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float).copy()
    bc_b = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float).copy()
    bc_t = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float).copy()

    if bc_x == 'D':
        for j in range(N):
            jj = j + 1
            # L_h contribution: edge + diagonal boundary neighbors
            val_L = 4.0 * bc_l[jj]
            if jj - 1 >= 0: val_L += bc_l[jj - 1]
            if jj + 1 < n: val_L += bc_l[jj + 1]
            bc_corr[0, j] -= coeff_L * val_L
            # sigma R_h contribution: edge boundary neighbor only (1/12 coeff)
            bc_corr[0, j] += sigma * coeff_R * bc_l[jj]

            val_L = 4.0 * bc_r[jj]
            if jj - 1 >= 0: val_L += bc_r[jj - 1]
            if jj + 1 < n: val_L += bc_r[jj + 1]
            bc_corr[N - 1, j] -= coeff_L * val_L
            bc_corr[N - 1, j] += sigma * coeff_R * bc_r[jj]

    if bc_y == 'D':
        for i in range(N):
            ii = i + 1
            # L_h contribution
            val_L = 4.0 * bc_b[ii]
            if ii - 1 >= 0: val_L += bc_b[ii - 1]
            if ii + 1 < n: val_L += bc_b[ii + 1]
            bc_corr[i, 0] -= coeff_L * val_L
            # sigma R_h contribution
            bc_corr[i, 0] += sigma * coeff_R * bc_b[ii]

            val_L = 4.0 * bc_t[ii]
            if ii - 1 >= 0: val_L += bc_t[ii - 1]
            if ii + 1 < n: val_L += bc_t[ii + 1]
            bc_corr[i, N - 1] -= coeff_L * val_L
            bc_corr[i, N - 1] += sigma * coeff_R * bc_t[ii]

    # Corner fix: diagonal boundary neighbors (coeff 1/(6h²) in L_h)
    # were double-counted by both x and y loops. Add back one copy.
    if bc_x == 'D' and bc_y == 'D':
        # Bottom-left: u(x[0], y[0]) affects interior (0,0)
        bc_corr[0, 0] += coeff_L * bc_l[0]
        # Top-left: u(x[0], y[-1]) affects interior (0,N-1)
        bc_corr[0, N - 1] += coeff_L * bc_l[n - 1]
        # Bottom-right: u(x[-1], y[0]) affects interior (N-1,0)
        bc_corr[N - 1, 0] += coeff_L * bc_r[0]
        # Top-right: u(x[-1], y[-1]) affects interior (N-1,N-1)
        bc_corr[N - 1, N - 1] += coeff_L * bc_r[n - 1]

    return bc_corr, bc_l, bc_r, bc_b, bc_t


def fft9_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet',
                    sx=1.0, sy=1.0, method='4th', sigma=None):
    """
    FFT9 solver for Helmholtz equation with 4th-order compact finite difference.

    Solves: (-L_h + sigma R_h) * u = R_h * f
    Equivalently: (L_h - sigma R_h) * u = -R_h * f

    Use sigma for the unified framework, or k2 for backward compatibility.

    Currently only supports Dirichlet BC (FFT9 + Neumann requires
    a different stencil treatment at boundaries).

    Parameters:
        method: '4th' for 4th-order compact, 'spectral' for exact eigenvalues
    """
    sigma = _resolve_sigma(sigma, k2)
    bc_x, bc_y = _normalize_bc_type(bc_type)

    if bc_x == 'N' or bc_y == 'N':
        # FFT9 with Neumann BC requires special treatment of the compact
        # stencil at boundaries. Fall back to FA solver for now.
        return fa_helmholtz(n, f_func, bc_func, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)

    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = sx / (n - 1)
    N = n - 2

    F = np.asarray(f_func(X, Y), dtype=float)

    # Apply R_h to g = -f
    G = -F
    Rg = apply_Rh_full(G)

    # BC correction for L_h and sigma R_h
    bc_corr, bc_l, bc_r, bc_b, bc_t = compute_bc_correction_9pt(
        n, bc_func, x, y, h, bc_x, bc_y, sigma)
    Rg += bc_corr

    # Get transforms
    fwd_x, inv_x, _ = _get_transform(bc_x, 'x')
    fwd_y, inv_y, _ = _get_transform(bc_y, 'y')

    # 2D transform
    Rg_hat = _transform_2d(Rg, fwd_x, fwd_y)

    if method == '4th':
        k_vals = np.arange(1, N + 1)
        cos_kx = np.cos(k_vals * np.pi / (N + 1))
        cos_ky = np.cos(k_vals * np.pi / (N + 1))
        CK, CL = np.meshgrid(cos_kx, cos_ky, indexing='ij')

        lam_L = (1.0 / (6.0 * h**2)) * (-20.0 + 8.0 * CK + 8.0 * CL + 4.0 * CK * CL)
        lam_R4 = 2.0 / 3.0 + (1.0 / 6.0) * (CK + CL)
        # Denom: (L_h - sigma R_h) in Fourier space → û = -ĝ / (λ̂_L - σ λ̂_R)
        denom = lam_L - sigma * lam_R4

        U_hat = np.zeros((N, N))
        mask = np.abs(denom) > 1e-14
        U_hat[mask] = Rg_hat[mask] / denom[mask]

    elif method == 'spectral':
        F_int = F[1:-1, 1:-1].copy()
        F_adj = F_int.copy()
        F_adj[0, :]  += bc_l[1:-1] / h**2
        F_adj[-1, :] += bc_r[1:-1] / h**2
        F_adj[:, 0]  += bc_b[1:-1] / h**2
        F_adj[:, -1] += bc_t[1:-1] / h**2
        F_hat = _transform_2d(F_adj, fwd_x, fwd_y)

        k_vals = np.arange(1, N + 1)
        P, Q = np.meshgrid(k_vals, k_vals, indexing='ij')
        lam_x_exact = (P * np.pi / sx)**2
        lam_y_exact = (Q * np.pi / sy)**2
        lam_exact = lam_x_exact + lam_y_exact + sigma

        U_hat = np.zeros((N, N))
        mask = np.abs(lam_exact) > 1e-14
        U_hat[mask] = F_hat[mask] / lam_exact[mask]

    # Inverse 2D transform
    U_int = _itransform_2d(U_hat, inv_x, inv_y)

    # Assemble (Dirichlet only)
    bc_values = {'bc_l': bc_l, 'bc_r': bc_r, 'bc_b': bc_b, 'bc_t': bc_t}
    grids = (x, y)
    d_info = {}
    U = _assemble_solution(U_int, n, bc_values, bc_x, bc_y, grids, d_info)
    return U


def fft9_oer_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0,
                        sigma=None):
    """
    FFT9 OER (Odd-Even Reduction) solver for Helmholtz equation.

    Solves: (-L_h + sigma R_h) * u = R_h * f

    Currently only supports Dirichlet BC.
    """
    sigma = _resolve_sigma(sigma, k2)
    bc_x, bc_y = _normalize_bc_type(bc_type)

    if bc_x == 'N' or bc_y == 'N':
        return fa_helmholtz(n, f_func, bc_func, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)

    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = sx / (n - 1)
    N = n - 2

    F = np.asarray(f_func(X, Y), dtype=float)

    # Apply R_h to g = -f
    G = -F
    Rg = apply_Rh_full(G)

    # BC correction
    bc_corr, bc_l, bc_r, bc_b, bc_t = compute_bc_correction_9pt(
        n, bc_func, x, y, h, bc_x, bc_y, sigma)
    Rg += bc_corr

    # 1D transform in x
    fwd_x, inv_x, _ = _get_transform(bc_x, 'x')
    Rg_hat = _transform_1d_x(Rg, fwd_x)

    # Eigenvalues of D and A_sub (9-point Laplacian blocks)
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    lam_D = (1.0 / (6.0 * h**2)) * (-20.0 + 8.0 * cos_k)
    lam_A = (1.0 / (6.0 * h**2)) * (4.0 + 2.0 * cos_k)

    # R_h block eigenvalues for Helmholtz correction
    lam_DR = 2.0 / 3.0 + cos_k / 6.0
    lam_AR = 1.0 / 12.0

    # Helmholtz: (L_h - sigma R_h) per mode
    lam_A_helm = lam_A - sigma * lam_AR
    lam_D_helm = lam_D - sigma * lam_DR

    # Thomas for each mode
    u_hat = np.zeros((N, N))
    for p in range(N):
        a = lam_A_helm[p]
        d = lam_D_helm[p]
        b = Rg_hat[p, :].copy()
        u_hat[p, :] = thomas(a, d, b)

    # Inverse transform in x
    U_int = _itransform_1d_x(u_hat, inv_x)

    # Assemble (Dirichlet only)
    bc_values = {'bc_l': bc_l, 'bc_r': bc_r, 'bc_b': bc_b, 'bc_t': bc_t}
    grids = (x, y)
    d_info = {}
    U = _assemble_solution(U_int, n, bc_values, bc_x, bc_y, grids, d_info)
    return U


# ============================================================================
# Method 5: 5-Point Baseline (2nd order)
# ============================================================================

def point5_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0,
                      sigma=None):
    """5-point FD + FFT for equation (-nabla^2 + sigma)u = f (2nd order)."""
    sigma = _resolve_sigma(sigma, k2)
    return fa_helmholtz(n, f_func, bc_func, sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)


# ============================================================================
# Solution Assembly Helper
# ============================================================================

def _assemble_solution(U_solve, n, bc_values, bc_x, bc_y, grids, d_info):
    """
    Assemble full n x n grid solution from solver output.

    For Dirichlet: U_solve has shape (n-2, n-2), add boundary values
    For Neumann:   U_solve has shape (n, n) on full grid with d_scale rescaling
    For Mixed:     appropriate combination
    """
    x_node, y_node = grids
    U = np.zeros((n, n))

    if bc_x == 'D' and bc_y == 'D':
        # Dirichlet: U_solve is (n-2, n-2) at interior nodes
        U[0, :] = bc_values['bc_l']
        U[-1, :] = bc_values['bc_r']
        U[:, 0] = bc_values['bc_b']
        U[:, -1] = bc_values['bc_t']
        U[1:-1, 1:-1] = U_solve

    elif bc_x == 'N' and bc_y == 'N':
        # Neumann: U_solve is already (n, n) on full grid
        # Apply d_scale rescaling: U = U_tilde * d_scale_x * d_scale_y
        d_scale_x = d_info['d_scale_x']
        d_scale_y = d_info['d_scale_y']
        U = U_solve * d_scale_x[:, None] * d_scale_y[None, :]

    elif bc_x == 'N' and bc_y == 'D':
        # x: Neumann (n nodes), y: Dirichlet (n-2 interior)
        d_scale_x = d_info['d_scale_x']
        U_temp = U_solve * d_scale_x[:, None]  # rescale x
        U[:, 0] = bc_values['bc_b']
        U[:, -1] = bc_values['bc_t']
        U[:, 1:-1] = U_temp

    elif bc_x == 'D' and bc_y == 'N':
        # x: Dirichlet (n-2 interior), y: Neumann (n nodes)
        d_scale_y = d_info['d_scale_y']
        U_temp = U_solve * d_scale_y[None, :]  # rescale y
        U[0, :] = bc_values['bc_l']
        U[-1, :] = bc_values['bc_r']
        U[1:-1, :] = U_temp

    return U


# ============================================================================
# Unified Interface
# ============================================================================

def solve_helmholtz(n, f_func, bc_func, k2=None, method='fa',
                     bc_type='dirichlet', sx=1.0, sy=1.0, sigma=None, **kwargs):
    """
    Unified interface for Helmholtz solvers.

    Parameters:
        n: grid size (n x n)
        f_func: RHS function f(x,y)
        bc_func: BC function (Dirichlet: u|dO, Neumann: du/dn|dO)
        k2: wavenumber squared (backward compat, implies sigma=+k^2)
        sigma: unified Helmholtz parameter in (-nabla^2 + sigma)u = f
               sigma > 0: modified Helmholtz, sigma < 0: true Helmholtz
        method: 'fa', 'cr', 'facr', 'fft9', 'fft9_oer', 'fft9_spectral', '5point'
        bc_type: 'dirichlet', 'neumann', or tuple ('D'/'N', 'D'/'N')
        sx, sy: domain size
        **kwargs: additional solver-specific parameters (e.g., l for FACR)

    Returns:
        U: solution on (n x n) grid
    """
    sigma = _resolve_sigma(sigma, k2)

    common_kw = dict(sigma=sigma, bc_type=bc_type, sx=sx, sy=sy)

    if method == 'fa':
        return fa_helmholtz(n, f_func, bc_func, **common_kw)
    elif method == 'cr':
        return cr_helmholtz(n, f_func, bc_func, **common_kw)
    elif method == 'facr':
        return facr_helmholtz(n, f_func, bc_func, l=kwargs.get('l', None), **common_kw)
    elif method == 'fft9':
        return fft9_helmholtz(n, f_func, bc_func, method='4th', **common_kw)
    elif method == 'fft9_oer':
        return fft9_oer_helmholtz(n, f_func, bc_func, **common_kw)
    elif method == 'fft9_spectral':
        return fft9_helmholtz(n, f_func, bc_func, method='spectral', **common_kw)
    elif method == '5point':
        return point5_helmholtz(n, f_func, bc_func, **common_kw)
    else:
        raise ValueError(f"Unknown method '{method}'. Available: ['fa','cr','facr','fft9','fft9_oer','fft9_spectral','5point']")


# ============================================================================
# Test Suite
# ============================================================================

def _helmholtz_test_problem_dirichlet(k2, sx=1.0, sy=1.0):
    """
    Test problem with known analytical solution for Dirichlet BC.

    u(x,y) = sin(pi*x/sx) * sin(pi*y/sy)
    f(x,y) = ((pi/sx)^2 + (pi/sy)^2 + k^2) * sin(pi*x/sx) * sin(pi*y/sy)
    BC: u = 0 on all boundaries (homogeneous Dirichlet)
    """
    def u_exact(x, y):
        return np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)

    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + k2) * \
               np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)

    def bc(x, y):
        return 0.0

    return u_exact, f_rhs, bc


def _helmholtz_test_problem_neumann(k2, sx=1.0, sy=1.0):
    """
    Test problem with known analytical solution for Neumann BC.

    u(x,y) = cos(pi*x/sx) * cos(pi*y/sy)
    f(x,y) = ((pi/sx)^2 + (pi/sy)^2 + k^2) * cos(pi*x/sx) * cos(pi*y/sy)
    BC: du/dn = 0 on all boundaries (homogeneous Neumann)

    With the ghost-point approach, f is evaluated at ALL grid nodes
    (including boundaries), and the solution is computed on the full grid.
    """
    def u_exact(x, y):
        return np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)

    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + k2) * \
               np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)

    def bc(x, y):
        return 0.0  # homogeneous Neumann: du/dn = 0

    return u_exact, f_rhs, bc


def _helmholtz_test_problem_general(k2, sx=1.0, sy=1.0):
    """
    General test problem with non-zero BC for Dirichlet.

    u(x,y) = sin(2*pi*x/sx) * sin(3*pi*y/sy)
    f(x,y) = ((2*pi/sx)^2 + (3*pi/sy)^2 + k^2) * sin(2*pi*x/sx) * sin(3*pi*y/sy)
    """
    def u_exact(x, y):
        return np.sin(2 * np.pi * x / sx) * np.sin(3 * np.pi * y / sy)

    def f_rhs(x, y):
        return ((2 * np.pi / sx)**2 + (3 * np.pi / sy)**2 + k2) * \
               np.sin(2 * np.pi * x / sx) * np.sin(3 * np.pi * y / sy)

    def bc(x, y):
        return u_exact(x, y)  # non-homogeneous Dirichlet

    return u_exact, f_rhs, bc


def test_correctness():
    """Test correctness of all Helmholtz solvers against analytical solutions."""
    print("=" * 80)
    print("HELMHOLTZ SOLVER -- CORRECTNESS TEST (Dirichlet)")
    print("=" * 80)

    k2_values = [0.0, 1.0, 10.0, 100.0]

    for k2 in k2_values:
        print(f"\n--- k^2 = {k2} (k = {np.sqrt(k2):.4f}) ---")
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

        ns = [9, 17, 33, 65]
        print(f"\n  {'n':>5s}|{'FA':>12s}|{'CR':>12s}|{'FACR':>12s}|{'FFT9':>12s}|{'5pt':>12s}")
        print("  " + "-" * 72)

        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U_fa = fa_helmholtz(n, f_rhs, bc, k2=k2)
                U_cr = cr_helmholtz(n, f_rhs, bc, k2=k2)
                U_fc = facr_helmholtz(n, f_rhs, bc, k2=k2)
                U_f9 = fft9_helmholtz(n, f_rhs, bc, k2=k2)
                U_5p = point5_helmholtz(n, f_rhs, bc, k2=k2)

                e_fa = np.max(np.abs(U_fa - Ue))
                e_cr = np.max(np.abs(U_cr - Ue))
                e_fc = np.max(np.abs(U_fc - Ue))
                e_f9 = np.max(np.abs(U_f9 - Ue))
                e_5p = np.max(np.abs(U_5p - Ue))

                print(f"  {n:5d}|{e_fa:12.2e}|{e_cr:12.2e}|{e_fc:12.2e}|{e_f9:12.2e}|{e_5p:12.2e}")
            except Exception as ex:
                print(f"  {n:5d}| ERROR: {ex}")


def test_convergence():
    """Test convergence rates of all solvers (Dirichlet)."""
    print("\n" + "=" * 80)
    print("CONVERGENCE RATE TEST (Dirichlet)")
    print("=" * 80)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

    ns = [9, 17, 33, 65, 129]
    prev_errors = {}

    print(f"\n  k^2 = {k2}")
    print(f"\n  {'n':>5s}|{'FA':>12s}|{'rate':>6s}|{'CR':>12s}|{'rate':>6s}|{'FFT9':>12s}|{'rate':>6s}")
    print("  " + "-" * 72)

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        e_fa = np.max(np.abs(fa_helmholtz(n, f_rhs, bc, k2=k2) - Ue))
        e_cr = np.max(np.abs(cr_helmholtz(n, f_rhs, bc, k2=k2) - Ue))
        e_f9 = np.max(np.abs(fft9_helmholtz(n, f_rhs, bc, k2=k2) - Ue))

        rates = {}
        for name, err in [('FA', e_fa), ('CR', e_cr), ('FFT9', e_f9)]:
            if name in prev_errors:
                rates[name] = np.log2(prev_errors[name] / err)
            else:
                rates[name] = np.nan
            prev_errors[name] = err

        r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
        r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
        r_f9 = f"{rates['FFT9']:6.2f}" if not np.isnan(rates['FFT9']) else "   --"

        print(f"  {n:5d}|{e_fa:12.2e}|{r_fa}|{e_cr:12.2e}|{r_cr}|{e_f9:12.2e}|{r_f9}")


def test_neumann_bc():
    """Test Neumann boundary condition solvers (ghost-point + symmetrization)."""
    print("\n" + "=" * 80)
    print("NEUMANN BC TEST (Ghost-Point + Symmetrized DCT-I)")
    print("=" * 80)

    k2_values = [1.0, 10.0]

    for k2 in k2_values:
        print(f"\n--- k^2 = {k2} ---")
        u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

        ns = [9, 17, 33, 65, 129]
        prev_errors = {}

        print(f"\n  {'n':>5s}|{'FA':>12s}|{'rate':>6s}|{'CR':>12s}|{'rate':>6s}|{'FACR':>12s}|{'rate':>6s}")
        print("  " + "-" * 72)

        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U_fa = fa_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
                U_cr = cr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
                U_fc = facr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')

                e_fa = np.max(np.abs(U_fa - Ue))
                e_cr = np.max(np.abs(U_cr - Ue))
                e_fc = np.max(np.abs(U_fc - Ue))

                rates = {}
                for name, err in [('FA', e_fa), ('CR', e_cr), ('FACR', e_fc)]:
                    if name in prev_errors:
                        rates[name] = np.log2(prev_errors[name] / err)
                    else:
                        rates[name] = np.nan
                    prev_errors[name] = err

                r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
                r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
                r_fc = f"{rates['FACR']:6.2f}" if not np.isnan(rates['FACR']) else "   --"

                print(f"  {n:5d}|{e_fa:12.2e}|{r_fa}|{e_cr:12.2e}|{r_cr}|{e_fc:12.2e}|{r_fc}")
            except Exception as ex:
                import traceback
                print(f"  {n:5d}| ERROR: {ex}")
                traceback.print_exc()


def test_neumann_poisson():
    """Test Neumann Poisson (k^2=0) solver."""
    print("\n" + "=" * 80)
    print("NEUMANN POISSON TEST (k^2=0)")
    print("=" * 80)

    k2 = 0.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_neumann(k2)

    ns = [9, 17, 33, 65, 129]
    prev_errors = {}

    print(f"\n  {'n':>5s}|{'FA':>12s}|{'rate':>6s}|{'CR':>12s}|{'rate':>6s}|{'FACR':>12s}|{'rate':>6s}")
    print("  " + "-" * 72)

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        try:
            U_fa = fa_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
            U_cr = cr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')
            U_fc = facr_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann')

            e_fa = np.max(np.abs(U_fa - Ue))
            e_cr = np.max(np.abs(U_cr - Ue))
            e_fc = np.max(np.abs(U_fc - Ue))

            rates = {}
            for name, err in [('FA', e_fa), ('CR', e_cr), ('FACR', e_fc)]:
                if name in prev_errors:
                    rates[name] = np.log2(prev_errors[name] / err)
                else:
                    rates[name] = np.nan
                prev_errors[name] = err

            r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
            r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
            r_fc = f"{rates['FACR']:6.2f}" if not np.isnan(rates['FACR']) else "   --"

            print(f"  {n:5d}|{e_fa:12.2e}|{r_fa}|{e_cr:12.2e}|{r_cr}|{e_fc:12.2e}|{r_fc}")
        except Exception as ex:
            import traceback
            print(f"  {n:5d}| ERROR: {ex}")
            traceback.print_exc()


def test_k_zero_regression():
    """Verify that k^2=0 gives same results as Poisson solvers."""
    print("\n" + "=" * 80)
    print("k^2=0 REGRESSION TEST (should match Poisson)")
    print("=" * 80)

    def u_exact(x, y): return np.sin(np.pi * x) * np.sin(np.pi * y)
    def f_rhs(x, y): return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    def bc(x, y): return 0.0

    n = 33
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Ue = u_exact(X, Y)

    U_helm = fa_helmholtz(n, f_rhs, bc, k2=0.0, bc_type='dirichlet')
    e_helm = np.max(np.abs(U_helm - Ue))

    print(f"\n  FA Helmholtz (k^2=0) error: {e_helm:.6e}")
    print(f"  Expected: ~1e-3 (2nd order, n=33)")

    # Also check that k^2=0 and k^2=0.001 give similar results for small k^2
    U_ksmall = fa_helmholtz(n, f_rhs, bc, k2=0.001, bc_type='dirichlet')
    e_diff = np.max(np.abs(U_helm - U_ksmall))
    print(f"  |U(k^2=0) - U(k^2=0.001)|: {e_diff:.6e} (should be small)")


def test_large_k():
    """Test with large wavenumber k."""
    print("\n" + "=" * 80)
    print("LARGE WAVENUMBER TEST")
    print("=" * 80)

    for k2 in [1.0, 10.0, 100.0, 1000.0]:
        u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
        n = 65
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        try:
            U = fa_helmholtz(n, f_rhs, bc, k2=k2)
            e = np.max(np.abs(U - Ue))
            print(f"  k^2={k2:8.1f}  (k={np.sqrt(k2):8.4f})  error={e:.6e}")
        except Exception as ex:
            print(f"  k^2={k2:8.1f}  ERROR: {ex}")


def test_performance():
    """Performance benchmark."""
    print("\n" + "=" * 80)
    print("PERFORMANCE BENCHMARK")
    print("=" * 80)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

    ns = [33, 65, 129, 257, 513]

    solvers = [
        ("FA", lambda n: fa_helmholtz(n, f_rhs, bc, k2=k2)),
        ("CR", lambda n: cr_helmholtz(n, f_rhs, bc, k2=k2)),
        ("FACR", lambda n: facr_helmholtz(n, f_rhs, bc, k2=k2)),
        ("FFT9", lambda n: fft9_helmholtz(n, f_rhs, bc, k2=k2)),
    ]

    print(f"\n  k^2={k2}")
    print(f"\n  {'n':>5s}", end="")
    for nm, _ in solvers:
        print(f"|{nm:>12s}", end="")
    print()
    print("  " + "-" * (6 + 14 * len(solvers)))

    for n in ns:
        print(f"  {n:5d}", end="")
        for _, fn in solvers:
            t0 = time.time()
            fn(n)
            dt = time.time() - t0
            print(f"|{dt:12.4f}s", end="")
        print()


def test_facr_l():
    """Test FACR(l) with different l values."""
    print("\n" + "=" * 80)
    print("FACR(l) VARIANTS FOR HELMHOLTZ")
    print("=" * 80)

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)
    n = 257
    N = n - 2
    lopt = max(0, int(np.round(np.log2(max(1, np.log2(N))))))

    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    Uex = u_exact(X, Y)

    print(f"\n  Grid n={n}, k^2={k2}, optimal l~{lopt}\n")
    print(f"  {'l':>3s}|{'Time(ms)':>10s}|{'Error':>12s}|Note")
    print("  " + "-" * 48)

    for lv in range(0, min(9, int(np.log2(N)))):
        t0 = time.time()
        U = facr_helmholtz(n, f_rhs, bc, k2=k2, l=lv)
        ms = (time.time() - t0) * 1000
        e = np.max(np.abs(U - Uex))
        note = "Pure FA" if lv == 0 else ("***OPT***" if lv == lopt else "")
        print(f"  {lv:3d}|{ms:10.2f}|{e:12.2e}|{note}")


def test_mixed_bc():
    """Test mixed boundary conditions."""
    print("\n" + "=" * 80)
    print("MIXED BC TEST (x=Neumann, y=Dirichlet)")
    print("=" * 80)

    k2 = 10.0

    # u(x,y) = cos(pi*x) * sin(pi*y)
    # du/dx = -pi*sin(pi*x) = 0 at x=0,1  (Neumann in x)
    # u = 0 at y=0,1  (Dirichlet in y)
    def u_exact(x, y):
        return np.cos(np.pi * x) * np.sin(np.pi * y)

    def f_rhs(x, y):
        return (2 * np.pi**2 + k2) * np.cos(np.pi * x) * np.sin(np.pi * y)

    def bc_n(x, y):
        return 0.0  # homogeneous Neumann in x

    def bc_d(x, y):
        return 0.0  # homogeneous Dirichlet in y

    # For mixed BC, we need a single bc_func. Since both are homogeneous,
    # this works. For non-homogeneous BC, we'd need separate bc functions.
    def bc_func(x, y):
        return 0.0

    ns = [9, 17, 33, 65, 129]
    prev_errors = {}

    print(f"\n  k^2 = {k2}")
    print(f"  u(x,y) = cos(pi*x) * sin(pi*y)")
    print(f"\n  {'n':>5s}|{'FA':>12s}|{'rate':>6s}|{'CR':>12s}|{'rate':>6s}|{'FACR':>12s}|{'rate':>6s}")
    print("  " + "-" * 72)

    for n in ns:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        try:
            U_fa = fa_helmholtz(n, f_rhs, bc_func, k2=k2, bc_type=('N', 'D'))
            U_cr = cr_helmholtz(n, f_rhs, bc_func, k2=k2, bc_type=('N', 'D'))
            U_fc = facr_helmholtz(n, f_rhs, bc_func, k2=k2, bc_type=('N', 'D'))

            e_fa = np.max(np.abs(U_fa - Ue))
            e_cr = np.max(np.abs(U_cr - Ue))
            e_fc = np.max(np.abs(U_fc - Ue))

            rates = {}
            for name, err in [('FA', e_fa), ('CR', e_cr), ('FACR', e_fc)]:
                if name in prev_errors:
                    rates[name] = np.log2(prev_errors[name] / err)
                else:
                    rates[name] = np.nan
                prev_errors[name] = err

            r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
            r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
            r_fc = f"{rates['FACR']:6.2f}" if not np.isnan(rates['FACR']) else "   --"

            print(f"  {n:5d}|{e_fa:12.2e}|{r_fa}|{e_cr:12.2e}|{r_cr}|{e_fc:12.2e}|{r_fc}")
        except Exception as ex:
            import traceback
            print(f"  {n:5d}| ERROR: {ex}")
            traceback.print_exc()


def test_nonhomogeneous_dirichlet():
    """Test non-homogeneous Dirichlet BC — validates the F += bc/h² sign fix.

    Uses u(x,y) = sin(2*pi*x) * sin(3*pi*y) which has non-zero BC values.
    f(x,y) = ((2*pi)^2 + (3*pi)^2 + k^2) * sin(2*pi*x) * sin(3*pi*y)
    BC: u|_{dO} = sin(2*pi*x) * sin(3*pi*y)  (non-homogeneous)
    """
    print("\n" + "=" * 80)
    print("NON-HOMOGENEOUS DIRICHLET BC TEST (Sign Bug Verification)")
    print("=" * 80)

    for k2 in [0.0, 10.0, 100.0]:
        print(f"\n--- k^2 = {k2} ---")

        def u_exact(x, y):
            return np.sin(2 * np.pi * x) * np.sin(3 * np.pi * y)

        def f_rhs(x, y):
            return ((2 * np.pi)**2 + (3 * np.pi)**2 + k2) * \
                   np.sin(2 * np.pi * x) * np.sin(3 * np.pi * y)

        def bc(x, y):
            return u_exact(x, y)

        ns = [9, 17, 33, 65]
        prev_errors = {}

        print(f"\n  {'n':>5s}|{'FA':>12s}|{'rate':>6s}|{'CR':>12s}|{'rate':>6s}|{'FFT9':>12s}|{'rate':>6s}")
        print("  " + "-" * 72)

        for n in ns:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            try:
                U_fa = fa_helmholtz(n, f_rhs, bc, k2=k2)
                U_cr = cr_helmholtz(n, f_rhs, bc, k2=k2)
                U_f9 = fft9_helmholtz(n, f_rhs, bc, k2=k2)

                e_fa = np.max(np.abs(U_fa - Ue))
                e_cr = np.max(np.abs(U_cr - Ue))
                e_f9 = np.max(np.abs(U_f9 - Ue))

                rates = {}
                for name, err in [('FA', e_fa), ('CR', e_cr), ('FFT9', e_f9)]:
                    if name in prev_errors:
                        rates[name] = np.log2(prev_errors[name] / err)
                    else:
                        rates[name] = np.nan
                    prev_errors[name] = err

                r_fa = f"{rates['FA']:6.2f}" if not np.isnan(rates['FA']) else "   --"
                r_cr = f"{rates['CR']:6.2f}" if not np.isnan(rates['CR']) else "   --"
                r_f9 = f"{rates['FFT9']:6.2f}" if not np.isnan(rates['FFT9']) else "   --"

                print(f"  {n:5d}|{e_fa:12.2e}|{r_fa}|{e_cr:12.2e}|{r_cr}|{e_f9:12.2e}|{r_f9}")
            except Exception as ex:
                import traceback
                print(f"  {n:5d}| ERROR: {ex}")
                traceback.print_exc()


def run_all_tests():
    """Run all tests."""
    test_correctness()
    test_convergence()
    test_neumann_bc()
    test_neumann_poisson()
    test_mixed_bc()
    test_nonhomogeneous_dirichlet()
    test_k_zero_regression()
    test_large_k()
    test_facr_l()
    test_performance()


if __name__ == "__main__":
    run_all_tests()
