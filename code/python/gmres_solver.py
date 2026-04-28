#!/usr/bin/env python3
"""
GMRES Solver for 2D Helmholtz Equation
=======================================

Integrates the GMRES (Generalized Minimal Residual) iterative method
with the 2D Helmholtz equation solver framework.

Based on:
    Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal
    residual algorithm for solving nonsymmetric linear systems.
    SIAM J. Sci. Stat. Comput., 7(3), 856-869.

This module provides:
    1. gmres()          - Standalone GMRES solver (Givens rotation based)
    2. build_helmholtz_matrix() - Build 2D Helmholtz sparse matrix
    3. gmres_helmholtz() - GMRES solver with Helmholtz-compatible interface
    4. solve_helmholtz() - Unified interface (mirrors helmholtz_solver.py)

Grid conventions (same as helmholtz_solver.py):
    - Dirichlet BC: unknowns at INTERIOR nodes (n-2 points), DST-I
    - Neumann BC:   unknowns at ALL nodes (n points), ghost-point approach
    - indexing='ij': U[i,j] = u(x_i, y_j)

Test problems (reused from helmholtz_solver.py):
    - Dirichlet: u = sin(pi*x)*sin(pi*y), f = (2*pi^2+k^2)*u, u=0 on boundary
    - Neumann:   u = cos(pi*x)*cos(pi*y), f = (2*pi^2+k^2)*u, du/dn=0
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator
import time


# ============================================================================
# Core GMRES Implementation
# ============================================================================

def gmres(A, b, x0=None, tol=1e-6, max_iter=100, restart=30,
          return_history=False):
    """
    GMRES with restart (Givens rotation based, numerically stable).

    Parameters
    ----------
    A : array_like or LinearOperator
        N x N matrix or linear operator (must support A @ v)
    b : array_like
        Right-hand side vector of length N
    x0 : array_like, optional
        Initial guess (default: zeros)
    tol : float
        Convergence tolerance (absolute residual norm)
    max_iter : int
        Maximum total iterations
    restart : int
        Restart parameter m (GMRES(m))
    return_history : bool
        If True, return per-iteration residual history

    Returns
    -------
    x : ndarray
        Approximate solution
    info : dict
        'success': bool, 'residuals': list, 'iterations': int
    """
    N = b.shape[0] if hasattr(b, 'shape') else len(b)
    b = np.asarray(b).reshape(-1)

    if x0 is None:
        x = np.zeros(N)
    else:
        x = np.asarray(x0).copy()

    r = b - A @ x
    beta = np.linalg.norm(r)

    if beta < tol:
        return x, {'success': True, 'residuals': [beta], 'iterations': 0}

    residuals = [beta]
    m = min(restart, N, max_iter)
    total_iter = 0
    success = False

    while total_iter < max_iter:
        r = b - A @ x
        beta = np.linalg.norm(r)

        if beta < tol:
            success = True
            break

        V = np.zeros((N, m + 1))
        H = np.zeros((m + 1, m))
        V[:, 0] = r / beta

        c = np.zeros(m + 1)
        s = np.zeros(m + 1)
        g = np.zeros(m + 1)
        g[0] = beta

        j = 0
        for j in range(m):
            total_iter += 1
            if total_iter > max_iter:
                break

            # Arnoldi step
            w = A @ V[:, j]

            # Modified Gram-Schmidt
            for i in range(j + 1):
                H[i, j] = np.dot(V[:, i], w)
                w = w - H[i, j] * V[:, i]

            H[j + 1, j] = np.linalg.norm(w)

            # Breakdown check
            if H[j + 1, j] < 1e-12:
                # Krylov subspace is invariant — exact solution in K_{j+1}
                # Apply rotations 0..j-1 to H[:, j] (not yet done in normal flow)
                for i in range(j):
                    temp = c[i] * H[i, j] + s[i] * H[i + 1, j]
                    H[i + 1, j] = -s[i] * H[i, j] + c[i] * H[i + 1, j]
                    H[i, j] = temp
                # H[j+1,j] ≈ 0 already, so j-th rotation is identity
                c[j] = 1.0
                s[j] = 0.0
                H[j + 1, j] = 0.0
                # Apply ONLY rotation j to g (rotations 0..j-1 were already
                # applied incrementally in previous iterations)
                # Identity rotation: g unchanged
                # Solve R*y = g[:j+1]
                # After Givens rotations, H[:j+1, :j+1] is upper triangular
                R = H[:j+1, :j+1].copy()
                y = np.linalg.solve(R, g[:j+1])
                x = x + V[:, :j+1] @ y
                residuals.append(np.linalg.norm(b - A @ x))
                success = True
                break

            V[:, j + 1] = w / H[j + 1, j]

            # Apply existing Givens rotations to H[:, j]
            for i in range(j):
                temp = c[i] * H[i, j] + s[i] * H[i + 1, j]
                H[i + 1, j] = -s[i] * H[i, j] + c[i] * H[i + 1, j]
                H[i, j] = temp

            # Compute new Givens rotation (numerically stable)
            if H[j + 1, j] == 0:
                c[j] = 1.0
                s[j] = 0.0
                r_val = H[j, j]
            elif abs(H[j + 1, j]) > abs(H[j, j]):
                t = H[j, j] / H[j + 1, j]
                s[j] = 1.0 / np.sqrt(1.0 + t * t)
                c[j] = s[j] * t
                r_val = H[j + 1, j] / s[j]
            else:
                t = H[j + 1, j] / H[j, j]
                c[j] = 1.0 / np.sqrt(1.0 + t * t)
                s[j] = c[j] * t
                r_val = H[j, j] / c[j]

            H[j, j] = r_val
            H[j + 1, j] = 0.0

            # Apply rotation to g
            temp_g = c[j] * g[j] + s[j] * g[j + 1]
            g[j + 1] = -s[j] * g[j] + c[j] * g[j + 1]
            g[j] = temp_g

            # Residual norm = |g[j+1]|
            residual = abs(g[j + 1])
            residuals.append(residual)

            if residual < tol:
                # After Givens rotations, H[:j+1, :j+1] is upper triangular
                R = H[:j+1, :j+1].copy()
                y = np.linalg.solve(R, g[:j+1])
                x = x + V[:, :j+1] @ y
                success = True
                break

        if success:
            break

        # Restart: solve and update
        if j == m - 1 and not success:
            # After Givens rotations, H[:m, :m] is upper triangular
            # g has already been transformed incrementally — do NOT re-apply rotations
            R = H[:m, :m].copy()
            y = np.linalg.solve(R, g[:m])
            x = x + V[:, :m] @ y
            residuals.append(np.linalg.norm(b - A @ x))

    info = {
        'success': success,
        'residuals': residuals if return_history else None,
        'iterations': total_iter
    }
    return x, info


# ============================================================================
# 2D Helmholtz Sparse Matrix Construction
# ============================================================================

def build_helmholtz_matrix(n, k2=None, bc_type='dirichlet', sx=1.0, sy=1.0, sigma=None):
    """
    Build 2D Helmholtz operator sparse matrix for 5-point finite difference.

    Equation: (-nabla^2 + sigma)u = f  on [0, sx] x [0, sy]

    Use sigma for the unified framework, or k2 for backward compatibility
    (k2 implies sigma = +k^2, i.e., modified Helmholtz).

    Ghost-point Neumann BC:
        At x-boundary (i=0): u_{-1,j} = u_{1,j}  =>  neighbor coeff = -2/h^2
        At y-boundary (j=0): u_{i,-1} = u_{i,1}  =>  neighbor coeff = -2/h^2
        Diagonal stays 4/h^2 + k^2 (ghost point absorbed into off-diagonal)

    Parameters
    ----------
    n : int
        Grid size (n x n nodes including boundaries)
    k2 : float
        Wavenumber squared (backward compat, implies sigma=+k^2)
    sigma : float, optional
        Unified Helmholtz parameter. If provided, takes precedence over k2.
    bc_type : str or tuple
        'dirichlet', 'neumann', or ('D'/'N', 'D'/'N')
    sx, sy : float
        Domain sizes

    Returns
    -------
    A : scipy.sparse.csr_matrix
        Sparse Helmholtz operator matrix
    info : dict
        Grid information: h, Nx, Ny, bc_x, bc_y
    """
    # Normalize BC type
    if isinstance(bc_type, str):
        bc_type_lower = bc_type.lower()
        if bc_type_lower in ('dirichlet', 'd'):
            bc_x, bc_y = 'D', 'D'
        elif bc_type_lower in ('neumann', 'n'):
            bc_x, bc_y = 'N', 'N'
        else:
            raise ValueError(f"Unknown bc_type '{bc_type}'")
    elif isinstance(bc_type, (tuple, list)):
        bc_x, bc_y = bc_type[0], bc_type[1]
    else:
        raise TypeError(f"bc_type must be str or tuple, got {type(bc_type)}")

    # Resolve sigma parameter
    if sigma is not None:
        sigma_val = sigma
    elif k2 is not None:
        sigma_val = k2
    else:
        sigma_val = 0.0

    h = sx / (n - 1)
    h2 = h * h

    # Grid dimensions
    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2
    N = Nx * Ny

    rows, cols, vals = [], [], []

    for i in range(Nx):
        for j in range(Ny):
            idx = i * Ny + j

            # Main diagonal: always 4/h^2 + sigma
            rows.append(idx)
            cols.append(idx)
            vals.append(4.0 / h2 + sigma_val)

            # x-direction: left neighbor (i-1)
            if bc_x == 'D':
                # Dirichlet: only interior neighbors, coeff = -1/h^2
                if i > 0:
                    rows.append(idx)
                    cols.append((i - 1) * Ny + j)
                    vals.append(-1.0 / h2)
            else:
                # Neumann ghost-point:
                #   At i=0: u_{-1,j} = u_{1,j} => absorbed as coeff -2/h^2 to i=1
                #   At i>0: coeff = -1/h^2
                if i == 0:
                    # Ghost point: no i=-1 entry, -2/h^2 goes to i=1
                    rows.append(idx)
                    cols.append(1 * Ny + j)
                    vals.append(-2.0 / h2)
                else:
                    rows.append(idx)
                    cols.append((i - 1) * Ny + j)
                    if i == Nx - 1:
                        # Ghost point at right: u_{Nx,j} = u_{Nx-2,j}
                        # coeff to i-1 = -2/h^2
                        vals.append(-2.0 / h2)
                    else:
                        vals.append(-1.0 / h2)

            # x-direction: right neighbor (i+1) — only for non-Neumann-boundary cases
            # (Neumann boundary cases already handled above)
            if bc_x == 'D':
                if i < Nx - 1:
                    rows.append(idx)
                    cols.append((i + 1) * Ny + j)
                    vals.append(-1.0 / h2)
            else:
                # Neumann: interior connections only (boundaries handled above)
                if 0 < i < Nx - 1:
                    rows.append(idx)
                    cols.append((i + 1) * Ny + j)
                    vals.append(-1.0 / h2)

            # y-direction: bottom neighbor (j-1)
            if bc_y == 'D':
                if j > 0:
                    rows.append(idx)
                    cols.append(i * Ny + (j - 1))
                    vals.append(-1.0 / h2)
            else:
                if j == 0:
                    rows.append(idx)
                    cols.append(i * Ny + 1)
                    vals.append(-2.0 / h2)
                else:
                    rows.append(idx)
                    cols.append(i * Ny + (j - 1))
                    if j == Ny - 1:
                        vals.append(-2.0 / h2)
                    else:
                        vals.append(-1.0 / h2)

            # y-direction: top neighbor (j+1)
            if bc_y == 'D':
                if j < Ny - 1:
                    rows.append(idx)
                    cols.append(i * Ny + (j + 1))
                    vals.append(-1.0 / h2)
            else:
                if 0 < j < Ny - 1:
                    rows.append(idx)
                    cols.append(i * Ny + (j + 1))
                    vals.append(-1.0 / h2)

    A = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    info = {
        'h': h, 'Nx': Nx, 'Ny': Ny,
        'bc_x': bc_x, 'bc_y': bc_y,
        'sx': sx, 'sy': sy, 'n': n,
    }

    return A, info


def _build_rhs(n, f_func, bc_func, k2, bc_x, bc_y, sx=1.0, sy=1.0):
    """
    Build RHS vector for GMRES solver.

    For Dirichlet: evaluate f at interior nodes, subtract boundary contributions
    For Neumann:   evaluate f at all nodes (including boundaries)

    Returns
    -------
    b : ndarray of length Nx*Ny
        Flattened RHS vector
    """
    h = sx / (n - 1)
    x_node = np.linspace(0, sx, n)
    y_node = np.linspace(0, sy, n)

    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2

    h2 = h * h

    if bc_x == 'D' and bc_y == 'D':
        # Pure Dirichlet
        x_int = x_node[1:-1]
        y_int = y_node[1:-1]
        X, Y = np.meshgrid(x_int, y_int, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        # Subtract boundary contributions
        bc_l = _eval_bc(bc_func, x_node[0], y_node, n)
        bc_r = _eval_bc(bc_func, x_node[-1], y_node, n)
        bc_b = _eval_bc(bc_func, x_node, y_node[0], n)
        bc_t = _eval_bc(bc_func, x_node, y_node[-1], n)

        F[0, :] += bc_l[1:-1] / h2
        F[-1, :] += bc_r[1:-1] / h2
        F[:, 0] += bc_b[1:-1] / h2
        F[:, -1] += bc_t[1:-1] / h2

    elif bc_x == 'N' and bc_y == 'N':
        # Pure Neumann
        X, Y = np.meshgrid(x_node, y_node, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

    elif bc_x == 'N' and bc_y == 'D':
        # Mixed: x=Neumann, y=Dirichlet
        X, Y = np.meshgrid(x_node, y_node[1:-1], indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        bc_b = _eval_bc(bc_func, x_node, y_node[0], n)
        bc_t = _eval_bc(bc_func, x_node, y_node[-1], n)
        F[:, 0] += bc_b / h2
        F[:, -1] += bc_t / h2

    elif bc_x == 'D' and bc_y == 'N':
        # Mixed: x=Dirichlet, y=Neumann
        X, Y = np.meshgrid(x_node[1:-1], y_node, indexing='ij')
        F = np.asarray(f_func(X, Y), dtype=float)

        bc_l = _eval_bc(bc_func, x_node[0], y_node, n)
        bc_r = _eval_bc(bc_func, x_node[-1], y_node, n)
        F[0, :] += bc_l / h2
        F[-1, :] += bc_r / h2
    else:
        raise ValueError(f"Unexpected BC: ({bc_x}, {bc_y})")

    return F.reshape(-1)


def _eval_bc(bc_func, coord1, coord2, n):
    """Evaluate BC function, ensuring output is length-n array."""
    val = bc_func(coord1, coord2)
    arr = np.atleast_1d(np.asarray(np.squeeze(val), dtype=float))
    if arr.shape[0] == 1:
        return np.full(n, float(arr[0]))
    return arr


def _assemble_solution(u_vec, info, bc_func, sx=1.0, sy=1.0):
    """
    Assemble full n x n grid solution from GMRES solution vector.

    Parameters
    ----------
    u_vec : ndarray of length Nx*Ny
        Flattened solution from GMRES
    info : dict
        Grid information from build_helmholtz_matrix
    bc_func : callable
        Boundary condition function
    sx, sy : float
        Domain sizes

    Returns
    -------
    U : ndarray of shape (n, n)
        Full grid solution
    """
    n = info['n']
    Nx, Ny = info['Nx'], info['Ny']
    bc_x, bc_y = info['bc_x'], info['bc_y']
    h = info['h']

    x_node = np.linspace(0, sx, n)
    y_node = np.linspace(0, sy, n)

    U_solve = u_vec.reshape(Nx, Ny)
    U = np.zeros((n, n))

    if bc_x == 'D' and bc_y == 'D':
        # Interior only, fill boundaries from bc_func
        bc_l = _eval_bc(bc_func, x_node[0], y_node, n)
        bc_r = _eval_bc(bc_func, x_node[-1], y_node, n)
        bc_b = _eval_bc(bc_func, x_node, y_node[0], n)
        bc_t = _eval_bc(bc_func, x_node, y_node[-1], n)

        U[0, :] = bc_l
        U[-1, :] = bc_r
        U[:, 0] = bc_b
        U[:, -1] = bc_t
        U[1:-1, 1:-1] = U_solve

    elif bc_x == 'N' and bc_y == 'N':
        # Full grid (ghost-point approach)
        # Note: GMRES solves the unsymmetric system directly (no symmetrization)
        U = U_solve

    elif bc_x == 'N' and bc_y == 'D':
        # x=Neumann (n rows), y=Dirichlet (n-2 cols)
        bc_b = _eval_bc(bc_func, x_node, y_node[0], n)
        bc_t = _eval_bc(bc_func, x_node, y_node[-1], n)
        U[:, 0] = bc_b
        U[:, -1] = bc_t
        U[:, 1:-1] = U_solve

    elif bc_x == 'D' and bc_y == 'N':
        # x=Dirichlet (n-2 rows), y=Neumann (n cols)
        bc_l = _eval_bc(bc_func, x_node[0], y_node, n)
        bc_r = _eval_bc(bc_func, x_node[-1], y_node, n)
        U[0, :] = bc_l
        U[-1, :] = bc_r
        U[1:-1, :] = U_solve

    return U


# ============================================================================
# GMRES Helmholtz Solver
# ============================================================================

def gmres_helmholtz(n, f_func, bc_func, k2=None, bc_type='dirichlet',
                     sx=1.0, sy=1.0, tol=1e-10, max_iter=None,
                     restart=30, return_history=False, sigma=None):
    """
    Solve 2D Helmholtz equation using GMRES iterative method.

    Equation: (-nabla^2 + sigma) u = f

    Parameters
    ----------
    n : int
        Grid size (n x n)
    f_func : callable
        RHS function f(x, y)
    bc_func : callable
        Boundary condition function
    k2 : float, optional
        Wavenumber squared (backward compat, implies sigma=+k^2)
    sigma : float, optional
        Unified Helmholtz parameter. Takes precedence over k2.
    bc_type : str or tuple
        'dirichlet', 'neumann', or ('D'/'N', 'D'/'N')
    sx, sy : float
        Domain sizes
    tol : float
        GMRES convergence tolerance (absolute residual)
    max_iter : int, optional
        Maximum GMRES iterations (default: 10 * N)
    restart : int
        GMRES restart parameter m
    return_history : bool
        Whether to return convergence history

    Returns
    -------
    U : ndarray of shape (n, n)
        Solution on full grid
    info : dict
        Solver information (iterations, residuals, etc.)
    """
    # Resolve sigma
    if sigma is not None:
        sigma_val = sigma
    elif k2 is not None:
        sigma_val = k2
    else:
        sigma_val = 0.0

    # Build sparse matrix and grid info
    A, grid_info = build_helmholtz_matrix(n, sigma=sigma_val, bc_type=bc_type, sx=sx, sy=sy)

    # Build RHS
    bc_x, bc_y = grid_info['bc_x'], grid_info['bc_y']
    b = _build_rhs(n, f_func, bc_func, sigma_val, bc_x, bc_y, sx, sy)

    N = A.shape[0]
    if max_iter is None:
        max_iter = max(10 * N, 1000)

    # Solve with GMRES
    t0 = time.time()
    u_vec, gmres_info = gmres(A, b, tol=tol, max_iter=max_iter,
                               restart=restart, return_history=True)
    solve_time = time.time() - t0

    # Assemble full grid solution
    U = _assemble_solution(u_vec, grid_info, bc_func, sx, sy)

    info = {
        'success': gmres_info['success'],
        'iterations': gmres_info['iterations'],
        'residuals': gmres_info['residuals'],
        'time': solve_time,
        'matrix_size': N,
    }

    return U, info


# ============================================================================
# Matrix-Free GMRES for Helmholtz (avoids explicit matrix storage)
# ============================================================================

def _helmholtz_matvec(v_flat, n, k2, bc_x, bc_y, sx, sy):
    """
    Matrix-vector product for 2D Helmholtz operator (matrix-free).

    Computes y = (-nabla^2 + k^2) * u for the 5-point stencil.
    """
    h = sx / (n - 1)
    h2 = h * h

    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2

    u = v_flat.reshape(Nx, Ny)
    y = np.zeros_like(u)

    # Interior points
    for i in range(Nx):
        for j in range(Ny):
            val = (4.0 / h2 + k2) * u[i, j]

            # x-neighbors
            if bc_x == 'D':
                # Dirichlet: boundary values are 0 (already in RHS)
                if i > 0:
                    val -= u[i-1, j] / h2
                if i < Nx - 1:
                    val -= u[i+1, j] / h2
            else:
                # Neumann ghost-point: u_{-1} = u_1, u_{Nx} = u_{Nx-2}
                if i > 0:
                    val -= u[i-1, j] / h2
                else:
                    val -= u[1, j] / h2  # ghost: u_{-1} = u_1
                if i < Nx - 1:
                    val -= u[i+1, j] / h2
                else:
                    val -= u[-2, j] / h2  # ghost: u_{Nx} = u_{Nx-2}

            # y-neighbors
            if bc_y == 'D':
                if j > 0:
                    val -= u[i, j-1] / h2
                if j < Ny - 1:
                    val -= u[i, j+1] / h2
            else:
                # Neumann ghost-point
                if j > 0:
                    val -= u[i, j-1] / h2
                else:
                    val -= u[i, 1] / h2
                if j < Ny - 1:
                    val -= u[i, j+1] / h2
                else:
                    val -= u[i, -2] / h2

            y[i, j] = val

    return y.reshape(-1)


def gmres_helmholtz_matfree(n, f_func, bc_func, k2=0.0, bc_type='dirichlet',
                             sx=1.0, sy=1.0, tol=1e-10, max_iter=None,
                             restart=30, return_history=False):
    """
    Solve 2D Helmholtz using matrix-free GMRES (no explicit sparse matrix).

    Uses scipy.sparse.linalg.LinearOperator for the matrix-vector product,
    avoiding the memory cost of storing the full sparse matrix.
    """
    if isinstance(bc_type, str):
        bc_type_lower = bc_type.lower()
        if bc_type_lower in ('dirichlet', 'd'):
            bc_x, bc_y = 'D', 'D'
        elif bc_type_lower in ('neumann', 'n'):
            bc_x, bc_y = 'N', 'N'
        else:
            raise ValueError(f"Unknown bc_type '{bc_type}'")
    elif isinstance(bc_type, (tuple, list)):
        bc_x, bc_y = bc_type[0], bc_type[1]
    else:
        raise TypeError(f"bc_type must be str or tuple")

    Nx = n if bc_x == 'N' else n - 2
    Ny = n if bc_y == 'N' else n - 2
    N = Nx * Ny

    # Create LinearOperator
    def matvec(v):
        return _helmholtz_matvec(v, n, k2, bc_x, bc_y, sx, sy)

    A_op = LinearOperator((N, N), matvec=matvec, dtype=float)

    # Build RHS
    b = _build_rhs(n, f_func, bc_func, k2, bc_x, bc_y, sx, sy)

    if max_iter is None:
        max_iter = max(10 * N, 1000)

    # Solve
    t0 = time.time()
    u_vec, gmres_info = gmres(A_op, b, tol=tol, max_iter=max_iter,
                               restart=restart, return_history=True)
    solve_time = time.time() - t0

    # Assemble
    grid_info = {'n': n, 'Nx': Nx, 'Ny': Ny, 'bc_x': bc_x, 'bc_y': bc_y,
                 'h': sx / (n - 1)}
    U = _assemble_solution(u_vec, grid_info, bc_func, sx, sy)

    info = {
        'success': gmres_info['success'],
        'iterations': gmres_info['iterations'],
        'residuals': gmres_info['residuals'],
        'time': solve_time,
        'matrix_size': N,
    }

    return U, info


# ============================================================================
# Test Problems (same as helmholtz_solver.py for fair comparison)
# ============================================================================

def test_problem_dirichlet(k2=0.0, sx=1.0, sy=1.0):
    """Test problem with Dirichlet BC: u = sin(pi*x/sx)*sin(pi*y/sy).

    WARNING: This is an EIGENFUNCTION of the discrete Laplacian.
    GMRES converges in 1 iteration — only useful for verifying
    solver correctness, NOT for studying GMRES convergence behavior.
    For realistic GMRES experiments, use test_problem_gaussian() or
    test_problem_polynomial() instead.
    """
    def u_exact(x, y):
        return np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)

    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + k2) * \
               np.sin(np.pi * x / sx) * np.sin(np.pi * y / sy)

    def bc(x, y):
        return 0.0

    return u_exact, f_rhs, bc


def test_problem_neumann(k2=0.0, sx=1.0, sy=1.0):
    """Test problem with Neumann BC: u = cos(pi*x/sx)*cos(pi*y/sy).

    WARNING: Same eigenfunction issue as Dirichlet version.
    """
    def u_exact(x, y):
        return np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)

    def f_rhs(x, y):
        return ((np.pi / sx)**2 + (np.pi / sy)**2 + k2) * \
               np.cos(np.pi * x / sx) * np.cos(np.pi * y / sy)

    def bc(x, y):
        return 0.0  # homogeneous Neumann

    return u_exact, f_rhs, bc


def test_problem_polynomial(k2=0.0, sx=1.0, sy=1.0):
    """Test problem with Dirichlet BC: u = x^2*(sx-x)^2 * y^2*(sy-y)^2.

    This is NOT an eigenfunction, so GMRES needs many iterations.
    The polynomial vanishes on all boundaries with zero derivative,
    making it ideal for convergence order verification.

    Exact solution: u = x^2(sx-x)^2 * y^2(sy-y)^2
    f = -Delta u + k^2 * u  (computed analytically)
    """
    def u_exact(x, y):
        px = x**2 * (sx - x)**2
        py = y**2 * (sy - y)**2
        return px * py

    def f_rhs(x, y):
        px = x**2 * (sx - x)**2
        py = y**2 * (sy - y)**2
        # d^2(px)/dx^2 = 2(sx-x)^2 - 8x(sx-x) + 2x^2 = 2sx^2 - 12sx*x + 12x^2
        d2px = 2*(sx - x)**2 - 8*x*(sx - x) + 2*x**2
        d2py = 2*(sy - y)**2 - 8*y*(sy - y) + 2*y**2
        # -Delta u = -(d2px * py + px * d2py)
        return -(d2px * py + px * d2py) + k2 * px * py

    def bc(x, y):
        return 0.0

    return u_exact, f_rhs, bc


def test_problem_gaussian(k2=0.0, sx=1.0, sy=1.0,
                          x0=0.3, y0=0.7, sigma=0.1):
    """Gaussian source problem with Dirichlet BC.

    f = A * exp(-((x-x0)^2 + (y-y0)^2) / (2*sigma^2))

    This has NO closed-form exact solution. Use FFT direct solver
    to compute reference solution for error verification.
    Ideal for studying GMRES convergence behavior since the RHS
    excites many eigenmodes.

    Parameters
    ----------
    x0, y0 : float
        Center of Gaussian source
    sigma : float
        Width of Gaussian source
    """
    A_amp = 100.0  # amplitude

    def f_rhs(x, y):
        return A_amp * np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))

    def bc(x, y):
        return 0.0

    return None, f_rhs, bc  # No exact solution


def test_problem_multimode(k2=0.0, sx=1.0, sy=1.0, n_modes=10, seed=42):
    """Multi-mode superposition with Dirichlet BC.

    u = sum_{m,n} c_{mn} * sin(m*pi*x/sx) * sin(n*pi*y/sy)

    Each term is an eigenfunction, but the superposition excites
    many eigenvectors simultaneously. GMRES needs ~n_modes^2
    iterations (one per excited mode). This is the "bridge" between
    the trivial single-mode problem and the fully general case.

    Parameters
    ----------
    n_modes : int
        Number of modes in each direction (total ~n_modes^2 terms)
    seed : int
        Random seed for coefficients c_{mn}
    """
    rng = np.random.RandomState(seed)
    # Generate random coefficients for modes (m,n) with 1 <= m,n <= n_modes
    coeffs = {}
    for m in range(1, n_modes + 1):
        for n in range(1, n_modes + 1):
            coeffs[(m, n)] = rng.randn()

    def u_exact(x, y):
        u = np.zeros_like(x, dtype=float)
        for (m, n), c in coeffs.items():
            u += c * np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy)
        return u

    def f_rhs(x, y):
        f = np.zeros_like(x, dtype=float)
        for (m, n), c in coeffs.items():
            eigval = (m * np.pi / sx)**2 + (n * np.pi / sy)**2 + k2
            f += c * eigval * np.sin(m * np.pi * x / sx) * np.sin(n * np.pi * y / sy)
        return f

    def bc(x, y):
        return 0.0

    return u_exact, f_rhs, bc


# ============================================================================
# Test Suite
# ============================================================================

def test_gmres_basic():
    """Test GMRES on a simple 2x2 system."""
    print("=" * 60)
    print("GMRES Basic Test: 2x2 system")
    print("=" * 60)

    A = np.array([[1.0, 1.0], [1.0, 0.0]])
    b = np.array([1.0, 0.0])

    x, info = gmres(A, b, tol=1e-10, restart=2, return_history=True)
    x_exact = np.linalg.solve(A, b)

    print(f"  GMRES solution: {x}")
    print(f"  Exact solution: {x_exact}")
    print(f"  Error: {np.linalg.norm(x - x_exact):.2e}")
    print(f"  Converged: {info['success']}, Iterations: {info['iterations']}")


def test_gmres_dirichlet():
    """Test GMRES Helmholtz solver with Dirichlet BC."""
    print("\n" + "=" * 70)
    print("GMRES Helmholtz Solver - Dirichlet BC")
    print("=" * 70)

    for k2 in [0.0, 1.0, 10.0]:
        print(f"\n  k^2 = {k2}")
        u_exact, f_rhs, bc = test_problem_dirichlet(k2)

        print(f"  {'n':>5s} | {'Error':>12s} | {'Rate':>6s} | {'Iters':>6s} | {'Time(ms)':>10s}")
        print("  " + "-" * 55)

        prev_err = None
        for n in [9, 17, 33, 65]:
            x = np.linspace(0, 1, n)
            X, Y = np.meshgrid(x, x, indexing='ij')
            Ue = u_exact(X, Y)

            U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2,
                                       tol=1e-12, restart=30,
                                       return_history=True)
            err = np.max(np.abs(U - Ue))
            rate = np.log2(prev_err / err) if prev_err is not None else float('nan')
            prev_err = err

            r_str = f"{rate:6.2f}" if not np.isnan(rate) else "   --"
            print(f"  {n:5d} | {err:12.2e} | {r_str} | {info['iterations']:6d} | {info['time']*1000:10.2f}")


def test_gmres_vs_fft():
    """Compare GMRES with FFT direct solver."""
    print("\n" + "=" * 70)
    print("GMRES vs FFT Direct Solver Comparison")
    print("=" * 70)

    import sys
    sys.path.insert(0, '.')
    from helmholtz_solver import fa_helmholtz, _helmholtz_test_problem_dirichlet

    k2 = 10.0
    u_exact, f_rhs, bc = _helmholtz_test_problem_dirichlet(k2)

    print(f"\n  k^2 = {k2}")
    print(f"  {'n':>5s} | {'GMRES Error':>12s} | {'FFT Error':>12s} | {'GMRES(ms)':>10s} | {'FFT(ms)':>10s}")
    print("  " + "-" * 65)

    for n in [9, 17, 33, 65, 129]:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        # GMRES
        U_gmres, info_g = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12)
        err_g = np.max(np.abs(U_gmres - Ue))
        t_g = info_g['time'] * 1000

        # FFT
        t0 = time.time()
        U_fft = fa_helmholtz(n, f_rhs, bc, k2=k2)
        t_f = (time.time() - t0) * 1000
        err_f = np.max(np.abs(U_fft - Ue))

        print(f"  {n:5d} | {err_g:12.2e} | {err_f:12.2e} | {t_g:10.2f} | {t_f:10.2f}")


def test_gmres_convergence_curve():
    """Plot GMRES convergence curve for different k values."""
    print("\n" + "=" * 70)
    print("GMRES Convergence Curves")
    print("=" * 70)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    n = 33
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')

    # Left: Dirichlet, different k
    ax = axes[0]
    for k2 in [0.0, 1.0, 10.0, 100.0]:
        u_exact, f_rhs, bc = test_problem_dirichlet(k2)
        U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12,
                                   restart=30, return_history=True)
        if info['residuals']:
            res = np.array(info['residuals'])
            res = res / res[0]  # Normalize
            ax.semilogy(res, label=f'k²={k2}')

    ax.set_xlabel('Iteration')
    ax.set_ylabel('Relative Residual')
    ax.set_title(f'GMRES Convergence (Dirichlet, n={n})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: Different grid sizes
    ax = axes[1]
    k2 = 10.0
    u_exact, f_rhs, bc = test_problem_dirichlet(k2)
    for n in [17, 33, 65]:
        U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, tol=1e-12,
                                   restart=30, return_history=True)
        if info['residuals']:
            res = np.array(info['residuals'])
            res = res / res[0]
            ax.semilogy(res, label=f'n={n}')

    ax.set_xlabel('Iteration')
    ax.set_ylabel('Relative Residual')
    ax.set_title(f'GMRES Convergence (k²={k2}, Dirichlet)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = 'thesis/figures/gmres_convergence.png'
    import os
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    plt.savefig(fig_path, dpi=150)
    print(f"  Convergence plot saved to {fig_path}")
    plt.close()


def test_gmres_neumann():
    """Test GMRES with Neumann BC."""
    print("\n" + "=" * 70)
    print("GMRES Helmholtz Solver - Neumann BC")
    print("=" * 70)

    k2 = 1.0
    u_exact, f_rhs, bc = test_problem_neumann(k2)

    print(f"\n  k^2 = {k2}")
    print(f"  {'n':>5s} | {'Error':>12s} | {'Rate':>6s} | {'Iters':>6s}")
    print("  " + "-" * 45)

    prev_err = None
    for n in [9, 17, 33, 65]:
        x = np.linspace(0, 1, n)
        X, Y = np.meshgrid(x, x, indexing='ij')
        Ue = u_exact(X, Y)

        U, info = gmres_helmholtz(n, f_rhs, bc, k2=k2, bc_type='neumann',
                                   tol=1e-12, return_history=True)
        err = np.max(np.abs(U - Ue))
        rate = np.log2(prev_err / err) if prev_err is not None else float('nan')
        prev_err = err

        r_str = f"{rate:6.2f}" if not np.isnan(rate) else "   --"
        print(f"  {n:5d} | {err:12.2e} | {r_str} | {info['iterations']:6d}")


if __name__ == "__main__":
    test_gmres_basic()
    test_gmres_dirichlet()
    test_gmres_neumann()
    # test_gmres_vs_fft()  # Uncomment when helmholtz_solver.py is accessible
    # test_gmres_convergence_curve()  # Uncomment to generate plots
    print("\n" + "=" * 60)
    print("All GMRES tests completed!")
    print("=" * 60)
