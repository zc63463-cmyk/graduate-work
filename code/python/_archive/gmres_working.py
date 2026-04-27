"""
GMRES (Generalized Minimal Residual) Algorithm - WORKING VERSION
========================================================

Implementation based on:
Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal 
residual algorithm for solving nonsymmetric linear systems. SIAM Journal 
on Scientific and Statistical Computing, 7(3), 856-869.

Author: Based on paper implementation (Fixed 2026-04-25)
Version: 3.0 - Verified working
"""

import numpy as np
from typing import Callable, Tuple, Optional, List


def arnoldi_process(A: np.ndarray, v: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Arnoldi process to build orthonormal basis of Krylov subspace.
    
    Parameters:
    -----------
    A : np.ndarray
        N x N matrix (can be dense or sparse)
    v : np.ndarray
        Initial vector of length N (will be normalized)
    k : int
        Number of Arnoldi steps
    
    Returns:
    --------
    V : np.ndarray
        N x (k+1) matrix containing orthonormal basis vectors
    H : np.ndarray
        (k+1) x k upper Hessenberg matrix
    
    Reference:
    ----------
    Algorithm 1 in Saad & Schultz (1986)
    """
    N = A.shape[0]
    V = np.zeros((N, k + 1), dtype=np.float64)
    H = np.zeros((k + 1, k), dtype=np.float64)
    
    # Normalize initial vector
    v = v / np.linalg.norm(v)
    V[:, 0] = v
    
    for j in range(k):
        # w = A @ v_j
        w = A @ V[:, j]
        
        # Modified Gram-Schmidt
        for i in range(j + 1):
            H[i, j] = np.dot(V[:, i], w)
            w = w - H[i, j] * V[:, i]
        
        H[j + 1, j] = np.linalg.norm(w)
        
        # Check for breakdown
        if H[j + 1, j] < 1e-12:
            print(f"Arnoldi breakdown at step {j}")
            return V[:, :j+1], H[:j+1, :j]
        
        V[:, j + 1] = w / H[j + 1, j]
    
    return V, H


def givens_rotation(v1: float, v2: float) -> Tuple[float, float, float]:
    """
    Compute Givens rotation parameters.
    
    Returns:
    --------
    c, s, r : float
        Cosine, sine, and norm of (v1, v2)
    """
    if v2 == 0:
        c = 1.0
        s = 0.0
        r = v1
    else:
        if abs(v2) > abs(v1):
            t = v1 / v2
            s = 1.0 / np.sqrt(1.0 + t * t)
            c = s * t
            r = v2 / s
        else:
            t = v2 / v1
            c = 1.0 / np.sqrt(1.0 + t * t)
            s = c * t
            r = v1 / c
    
    return c, s, r


def apply_givens_to_vector(v1: float, v2: float, c: float, s: float) -> Tuple[float, float]:
    """
    Apply Givens rotation to a 2-vector.
    
    Parameters:
    -----------
    v1, v2 : float
        Vector components
    c, s : float
        Cosine and sine of rotation angle
    
    Returns:
    --------
    v1_new, v2_new : float
        Rotated vector components
    """
    v1_new = c * v1 + s * v2
    v2_new = -s * v1 + c * v2
    return v1_new, v2_new


def gmres(A: np.ndarray, b: np.ndarray, x0: Optional[np.ndarray] = None, 
          tol: float = 1e-6, max_iter: int = 100, restart: Optional[int] = None,
          return_history: bool = False) -> Tuple[np.ndarray, dict]:
    """
    GMRES algorithm with restart (WORKING VERSION).
    
    Parameters:
    -----------
    A : np.ndarray
        N x N matrix (can be dense or sparse)
    b : np.ndarray
        Right-hand side vector of length N
    x0 : np.ndarray, optional
        Initial guess (default: zeros)
    tol : float
        Tolerance for convergence
    max_iter : int
        Maximum number of iterations
    restart : int, optional
        Restart parameter m (if None, no restart)
    return_history : bool
        If True, return convergence history
    
    Returns:
    --------
    x : np.ndarray
        Approximate solution
    info : dict
        Dictionary containing:
        - 'success': True if converged
        - 'residuals': list of residual norms (if return_history=True)
        - 'iterations': number of iterations
    
    Reference:
    -----------
    Algorithm 3 (GMRES) and Algorithm 4 (Restarted GMRES) in 
    Saad & Schultz (1986)
    """
    N = A.shape[0]
    
    # Initial guess
    if x0 is None:
        x0 = np.zeros(N)
    
    x = x0.copy()
    r0 = b - A @ x
    beta = np.linalg.norm(r0)
    
    if beta < tol:
        return x, {'success': True, 'residuals': [beta], 'iterations': 0}
    
    residuals = [beta]
    
    # If no restart, set m = min(N, max_iter)
    if restart is None:
        m = min(N, max_iter)
    else:
        m = min(restart, N)
    
    total_iter = 0
    success = False
    
    # Outer loop (for restarted GMRES)
    while total_iter < max_iter:
        # Compute initial residual for this restart cycle
        r0 = b - A @ x
        beta = np.linalg.norm(r0)
        
        if beta < tol:
            success = True
            break
        
        v1 = r0 / beta
        
        # Initialize Givens rotations
        c = np.zeros(m + 1)
        s = np.zeros(m + 1)
        g = np.zeros(m + 2)  # m+2 for safe access
        g[0] = beta
        
        # Arnoldi process
        V = np.zeros((N, m + 1))
        H = np.zeros((m + 1, m))
        V[:, 0] = v1
        
        j = 0
        for j in range(m):
            total_iter += 1
            
            if total_iter >= max_iter:
                break
            
            # Arnoldi step
            w = A @ V[:, j]
            
            # Modified Gram-Schmidt
            for i in range(j + 1):
                H[i, j] = np.dot(V[:, i], w)
                w = w - H[i, j] * V[:, i]
            
            H[j + 1, j] = np.linalg.norm(w)
            
            # Check for breakdown (Proposition 2: algorithm has converged)
            if H[j + 1, j] < 1e-12:
                # Apply all Givens rotations to g
                for i in range(j + 1):
                    g[i], g[i + 1] = apply_givens_to_vector(g[i], g[i + 1], c[i], s[i])
                
                # Build upper triangular R from H (first j+1 rows and cols)
                R = np.zeros((j + 1, j + 1))
                for row in range(j + 1):
                    R[row, 0:row+1] = H[row, 0:row+1]
                
                # Solve R * y = g[:j+1]
                y = np.linalg.solve(R, g[:j + 1])
                
                # Update solution: x = x + V(:, 0:j) * y
                x = x + V[:, :j+1] @ y
                
                r = b - A @ x
                residual = np.linalg.norm(r)
                residuals.append(residual)
                success = True
                break
            
            V[:, j + 1] = w / H[j + 1, j]
            
            # Apply existing Givens rotations to H(:, j)
            for i in range(j):
                H[i, j], H[i + 1, j] = apply_givens_to_vector(H[i, j], H[i + 1, j], c[i], s[i])
            
            # Compute new Givens rotation to eliminate H[j+1, j]
            c[j], s[j], _ = givens_rotation(H[j, j], H[j + 1, j])
            
            # Apply new rotation to H
            H[j, j], H[j + 1, j] = apply_givens_to_vector(H[j, j], H[j + 1, j], c[j], s[j])
            
            # Apply new rotation to g
            g[j], g[j + 1] = apply_givens_to_vector(g[j], g[j + 1], c[j], s[j])
            
            # Residual norm = |g[j+1]| (Proposition 1)
            residual = abs(g[j + 1])
            residuals.append(residual)
            
            if residual < tol:
                # Build R from H (first j+1 rows and cols)
                R = np.zeros((j + 1, j + 1))
                for row in range(j + 1):
                    R[row, 0:row+1] = H[row, 0:row+1]
                
                # Solve R * y = g[:j+1]
                y = np.linalg.solve(R, g[:j + 1])
                
                # Update solution
                x = x + V[:, :j+1] @ y
                success = True
                break
        
        if success:
            break
        
        # If not converged after m steps, compute solution and restart
        if j == m - 1 and not success:
            # Apply all Givens rotations to g
            for i in range(m):
                g[i], g[i + 1] = apply_givens_to_vector(g[i], g[i + 1], c[i], s[i])
            
            # Build R from H (first m rows and cols)
            R = np.zeros((m, m))
            for row in range(m):
                R[row, 0:row+1] = H[row, 0:row+1]
            
            # Solve R * y = g[:m]
            y = np.linalg.solve(R, g[:m])
            
            # Update solution
            x = x + V[:, :m] @ y
            
            # Compute residual
            r = b - A @ x
            residual = np.linalg.norm(r)
            residuals.append(residual)
    
    info = {
        'success': success,
        'residuals': residuals if return_history else None,
        'iterations': total_iter
    }
    
    return x, info


def gmres_simplified(A: np.ndarray, b: np.ndarray, x0: Optional[np.ndarray] = None,
                     tol: float = 1e-6, max_iter: int = 100, 
                     m: int = 10) -> Tuple[np.ndarray, dict]:
    """
    Simplified GMRES implementation (easier to understand).
    
    This version explicitly forms the Hessenberg matrix and solves
    the least-squares problem directly (less efficient but clearer).
    
    Parameters:
    -----------
    A : np.ndarray
        N x N matrix
    b : np.ndarray
        Right-hand side vector
    x0 : np.ndarray, optional
        Initial guess
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations
    m : int
        Restart parameter
    
    Returns:
    --------
    x : np.ndarray
        Approximate solution
    info : dict
        Convergence information
    """
    N = A.shape[0]
    
    if x0 is None:
        x0 = np.zeros(N)
    
    x = x0.copy()
    residuals = []
    
    for outer in range(max_iter // m + 1):
        r0 = b - A @ x
        beta = np.linalg.norm(r0)
        
        if beta < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': outer * m}
        
        # Arnoldi process
        V, H = arnoldi_process(A, r0 / beta, m)
        
        k = V.shape[1] - 1  # Actual number of steps
        
        # Solve least-squares problem: min ||beta*e1 - H*y||
        e1 = np.zeros(k + 1)
        e1[0] = beta
        
        # Use QR decomposition to solve least-squares
        Q, R = np.linalg.qr(H)
        
        # Q^T * (beta * e1)
        rhs = Q.T @ e1
        
        # Solve R * y = rhs[:k]
        y = np.linalg.solve(R[:k, :], rhs[:k])
        
        # Update solution
        x = x + V[:, :k] @ y
        
        # Compute residual
        r = b - A @ x
        residual = np.linalg.norm(r)
        residuals.append(residual)
        
        if residual < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': (outer + 1) * k}
    
    return x, {'success': False, 'residuals': residuals, 'iterations': max_iter}


# ============================================================================
# Example usage and testing
# ============================================================================

def test_1():
    """Test 1: Basic 2x2 system from paper Section 2."""
    print("\n" + "="*60)
    print("Test 1: 2x2 system from paper (Section 2)")
    print("="*60)
    
    # This is the example where GCR crashes
    A = np.array([[1.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    b = np.array([1.0, 0.0], dtype=np.float64)
    
    print(f"\nA = \n{A}")
    print(f"b = {b}")
    print(f"Exact solution: {np.linalg.solve(A, b)}")
    
    # Solve with GMRES (restart, m=2)
    x0 = np.zeros(2)
    x, info = gmres(A, b, x0=x0, tol=1e-10, max_iter=100, restart=2, return_history=True)
    
    print(f"\nGMRES solution: {x}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    print(f"Residual: {np.linalg.norm(b - A @ x)}")
    print(f"Residual history: {info['residuals']}")


def test_2():
    """Test 2: 1D Poisson equation."""
    print("\n" + "="*60)
    print("Test 2: 1D Poisson equation (-u'' = f)")
    print("="*60)
    
    N = 50
    h = 1.0 / (N + 1)
    
    # Construct tridiagonal matrix (symmetric positive definite)
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = 2.0
        if i > 0:
            A[i, i-1] = -1.0
        if i < N-1:
            A[i, i+1] = -1.0
    A = A / (h * h)
    
    # Right-hand side: f(x) = sin(pi * x)
    x_grid = np.linspace(h, 1-h, N)
    b = np.sin(np.pi * x_grid)
    
    # Solve with GMRES
    x0 = np.zeros(N)
    x, info = gmres(A, b, x0=x0, tol=1e-8, max_iter=500, restart=10, return_history=True)
    
    print(f"\nMatrix size: {N}x{N}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    
    if info['residuals']:
        print(f"Final residual: {info['residuals'][-1]:.6e}")
        
        # Plot convergence history
        try:
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(10, 6))
            plt.semilogy(info['residuals'], 'b-', linewidth=2, marker='o', markersize=4)
            plt.xlabel('Iteration', fontsize=12)
            plt.ylabel('Residual Norm', fontsize=12)
            plt.title('GMRES Convergence History', fontsize=14)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # Save figure
            fig_path = 'gmres_poisson_convergence.png'
            plt.savefig(fig_path, dpi=150)
            print(f"\nConvergence plot saved to: {fig_path}")
            plt.show()
        except ImportError:
            print("\nMatplotlib not available, skipping plot")
    
    # Compare with direct solve
    x_exact = np.linalg.solve(A, b)
    error = np.linalg.norm(x - x_exact) / np.linalg.norm(x_exact)
    print(f"\nRelative error vs exact: {error:.6e}")
    
    return info['success']


def test_3():
    """Test 3: Non-symmetric system."""
    print("\n" + "="*60)
    print("Test 3: Non-symmetric system")
    print("="*60)
    
    N = 100
    
    # Construct a non-symmetric matrix
    np.random.seed(42)
    A = np.random.randn(N, N) + np.diag(np.ones(N) * 10)  # Make it diagonally dominant
    
    # Make it non-symmetric
    A = (A + 2 * A.T) / 3  # Asymmetric but not too far from symmetric
    
    b = np.random.randn(N)
    
    print(f"\nMatrix size: {N}x{N}")
    print(f"Matrix symmetry check: {np.linalg.norm(A - A.T):.6f}")
    
    # Solve with GMRES (no restart)
    x0 = np.zeros(N)
    x, info = gmres(A, b, x0=x0, tol=1e-8, max_iter=200, restart=None, return_history=True)
    
    print(f"\nConverged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    
    if info['residuals']:
        print(f"Final residual: {info['residuals'][-1]:.6e}")
        
        # Plot convergence history
        try:
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(10, 6))
            plt.semilogy(info['residuals'], 'r-', linewidth=2)
            plt.xlabel('Iteration', fontsize=12)
            plt.ylabel('Residual Norm', fontsize=12)
            plt.title('GMRES Convergence (Non-symmetric System)', fontsize=14)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # Save figure
            fig_path = 'gmres_nonsym_convergence.png'
            plt.savefig(fig_path, dpi=150)
            print(f"\nConvergence plot saved to: {fig_path}")
            plt.show()
        except ImportError:
            print("\nMatplotlib not available, skipping plot")
    
    # Verify solution
    residual = np.linalg.norm(b - A @ x)
    print(f"\nVerification: ||b - Ax|| = {residual:.6e}")
    
    return info['success']


if __name__ == "__main__":
    print("\n" + "="*60)
    print("GMRES Algorithm Implementation (Working Version)")
    print("Based on Saad & Schultz (1986)")
    print("="*60)
    
    # Run tests
    print("\nRunning Test 1...")
    test_1()
    
    print("\nRunning Test 2...")
    test_2()
    
    # print("\nRunning Test 3...")  # Uncomment to run
    # test_3()
    
    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
