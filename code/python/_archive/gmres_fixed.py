"""
GMRES (Generalized Minimal Residual) Algorithm - FIXED VERSION
=====================================================

Implementation based on:
Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal 
residual algorithm for solving nonsymmetric linear systems. SIAM Journal 
on Scientific and Statistical Computing, 7(3), 856-869.

Author: Based on paper implementation (Fixed 2026-04-25)
"""

import numpy as np
from scipy.linalg import solve_triangular
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


def gmres(A: np.ndarray, b: np.ndarray, x0: Optional[np.ndarray] = None, 
          tol: float = 1e-6, max_iter: int = 100, restart: Optional[int] = None,
          return_history: bool = False) -> Tuple[np.ndarray, dict]:
    """
    GMRES algorithm with restart (FIXED VERSION).
    
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
        
        # Arnoldi process with incremental QR
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
                # Compute solution using Eq. (8): x = x0 + V_k * y_k
                # Here k = j+1 (1-indexed), so we have j+1 columns in V
                # Build upper triangular R from H (first j+1 rows and cols)
                R = np.zeros((j + 1, j + 1))
                for row in range(j + 1):
                    for col in range(row + 1):  # Only upper triangular + diagonal
                        R[row, col] = H[row, col]
                
                # Apply Givens rotations to g (first j+2 elements)
                for i in range(j + 1):
                    temp = c[i] * g[i] + s[i] * g[i + 1]
                    g[i + 1] = -s[i] * g[i] + c[i] * g[i + 1]
                    g[i] = temp
                
                # Solve R * y = g[:j+1]
                y = solve_triangular(R, g[:j + 1], lower=False)
                
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
                temp = c[i] * H[i, j] + s[i] * H[i + 1, j]
                H[i + 1, j] = -s[i] * H[i, j] + c[i] * H[i + 1, j]
                H[i, j] = temp
            
            # Compute new Givens rotation to eliminate H[j+1, j]
            c[j], s[j], r_temp = givens_rotation(H[j, j], H[j + 1, j])
            H[j, j] = r_temp
            H[j + 1, j] = 0.0
            
            # Apply new rotation to g
            temp = c[j] * g[j] + s[j] * g[j + 1]
            g[j + 1] = -s[j] * g[j] + c[j] * g[j + 1]
            g[j] = temp
            
            # Residual norm = |g[j+1]| (Proposition 1)
            residual = abs(g[j + 1])
            residuals.append(residual)
            
            if residual < tol:
                # Compute solution
                R = np.zeros((j + 1, j + 1))
                for row in range(j + 1):
                    for col in range(min(row, j) + 1):
                        R[row, col] = H[row, col]
                
                g_small = g[:j + 1]
                y = solve_triangular(R, g_small, lower=False)
                x = x + V[:, :j+1] @ y
                success = True
                break
        
        if success:
            break
        
        # If not converged after m steps, compute solution and restart
        if j == m - 1 and not success:
            # Apply all Givens rotations to g
            for i in range(m):
                temp = c[i] * g[i] + s[i] * g[i + 1]
                g[i + 1] = -s[i] * g[i] + c[i] * g[i + 1]
                g[i] = temp
            
            # Compute solution after m steps
            R = np.zeros((m, m))
            for row in range(m):
                for col in range(row + 1):
                    R[row, col] = H[row, col]
            
            g_small = g[:m]
            y = solve_triangular(R, g_small, lower=False)
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
        y = solve_triangular(R[:k, :], rhs[:k], lower=False)
        
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

def example_1():
    """Simple example: 2x2 system from paper (Section 2)."""
    print("=" * 60)
    print("Example 1: 2x2 system from paper (Section 2)")
    print("=" * 60)
    
    # This is the example where GCR crashes
    A = np.array([[1.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    b = np.array([1.0, 0.0], dtype=np.float64)
    
    print(f"\nA = \n{A}")
    print(f"b = {b}")
    
    # Solve with GMRES
    x0 = np.zeros(2)
    x, info = gmres(A, b, x0=x0, tol=1e-10, max_iter=100, restart=2, return_history=True)
    
    print(f"\nSolution: {x}")
    print(f"Exact solution: {np.linalg.solve(A, b)}")
    print(f"Residual: {np.linalg.norm(b - A @ x)}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    print(f"Residual history: {info['residuals']}")


def example_2():
    """Larger example: 1D Poisson equation."""
    print("\n" + "=" * 60)
    print("Example 2: 1D Poisson equation (-u'' = f)")
    print("=" * 60)
    
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
    
    # Right-hand side
    x_true = np.linspace(h, 1-h, N)
    b = A @ np.sin(np.pi * x_true)  # Example RHS
    
    # Solve with GMRES
    x0 = np.zeros(N)
    x, info = gmres(A, b, x0=x0, tol=1e-8, max_iter=500, restart=10, return_history=True)
    
    print(f"\nMatrix size: {N}x{N}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    if info['residuals']:
        print(f"Final residual: {info['residuals'][-1]:.6e}")
        
        # Plot convergence history
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.semilogy(info['residuals'], 'b-', linewidth=2, marker='o', markersize=4)
        plt.xlabel('Iteration', fontsize=12)
        plt.ylabel('Residual Norm', fontsize=12)
        plt.title('GMRES Convergence History', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('gmres_convergence.png', dpi=150)
        print("\nConvergence plot saved to 'gmres_convergence.png'")
        plt.show()


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("GMRES Algorithm Implementation (Fixed Version)")
    print("Based on Saad & Schultz (1986)")
    print("=" * 60)
    
    # Run examples
    example_1()
    example_2()
    
    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)
