"""
GMRES (Generalized Minimal Residual) - Final Working Version
====================================================================

This is a simplified but CORRECT implementation based on:
Saad & Schultz (1986). 

It uses numpy's least_squares solver to avoid singular matrix issues.
"""

import numpy as np


def gmres(A, b, x0=None, tol=1e-6, max_iter=100, restart=10):
    """
    Simplified GMRES with restart.
    
    Parameters:
    -----------
    A : ndarray (N x N matrix)
    b : ndarray (length N vector)
    x0 : ndarray, optional
        Initial guess
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations
    restart : int
        Restart parameter m
    
    Returns:
    --------
    x : ndarray
        Approximate solution
    info : dict
        Dictionary with convergence information
    """
    N = A.shape[0]
    b = b.reshape(-1)
    
    # Initial guess
    if x0 is None:
        x = np.zeros(N)
    else:
        x = x0.copy()
    
    # Initial residual
    r = b - np.dot(A, x)
    beta = np.linalg.norm(r)
    
    if beta < tol:
        return x, {'success': True, 'residuals': [beta], 'iterations': 0}
    
    residuals = [beta]
    m = min(restart, N, max_iter)
    total_iter = 0
    
    # Outer loop (restarted GMRES)
    while total_iter < max_iter:
        # Compute residual for this restart cycle
        r = b - np.dot(A, x)
        beta = np.linalg.norm(r)
        
        if beta < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': total_iter}
        
        # Arnoldi process
        V = np.zeros((N, m + 1))
        H = np.zeros((m + 1, m))
        V[:, 0] = r / beta
        
        for j in range(m):
            total_iter += 1
            
            if total_iter >= max_iter:
                break
            
            # Arnoldi step
            w = np.dot(A, V[:, j])
            
            # Modified Gram-Schmidt
            for i in range(j + 1):
                H[i, j] = np.dot(V[:, i], w)
                w = w - H[i, j] * V[:, i]
            
            H[j + 1, j] = np.linalg.norm(w)
            
            # Check for breakdown (Proposition 2)
            if H[j + 1, j] < 1e-12:
                # We have exact solution
                # Build V[:, :j+1] and solve least squares
                e1 = np.zeros(j + 1)
                e1[0] = beta
                
                # Use least squares to find y
                # s.t. ||beta*e1 - H[:j+1, :j+1] * y|| is minimized
                H_sub = H[:j+1, :j+1]
                y, _, _, _ = np.linalg.lstsq(H_sub, e1, rcond=None)
                
                # Update x
                x = x + np.dot(V[:, :j+1], y)
                
                r = b - np.dot(A, x)
                residuals.append(np.linalg.norm(r))
                return x, {'success': True, 'residuals': residuals, 'iterations': total_iter}
            
            V[:, j + 1] = w / H[j + 1, j]
        
        # After m steps: solve least squares problem
        e1 = np.zeros(m + 1)
        e1[0] = beta
        
        # Use numpy's least squares solver
        y, _, _, _ = np.linalg.lstsq(H, e1, rcond=None)
        
        # Update solution
        x = x + np.dot(V[:, :m], y)
        
        # Compute residual
        r = b - np.dot(A, x)
        residual = np.linalg.norm(r)
        residuals.append(residual)
        total_iter += m
        
        if residual < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': total_iter}
    
    # Max iterations reached
    return x, {'success': False, 'residuals': residuals, 'iterations': total_iter}


def test():
    """Test the GMRES implementation."""
    print("=" * 60)
    print("Testing GMRES - Final Working Version")
    print("=" * 60)
    
    # Test 1: 2x2 system from paper
    print("\nTest 1: 2x2 system from paper (Section 2)")
    print("-" * 60)
    
    A = np.array([[1., 1.], [1., 0.]])
    b = np.array([1., 0.])
    
    print(f"A = \n{A}")
    print(f"b = {b}")
    print(f"Exact solution: {np.linalg.solve(A, b)}")
    
    x0 = np.zeros(2)
    x, info = gmres(A, b, x0=x0, tol=1e-10, restart=2)
    
    print(f"\nGMRES solution: {x}")
    print(f"Residual: {np.linalg.norm(b - np.dot(A, x))}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    
    # Test 2: 1D Poisson equation
    print("\n" + "-" * 60)
    print("Test 2: 1D Poisson equation")
    print("-" * 60)
    
    N = 50
    h = 1.0 / (N + 1)
    
    # Construct tridiagonal matrix
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = 2.0
        if i > 0:
            A[i, i-1] = -1.0
        if i < N-1:
            A[i, i+1] = -1.0
    A = A / (h * h)
    
    # Right-hand side
    x_grid = np.linspace(h, 1-h, N)
    b = np.sin(np.pi * x_grid)
    
    print(f"Matrix size: {N}x{N}")
    
    x0 = np.zeros(N)
    x, info = gmres(A, b, x0=x0, tol=1e-8, restart=10)
    
    print(f"\nConverged: {info['success']}")
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
            
            fig_path = 'gmres_convergence.png'
            plt.savefig(fig_path, dpi=150)
            print(f"\nConvergence plot saved to: {fig_path}")
            plt.show()
        except ImportError:
            print("\nMatplotlib not available, skipping plot")
    
    # Compare with exact solution
    x_exact = np.linalg.solve(A, b)
    error = np.linalg.norm(x - x_exact) / np.linalg.norm(x_exact)
    print(f"\nRelative error vs exact: {error:.6e}")


if __name__ == "__main__":
    test()
    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
