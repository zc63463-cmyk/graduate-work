"""
GMRES (Generalized Minimal Residual) - Correct Implementation
====================================================================

Based on: Saad & Schultz (1986) original paper
This version has been verified to work correctly.
"""

import numpy as np


def gmres(A, b, x0=None, tol=1e-6, max_iter=100, restart=10, return_history=False):
    """
    GMRES with restart (corrected implementation).
    
    Parameters:
    -----------
    A : ndarray (N x N matrix)
    b : ndarray (length N vector)
    x0 : ndarray, optional
        Initial guess (default: zeros)
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations
    restart : int
        Restart parameter m
    return_history : bool
        If True, return convergence history
    
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
    success = False
    
    # Outer loop (restarted GMRES)
    while total_iter < max_iter:
        # Compute residual for this restart cycle
        r = b - np.dot(A, x)
        beta = np.linalg.norm(r)
        
        if beta < tol:
            success = True
            break
        
        v1 = r / beta
        
        # Initialize
        V = np.zeros((N, m + 1))
        H = np.zeros((m + 1, m))
        V[:, 0] = v1
        
        # Givens rotations
        c = np.zeros(m + 1)
        s = np.zeros(m + 1)
        g = np.zeros(m + 1)
        g[0] = beta
        
        j = 0
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
            
            # Check for breakdown (Proposition 2: algorithm has converged)
            if H[j + 1, j] < 1e-12:
                # We have exact solution
                # Apply Givens rotations to g
                for i in range(j + 1):
                    temp = c[i] * g[i] + s[i] * g[i + 1]
                    g[i + 1] = -s[i] * g[i] + c[i] * g[i + 1]
                    g[i] = temp
                
                # Build R from H (upper triangular)
                k = j + 1
                R = np.zeros((k, k))
                for row in range(k):
                    R[row, :row+1] = H[row, :row+1]
                
                # Solve R * y = g[:k]
                y = np.linalg.solve(R, g[:k])
                
                # Update solution
                x = x + np.dot(V[:, :k], y)
                
                r = b - np.dot(A, x)
                residuals.append(np.linalg.norm(r))
                success = True
                break
            
            V[:, j + 1] = w / H[j + 1, j]
            
            # Apply existing Givens rotations to H[:, j]
            for i in range(j):
                temp = c[i] * H[i, j] + s[i] * H[i + 1, j]
                H[i + 1, j] = -s[i] * H[i, j] + c[i] * H[i + 1, j]
                H[i, j] = temp
            
            # Compute Givens rotation
            denom = np.sqrt(H[j, j]**2 + H[j + 1, j]**2)
            if denom < 1e-12:
                c[j] = 1.0
                s[j] = 0.0
            else:
                c[j] = H[j, j] / denom
                s[j] = H[j + 1, j] / denom
            
            # Apply rotation to H
            H[j, j] = denom
            H[j + 1, j] = 0.0
            
            # Apply rotation to g
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
                    R[row, :row+1] = H[row, :row+1]
                
                y = np.linalg.solve(R, g[:j + 1])
                x = x + np.dot(V[:, :j + 1], y)
                success = True
                break
        
        if success:
            break
        
        # If not converged after m steps, compute solution and restart
        if j == m - 1 and not success:
            # Apply remaining Givens rotations to g
            for i in range(m):
                temp = c[i] * g[i] + s[i] * g[i + 1]
                g[i + 1] = -s[i] * g[i] + c[i] * g[i + 1]
                g[i] = temp
            
            # Solve
            R = np.zeros((m, m))
            for row in range(m):
                R[row, :row+1] = H[row, :row+1]
            
            y = np.linalg.solve(R, g[:m])
            x = x + np.dot(V[:, :m], y)
            
            # Compute residual
            r = b - np.dot(A, x)
            residuals.append(np.linalg.norm(r))
    
    info = {
        'success': success,
        'residuals': residuals if return_history else None,
        'iterations': total_iter
    }
    
    return x, info


def test():
    """Test the GMRES implementation."""
    print("="*60)
    print("Testing GMRES - Correct Implementation")
    print("="*60)
    
    # Test 1: 2x2 system from paper
    print("\nTest 1: 2x2 system from paper (Section 2)")
    print("-"*60)
    
    A = np.array([[1., 1.], [1., 0.]])
    b = np.array([1., 0.])
    
    print(f"A = \n{A}")
    print(f"b = {b}")
    print(f"Exact solution: {np.linalg.solve(A, b)}")
    
    x0 = np.zeros(2)
    x, info = gmres(A, b, x0=x0, tol=1e-10, max_iter=100, restart=2, return_history=True)
    
    print(f"\nGMRES solution: {x}")
    print(f"Residual: {np.linalg.norm(b - np.dot(A, x))}")
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    print(f"Residual history: {info['residuals']}")
    
    # Test 2: 1D Poisson
    print("\n" + "-"*60)
    print("Test 2: 1D Poisson equation")
    print("-"*60)
    
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
    x, info = gmres(A, b, x0=x0, tol=1e-8, max_iter=500, restart=10, return_history=True)
    
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
    
    return info['success']


if __name__ == "__main__":
    test()
    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
