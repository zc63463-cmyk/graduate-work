#!/usr/bin/env python
"""
Test GMRES implementation - focusing on 2x2 example first
"""

import numpy as np

def gmres_simple(A, b, x0=None, tol=1e-6, max_iter=100, m=10):
    """
    Super simple GMRES - just use least squares directly.
    This avoids all the Givens rotation complexities.
    """
    N = A.shape[0]
    b = b.reshape(-1)
    
    if x0 is None:
        x = np.zeros(N)
    else:
        x = x0.copy()
    
    r = b - np.dot(A, x)
    beta = np.linalg.norm(r)
    
    if beta < tol:
        return x, {'success': True, 'residuals': [beta], 'iterations': 0}
    
    residuals = [beta]
    total_iter = 0
    
    while total_iter < max_iter:
        r = b - np.dot(A, x)
        beta = np.linalg.norm(r)
        
        if beta < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': total_iter}
        
        # Arnoldi
        V = np.zeros((N, m +1))
        H = np.zeros((m +1, m))
        V[:, 0] = r / beta
        
        for j in range(m):
            w = np.dot(A, V[:, j])
            
            for i in range(j +1):
                H[i, j] = np.dot(V[:, i], w)
                w = w - H[i, j] * V[:, i]
            
            H[j +1, j] = np.linalg.norm(w)
            
            if H[j +1, j] < 1e-12:
                # Breakdown - exact solution found
                # Build R and solve
                k = j + 1
                R = np.zeros((k, k))
                for row in range(k):
                    R[row, :row+1] = H[row, :row+1]
                
                # RHS
                rhs = np.zeros(k)
                rhs[0] = beta
                
                # Apply Givens rotations to rhs (simplified: skip for now)
                # Just solve
                y = np.linalg.solve(R, rhs)
                
                # Update x
                x = x + np.dot(V[:, :k], y)
                
                r = b - np.dot(A, x)
                residuals.append(np.linalg.norm(r))
                return x, {'success': True, 'residuals': residuals, 'iterations': total_iter + 1}
            
            V[:, j +1] = w / H[j +1, j]
        
        # After m steps: solve least squares
        e1 = np.zeros(m +1)
        e1[0] = beta
        
        # Use QR
        Q, R = np.linalg.qr(H)
        rhs = np.dot(Q.T, e1)
        y = np.linalg.solve(R[:m, :], rhs[:m])
        
        x = x + np.dot(V[:, :m], y)
        
        r = b - np.dot(A, x)
        residual = np.linalg.norm(r)
        residuals.append(residual)
        total_iter += m
        
        if residual < tol:
            return x, {'success': True, 'residuals': residuals, 'iterations': total_iter}
    
    return x, {'success': False, 'residuals': residuals, 'iterations': total_iter}


def test_1():
    """Test 1: 2x2 system from paper."""
    print("Test 1: 2x2 system")
    print("-"*40)
    
    A = np.array([[1., 1.], [1., 0.]])
    b = np.array([1., 0.])
    
    print(f"A = \n{A}")
    print(f"b = {b}")
    print(f"Exact: {np.linalg.solve(A, b)}")
    
    x0 = np.zeros(2)
    x, info = gmres_simple(A, b, x0=x0, tol=1e-10, m=2)
    
    print(f"GMRES: {x}")
    print(f"Residual: {np.linalg.norm(b - np.dot(A, x))}")
    print(f"Success: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    print(f"History: {info['residuals']}")


def test_2():
    """Test 2: 1D Poisson."""
    print("\nTest 2: 1D Poisson")
    print("-"*40)
    
    N = 50
    h = 1.0 / (N +1)
    
    A = np.zeros((N, N))
    for i in range(N):
        A[i, i] = 2.0
        if i > 0:
            A[i, i-1] = -1.0
        if i < N-1:
            A[i, i+1] = -1.0
    A = A / (h * h)
    
    x_grid = np.linspace(h, 1-h, N)
    b = np.sin(np.pi * x_grid)
    
    x0 = np.zeros(N)
    x, info = gmres_simple(A, b, x0=x0, tol=1e-8, m=10)
    
    print(f"Converged: {info['success']}")
    print(f"Iterations: {info['iterations']}")
    if info['residuals']:
        print(f"Final residual: {info['residuals'][-1]:.6e}")
    
    # Plot
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.semilogy(info['residuals'], 'b-', linewidth=2, marker='o', markersize=4)
        plt.xlabel('Outer Iteration')
        plt.ylabel('Residual Norm')
        plt.title('GMRES Convergence')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('gmres_convergence.png', dpi=150)
        print("Plot saved to 'gmres_convergence.png'")
        plt.show()
    except ImportError:
        pass


if __name__ == "__main__":
    test_1()
    test_2()
    print("\nDone!")
