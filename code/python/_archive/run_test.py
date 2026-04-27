#!/usr/bin/env python
"""
Quick test for GMRES implementation.
"""

import sys
sys.path.insert(0, '.')

from gmres import gmres
import numpy as np

def test_1():
    """Test 1: 2x2 system from paper."""
    print("\nTest 1: 2x2 system from paper")
    print("-" * 40)
    
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


def test_2():
    """Test 2: 1D Poisson equation."""
    print("\nTest 2: 1D Poisson equation")
    print("-" * 40)
    
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
        plt.xlabel('Outer Iteration', fontsize=12)
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
    print("="*60)
    print("Testing GMRES - Simplified Version")
    print("="*60)
    
    test_1()
    test_2()
    
    print("\n" + "="*60)
    print("All tests completed!")
    print("="*60)
