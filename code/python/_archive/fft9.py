#!/usr/bin/env python3
"""
FFT9 Algorithm - Python Implementation
Based on: 
  - Houstis, E.N. and Papatheodorou, T.S. (1979)
    "High-Order Fast Elliptic Equation Solvers"
    ACM TOMS, Vol.5, No.4, pp.431-441
  - ACM Algorithm 543: FFT9

This is a simplified Python implementation for educational purposes.
For production use, consider using FFTW or other optimized libraries.
"""

import numpy as np
from scipy.fft import fft, ifft
import time


class FFT9Solver:
    """
    FFT9 Solver for Helmholtz-type and Poisson equations.
    
    Solves: alpha * u_xx + beta * u_yy + gamma * u = f(x,y)
    on rectangular domain [0, SX] x [0, SY]
    with Dirichlet boundary conditions.
    """
    
    def __init__(self, sx, sy, nx, ny, order=4):
        """
        Initialize FFT9 solver.
        
        Parameters:
        -----------
        sx, sy : float
            Domain dimensions in x and y directions
        nx, ny : int
            Number of grid lines (should be 2^k + 1)
        order : int
            Difference order (4 for Helmholtz, 6 for Poisson)
        """
        self.sx = sx
        self.sy = sy
        self.nx = nx
        self.ny = ny
        self.order = order
        
        # Grid spacing
        self.hx = sx / (nx - 1)
        self.hy = sy / (ny - 1)
        
        # Interior points
        self.imax = nx - 2
        self.jmax = ny - 2
        
        # Grid coordinates
        self.x = np.linspace(0, sx, nx)
        self.y = np.linspace(0, sy, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        
        print(f"FFT9 Solver initialized:")
        print(f"  Domain: [{0}, {sx}] x [{0}, {sy}]")
        print(f"  Grid: {nx} x {ny} ({self.imax} x {self.jmax} interior)")
        print(f"  Order: {order}")
        print(f"  Grid spacing: hx={self.hx:.6f}, hy={self.hy:.6f}")
    
    def set_pde_coefficients(self, alpha, beta, gamma):
        """
        Set PDE coefficients.
        
        PDE: alpha * u_xx + beta * u_yy + gamma * u = f
        
        For Poisson: alpha = beta = 1, gamma = 0
        For Helmholtz: alpha = beta = 1, gamma = -k^2
        """
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        
        # Compute difference scheme coefficients (Eq.2.1 in Paper 1)
        p = self.hx**2 / alpha
        q = self.hy**2 / beta
        
        self.a = 12 / (p**2 + q)
        self.b = 12 * p**2 / (p**2 + q) - 2
        self.c = (self.b + 2) * (1 - 0 * p / alpha) * p**2 / alpha - 2 * (alpha + self.beta) - 4
        
        print(f"  PDE coefficients: alpha={alpha}, beta={beta}, gamma={gamma}")
        print(f"  Difference coefficients: a={self.a:.6f}, b={self.b:.6f}, c={self.c:.6f}")
    
    def set_rhs_function(self, f_func):
        """Set right-hand side function f(x,y)"""
        self.f_func = f_func
    
    def set_boundary_conditions(self, bc_func):
        """Set boundary condition function bc(x,y)"""
        self.bc_func = bc_func
    
    def discretize(self):
        """
        Discretize PDE using 9-point high-order difference scheme.
        
        Returns:
        --------
        rhs : ndarray (nx-2, ny-2)
            Right-hand side vector after discretization
        """
        nx, ny = self.nx, self.ny
        hx, hy = self.hx, self.hy
        
        # Initialize RHS
        rhs = np.zeros((nx, ny))
        
        # Interior points
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                x, y = self.x[i], self.y[j]
                rhs[i, j] = self.f_func(x, y)
        
        # Apply 9-point stencil weights
        if self.order == 4:
            # Fourth-order 9-point stencil (Eq.2.1)
            # Stencil for Helmholtz-type equation
            weight_center = self.c
            weight_side = self.b
            weight_corner = self.a
            
            # Note: This is simplified. Full implementation needs
            # to construct the block tridiagonal system as in Paper 1
            print("  Using 4th-order 9-point stencil")
            
        elif self.order == 6:
            # Sixth-order stencil for Poisson equation
            # Stencil from Lynch (1977)
            print("  Using 6th-order 9-point stencil for Poisson")
            
        return rhs[1:-1, 1:-1]
    
    def solve_poisson_fft(self, f, bc):
        """
        Solve Poisson equation using FFT method.
        
        This is a simplified version. The full FFT9 algorithm uses:
        1. Odd-even reduction
        2. Fourier analysis
        3. Fourier synthesis
        
        Here we implement the basic FFT-based Poisson solver.
        """
        nx, ny = self.nx, self.ny
        hx, hy = self.hx, self.hy
        
        # Create grid
        x = np.linspace(0, self.sx, nx)
        y = np.linspace(0, self.sy, ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # Evaluate RHS
        F = f(X, Y)
        
        # Apply boundary conditions
        U = np.zeros((nx, ny))
        U[0, :] = bc(x[0], y)  # Left
        U[-1, :] = bc(x[-1], y)  # Right
        U[:, 0] = bc(x, y[0])  # Bottom
        U[:, -1] = bc(x, y[-1])  # Top
        
        # For Dirichlet BC, adjust RHS
        # (Simplified - full implementation needed for exact adjustment)
        
        # Fourier method for Poisson equation
        # Eigenvalues: lambda_k = 2/hx^2 * (cos(k*pi/(nx-1)) - 1) 
        #             + 2/hy^2 * (cos(l*pi/(ny-1)) - 1)
        
        # FFT-based solution (simplified)
        # In practice, use sine transform for Dirichlet BC
        
        # Using 2D FFT (for periodic BC - need adjustment for Dirichlet)
        # This is a placeholder for the actual FFT9 algorithm
        
        print("  Solving using FFT method...")
        print("  (Note: This is simplified. Full FFT9 implementation needed.)")
        
        # Placeholder: Return zeros
        # Actual implementation would:
        # 1. Perform odd-even reduction
        # 2. Apply FFT to decouple equations
        # 3. Solve tridiagonal systems
        # 4. Apply inverse FFT
        
        return U
    
    def solve(self):
        """
        Solve the PDE using FFT9 algorithm.
        
        Returns:
        --------
        U : ndarray (nx, ny)
            Solution matrix
        """
        print("\n" + "="*60)
        print("Solving with FFT9 algorithm...")
        print("="*60)
        
        start_time = time.time()
        
        if self.order == 6 and self.gamma == 0:
            print("Poisson equation with 6th-order discretization")
            # Use FFT method for Poisson
            U = self.solve_poisson_fft(self.f_func, self.bc_func)
        else:
            print(f"Helmholtz-type equation with {self.order}th-order discretization")
            # Would implement full FFT9 here
            print("Full FFT9 implementation needed (see Fortran code in Paper 2)")
            U = np.zeros((self.nx, self.ny))
        
        end_time = time.time()
        self.computation_time = end_time - start_time
        
        print(f"Solution computed in {self.computation_time:.4f} seconds")
        
        return U


def test_poisson_smooth():
    """
    Test 1: Poisson equation with smooth solution.
    
    Equation: u_xx + u_yy = f
    Domain: [0,1] x [0,1]
    Exact solution: u = 3xy(x-x^2)(y-y^2)
    """
    print("\n" + "="*60)
    print("Test 1: Poisson equation with smooth solution")
    print("="*60)
    
    # Exact solution
    def u_exact(x, y):
        return 3 * x * y * (x - x**2) * (y - y**2)
    
    # RHS
    def f_rhs(x, y):
        # u_xx + u_yy for u = 3xy(x-x^2)(y-y^2)
        # This is a simplified version
        return 6 * (y - y**2) * (1 - 6*x + 6*x**2) + 6 * (x - x**2) * (1 - 6*y + 6*y**2)
    
    # Boundary condition
    def bc(x, y):
        return u_exact(x, y)
    
    # Create solver
    sx, sy = 1.0, 1.0
    nx, ny = 33, 33  # 2^5 + 1
    order = 6  # Sixth-order for Poisson
    
    solver = FFT9Solver(sx, sy, nx, ny, order)
    solver.set_pde_coefficients(alpha=1.0, beta=1.0, gamma=0.0)
    solver.set_rhs_function(f_rhs)
    solver.set_boundary_conditions(bc)
    
    # Solve
    U = solver.solve()
    
    # Compute error at interior points
    X, Y = np.meshgrid(solver.x, solver.y, indexing='ij')
    U_exact = u_exact(X, Y)
    error = np.abs(U - U_exact)
    max_error = np.max(error[1:-1, 1:-1])
    
    print(f"\nResults:")
    print(f"  Max error (interior): {max_error:.6e}")
    print(f"  Grid: {nx} x {ny}")
    
    return max_error, solver.computation_time


def test_helmholtz():
    """
    Test 2: Helmholtz equation.
    
    Equation: u_xx + u_yy - k^2 * u = f
    """
    print("\n" + "="*60)
    print("Test 2: Helmholtz equation")
    print("="*60)
    
    k = 10.0  # Wave number
    
    # Exact solution: u = cosh(k*x) * cosh(k*y) / cosh(k)^2
    def u_exact(x, y):
        return np.cosh(k*x) * np.cosh(k*y) / (np.cosh(k)**2)
    
    # RHS (simplified)
    def f_rhs(x, y):
        return 0.0  # Simplified
    
    # Boundary condition
    def bc(x, y):
        return u_exact(x, y)
    
    # Create solver
    sx, sy = 1.0, 1.0
    nx, ny = 17, 17  # 2^4 + 1
    order = 4  # Fourth-order for Helmholtz
    
    solver = FFT9Solver(sx, sy, nx, ny, order)
    solver.set_pde_coefficients(alpha=1.0, beta=1.0, gamma=-k**2)
    solver.set_rhs_function(f_rhs)
    solver.set_boundary_conditions(bc)
    
    # Solve
    U = solver.solve()
    
    print(f"\nResults:")
    print(f"  Helmholtz equation with k={k}")
    print(f"  Grid: {nx} x {ny}")
    
    return solver.computation_time


def compare_with_paper():
    """
    Compare results with Paper 1 Table I.
    """
    print("\n" + "="*60)
    print("Comparison with Paper 1 results (Table I)")
    print("="*60)
    
    # Paper 1, Table I:
    # N = 4, 8, 16, 32, 64, 128
    # FFT9 (4th order): Max error and Time
    # FFT9 (6th order): Max error and Time
    
    paper_results = {
        'N': [4, 8, 16, 32, 64, 128],
        'FFT9_4_err': [6.82e-05, 4.27e-06, 2.68e-07, 1.68e-08, 1.04e-09, None],
        'FFT9_4_time': [0.05, 0.19, 0.76, 3.13, 13.03, None],
        'FFT9_6_err': [1.34e-07, 2.10e-09, 3.30e-11, 9.49e-13, 8.42e-13, None],
        'FFT9_6_time': [0.07, 0.28, 1.14, 4.66, 19.16, None],
    }
    
    print("\nPaper 1, Table I results (Poisson equation):")
    print("N\t4th-order Error\tTime\t6th-order Error\tTime")
    print("-" * 60)
    
    for i, N in enumerate(paper_results['N']):
        if i < 5:  # Only data for N=4,8,16,32,64
            err4 = paper_results['FFT9_4_err'][i]
            time4 = paper_results['FFT9_4_time'][i]
            err6 = paper_results['FFT9_6_err'][i]
            time6 = paper_results['FFT9_6_time'][i]
            
            print(f"{N}\t{err4:.2e}\t\t{time4:.2f}\t{err6:.2e}\t{time6:.2f}")
    
    print("\n(Note: Paper used CDC 6500, single precision)")
    print("Time comparison may not be direct due to hardware differences.")


def main():
    """Main test function"""
    print("="*60)
    print("FFT9 Algorithm - Python Implementation")
    print("Based on Houstis & Papatheodorou (1979)")
    print("="*60)
    
    # Run tests
    test_poisson_smooth()
    test_helmholtz()
    
    # Show paper results
    compare_with_paper()
    
    print("\n" + "="*60)
    print("Note: This is a simplified implementation.")
    print("For full FFT9 algorithm, see Paper 2 (Fortran code).")
    print("="*60)


if __name__ == "__main__":
    main()
