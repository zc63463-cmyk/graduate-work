///
// GMRES Algorithm - Example Usage
// ====================================
//
// This file demonstrates how to use the GMRES implementation
// provided in gmres.hpp.
//
// Compile with:
//   g++ -std=c++11 -o gmres_example main.cpp -I.
//
// Author: Based on Saad & Schultz (1986)
// Date: 2026-04-25
///

#include <iostream>
#include <vector>
#include <cmath>
#include "gmres.hpp"

using namespace gmres;

///
// @brief Create a 2x2 test matrix (from paper Section 2)
//
// A = [1 1; 1 0]
// This is the example where GCR crashes but GMRES works.
///
DenseMatrix create_2x2_test() {
    DenseMatrix A(2, 2);
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    A(1, 1) = 0.0;
    return A;
}

///
// @brief Create a simple tridiagonal matrix (1D Poisson)
///
DenseMatrix create_poisson_matrix(int N) {
    DenseMatrix A(N, N);
    
    for (int i = 0; i < N; ++i) {
        A(i, i) = 2.0;
        if (i > 0) {
            A(i, i-1) = -1.0;
        }
        if (i < N-1) {
            A(i, i+1) = -1.0;
        }
    }
    
    return A;
}

///
// @brief Print a vector
///
void print_vector(const std::string& name, const Vector& v) {
    std::cout << name << " = [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i];
        if (i < v.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

///
// @brief Example 1: 2x2 system from paper
///
void example_1() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Example 1: 2x2 system from paper (Section 2)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Create matrix A
    DenseMatrix A = create_2x2_test();
    std::cout << "\nMatrix A:" << std::endl;
    std::cout << "[1, 1; 1, 0]" << std::endl;
    
    // Right-hand side
    Vector b = {1.0, 0.0};
    print_vector("b", b);
    
    // Initial guess
    Vector x = {0.0, 0.0};
    
    // Define matrix-vector product function
    auto A_op = [&A](const Vector& v) -> Vector {
        return A * v;
    };
    
    // Solve with GMRES (no restart, m=2)
    std::cout << "\nSolving with GMRES (m=2)..." << std::endl;
    std::vector<double> residuals;
    bool converged = gmres::gmres(A_op, b, x, 1e-10, 100, 2, &residuals);
    
    print_vector("\nSolution x", x);
    
    // Exact solution
    Vector x_exact = {0.0, 1.0};
    print_vector("Exact solution", x_exact);
    
    // Compute residual
    Vector r = gmres::subtract(b, A_op(x));
    double residual_norm = gmres::norm(r);
    std::cout << "Residual norm: " << residual_norm << std::endl;
    std::cout << "Converged: " << (converged ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << residuals.size() << std::endl;
}

///
// @brief Example 2: 1D Poisson equation
///
void example_2() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Example 2: 1D Poisson equation (-u'' = f)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int N = 50;
    std::cout << "\nMatrix size: " << N << "x" << N << std::endl;
    
    // Create Poisson matrix
    DenseMatrix A = create_poisson_matrix(N);
    
    // Right-hand side: f(x) = sin(pi * x)
    Vector b(N, 0.0);
    double h = 1.0 / (N + 1);
    for (int i = 0; i < N; ++i) {
        double x = (i + 1) * h;
        b[i] = std::sin(M_PI * x);
    }
    
    // Initial guess
    Vector x(N, 0.0);
    
    // Define matrix-vector product
    auto A_op = [&A](const Vector& v) -> Vector {
        return A * v;
    };
    
    // Solve with GMRES (restart m=10)
    std::cout << "\nSolving with GMRES (m=10)..." << std::endl;
    std::vector<double> residuals;
    bool converged = gmres::gmres(A_op, b, x, 1e-8, 500, 10, &residuals);
    
    std::cout << "\nConverged: " << (converged ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << residuals.size() << std::endl;
    std::cout << "Final residual: " << residuals.back() << std::endl;
    
    // Print first few values of solution
    std::cout << "\nFirst 10 values of solution:" << std::endl;
    for (int i = 0; i < std::min(10, N); ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    
    // Print convergence history
    std::cout << "\nConvergence history (log scale):" << std::endl;
    for (size_t i = 0; i < residuals.size(); ++i) {
        std::cout << "Iteration " << i << ": residual = " << residuals[i] << std::endl;
    }
}

///
// @brief Example 3: Non-symmetric system
///
void example_3() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Example 3: Non-symmetric system" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int N = 100;
    std::cout << "\nMatrix size: " << N << "x" << N << std::endl;
    
    // Create a non-symmetric matrix (random + diagonal dominance)
    DenseMatrix A(N, N);
    std::srand(42);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A(i, j) = static_cast<double>(std::rand()) / RAND_MAX;
        }
        A(i, i) += 10.0;  // Make diagonally dominant
    }
    
    // Make it non-symmetric
    std::cout << "Matrix symmetry check: ";
    double sym_error = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            sym_error += std::abs(A(i, j) - A(j, i));
        }
    }
    std::cout << sym_error << std::endl;
    
    // Right-hand side
    Vector b(N, 1.0);
    
    // Initial guess
    Vector x(N, 0.0);
    
    // Define matrix-vector product
    auto A_op = [&A](const Vector& v) -> Vector {
        return A * v;
    };
    
    // Solve with GMRES (no restart)
    std::cout << "\nSolving with GMRES (no restart)..." << std::endl;
    std::vector<double> residuals;
    bool converged = gmres::gmres(A_op, b, x, 1e-8, 200, 0, &residuals);
    
    std::cout << "\nConverged: " << (converged ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << residuals.size() << std::endl;
    if (!residuals.empty()) {
        std::cout << "Final residual: " << residuals.back() << std::endl;
    }
}

///
// @brief Helper function to save convergence history to file
///
void save_convergence_history(const std::vector<double>& residuals, 
                             const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    file << "# GMRES Convergence History" << std::endl;
    file << "# Iteration\tResidual" << std::endl;
    for (size_t i = 0; i < residuals.size(); ++i) {
        file << i << "\t" << residuals[i] << std::endl;
    }
    
    file.close();
    std::cout << "Convergence history saved to " << filename << std::endl;
}

int main() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "GMRES Algorithm Implementation" << std::endl;
    std::cout << "Based on Saad & Schultz (1986)" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    try {
        // Run examples
        example_1();
        example_2();
        // example_3();  // Uncomment to run
        
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "All examples completed successfully!" << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
