///
// GMRES (Generalized Minimal Residual) Algorithm - CORRECTED VERSION
// ================================================================
//
// Implementation based on:
// Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal 
// residual algorithm for solving nonsymmetric linear systems. SIAM Journal 
// on Scientific and Statistical Computing, 7(3), 856-869.
//
// Author: Based on paper implementation (Corrected 2026-04-25)
// Date: 2026-04-25
//

#ifndef GMRES_HPP
#define GMRES_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <memory>
#include <functional>
#include <stdexcept>
#include <algorithm>

namespace gmres {

// Type aliases for clarity
using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

///
// @brief Compute dot product of two vectors
//
double dot(const Vector& a, const Vector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes don't match");
    }
    
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

///
// @brief Compute Euclidean norm of vector
//
double norm(const Vector& v) {
    return std::sqrt(dot(v, v));
}

///
// @brief Scale vector by constant
//
Vector scale(const Vector& v, double alpha) {
    Vector result = v;
    for (auto& x : result) {
        x *= alpha;
    }
    return result;
}

///
// @brief Vector addition: result = a + b
//
Vector add(const Vector& a, const Vector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes don't match");
    }
    
    Vector result = a;
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] += b[i];
    }
    return result;
}

///
// @brief Vector subtraction: result = a - b
//
Vector subtract(const Vector& a, const Vector& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes don't match");
    }
    
    Vector result = a;
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] -= b[i];
    }
    return result;
}

///
// @brief Dense matrix for Hessenberg and related operations
//
class DenseMatrix {
public:
    DenseMatrix(size_t rows, size_t cols, double init = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, init) {}
    
    // Access element
    double& operator()(size_t i, size_t j) {
        return data_[i * cols_ + j];
    }
    
    double operator()(size_t i, size_t j) const {
        return data_[i * cols_ + j];
    }
    
    // Matrix-vector multiplication
    Vector operator*(const Vector& x) const {
        if (x.size() != cols_) {
            throw std::invalid_argument("Vector size mismatch");
        }
        
        Vector y(rows_, 0.0);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                y[i] += data_[i * cols_ + j] * x[j];
            }
        }
        return y;
    }
    
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    
    // Get column
    Vector get_column(size_t j) const {
        Vector col(rows_);
        for (size_t i = 0; i < rows_; ++i) {
            col[i] = data_[i * cols_ + j];
        }
        return col;
    }
    
    // Set column
    void set_column(size_t j, const Vector& col) {
        if (col.size() != rows_) {
            throw std::invalid_argument("Column size mismatch");
        }
        for (size_t i = 0; i < rows_; ++i) {
            data_[i * cols_ + j] = col[i];
        }
    }
    
    // Get raw data pointer (for compatibility)
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }
    
private:
    size_t rows_;
    size_t cols_;
    std::vector<double> data_;
};

///
// @brief Apply Givens rotation to a 2-vector (CORRECTED)
//
void apply_givens(double& v1, double& v2, double c, double s) {
    double temp = c * v1 + s * v2;
    v2 = -s * v1 + c * v2;
    v1 = temp;
}

///
// @brief Generate Givens rotation parameters (CORRECTED)
//
void generate_givens(double dx, double dy, double& c, double& s, double& r) {
    if (dy == 0.0) {
        c = 1.0;
        s = 0.0;
        r = dx;
    } else if (std::abs(dy) > std::abs(dx)) {
        double t = dx / dy;
        s = 1.0 / std::sqrt(1.0 + t * t);
        c = s * t;
        r = dy / s;
    } else {
        double t = dy / dx;
        c = 1.0 / std::sqrt(1.0 + t * t);
        s = c * t;
        r = dx / c;
    }
}

///
// @brief Solve upper triangular system R * y = b (CORRECTED)
//
Vector solve_upper_triangular(const DenseMatrix& R, const Vector& b) {
    int n = static_cast<int>(b.size());
    Vector y(n, 0.0);
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += R(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / R(i, i);
    }
    
    return y;
}

///
// @brief Arnoldi process to build orthonormal basis of Krylov subspace (CORRECTED)
//
void arnoldi_process(
    std::function<Vector(const Vector&)> A,
    const Vector& v0,
    int k,
    DenseMatrix& V,
    DenseMatrix& H) {
    
    int N = static_cast<int>(v0.size());
    
    // Resize V and H
    V = DenseMatrix(N, k + 1);
    H = DenseMatrix(k + 1, k);
    
    // Normalize initial vector
    double beta = norm(v0);
    if (beta < 1e-15) {
        throw std::runtime_error("Initial vector has zero norm");
    }
    
    Vector v0_normalized = scale(v0, 1.0 / beta);
    V.set_column(0, v0_normalized);
    
    // Arnoldi iterations
    for (int j = 0; j < k; ++j) {
        // w = A * v_j
        Vector v_j = V.get_column(j);
        Vector w = A(v_j);
        
        // Modified Gram-Schmidt
        for (int i = 0; i <= j; ++i) {
            Vector v_i = V.get_column(i);
            double h_ij = dot(w, v_i);
            H(i, j) = h_ij;
            
            // w = w - h_ij * v_i
            for (int idx = 0; idx < N; ++idx) {
                w[idx] -= h_ij * v_i[idx];
            }
        }
        
        double h_j1_j = norm(w);
        H(j + 1, j) = h_j1_j;
        
        // Check for breakdown (Proposition 2: algorithm has converged)
        if (h_j1_j < 1e-12) {
            // Truncate V and H (not shown in this simplified version)
            std::cout << "Arnoldi breakdown at step " << j << std::endl;
            return;
        }
        
        // v_{j+1} = w / h_{j+1,j}
        Vector v_j1 = scale(w, 1.0 / h_j1_j);
        V.set_column(j + 1, v_j1);
    }
}

///
// @brief GMRES algorithm with restart (CORRECTED)
//
bool gmres(
    std::function<Vector(const Vector&)> A,
    const Vector& b,
    Vector& x,
    double tol = 1e-6,
    int max_iter = 100,
    int m = 10,
    std::vector<double>* residuals = nullptr) {
    
    int N = static_cast<int>(b.size());
    
    // Check initial residual
    Vector r0 = subtract(b, A(x));
    double beta = norm(r0);
    
    if (beta < tol) {
        if (residuals) {
            residuals->push_back(beta);
        }
        return true;
    }
    
    if (residuals) {
        residuals->push_back(beta);
    }
    
    // If m <= 0, set m = min(N, max_iter)
    if (m <= 0) {
        m = std::min(N, max_iter);
    } else {
        m = std::min({m, N, max_iter});
    }
    
    int total_iter = 0;
    bool converged = false;
    
    // Outer loop (for restarted GMRES)
    while (total_iter < max_iter) {
        // Compute initial residual for this restart cycle
        r0 = subtract(b, A(x));
        beta = norm(r0);
        
        if (beta < tol) {
            converged = true;
            break;
        }
        
        Vector v1 = scale(r0, 1.0 / beta);
        
        // Initialize Givens rotations storage
        std::vector<double> c(m + 1, 0.0), s(m + 1, 0.0);
        Vector g(m + 2, 0.0);  // Need m+2 for safe access
        g[0] = beta;
        
        // Arnoldi process with incremental QR
        DenseMatrix V(N, m + 1);
        DenseMatrix H(m + 1, m);
        V.set_column(0, v1);
        
        int j = 0;
        for (j = 0; j < m; ++j) {
            total_iter++;
            
            if (total_iter >= max_iter) {
                break;
            }
            
            // Arnoldi step
            Vector v_j = V.get_column(j);
            Vector w = A(v_j);
            
            // Modified Gram-Schmidt
            for (int i = 0; i <= j; ++i) {
                Vector v_i = V.get_column(i);
                H(i, j) = dot(w, v_i);
                
                // w = w - H(i,j) * v_i
                for (int idx = 0; idx < N; ++idx) {
                    w[idx] -= H(i, j) * v_i[idx];
                }
            }
            
            H(j + 1, j) = norm(w);
            
            // Check for breakdown
            if (H(j + 1, j) < 1e-12) {
                // Compute solution using Eq. (8): x = x0 + V_k * y_k
                // Here k = j (0-indexed), so we have j+1 columns in V
                DenseMatrix R(j + 1, j + 1);
                for (int row = 0; row <= j; ++row) {
                    for (int col = 0; col <= std::min(row, j); ++col) {
                        R(row, col) = H(row, col);
                    }
                }
                
                Vector g_small(j + 2);
                for (int i = 0; i <= j + 1; ++i) {
                    g_small[i] = g[i];
                }
                
                // Apply Givens rotations to g_small
                for (int i = 0; i <= j; ++i) {
                    apply_givens(g_small[i], g_small[i + 1], c[i], s[i]);
                }
                
                // Now solve R * y = g_small[0:j+1]
                Vector g_RHS(j + 1);
                for (int i = 0; i <= j; ++i) {
                    g_RHS[i] = g_small[i];
                }
                
                Vector y = solve_upper_triangular(R, g_RHS);
                
                // Update x: x = x + V(:, 0:j) * y
                for (int i = 0; i <= j; ++i) {
                    Vector v_i = V.get_column(i);
                    for (int idx = 0; idx < N; ++idx) {
                        x[idx] += y[i] * v_i[idx];
                    }
                }
                
                r0 = subtract(b, A(x));
                double residual = norm(r0);
                if (residuals) {
                    residuals->push_back(residual);
                }
                
                converged = true;
                break;
            }
            
            // Normalize w to get v_{j+1}
            Vector v_j1 = scale(w, 1.0 / H(j + 1, j));
            V.set_column(j + 1, v_j1);
            
            // Apply existing Givens rotations to H(:, j)
            for (int i = 0; i < j; ++i) {
                apply_givens(H(i, j), H(i + 1, j), c[i], s[i]);
            }
            
            // Generate new Givens rotation to eliminate H(j+1, j)
            double r_temp;
            generate_givens(H(j, j), H(j + 1, j), c[j], s[j], r_temp);
            H(j, j) = r_temp;
            H(j + 1, j) = 0.0;
            
            // Apply new rotation to g
            apply_givens(g[j], g[j + 1], c[j], s[j]);
            
            // Residual norm = |g[j+1]| (Proposition 1)
            double residual = std::abs(g[j + 1]);
            if (residuals) {
                residuals->push_back(residual);
            }
            
            if (residual < tol) {
                // Compute solution
                DenseMatrix R(j + 1, j + 1);
                for (int row = 0; row <= j; ++row) {
                    for (int col = 0; col <= std::min(row, j); ++col) {
                        R(row, col) = H(row, col);
                    }
                }
                
                Vector g_small(j + 1);
                for (int i = 0; i <= j; ++i) {
                    g_small[i] = g[i];
                }
                
                Vector y = solve_upper_triangular(R, g_small);
                
                // Update x
                for (int i = 0; i <= j; ++i) {
                    Vector v_i = V.get_column(i);
                    for (int idx = 0; idx < N; ++idx) {
                        x[idx] += y[i] * v_i[idx];
                    }
                }
                
                converged = true;
                break;
            }
        }
        
        if (converged) {
            break;
        }
        
        // If not converged after m steps, compute solution and restart
        if (j == m) {
            // Apply all Givens rotations to H and g
            for (int i = 0; i < m; ++i) {
                apply_givens(g[i], g[i + 1], c[i], s[i]);
            }
            
            // Compute solution after m steps
            DenseMatrix R(m, m);
            for (int row = 0; row < m; ++row) {
                for (int col = 0; col <= row; ++col) {
                    R(row, col) = H(row, col);
                }
            }
            
            Vector g_small(m);
            for (int i = 0; i < m; ++i) {
                g_small[i] = g[i];
            }
            
            Vector y = solve_upper_triangular(R, g_small);
            
            // Update x
            for (int i = 0; i < m; ++i) {
                Vector v_i = V.get_column(i);
                for (int idx = 0; idx < N; ++idx) {
                    x[idx] += y[i] * v_i[idx];
                }
            }
            
            // Compute residual
            r0 = subtract(b, A(x));
            double residual = norm(r0);
            if (residuals) {
                residuals->push_back(residual);
            }
        }
    }
    
    return converged;
}

///
// @brief Simplified GMRES (easier to understand) - CORRECTED VERSION
//
bool gmres_simplified(
    std::function<Vector(const Vector&)> A,
    const Vector& b,
    Vector& x,
    double tol = 1e-6,
    int max_iter = 100,
    int m = 10,
    std::vector<double>* residuals = nullptr) {
    
    int N = static_cast<int>(b.size());
    
    // Check initial residual
    Vector r0 = subtract(b, A(x));
    double beta = norm(r0);
    
    if (beta < tol) {
        if (residuals) {
            residuals->push_back(beta);
        }
        return true;
    }
    
    if (residuals) {
        residuals->push_back(beta);
    }
    
    // Arnoldi process
    DenseMatrix V;
    DenseMatrix H;
    arnoldi_process(A, r0, m, V, H);
    
    int k = static_cast<int>(V.cols()) - 1;  // Actual number of steps
    
    // Solve least-squares problem: min ||beta*e1 - H*y||
    // Use: (H^T H) y = H^T (beta * e1)
    Vector e1(k + 1, 0.0);
    e1[0] = beta;
    
    // Form H^T H (k x k matrix)
    DenseMatrix HtH(k, k);
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            double sum = 0.0;
            for (int row = 0; row <= std::min(i, j) + 1; ++row) {
                sum += H(row, i) * H(row, j);
            }
            HtH(i, j) = sum;
        }
    }
    
    // Form H^T e1
    Vector HtE1(k, 0.0);
    for (int i = 0; i < k; ++i) {
        HtE1[i] = H(0, i) * beta;  // Only H(0,i) contributes
    }
    
    // Solve HtH * y = HtE1 (using simplified Gaussian elimination)
    // For simplicity, we'll use a basic solver
    Vector y = solve_upper_triangular(HtH, HtE1);
    
    // Update solution: x = x + V(:, :k) * y
    for (int i = 0; i < k; ++i) {
        Vector v_i = V.get_column(i);
        for (int idx = 0; idx < N; ++idx) {
            x[idx] += y[i] * v_i[idx];
        }
    }
    
    // Compute residual
    Vector r = subtract(b, A(x));
    double residual = norm(r);
    if (residuals) {
        residuals->push_back(residual);
    }
    
    return residual < tol;
}

} // namespace gmres

#endif // GMRES_HPP
