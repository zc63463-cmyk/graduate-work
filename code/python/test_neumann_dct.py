#!/usr/bin/env python3
"""Test Neumann solver with symmetrized ghost-point + DCT-I approach."""
import numpy as np
from scipy.fft import dct, idct

def solve_neumann_1d(n, f_func, k2=0.0):
    h = 1.0 / (n - 1)
    N = n
    x = np.linspace(0, 1, n)
    f = f_func(x)
    
    d_scale = np.ones(N)
    d_scale[0] = np.sqrt(2)
    d_scale[-1] = np.sqrt(2)
    
    b_tilde = h**2 * f / d_scale
    b_hat = dct(b_tilde, type=1, norm='ortho')
    
    k_vals = np.arange(N)
    lam = 2 * (1 - np.cos(np.pi * k_vals / (N - 1)))
    
    denom = lam + h**2 * k2
    u_hat = np.zeros(N)
    mask = np.abs(denom) > 1e-14
    u_hat[mask] = b_hat[mask] / denom[mask]
    
    u_tilde = idct(u_hat, type=1, norm='ortho')
    u = u_tilde * d_scale
    return u

def solve_neumann_2d_fa(n, f_func, k2=0.0):
    h = 1.0 / (n - 1)
    N = n
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    F = f_func(X, Y)
    
    d_scale = np.ones(N)
    d_scale[0] = np.sqrt(2)
    d_scale[-1] = np.sqrt(2)
    
    F_tilde = h**2 * F / d_scale[:, None] / d_scale[None, :]
    
    Fh = np.zeros((N, N))
    for i in range(N):
        Fh[i, :] = dct(F_tilde[i, :], type=1, norm='ortho')
    for j in range(N):
        Fh[:, j] = dct(Fh[:, j], type=1, norm='ortho')
    
    k_vals = np.arange(N)
    lam_1d = 2 * (1 - np.cos(np.pi * k_vals / (N - 1)))
    LXX, LYY = np.meshgrid(lam_1d, lam_1d, indexing='ij')
    denom = LXX + LYY + h**2 * k2
    
    Uh = np.zeros((N, N))
    mask = np.abs(denom) > 1e-14
    Uh[mask] = Fh[mask] / denom[mask]
    
    Ui = np.zeros((N, N))
    for i in range(N):
        Ui[i, :] = idct(Uh[i, :], type=1, norm='ortho')
    for j in range(N):
        Ui[:, j] = idct(Ui[:, j], type=1, norm='ortho')
    
    U = Ui * d_scale[:, None] * d_scale[None, :]
    return U

# Test 1D
k2 = 10.0
def u1(x): return np.cos(np.pi * x)
def f1(x): return (np.pi**2 + k2) * np.cos(np.pi * x)

print('1D Neumann Helmholtz (k^2=10) convergence:')
prev_err = None
for n in [9, 17, 33, 65, 129, 257]:
    x = np.linspace(0, 1, n)
    u = solve_neumann_1d(n, f1, k2)
    err = np.max(np.abs(u - u1(x)))
    rate = np.log2(prev_err/err) if prev_err is not None else None
    if rate is not None:
        print(f'  n={n:5d}  error={err:.6e}  rate={rate:.2f}')
    else:
        print(f'  n={n:5d}  error={err:.6e}')
    prev_err = err

# Test 2D
print()
print('2D Neumann Helmholtz (k^2=10) convergence:')

def u2(x, y): return np.cos(np.pi * x) * np.cos(np.pi * y)
def f2(x, y): return (2 * np.pi**2 + k2) * np.cos(np.pi * x) * np.cos(np.pi * y)

prev_err = None
for n in [9, 17, 33, 65, 129]:
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    U = solve_neumann_2d_fa(n, f2, k2)
    err = np.max(np.abs(U - u2(X, Y)))
    rate = np.log2(prev_err/err) if prev_err is not None else None
    if rate is not None:
        print(f'  n={n:5d}  error={err:.6e}  rate={rate:.2f}')
    else:
        print(f'  n={n:5d}  error={err:.6e}')
    prev_err = err

# Test Poisson (k^2=0)
print()
print('2D Neumann Poisson (k^2=0) convergence:')

def u2p(x, y): return np.cos(np.pi * x) * np.cos(np.pi * y)
def f2p(x, y): return 2 * np.pi**2 * np.cos(np.pi * x) * np.cos(np.pi * y)

prev_err = None
for n in [9, 17, 33, 65, 129]:
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    U = solve_neumann_2d_fa(n, f2p, k2=0.0)
    err = np.max(np.abs(U - u2p(X, Y)))
    rate = np.log2(prev_err/err) if prev_err is not None else None
    if rate is not None:
        print(f'  n={n:5d}  error={err:.6e}  rate={rate:.2f}')
    else:
        print(f'  n={n:5d}  error={err:.6e}')
    prev_err = err
