#!/usr/bin/env python3
"""
FFT9 6th-Order Scheme - Eigenvalue Analysis

Instead of guessing the correction term, let me compute the EXACT eigenvalue
correction needed for 6th order and find the right functional form.

For -Laplacian(u) = f on [0,1]^2 with Dirichlet BC:
  g = Laplacian(u) = -f

4th-order scheme: L_h u = R_h^(4) g
  In Fourier: lambda_L * u_hat = lambda_R * g_hat
  => u_hat = (lambda_R / lambda_L) * g_hat

Exact solution: u_hat = g_hat / (-(k^2 + l^2))
  where k = p*pi, l = q*pi for DST modes (p,q)

So the exact eigenvalue ratio is:
  lambda_R^(exact) / lambda_L = 1 / (-(k^2 + l^2))

The correction needed:
  Delta_lambda_R = lambda_L / (-(k^2+l^2)) - lambda_R^(4)
"""

import numpy as np
from scipy.fft import dst, idst


def compute_eigenvalue_corrections(N_max=32):
    """Compute exact eigenvalue corrections for different modes"""
    h = 1.0 / (N_max + 1)
    k_vals = np.arange(1, N_max + 1)
    cos_k = np.cos(k_vals * np.pi / (N_max + 1))

    print("Eigenvalue Analysis for 6th-Order Correction")
    print("=" * 80)

    # For selected modes, compute the exact correction and compare with
    # candidate functional forms
    modes = [(1,1), (1,2), (2,1), (2,2), (1,3), (3,1), (2,3), (3,3)]

    print(f"\nN = {N_max}, h = {h:.6f}")
    print(f"{'Mode':>8s}  {'Delta_R':>12s}  {'h4*lam5^2':>12s}  "
          f"{'h4*lam5*lamL':>12s}  {'h4*lamL^2':>12s}  {'Best_alpha':>12s}")

    for p, q in modes:
        # Eigenvalues
        lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[p-1] + 8.0*cos_k[q-1]
                                      + 4.0*cos_k[p-1]*cos_k[q-1])
        lam_R4 = 2.0/3.0 + (1.0/6.0)*(cos_k[p-1] + cos_k[q-1])
        lam_5 = (2.0/h**2)*(cos_k[p-1] + cos_k[q-1] - 2.0)

        # Exact wavenumbers
        kx = p * np.pi
        ky = q * np.pi
        lam_exact = -(kx**2 + ky**2)

        # Exact correction
        delta_R = lam_L / lam_exact - lam_R4

        # Candidate corrections
        c1 = h**4 * lam_5**2          # alpha * h^4 * lambda_5^2
        c2 = h**4 * lam_5 * lam_L     # alpha * h^4 * lambda_5 * lambda_L
        c3 = h**4 * lam_L**2          # alpha * h^4 * lambda_L^2

        # Best alpha for each form
        if abs(c1) > 1e-30: a1 = delta_R / c1
        else: a1 = 0
        if abs(c2) > 1e-30: a2 = delta_R / c2
        else: a2 = 0
        if abs(c3) > 1e-30: a3 = delta_R / c3
        else: a3 = 0

        print(f"  ({p},{q})  {delta_R:12.6e}  {c1:12.6e}  "
              f"{c2:12.6e}  {c3:12.6e}  {a3:12.6e}")

    return modes, cos_k, h


def search_correction_form(N=16):
    """
    Search for the correct functional form of the 6th-order correction.

    Try: lambda_R^(6) = lambda_R^(4) + alpha * correction(k,l)
    where correction is a function of the DST mode eigenvalues.
    """
    h = 1.0 / (N + 1)
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    print("\n\nSearching for correct 6th-order correction form")
    print("=" * 80)

    # Compute exact corrections for ALL modes
    alphas_by_form = {}

    # Form 1: alpha * h^4 * lambda_5^2
    # Form 2: alpha * h^4 * lambda_5 * lambda_L
    # Form 3: alpha * h^4 * lambda_L^2
    # Form 4: alpha * h^4 * (lambda_5^2 - beta * lambda_L * lambda_5)
    # Form 5: alpha * h^2 * lambda_5 (simple correction)

    a1_list, a2_list, a3_list, a5_list = [], [], [], []

    for p in range(1, N+1):
        for q in range(1, N+1):
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[p-1] + 8.0*cos_k[q-1]
                                          + 4.0*cos_k[p-1]*cos_k[q-1])
            lam_R4 = 2.0/3.0 + (1.0/6.0)*(cos_k[p-1] + cos_k[q-1])
            lam_5 = (2.0/h**2)*(cos_k[p-1] + cos_k[q-1] - 2.0)

            kx = p * np.pi
            ky = q * np.pi
            lam_exact = -(kx**2 + ky**2)

            delta_R = lam_L / lam_exact - lam_R4

            c1 = h**4 * lam_5**2
            c2 = h**4 * lam_5 * lam_L
            c3 = h**4 * lam_L**2
            c5 = h**2 * lam_5

            if abs(c1) > 1e-30: a1_list.append(delta_R / c1)
            if abs(c2) > 1e-30: a2_list.append(delta_R / c2)
            if abs(c3) > 1e-30: a3_list.append(delta_R / c3)
            if abs(c5) > 1e-30: a5_list.append(delta_R / c5)

    print(f"\nForm 1: alpha * h^4 * lambda_5^2")
    print(f"  alpha range: [{min(a1_list):.6e}, {max(a1_list):.6e}]")
    print(f"  alpha std/mean: {np.std(a1_list)/abs(np.mean(a1_list)):.4f}")
    print(f"  alpha mean: {np.mean(a1_list):.6e}")

    print(f"\nForm 2: alpha * h^4 * lambda_5 * lambda_L")
    print(f"  alpha range: [{min(a2_list):.6e}, {max(a2_list):.6e}]")
    print(f"  alpha std/mean: {np.std(a2_list)/abs(np.mean(a2_list)):.4f}")
    print(f"  alpha mean: {np.mean(a2_list):.6e}")

    print(f"\nForm 3: alpha * h^4 * lambda_L^2")
    print(f"  alpha range: [{min(a3_list):.6e}, {max(a3_list):.6e}]")
    print(f"  alpha std/mean: {np.std(a3_list)/abs(np.mean(a3_list)):.4f}")
    print(f"  alpha mean: {np.mean(a3_list):.6e}")

    print(f"\nForm 5: alpha * h^2 * lambda_5")
    print(f"  alpha range: [{min(a5_list):.6e}, {max(a5_list):.6e}]")
    print(f"  alpha std/mean: {np.std(a5_list)/abs(np.mean(a5_list)):.4f}")

    # The best form is the one with the smallest relative std
    forms = [
        ("h^4 * lambda_5^2", a1_list),
        ("h^4 * lambda_5 * lambda_L", a2_list),
        ("h^4 * lambda_L^2", a3_list),
    ]

    best_form = min(forms, key=lambda x: np.std(x[1])/abs(np.mean(x[1])))
    print(f"\nBest form: {best_form[0]}")
    print(f"  Best alpha: {np.mean(best_form[1]):.6e}")


def try_exact_eigenvalue_solver(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    Use EXACT eigenvalues in Fourier space (spectral method).
    This gives exponential convergence for smooth solutions.
    """
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    hx = sx / (n - 1)
    h = hx
    N = n - 2

    F = f_func(X, Y)
    G = -F  # g = Laplacian(u) = -f

    # Apply R_h^(4) to g (interior only)
    G_int = G[1:-1, 1:-1]
    Rg = (2.0/3.0) * G_int.copy()
    Rg[1:, :]  += (1.0/12.0) * G_int[:-1, :]
    Rg[:-1, :] += (1.0/12.0) * G_int[1:, :]
    Rg[:, 1:]  += (1.0/12.0) * G_int[:, :-1]
    Rg[:, :-1] += (1.0/12.0) * G_int[:, 1:]

    # BC correction for L_h
    bc_left = np.broadcast_to(bc_func(x[0], y), y.shape).astype(float)
    bc_right = np.broadcast_to(bc_func(x[-1], y), y.shape).astype(float)
    bc_bottom = np.broadcast_to(bc_func(x, y[0]), x.shape).astype(float)
    bc_top = np.broadcast_to(bc_func(x, y[-1]), x.shape).astype(float)

    bc_corr = np.zeros((N, N))
    for j in range(N):
        jj = j + 1
        val_l = 4.0 * bc_left[jj]
        if jj-1 >= 0: val_l += bc_left[jj-1]
        if jj+1 < n: val_l += bc_left[jj+1]
        bc_corr[0, j] -= (1.0/(6.0*h**2)) * val_l
        val_r = 4.0 * bc_right[jj]
        if jj-1 >= 0: val_r += bc_right[jj-1]
        if jj+1 < n: val_r += bc_right[jj+1]
        bc_corr[N-1, j] -= (1.0/(6.0*h**2)) * val_r
    for i in range(N):
        ii = i + 1
        val_b = 4.0 * bc_bottom[ii]
        if ii-1 >= 0: val_b += bc_bottom[ii-1]
        if ii+1 < n: val_b += bc_bottom[ii+1]
        bc_corr[i, 0] -= (1.0/(6.0*h**2)) * val_b
        val_t = 4.0 * bc_top[ii]
        if ii-1 >= 0: val_t += bc_top[ii-1]
        if ii+1 < n: val_t += bc_top[ii+1]
        bc_corr[i, N-1] -= (1.0/(6.0*h**2)) * val_t

    Rg += bc_corr

    # 2D DST
    Rg_hat = np.zeros((N, N))
    for i in range(N):
        Rg_hat[i, :] = dst(Rg[i, :], type=1, norm='ortho')
    for j in range(N):
        Rg_hat[:, j] = dst(Rg_hat[:, j], type=1, norm='ortho')

    # Solve using EXACT eigenvalues
    k_vals = np.arange(1, N + 1)
    cos_k = np.cos(k_vals * np.pi / (N + 1))

    U_hat = np.zeros((N, N))
    for ki in range(N):
        for li in range(N):
            lam_L = (1.0/(6.0*h**2)) * (-20.0 + 8.0*cos_k[ki] + 8.0*cos_k[li]
                                          + 4.0*cos_k[ki]*cos_k[li])
            # Use exact eigenvalue: lambda_exact = -((ki*pi)^2 + (li*pi)^2)
            # But we need lambda_R^(exact) = lambda_L / lambda_exact
            # And u_hat = lambda_R^(exact) / lambda_L * g_hat
            #           = g_hat / lambda_exact
            # But we have R_h * g, not g. So:
            # u_hat = Rg_hat / lambda_L  (4th order)
            # For exact: u_hat = g_hat / lambda_exact
            # But Rg_hat = lambda_R * g_hat (in Fourier)
            # So g_hat = Rg_hat / lambda_R
            # And u_hat = Rg_hat / (lambda_R * lambda_exact)
            #           = Rg_hat / (lambda_R * (-(kx^2 + ky^2)))

            # Compute 4th-order R_h eigenvalue
            lam_R4 = 2.0/3.0 + (1.0/6.0)*(cos_k[ki] + cos_k[li])

            # Exact continuous eigenvalue
            kx = (ki + 1) * np.pi / sx  # Wait, for DST type 1, the wavenumber is k*pi/(N+1)?
            # Actually, for DST-I on [0,1] with N interior points:
            # The eigenvalues of the 5-point Laplacian are (2/h^2)(cos(p*pi/(N+1)) - 1)
            # The corresponding continuous wavenumber is k_p = p*pi
            kx_cont = (ki + 1) * np.pi  # Hmm, this doesn't seem right

            # Actually, the DST-I eigenvectors are sin(p*j*pi/(N+1))
            # These correspond to the function sin(p*pi*x) on [0,1]
            # So the continuous wavenumber is p*pi

            # But wait, for the DST type 1 with norm='ortho':
            # The transform pair is defined with specific normalization

            # The key point: the DST mode with index p corresponds to
            # the continuous function sin(p*pi*x)
            # So k_x = p*pi, k_y = q*pi
            # lambda_exact = -(p^2*pi^2 + q^2*pi^2) = -(p^2+q^2)*pi^2

            lam_exact = -((ki+1)**2 + (li+1)**2) * np.pi**2

            if abs(lam_R4) > 1e-14 and abs(lam_exact) > 1e-14:
                # 6th order: use exact eigenvalue ratio
                # u_hat = Rg_hat * lam_exact_ratio / lam_L
                # where lam_exact_ratio = lam_R4 * lam_exact / lam_L ... no

                # Actually, the simplest approach:
                # Rg_hat = lambda_R4 * g_hat (since Rg = R_h * g)
                # u_hat = g_hat / lambda_exact = Rg_hat / (lambda_R4 * lambda_exact)
                # But lambda_exact is for -Laplacian, so:
                # -Laplacian(u) = f => u_hat = f_hat / (k^2+l^2)
                # g = -f => g_hat = -f_hat
                # u_hat = -g_hat / (k^2+l^2) = g_hat / (-(k^2+l^2)) = g_hat / lambda_exact
                # Rg_hat = lambda_R4 * g_hat
                # u_hat = Rg_hat / (lambda_R4 * lambda_exact)
                # But lambda_exact is NEGATIVE (it's -(k^2+l^2))
                # So u_hat = Rg_hat / (lambda_R4 * (-(k^2+l^2)))

                U_hat[ki, li] = Rg_hat[ki, li] / (lam_R4 * lam_exact)

    # Inverse 2D DST
    U_int = np.zeros((N, N))
    for ki in range(N):
        U_int[ki, :] = idst(U_hat[ki, :], type=1, norm='ortho')
    for li in range(N):
        U_int[:, li] = idst(U_int[:, li], type=1, norm='ortho')

    U = np.zeros((n, n))
    U[0, :] = bc_left; U[-1, :] = bc_right
    U[:, 0] = bc_bottom; U[:, -1] = bc_top
    U[1:-1, 1:-1] = U_int
    return U


def test_exact_eigenvalue():
    """Test the exact eigenvalue approach"""
    print("\n" + "=" * 80)
    print("Exact Eigenvalue (Spectral) Approach")
    print("=" * 80)

    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bc(x, y): return 0.0

    ns = [9, 17, 33, 65, 129]

    print("\n--- Exact Eigenvalue Method ---")
    errors = []
    for n in ns:
        U = try_exact_eigenvalue_solver(n, f_rhs, bc)
        x = np.linspace(0, 1, n); X, Y = np.meshgrid(x, x, indexing='ij')
        err = np.max(np.abs(U - u_exact(X, Y)))
        errors.append(err)
        print(f"  n={n:3d}: error = {err:.6e}")

    for i in range(1, len(ns)):
        rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
        print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    compute_eigenvalue_corrections(16)
    search_correction_form(16)
    test_exact_eigenvalue()
