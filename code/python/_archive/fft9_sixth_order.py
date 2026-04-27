#!/usr/bin/env python3
"""
FFT9 Sixth-Order Scheme - Theoretical Derivation & Implementation

For the Poisson equation: Lap(u) = g, or equivalently -Lap(u) = f where g = -f

The 4th-order 9-point compact scheme:
  L_h * u = R_h^(4) * g

  L_h = (1/(6h^2)) * [1  4  1; 4 -20  4; 1  4  1]
  R_h^(4) = I + (h^2/12)*Lap_h = [0  1/12  0; 1/12  2/3  1/12; 0  1/12  0]

For 6th order, we need R_h^(6) that includes h^4 correction terms.

The key insight: The 4th-order scheme has truncation error O(h^4).
For 6th order, we need to cancel the O(h^4) terms.

By Taylor expansion analysis:
  L_h * u = (Lap u) + (h^2/12)(Lap^2 u) + (h^4/360)(Lap^3 u) + O(h^6)
  R_h^(4) * g = g + (h^2/12)*Lap(g) = g + (h^2/12)*(Lap^2 u) 

So L_h * u - R_h^(4) * g = (h^4/360)*(Lap^3 u) + O(h^6)

To cancel the h^4 term, we need:
  R_h^(6) = R_h^(4) + (h^4/360)*Lap_h^2

where Lap_h^2 is the 9-point stencil approximation of Lap^2.

The Lap_h^2 stencil for a function f:
  Lap_h f = (1/h^2) * [0 1 0; 1 -4 1; 0 1 0] * f
  Lap_h^2 f = Lap_h(Lap_h f)

The 9-point stencil for Lap^2:
  (1/h^4) * [0  0   1   0  0;
              0  2  -8   2  0;
              1 -8  20  -8  1;
              0  2  -8   2  0;
              0  0   1   0  0]

Wait - this is a 5x5 stencil, not 3x3. For a compact 3x3 scheme, we need
a different approach.

Actually, let me reconsider. The paper's formula:
  (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] * u = R_h * f

For 6th order, R_h must be chosen so that:
  L_h * u - R_h * g = O(h^6)

We know:
  L_h * u = Lap(u) + (h^2/12)*Lap^2(u) + (h^4/360)*Lap^3(u) + O(h^6)

For R_h^(4) * g:
  R_h^(4) * g = g + (h^2/12)*Lap(g) + O(h^4)

So L_h*u - R_h^(4)*g = (h^4/360)*Lap^3(u) + O(h^6)

To get 6th order, we need R_h that satisfies:
  R_h * g = g + (h^2/12)*Lap(g) + (h^4/360)*Lap^2(g) + O(h^6)

This means R_h = I + (h^2/12)*Lap_h + (h^4/360)*Lap_h^2

But Lap_h^2 requires a 5x5 stencil... unless we use the 9-point Lap_h
instead of the 5-point Lap_h.

The 9-point Lap_h is L_h (scaled): (6h^2)*L_h = [1 4 1; 4 -20 4; 1 4 1]

So Lap_approx = 6*L_h (unnormalized). Then:
  Lap_h^2 ≈ (6*L_h)^2

But this gives a 5x5 stencil again. For a compact 3x3 method,
we need to find a 3x3 R_h that approximates the action of 
I + (h^2/12)*Lap + (h^4/360)*Lap^2 on g.

Alternative approach: Use the relationship between L_h and Lap_h.
Since L_h = Lap_h + (h^2/12)*Lap_h^2 + O(h^4), we can write:
  (h^2/12)*Lap_h^2 ≈ L_h - Lap_h

But this doesn't directly help for a 3x3 stencil.

Let me try a completely different approach: just use the general form
R_h = [a  b  a; b  c  b; a  b  a] and determine a, b, c by requiring
6th-order accuracy.

For u = sin(kx)*sin(ly), the eigenvalue of R_h acting on g:
  eig_R = c + 2b*(cos(kh) + cos(lh)) + 2a*cos(kh)*cos(lh) + 2a

Wait, let me be more careful. R_h has the form:
  [a  b  a; b  c  b; a  b  a]

For a mode u_{ij} = sin(ikθ)*sin(jlφ), the eigenvalue of R_h is:
  λ_R = c + 2b*(cos(kθ) + cos(lφ)) + 4a*cos(kθ)*cos(lφ)

Hmm wait, that's not right either. Let me think about this more carefully.

R_h applied to g at point (i,j):
  (R_h*g)_{i,j} = a*g_{i-1,j-1} + b*g_{i-1,j} + a*g_{i-1,j+1}
                 + b*g_{i,j-1} + c*g_{i,j} + b*g_{i,j+1}
                 + a*g_{i+1,j-1} + b*g_{i+1,j} + a*g_{i+1,j+1}

For g_{ij} = exp(i(kx_i + ly_j)):
  λ_R = a*exp(-i(kh+lh)) + b*exp(-ikh) + a*exp(-i(kh-lh))
       + b*exp(-ilh) + c + b*exp(+ilh)
       + a*exp(+i(kh-lh)) + b*exp(+ikh) + a*exp(+i(kh+lh))

  = c + 2b*(cos(kh) + cos(lh)) + 2a*(cos(kh+lh) + cos(kh-lh))
  = c + 2b*(cos(kh) + cos(lh)) + 4a*cos(kh)*cos(lh)

Similarly, L_h eigenvalue:
  λ_L = (-20 + 8*(cos(kh) + cos(lh)) + 4*cos(kh)*cos(lh)) / (6h^2)

For the equation L_h * u = R_h * g:
  u = λ_R/λ_L * g

For exact solution: u = g/λ(Lap) where λ(Lap) = -(k^2 + l^2)

So we need: λ_R/λ_L = 1/(-(k^2+l^2)) for all k, l up to O(h^6)

This gives us: λ_R * (-(k^2+l^2)) = λ_L for all k, l

Let me expand in powers of h (with k and l fixed):
  cos(kh) = 1 - (kh)^2/2 + (kh)^4/24 - (kh)^6/720 + ...
  cos(lh) = 1 - (lh)^2/2 + (lh)^4/24 - (lh)^6/720 + ...

  cos(kh)*cos(lh) = 1 - (k^2+l^2)h^2/2 + ((k^4+l^4)h^4/24 + k^2l^2h^4/4) - ...

Let s = k^2+l^2, p = k^4+l^4, q = k^2l^2

  cos(kh) + cos(lh) = 2 - sh^2/2 + ph^4/24 - ...
  cos(kh)*cos(lh) = 1 - sh^2/2 + (ph/24 + q/4 - k^2l^2/4)h^4 ...

Hmm, this is getting complicated. Let me just do it numerically.

For 6th order, we need 3 parameters (a, b, c) to satisfy conditions
at multiple (k, l) values. With 3 unknowns, we can satisfy up to 3 
conditions. But 4th order already uses all 3 to get O(h^4).

For 6th order with the SAME L_h, we can't do better than 4th order
with a 3x3 R_h! The 6th order scheme must use a DIFFERENT L_h as well,
or use additional terms.

Actually, looking at this more carefully, I think the paper's "6th order"
scheme for Poisson might use:
  L_h * u = R_h^(6) * f
where f = g (the RHS of Lap(u) = f), and R_h^(6) involves h^4 terms
that also require computing Lap(f) or similar.

This means the 6th-order scheme modifies the RHS but keeps L_h the same!
The RHS correction involves computing Lap(f) using the same grid,
which adds more stencil points.

So the paper's formula:
  L_h * u = (1/360)*[0 48 0; 48 0 48; 0 48 0] * f

This has NO center term and only 4 off-center terms. Let's check:

R_h = (1/360)*[0 48 0; 48 0 48; 0 48 0] = [0 2/15 0; 2/15 0 2/15; 0 2/15 0]

λ_R = 0 + 2*(2/15)*(cos(kh) + cos(lh)) + 4*0*cos(kh)*cos(lh)
    = (4/15)*(cos(kh) + cos(lh))

For the equation to be 6th order, we need:
  (4/15)*(cos(kh) + cos(lh)) * (-(k^2+l^2)) ≈ λ_L

Let's check: λ_L = (-20 + 8*(cos(kh)+cos(lh)) + 4*cos(kh)*cos(lh))/(6h^2)

If R_h * (-s) = λ_L, then:
  -(4/15)*s*(cos(kh) + cos(lh)) = (-20 + 8*(cos(kh)+cos(lh)) + 4*cos(kh)*cos(lh))/(6h^2)

This can't be right because the LHS depends on k,l but RHS has 1/h^2.

I think the issue is that the paper's formula is for Lap(u) = f (NOT -Lap(u) = f),
and the coefficients in the paper already include the 1/(6h^2) factor on the left.

So the full equation is:
  (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] * u = (1/360) * [0 48 0; 48 0 48; 0 48 0] * f

This means: L_h * u = R_h * f  where f = g = Lap(u)

So: λ_L * u_hat = λ_R * f_hat = λ_R * (-s) * u_hat
=> λ_L = -s * λ_R
=> λ_R = -λ_L / s

Let's verify for the 4th order scheme:
  λ_R^(4) = 2/3 + (2/12)*(cos(kh)+cos(lh))
          = 2/3 + (1/6)*(cos(kh)+cos(lh))

  -s * λ_R^(4) should equal λ_L up to O(h^4)

For k=l=1, h small:
  λ_L ≈ (-20 + 8*(2-sh^2/2+...) + 4*(1-sh^2/2+...))/(6h^2)
      = (-20 + 16 + 4 + O(h^2))/(6h^2)
      = 0/(6h^2) ... 

Wait, that's zero for kh→0! That's the zero eigenvalue for the constant mode.

Let me use specific values. k=1, l=1, h=0.1, s=2:
  cos(0.1) = 0.99500
  cos(0.1)*cos(0.1) = 0.99002

  λ_L = (-20 + 8*1.99001 + 4*0.99002)/(6*0.01)
      = (-20 + 15.92008 + 3.96008)/0.06
      = -0.11984/0.06 = -1.99733

  Exact: -s = -2

  λ_R^(4) = 2/3 + (1/6)*1.99001 = 0.66667 + 0.33167 = 0.99834

  -s * λ_R^(4) = -2 * 0.99834 = -1.99668

  Error: -1.99733 - (-1.99668) = -0.00065, which is O(h^2) ≈ 0.01

Hmm, the error is O(h^2) not O(h^4). Let me reconsider...

Actually wait - for the 4th order scheme, the LOCAL truncation error is O(h^4),
but the GLOBAL error is also O(h^4) (because the scheme is consistent and stable).

The relationship λ_L = -s * λ_R holds up to the truncation error level.
Let me be more precise.

For u = sin(πx)*sin(πy) on [0,1]^2 with N interior points:
  k = π, l = π, h = 1/(N+1)
  θ = kh = π/(N+1)

  λ_L(θ) = (-20 + 16cos(θ) + 4cos²(θ))/(6h²)

  For the 4th order R_h:
  λ_R(θ) = 2/3 + (1/3)*cos(θ)  [since cos(kh)=cos(lh)=cos(θ)]

  We need: λ_L(θ) / λ_R(θ) = -(k²+l²) = -2π²/h₀² where h₀=1

  Actually no, k and l are the continuous wavenumbers, and the grid 
  wavenumber is θ = kh.

Let me just do numerical verification with a matrix approach.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.fft import dst, idst


def build_Lh_matrix(N, h):
    """Build L_h = (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] as sparse matrix"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # Center: -20/(6h^2)
            rows.append(idx); cols.append(idx); vals.append(-20.0/(6.0*h**2))
            # Edge neighbors: 4/(6h^2)
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(4.0/(6.0*h**2))
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(4.0/(6.0*h**2))
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(4.0/(6.0*h**2))
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(4.0/(6.0*h**2))
            # Corner neighbors: 1/(6h^2)
            if i>0 and j>0: rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j>0: rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i>0 and j<N-1: rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j<N-1: rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def apply_Rh_4th(g, N, h):
    """Apply 4th-order R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0] to g (N x N grid)"""
    G = g.reshape((N, N))
    Rg = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            Rg[i,j] = (2.0/3.0) * G[i,j]
            if i > 0: Rg[i,j] += (1.0/12.0) * G[i-1,j]
            if i < N-1: Rg[i,j] += (1.0/12.0) * G[i+1,j]
            if j > 0: Rg[i,j] += (1.0/12.0) * G[i,j-1]
            if j < N-1: Rg[i,j] += (1.0/12.0) * G[i,j+1]
    return Rg.flatten()


def apply_Rh_6th_candidate(g, N, h):
    """
    Apply a candidate 6th-order R_h to g.
    
    For Lap(u) = g, the 6th order scheme needs:
    R_h * g = g + (h^2/12)*Lap(g) + (h^4/360)*Lap^2(g) + O(h^6)
    
    We compute Lap(g) using 5-point stencil, then Lap^2(g) using
    Lap(Lap(g)).
    """
    G = g.reshape((N, N))
    # Pad with zeros for boundary
    G_pad = np.zeros((N+2, N+2))
    G_pad[1:-1, 1:-1] = G
    
    # Lap(g) using 5-point stencil
    Lap_g = (G_pad[2:, 1:-1] + G_pad[:-2, 1:-1] + 
             G_pad[1:-1, 2:] + G_pad[1:-1, :-2] - 4*G) / h**2
    
    # Lap^2(g) = Lap(Lap(g))
    Lap_g_pad = np.zeros((N+2, N+2))
    Lap_g_pad[1:-1, 1:-1] = Lap_g
    Lap2_g = (Lap_g_pad[2:, 1:-1] + Lap_g_pad[:-2, 1:-1] + 
              Lap_g_pad[1:-1, 2:] + Lap_g_pad[1:-1, :-2] - 4*Lap_g) / h**2
    
    Rg = G + (h**2/12.0)*Lap_g + (h**4/360.0)*Lap2_g
    return Rg.flatten()


def apply_Rh_6th_9pt(g, N, h):
    """
    Apply 6th-order R_h using 9-point Laplacian for the correction terms.
    
    R_h * g = g + (h^2/12)*Lap_9pt(g) + (h^4/360)*Lap_9pt(Lap_9pt(g))
    
    where Lap_9pt uses the L_h operator (6h^2 * L_h).
    """
    G = g.reshape((N, N))
    G_pad = np.zeros((N+2, N+2))
    G_pad[1:-1, 1:-1] = G
    
    # 9-point Laplacian of g
    Lap9_g = (1.0/(6.0*h**2)) * (
        -20*G +
        4*(G_pad[2:,1:-1] + G_pad[:-2,1:-1] + G_pad[1:-1,2:] + G_pad[1:-1,:-2]) +
        1*(G_pad[2:,2:] + G_pad[2:,:-2] + G_pad[:-2,2:] + G_pad[:-2,:-2])
    )
    
    # 9-point Laplacian of Lap9_g
    Lap9_g_pad = np.zeros((N+2, N+2))
    Lap9_g_pad[1:-1, 1:-1] = Lap9_g
    Lap9_2g = (1.0/(6.0*h**2)) * (
        -20*Lap9_g +
        4*(Lap9_g_pad[2:,1:-1] + Lap9_g_pad[:-2,1:-1] + Lap9_g_pad[1:-1,2:] + Lap9_g_pad[1:-1,:-2]) +
        1*(Lap9_g_pad[2:,2:] + Lap9_g_pad[2:,:-2] + Lap9_g_pad[:-2,2:] + Lap9_g_pad[:-2,:-2])
    )
    
    Rg = G + (h**2/12.0)*Lap9_g + (h**4/360.0)*Lap9_2g
    return Rg.flatten()


def apply_Rh_6th_compact(g, N, h):
    """
    6th-order R_h using a compact 3x3 stencil derived from
    combining the 4th-order R_h with h^4 correction.
    
    From Taylor expansion:
    R_h^(6) = I + (h^2/12)*Lap_h + (h^4/360)*Lap_h^2
    
    But Lap_h^2 can be approximated using the 9-point Laplacian
    applied twice. However, for a compact stencil, we note that
    for Lap(u) = g:
    
    Lap^2(g) = Lap(Lap(u)) = Lap^3(u)
    
    And we can approximate Lap^3 using L_h^3 / (6h^2)^3 which
    would give a 7x7 stencil - not compact.
    
    Alternative: Use the Mehrstellenverfahren (Hermitian) approach.
    The 6th-order scheme for -Lap(u) = f on a uniform grid is:
    
    L_h * u = R_h * f
    
    where R_h has the form:
    R_h = alpha*I + beta*Lap_h (5pt) + gamma*Lap_h (9pt)
    
    For 6th order with h_x = h_y = h:
    R_h = [a  b  a; b  c  b; a  b  a]
    
    We need 3 conditions from matching the Fourier symbol:
    lambda_R(k,l) / lambda_L(k,l) = 1/(k^2+l^2) + O(h^6)
    
    i.e., -(k^2+l^2) * lambda_R = lambda_L + O(h^8)
    
    Expanding both sides in powers of h for general (k,l):
    """
    # Use the approach: R_h = c1*I + c2*Lap5_h + c3*Lap9_h
    # where Lap5_h is 5-point Laplacian and Lap9_h is 9-point Laplacian
    # These three operators are linearly independent on the 3x3 stencil
    
    G = g.reshape((N, N))
    G_pad = np.zeros((N+2, N+2))
    G_pad[1:-1, 1:-1] = G
    
    # 5-point Laplacian of g (stencil: [0 1 0; 1 -4 1; 0 1 0]/h^2)
    Lap5_g = (G_pad[2:,1:-1] + G_pad[:-2,1:-1] + 
              G_pad[1:-1,2:] + G_pad[1:-1,:-2] - 4*G) / h**2
    
    # 9-point Laplacian of g (stencil: [1 4 1; 4 -20 4; 1 4 1]/(6h^2))
    Lap9_g = (1.0/(6.0*h**2)) * (
        -20*G +
        4*(G_pad[2:,1:-1] + G_pad[:-2,1:-1] + G_pad[1:-1,2:] + G_pad[1:-1,:-2]) +
        1*(G_pad[2:,2:] + G_pad[2:,:-2] + G_pad[:-2,2:] + G_pad[:-2,:-2])
    )
    
    # For 4th order: R_h = I + (h^2/12)*Lap5_h  (already gives 4th order)
    # For 6th order, we need additional terms
    # R_h^(6) = I + (h^2/12)*Lap5_h + c*Lap9_h  
    # But this overparameterizes...
    
    # Let me try: R_h^(6) = alpha*I + beta*h^2*Lap5_h + gamma*h^2*Lap9_h
    # With alpha, beta, gamma chosen for 6th order
    
    # For now, use the sequential correction approach
    Rg = G + (h**2/12.0)*Lap5_g + (h**4/360.0)*(4.0*Lap5_g + Lap9_g)/(5.0) * h**(-2)
    
    # This is getting complicated. Let me just try the simplest approach:
    # R_h^(6) = R_h^(4) + (h^4/360)*Lap_h(Lap_h(g))
    # where Lap_h is the 5-point stencil
    
    Lap5_g_pad = np.zeros((N+2, N+2))
    Lap5_g_pad[1:-1, 1:-1] = Lap5_g
    Lap5_2g = (Lap5_g_pad[2:,1:-1] + Lap5_g_pad[:-2,1:-1] + 
               Lap5_g_pad[1:-1,2:] + Lap5_g_pad[1:-1,:-2] - 4*Lap5_g) / h**2
    
    Rg = G + (h**2/12.0)*Lap5_g + (h**4/360.0)*Lap5_2g
    return Rg.flatten()


def test_sixth_order_matrix():
    """Test various 6th-order R_h candidates using matrix method"""
    print("=" * 80)
    print("Sixth-Order R_h Candidates - Matrix Method Test")
    print("u = sin(pi*x)*sin(pi*y), -Lap(u) = f, Lap(u) = g = -f")
    print("=" * 80)
    
    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    
    ns = [9, 17, 33, 65, 129]
    
    schemes = [
        ("4th-order R_h (baseline)", "4th"),
        ("6th: R_h = I + h^2/12*Lap5 + h^4/360*Lap5^2", "6th_5pt"),
        ("6th: R_h = I + h^2/12*Lap5 + h^4/360*Lap9(Lap5)", "6th_9pt"),
        ("6th: R_h = I + h^2/12*Lap9 + h^4/360*Lap9^2", "6th_9pt2"),
    ]
    
    all_errors = {}
    
    for name, scheme in schemes:
        print(f"\n--- {name} ---")
        errors = []
        for n in ns:
            N = n - 2
            h = 1.0/(n-1)
            x = np.linspace(0,1,n); y = np.linspace(0,1,n)
            X, Y = np.meshgrid(x, y, indexing='ij')
            U_ex = u_exact(X[1:-1,1:-1], Y[1:-1,1:-1])
            
            # g = Lap(u) = -f
            g_int = -f_rhs(X[1:-1,1:-1], Y[1:-1,1:-1]).flatten()
            
            Lh = build_Lh_matrix(N, h)
            
            if scheme == "4th":
                rhs = apply_Rh_4th(g_int, N, h)
            elif scheme == "6th_5pt":
                rhs = apply_Rh_6th_candidate(g_int, N, h)
            elif scheme == "6th_9pt":
                rhs = apply_Rh_6th_9pt(g_int, N, h)
            elif scheme == "6th_9pt2":
                rhs = apply_Rh_6th_compact(g_int, N, h)
            
            u_num = spsolve(Lh, rhs).reshape((N, N))
            err = np.max(np.abs(u_num - U_ex))
            errors.append(err)
            print(f"  n={n:3d}: error = {err:.6e}")
        
        all_errors[scheme] = errors
    
    # Convergence rates
    print("\n" + "=" * 80)
    print("Convergence Rates")
    print("=" * 80)
    for name, scheme in schemes:
        errors = all_errors[scheme]
        print(f"\n{name}:")
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"  n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    test_sixth_order_matrix()
