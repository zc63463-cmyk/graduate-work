#!/usr/bin/env python3
"""
FFT9 Sixth-Order Scheme - Correct Approach

The key insight I was missing: For the POISSON equation specifically,
the 6th-order scheme uses the SAME L_h but with a DIFFERENT R_h that
exploits the fact that g = Lap(u) and thus Lap(g) = Lap^2(u) can be
expressed in terms of g using the PDE itself.

For -Lap(u) = f (Poisson equation):
  Lap(u) = -f = g

The 4th-order scheme: L_h * u = R_h^(4) * g
  L_h = (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1]
  R_h^(4) = I + (h^2/12)*Lap5_h

This has truncation error T = (h^4/360)*Lap^3(u) + O(h^6)

Now, Lap^3(u) = Lap^2(Lap(u)) = Lap^2(g) = Lap^2(-f)

For the 6th-order scheme, we need to cancel this h^4 term:
  R_h^(6) = R_h^(4) + (h^4/360)*Lap^2_h

But Lap^2_h(g) requires a 5x5 stencil. HOWEVER, we can compute
Lap^2(g) numerically by applying Lap5_h twice, which DOES use a 5x5 
stencil but can be computed on the same grid.

The resulting scheme is NOT a compact 3x3 scheme, but it IS 6th-order
accurate. The paper calls it "compact" in the sense that L_h is compact
(3x3), even though R_h^(6) effectively uses a wider stencil.

Wait, but the paper writes R_h as (1/360)*[0 48 0; 48 0 48; 0 48 0],
which IS a 3x3 stencil (5-point). How can a 5-point R_h give 6th order?

Let me re-examine the paper's formula more carefully.

The paper states for the Poisson case:
  L_h * u = R_h * f   (where f = g = Lap(u))

with R_h = (1/360) * [0 48 0; 48 0 48; 0 48 0]

This R_h has NO center coefficient! That means:
  R_h * f at (i,j) = (48/360)*(f_{i-1,j} + f_{i+1,j} + f_{i,j-1} + f_{i,j+1})
                    = (2/15)*(f_{i-1,j} + f_{i+1,j} + f_{i,j-1} + f_{i,j+1})

For the eigenvalue on mode (k,l):
  lambda_R = (2/15)*(2*cos(kh) + 2*cos(lh))
           = (4/15)*(cos(kh) + cos(lh))

Let's check: For -s * lambda_R = lambda_L (where s = k^2+l^2):
  -(k^2+l^2) * (4/15)*(cos(kh) + cos(lh)) 
  = (-20 + 8*(cos(kh)+cos(lh)) + 4*cos(kh)*cos(lh))/(6h^2)

The left side is O(1) while right is O(1/h^2). This can't work unless
there's a factor of h^2 somewhere I'm missing.

I think the issue is that the paper's notation uses DIFFERENT normalization.
Let me reconsider: maybe the paper writes:
  L_h * u = R_h * f

where L_h and R_h are the OPERATORS (without the 1/(6h^2) factor), and
the 1/(6h^2) is applied to the left side only:

  (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1] * u = (1/360) * [0 48 0; 48 0 48; 0 48 0] * f

So: lambda_L = (-20 + 8*c_k + 8*c_l + 4*c_k*c_l) / (6h^2)
    lambda_R = (48*c_k + 48*c_l) / 360 = (4/15)*(c_k + c_l)

Wait, that's the same thing. Let me check if maybe the paper's R_h 
is for a DIFFERENT equation form.

Actually, I just realized: maybe the paper uses the form
  L_h * u = R_h * f
where L_h and R_h both include the 1/(6h^2) factor.

Then: L_h = (1/(6h^2)) * [-20, 8, 8, 4] (in eigenvalue terms)
      R_h * f means applying R_h stencil and then the equation reads
      L_h u = R_h f where both sides are on the grid.

Let me try a completely different interpretation. The paper's formula might be:
  Sum over 9-point stencil of L_h * u = Sum over 9-point stencil of R_h * f

where both sums are at the same grid point, and the equation is for
Lap(u) = f (positive Laplacian).

With h_x = h_y = h, the stencil coefficients for L_h are divided by 6h^2,
and the R_h coefficients are divided by 360.

For the eigenvalue equation:
  [(-20 + 8c_k + 8c_l + 4c_k*c_l)/(6h^2)] * u_hat = [(48c_k + 48c_l)/360] * f_hat

Since f = Lap(u), f_hat = -(k^2+l^2)*u_hat, so:
  (-20 + 8c_k + 8c_l + 4c_k*c_l)/(6h^2) = -(k^2+l^2) * (48c_k + 48c_l)/360

Let's check for k=l=pi, h=1/10:
  c_k = c_l = cos(pi/10) = 0.95106
  s = 2*pi^2 = 19.739

  LHS = (-20 + 8*0.95106 + 8*0.95106 + 4*0.95106^2)/(6*0.01)
      = (-20 + 15.217 + 15.217 + 3.614)/0.06
      = 14.048/0.06 = 234.13

  RHS = -19.739 * (48*0.95106 + 48*0.95106)/360
      = -19.739 * 91.302/360
      = -19.739 * 0.25362
      = -5.006

LHS = 234.13, RHS = -5.006. These are completely different!

Something is fundamentally wrong with my interpretation. Let me try yet
another interpretation.

Maybe the paper means:
  (1/(6h^2)) * L_9 * u = (1/(6h^2)) * R_9 * f

where L_9 = [1 4 1; 4 -20 4; 1 4 1] and R_9 = (6h^2/360) * [0 48 0; 48 0 48; 0 48 0]
                                                          = (h^2/60) * [0 48 0; 48 0 48; 0 48 0]

No that doesn't make sense either.

OK let me just try ALL possible sign/factor combinations numerically.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def build_Lh_matrix(N, h):
    """Build L_h = (1/(6h^2)) * [1 4 1; 4 -20 4; 1 4 1]"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            rows.append(idx); cols.append(idx); vals.append(-20.0/(6.0*h**2))
            if i > 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(4.0/(6.0*h**2))
            if i < N-1: rows.append(idx); cols.append((i+1)*N+j); vals.append(4.0/(6.0*h**2))
            if j > 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(4.0/(6.0*h**2))
            if j < N-1: rows.append(idx); cols.append(i*N+(j+1)); vals.append(4.0/(6.0*h**2))
            if i>0 and j>0: rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j>0: rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(1.0/(6.0*h**2))
            if i>0 and j<N-1: rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
            if i<N-1 and j<N-1: rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(1.0/(6.0*h**2))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def build_Rh_matrix(N, h, center, edge, corner=0.0):
    """Build R_h = [corner edge corner; edge center edge; corner edge corner]"""
    total = N * N
    rows, cols, vals = [], [], []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            if center != 0: rows.append(idx); cols.append(idx); vals.append(center)
            if i > 0 and edge != 0: rows.append(idx); cols.append((i-1)*N+j); vals.append(edge)
            if i < N-1 and edge != 0: rows.append(idx); cols.append((i+1)*N+j); vals.append(edge)
            if j > 0 and edge != 0: rows.append(idx); cols.append(i*N+(j-1)); vals.append(edge)
            if j < N-1 and edge != 0: rows.append(idx); cols.append(i*N+(j+1)); vals.append(edge)
            if i>0 and j>0 and corner != 0: rows.append(idx); cols.append((i-1)*N+(j-1)); vals.append(corner)
            if i<N-1 and j>0 and corner != 0: rows.append(idx); cols.append((i+1)*N+(j-1)); vals.append(corner)
            if i>0 and j<N-1 and corner != 0: rows.append(idx); cols.append((i-1)*N+(j+1)); vals.append(corner)
            if i<N-1 and j<N-1 and corner != 0: rows.append(idx); cols.append((i+1)*N+(j+1)); vals.append(corner)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(total, total))


def exhaustive_6th_order_search():
    """
    Try many combinations of R_h coefficients to find 6th order.
    
    For Lap(u) = g where g = -f (with -Lap(u) = f):
    L_h * u = R_h * g
    
    We know 4th order works with R_h = [0 1/12 0; 1/12 2/3 1/12; 0 1/12 0]
    
    For 6th order, try R_h with additional h^4-dependent terms.
    General form: R_h = alpha*I + beta*h^2*Lap5 + gamma*h^4*Lap5^2
    
    where Lap5 = [0 1 0; 1 -4 1; 0 1 0]/h^2
    
    R_h stencil: 
    corner = gamma*h^4/h^4 = gamma  (from Lap5^2)
    edge = beta*h^2/h^2 + 2*gamma*h^4/h^4 = beta + 2*gamma  (from Lap5 + Lap5^2)
    center = alpha - 4*beta - 4*gamma  (wait this isn't right)
    
    Actually, let me just expand:
    R_h = alpha*I + beta*(h^2)*Lap5 + gamma*(h^4)*Lap5^2
    
    Lap5 stencil: [0 1/h^2 0; 1/h^2 -4/h^2 1/h^2; 0 1/h^2 0]
    h^2*Lap5 stencil: [0 1 0; 1 -4 1; 0 1 0]
    
    Lap5^2 stencil (5x5, applied twice):
    The Lap5^2 is NOT a 3x3 stencil. So this approach gives a non-compact R_h.
    
    Let me instead try the general 3x3 symmetric R_h:
    R_h = [a b a; b c b; a b a]
    
    and determine a, b, c from the Fourier symbol matching condition.
    
    For the equation L_h * u = R_h * g (where g = Lap(u)):
    lambda_L / lambda_R should equal -1/(k^2+l^2) for 4th order (3 conditions)
    and also match the h^4 terms for 6th order.
    
    But with only 3 parameters (a,b,c), we can only satisfy 3 conditions,
    which gives 4th order. For 6th order we need more parameters.
    
    UNLESS the Poisson equation provides additional constraints.
    Since g = Lap(u), we have additional relationships between g and u.
    
    The Mehrstellenverfahren (Hermitian method) for 6th order uses:
    L_h * u = R_h * f + S_h * (Lap(f))
    
    where S_h is another operator. But this requires computing Lap(f),
    which is a second derivative of the RHS.
    
    For the specific Poisson case (-Lap(u) = f), we have:
    Lap(f) = Lap(-Lap(u)) = -Lap^2(u)
    
    And from the 4th order scheme:
    L_h * u = R_h^(4) * (-f) + O(h^4)
    
    The O(h^4) term involves Lap^2(f), which for Poisson equals Lap^3(u).
    By computing Lap(f) on the grid, we can correct this term.
    
    So the 6th order scheme is:
    L_h * u = R_h^(4) * (-f) + correction
    
    where correction = (h^4/360) * Lap_h(f) (or something similar)
    
    Let me test this numerically.
    """
    print("=" * 80)
    print("Sixth-Order Search: Using Lap(f) correction")
    print("=" * 80)
    
    def u_exact(x, y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f_rhs(x, y): return 2.0*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    
    ns = [9, 17, 33, 65, 129]
    
    # Test: L_h * u = R_h^(4) * (-f) + c * h^4 * Lap5(f)
    print("\n--- Scheme: R_h^(4)*(-f) + c*h^4*Lap5(f) ---")
    for c_val in [1.0/360.0, 1.0/180.0, 1.0/90.0, 1.0/720.0, -1.0/360.0]:
        print(f"\n  c = {c_val:.6f}")
        errors = []
        for n in ns:
            N = n - 2; h = 1.0/(n-1)
            x = np.linspace(0,1,n); y = np.linspace(0,1,n)
            X, Y = np.meshgrid(x, y, indexing='ij')
            U_ex = u_exact(X[1:-1,1:-1], Y[1:-1,1:-1])
            F = f_rhs(X[1:-1,1:-1], Y[1:-1,1:-1])
            
            Lh = build_Lh_matrix(N, h)
            
            # R_h^(4) * (-f) 
            G = -F  # g = -f
            rhs_4th = np.zeros((N, N))
            for i in range(N):
                for j in range(N):
                    rhs_4th[i,j] = (2.0/3.0)*G[i,j]
                    if i > 0: rhs_4th[i,j] += (1.0/12.0)*G[i-1,j]
                    if i < N-1: rhs_4th[i,j] += (1.0/12.0)*G[i+1,j]
                    if j > 0: rhs_4th[i,j] += (1.0/12.0)*G[i,j-1]
                    if j < N-1: rhs_4th[i,j] += (1.0/12.0)*G[i,j+1]
            
            # Lap5(f) with zero BC
            F_pad = np.zeros((N+2, N+2))
            F_pad[1:-1, 1:-1] = F
            Lap5_f = (F_pad[2:,1:-1] + F_pad[:-2,1:-1] + 
                      F_pad[1:-1,2:] + F_pad[1:-1,:-2] - 4*F) / h**2
            
            rhs = rhs_4th.flatten() + c_val * h**4 * Lap5_f.flatten()
            
            u_num = spsolve(Lh, rhs).reshape((N, N))
            err = np.max(np.abs(u_num - U_ex))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")
    
    # Test: L_h * u = R_h^(4) * (-f) + c * h^4 * Lap9(f)
    print("\n\n--- Scheme: R_h^(4)*(-f) + c*h^4*Lap9(f) ---")
    for c_val in [1.0/360.0, 1.0/180.0, 1.0/90.0, 1.0/720.0, -1.0/360.0]:
        print(f"\n  c = {c_val:.6f}")
        errors = []
        for n in ns:
            N = n - 2; h = 1.0/(n-1)
            x = np.linspace(0,1,n); y = np.linspace(0,1,n)
            X, Y = np.meshgrid(x, y, indexing='ij')
            U_ex = u_exact(X[1:-1,1:-1], Y[1:-1,1:-1])
            F = f_rhs(X[1:-1,1:-1], Y[1:-1,1:-1])
            
            Lh = build_Lh_matrix(N, h)
            
            G = -F
            rhs_4th = np.zeros((N, N))
            for i in range(N):
                for j in range(N):
                    rhs_4th[i,j] = (2.0/3.0)*G[i,j]
                    if i > 0: rhs_4th[i,j] += (1.0/12.0)*G[i-1,j]
                    if i < N-1: rhs_4th[i,j] += (1.0/12.0)*G[i+1,j]
                    if j > 0: rhs_4th[i,j] += (1.0/12.0)*G[i,j-1]
                    if j < N-1: rhs_4th[i,j] += (1.0/12.0)*G[i,j+1]
            
            # Lap9(f) with zero BC
            F_pad = np.zeros((N+2, N+2))
            F_pad[1:-1, 1:-1] = F
            Lap9_f = (1.0/(6.0*h**2)) * (
                -20*F +
                4*(F_pad[2:,1:-1] + F_pad[:-2,1:-1] + F_pad[1:-1,2:] + F_pad[1:-1,:-2]) +
                1*(F_pad[2:,2:] + F_pad[2:,:-2] + F_pad[:-2,2:] + F_pad[:-2,:-2])
            )
            
            rhs = rhs_4th.flatten() + c_val * h**4 * Lap9_f.flatten()
            
            u_num = spsolve(Lh, rhs).reshape((N, N))
            err = np.max(np.abs(u_num - U_ex))
            errors.append(err)
            print(f"    n={n:3d}: error = {err:.6e}")
        
        for i in range(1, len(ns)):
            rate = np.log(errors[i-1]/errors[i])/np.log(ns[i]/ns[i-1])
            print(f"    n={ns[i]:3d}: rate = {rate:.2f}")


if __name__ == "__main__":
    exhaustive_6th_order_search()
