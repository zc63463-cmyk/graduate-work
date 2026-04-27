#!/usr/bin/env python3
"""
Cyclic Reduction & FACR Algorithm Implementation
==============================================

Clean implementation of:
  1. Fourier Analysis (2D-DST direct solver)     — O(N^2 log N)
  2. Cyclic Reduction (DST-x + Thomas-y)        — O(N^2 log N), same as FA
  3. FACR(l) hybrid (CR + FA)                   — O(N^2 log log N), optimal

For solving: -∇²u = f on [0,1]^2, u|∂Ω = 0 (Dirichlet)
Discretization: standard 5-point stencil, O(h^2) accuracy
"""

import numpy as np
from scipy.fft import dst, idst
import time


# ============================================================================
# Utility
# ============================================================================

def _make_bc_array(bc_func, coord_arr1, coord_arr2, n):
    """Evaluate BC function and ensure it returns a length-n array."""
    val = bc_func(coord_arr1, coord_arr2)
    arr = np.atleast_1d(np.asarray(np.squeeze(val), dtype=float))
    if arr.shape[0] == 1:
        return np.full(n, float(arr[0]))
    return arr


def build_rhs(n, f_func, bc_func, sx=1.0, sy=1.0):
    """Build adjusted RHS for 5-point Poisson equation."""
    x = np.linspace(0, sx, n)
    y = np.linspace(0, sy, n)
    X, Y = np.meshgrid(x, y, indexing='ij')
    h = sx / (n - 1)
    N = n - 2
    
    F = np.asarray(f_func(X, Y), dtype=float)
    
    bc_l = _make_bc_array(bc_func, x[0], y, n)
    bc_r = _make_bc_array(bc_func, x[-1], y, n)
    bc_b = _make_bc_array(bc_func, x, y[0], n)
    bc_t = _make_bc_array(bc_func, x, y[-1], n)
    
    F_adj = F[1:-1, 1:-1].copy()
    F_adj[0, :]  -= bc_l[1:-1] / h**2
    F_adj[-1, :] -= bc_r[1:-1] / h**2
    F_adj[:, 0]  -= bc_b[1:-1] / h**2
    F_adj[:, -1] -= bc_t[1:-1] / h**2
    
    return x, y, F_adj, (bc_l, bc_r, bc_b, bc_t), h


# ============================================================================
# Thomas Algorithm
# ============================================================================

def thomas(a, d, b_vec):
    """
    Thomas algorithm for: a*u_{j-1} + d*u_j + a*u_{j+1} = b_j
    Constant coefficients, Dirichlet BC (u_{-1}=u_M=0). O(M).
    """
    M = len(b_vec)
    if M == 0: return np.array([])
    if M == 1: return np.array([b_vec[0]/d]) if abs(d)>1e-30 else np.zeros(1)
    
    cp, dp = np.zeros(M), np.zeros(M)
    cp[0] = a/d if abs(d)>1e-30 else 0.0
    dp[0] = b_vec[0]/d if abs(d)>1e-30 else 0.0
    
    for j in range(1, M):
        denom = d - a*cp[j-1]
        if abs(denom) < 1e-30: denom = 1e-30*(1 if denom>=0 else -1)
        cp[j] = a/denom if j < M-1 else 0.0
        dp[j] = (b_vec[j] - a*dp[j-1])/denom
    
    u = np.zeros(M)
    u[M-1] = dp[M-1]
    for j in range(M-2, -1, -1):
        u[j] = dp[j] - cp[j]*u[j+1]
    return u


def thomas_var(a_vec, d_vec, b_vec):
    """Variable-coefficient Thomas: a[j]*u_{j-1}+d[j]*u_j+a[j]*u_{j+1}=b[j]"""
    M = len(b_vec)
    if M <= 1: 
        if M==1: return np.array([b_vec[0]/d_vec[0]]) if abs(d_vec[0])>1e-30 else np.zeros(1)
        return np.array([])
    cp, dp = np.zeros(M), np.zeros(M)
    cp[0]=a_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0
    dp[0]=b_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0
    for j in range(1,M):
        denom=d_vec[j]-a_vec[j]*cp[j-1]
        if abs(denom)<1e-30: denom=1e-30*(1 if denom>=0 else -1)
        cp[j]=a_vec[j]/denom if j<M-1 else 0.0
        dp[j]=(b_vec[j]-a_vec[j]*dp[j-1])/denom
    u=np.zeros(M); u[M-1]=dp[M-1]
    for j in range(M-2,-1,-1): u[j]=dp[j]-cp[j]*u[j+1]
    return u


# ============================================================================
# Method 1: Pure Fourier Analysis (2D-DST)
# ============================================================================

def fa_solver(n, f_func, bc_func, sx=1.0, sy=1.0):
    """2D-DST direct Poisson solver. O(N^2 log N)."""
    x, y, F_adj, bcs, h = build_rhs(n, f_func, bc_func, sx, sy)
    N = n-2; bc_l,bc_r,bc_b,bc_t = bcs
    
    # 2D DST
    Fh = np.zeros((N,N))
    for i in range(N): Fh[i,:] = dst(F_adj[i,:], type=1, norm='ortho')
    for j in range(N): Fh[:,j] = dst(Fh[:,j], type=1, norm='ortho')
    
    # Eigenvalues
    kv = np.arange(1,N+1)
    lam = 2.0*(1.0-np.cos(kv*np.pi/(N+1)))/h**2
    lxx,lyy = np.meshgrid(lam,lam,indexing='ij')
    ltot = lxx+lyy
    
    # Solve & inverse DST
    Uh = np.zeros((N,N))
    mask = np.abs(ltot)>1e-14; Uh[mask] = Fh[mask]/ltot[mask]
    
    Ui = np.zeros((N,N))
    for i in range(N): Ui[i,:] = idst(Uh[i,:],type=1,norm='ortho')
    for j in range(N): Ui[:,j] = idst(Ui[:,j],type=1,norm='ortho')
    
    U = np.zeros((n,n)); U[0,:]=bc_l;U[-1,:]=bc_r;U[:,0]=bc_b;U[:,-1]=bc_t
    U[1:-1,1:-1] = Ui
    return U


# ============================================================================
# Method 2: Cyclic Reduction (DST-x + Thomas-y)
# ============================================================================

def cr_solver(n, f_func, bc_func, sx=1.0, sy=1.0):
    """
    CR via diagonalization: 1D-DST in x → N scalar tridiagonal systems → Thomas.
    Mathematically identical to FA solver.
    """
    x,y,F_adj,bcs,h = build_rhs(n,f_func,bc_func,sx,sy); N=n-2
    bc_l,bc_r,bc_b,bc_t=bcs
    
    # DST in x
    Fh = np.zeros((N,N))
    for j in range(N): Fh[:,j] = dst(F_adj[:,j], type=1, norm='ortho')
    
    # Per-mode eigenvalues of D_block = tridiag(-1,4,-1)
    kv=np.arange(1,N+1); dp=4.0-2.0*np.cos(kv*np.pi/(N+1)); av=-1.0
    
    # Solve each mode's tridiagonal system
    Uh=np.zeros((N,N)); h2=h*h
    for p in range(N): Uh[p,:]=thomas(av, dp[p], h2*Fh[p,:])
    
    # Inverse DST
    Ui=np.zeros((N,N))
    for j in range(N): Ui[:,j]=idst(Uh[:,j],type=1,norm='ortho')
    
    U=np.zeros((n,n)); U[0,:]=bc_l;U[-1,:]=bc_r;U[:,0]=bc_b;U[:,-1]=bc_t
    U[1:-1,1:-1]=Ui
    return U


# ============================================================================
# Method 3: FACR(l) — Corrected implementation
# ============================================================================

def thomas_sym(d_vec, e_vec, b_vec):
    """
    Thomas algorithm for symmetric tridiagonal:
      e[j-1]*u_{j-1} + d[j]*u_j + e[j]*u_{j+1} = b[j]
    
    d_vec: length M (diagonal)
    e_vec: length M-1 (off-diagonal; e[k] couples row k and row k+1)
    b_vec: length M (RHS)
    """
    M = len(b_vec)
    if M == 0: return np.array([])
    if M == 1:
        return np.array([b_vec[0]/d_vec[0]]) if abs(d_vec[0])>1e-30 else np.zeros(1)
    
    cp = np.zeros(M)
    dp = np.zeros(M)
    
    cp[0] = e_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0
    dp[0] = b_vec[0]/d_vec[0] if abs(d_vec[0])>1e-30 else 0.0
    
    for j in range(1, M):
        c_j = e_vec[j-1]  # sub-diagonal at row j
        denom = d_vec[j] - c_j*cp[j-1]
        if abs(denom) < 1e-30: denom = 1e-30*(1 if denom>=0 else -1)
        cp[j] = (e_vec[j]/denom) if j < M-1 else 0.0
        dp[j] = (b_vec[j] - c_j*dp[j-1])/denom
    
    u = np.zeros(M)
    u[M-1] = dp[M-1]
    for j in range(M-2, -1, -1):
        u[j] = dp[j] - cp[j]*u[j+1]
    return u


def facr_solver(n, f_func, bc_func, sx=1.0, sy=1.0, l=None):
    """
    FACR(l) hybrid algorithm — vectorized across all Fourier modes.
    
    Algorithm:
      1. 1D-DST in x → decouple into N independent tridiagonal systems
      2. l steps of CR on each y-line system (vectorized across modes)
      3. Solve small reduced systems (general symmetric tridiagonal)
      4. Back-substitute to recover all y-line values (vectorized)
      5. Inverse DST in x
    
    Key insight: after each CR step, the reduced system is a symmetric
    tridiagonal with VARIABLE diagonal d[] and VARIABLE off-diagonal e[].
    We track the full (d[], e[]) structure per level for correct multi-step CR.
    
    Complexity: O(N^2 log log N) with optimal l ~ log2(log2(N)).
    """
    x,y,F_adj,bcs,h = build_rhs(n,f_func,bc_func,sx,sy); N=n-2
    bc_l,bc_r,bc_b,bc_t=bcs
    
    if N <= 2: return fa_solver(n,f_func,bc_func,sx,sy)
    
    if l is None:
        l=max(0,int(np.round(np.log2(max(1,np.log2(N))))))
        l=min(l,int(np.log2(N))-1)
    if l <= 0: return fa_solver(n,f_func,bc_func,sx,sy)
    
    h2=h*h
    
    # Phase 1: DST in x
    Fh = np.zeros((N,N))
    for j in range(N): Fh[:,j]=dst(F_adj[:,j],type=1,norm='ortho')
    
    kv=np.arange(1,N+1); d_eig=4.0-2.0*np.cos(kv*np.pi/(N+1)); a0=-1.0
    
    # Phase 2: Vectorized CR forward reduction across all modes
    # State shapes: (N_modes, M) for d and b, (N_modes, M-1) for e
    cur_d = np.tile(d_eig.reshape(-1,1), (1, N))   # (N, N)
    cur_e = np.full((N, N-1), a0)                    # (N, N-1)
    cur_b = h2 * Fh                                   # (N, N)
    
    levels = []  # save (d, e, b) at each level for back-substitution
    
    for step in range(l):
        cur_M = cur_d.shape[1]
        if cur_M <= 2:
            break
        
        # Save current level
        levels.append((cur_d.copy(), cur_e.copy(), cur_b.copy()))
        
        M_new = (cur_M + 1) // 2
        even_idx = np.arange(0, cur_M, 2)  # [0, 2, 4, ...]
        
        # Start with even-row values
        d_new = cur_d[:, even_idx].copy()   # (N, M_new)
        b_new = cur_b[:, even_idx].copy()
        
        # Left odd neighbor contribution (j-1 is odd, exists if j > 0)
        has_left = even_idx > 0
        left_odd = even_idx[has_left] - 1
        
        d_new[:, has_left] -= cur_e[:, left_odd]**2 / cur_d[:, left_odd]
        b_new[:, has_left] -= (cur_e[:, left_odd] / cur_d[:, left_odd]) * cur_b[:, left_odd]
        
        # Right odd neighbor contribution (j+1 is odd, exists if j < cur_M-1)
        has_right = even_idx < cur_M - 1
        right_odd = even_idx[has_right] + 1
        e_at_even = cur_e[:, even_idx[has_right]]  # e values at even positions
        
        d_new[:, has_right] -= e_at_even**2 / cur_d[:, right_odd]
        b_new[:, has_right] -= (e_at_even / cur_d[:, right_odd]) * cur_b[:, right_odd]
        
        # Reduced off-diagonal: e_new[k] = -e[j-1]*e[j]/d[j] where j = 2k+1
        if M_new > 1:
            odd_between = np.arange(1, cur_M, 2)[:M_new-1]  # odd indices between even rows
            e_new = -cur_e[:, odd_between-1] * cur_e[:, odd_between] / cur_d[:, odd_between]
        else:
            e_new = np.empty((N, 0))
        
        cur_d = d_new
        cur_e = e_new
        cur_b = b_new
    
    # Phase 3: Solve reduced systems (per-mode Thomas for small systems)
    M_red = cur_d.shape[1]
    u_cur = np.zeros((N, M_red))
    
    for p in range(N):
        if M_red == 1:
            u_cur[p, 0] = cur_b[p, 0] / cur_d[p, 0] if abs(cur_d[p, 0]) > 1e-30 else 0.0
        elif cur_e.shape[1] > 0:
            u_cur[p, :] = thomas_sym(cur_d[p], cur_e[p], cur_b[p])
        else:
            u_cur[p, :M_red] = cur_b[p, :M_red] / cur_d[p, :M_red]
    
    # Phase 4: Vectorized back-substitution
    for rev in reversed(range(len(levels))):
        d_lv, e_lv, b_lv = levels[rev]
        M_lv = d_lv.shape[1]
        M_even = u_cur.shape[1]
        
        u_full = np.zeros((N, M_lv))
        
        # Place even-indexed values: u_full[:, 0::2] = u_cur[:, :M_even]
        u_full[:, :2*M_even:2] = u_cur
        
        # Recover odd-indexed values (vectorized across modes)
        # Interior odd rows: j = 2k+1 with both neighbors
        n_interior_odd = min(M_even - 1, M_lv // 2)
        if n_interior_odd > 0:
            odd_int = np.arange(1, 2*n_interior_odd + 1, 2)  # [1, 3, 5, ...]
            u_full[:, odd_int] = (b_lv[:, odd_int] 
                - e_lv[:, odd_int-1] * u_full[:, odd_int-1] 
                - e_lv[:, odd_int] * u_full[:, odd_int+1]) / d_lv[:, odd_int]
        
        # Last odd row (M_lv even): only left neighbor
        if M_lv % 2 == 0:
            j_last = M_lv - 1
            u_full[:, j_last] = (b_lv[:, j_last] - e_lv[:, j_last-1] * u_full[:, j_last-1]) / d_lv[:, j_last]
        
        u_cur = u_full
    
    Uh = u_cur  # (N, N) — solution in Fourier space
    
    # Phase 5: Inverse DST
    Ui=np.zeros((N,N))
    for j in range(N): Ui[:,j]=idst(Uh[:,j],type=1,norm='ortho')
    
    U=np.zeros((n,n)); U[0,:]=bc_l;U[-1,:]=bc_r;U[:,0]=bc_b;U[:,-1]=bc_t
    U[1:-1,1:-1]=Ui
    return U


def _facr_mode(a0, d0, b_orig, l_steps):
    """
    FACR for a single Fourier mode — scalar version for reference/testing.
    
    Tracks full symmetric tridiagonal structure (d[], e[]) through all
    CR levels. See facr_solver docstring for algorithm details.
    """
    M = len(b_orig)
    if M <= 2 or l_steps <= 0:
        return thomas(a0, d0, b_orig)
    
    levels = []
    cur_d = np.full(M, d0)
    cur_e = np.full(M - 1, a0)
    cur_b = b_orig.copy()
    
    for step in range(l_steps):
        cur_M = len(cur_d)
        if cur_M <= 2:
            break
        levels.append({'d': cur_d.copy(), 'e': cur_e.copy(), 'b': cur_b.copy()})
        
        M_new = (cur_M + 1) // 2
        even_idx = np.arange(0, cur_M, 2)
        
        d_new = cur_d[even_idx].copy()
        b_new = cur_b[even_idx].copy()
        
        has_left = even_idx > 0
        left_odd = even_idx[has_left] - 1
        d_new[has_left] -= cur_e[left_odd]**2 / cur_d[left_odd]
        b_new[has_left] -= (cur_e[left_odd] / cur_d[left_odd]) * cur_b[left_odd]
        
        has_right = even_idx < cur_M - 1
        right_odd = even_idx[has_right] + 1
        e_at_even = cur_e[even_idx[has_right]]
        d_new[has_right] -= e_at_even**2 / cur_d[right_odd]
        b_new[has_right] -= (e_at_even / cur_d[right_odd]) * cur_b[right_odd]
        
        if M_new > 1:
            odd_between = np.arange(1, cur_M, 2)[:M_new-1]
            e_new = -cur_e[odd_between-1] * cur_e[odd_between] / cur_d[odd_between]
        else:
            e_new = np.array([])
        
        cur_d = d_new
        cur_e = e_new
        cur_b = b_new
    
    if len(cur_d) == 1:
        u_red = np.array([cur_b[0] / cur_d[0]] if abs(cur_d[0]) > 1e-30 else [0.0])
    elif len(cur_e) > 0:
        u_red = thomas_sym(cur_d, cur_e, cur_b)
    else:
        u_red = cur_b / cur_d
    
    u_cur = u_red
    for rev in reversed(range(len(levels))):
        ld = levels[rev]
        d_lv, e_lv, b_lv = ld['d'], ld['e'], ld['b']
        M_lv = len(d_lv)
        M_even = len(u_cur)
        u_full = np.zeros(M_lv)
        
        u_full[0::2] = u_cur[:M_even]
        
        n_int = min(M_even - 1, M_lv // 2)
        if n_int > 0:
            odd_int = np.arange(1, 2*n_int+1, 2)
            u_full[odd_int] = (b_lv[odd_int] - e_lv[odd_int-1]*u_full[odd_int-1]
                              - e_lv[odd_int]*u_full[odd_int+1]) / d_lv[odd_int]
        
        if M_lv % 2 == 0:
            j_last = M_lv - 1
            u_full[j_last] = (b_lv[j_last] - e_lv[j_last-1]*u_full[j_last-1]) / d_lv[j_last]
        
        u_cur = u_full
    
    return u_cur


# ============================================================================
# Tests
# ============================================================================

def test_correctness():
    print("=" * 75)
    print("CORRECTNESS TEST")
    print("=" * 75)
    
    def u1(x,y): return np.sin(np.pi*x)*np.sin(np.pi*y)
    def f1(x,y): return 2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bz(x,y): return 0.0
    
    ns=[9,17,33,65,129,257]
    print(f"\n{'n':>5s}|{'FA':>12s}|{'CR':>12s}|{'FACR':>12s}|{'FA=CR?':>10s}")
    print("-"*56)
    for n in ns:
        X,Y=np.meshgrid(np.linspace(0,1,n),np.linspace(0,1,n),indexing='ij'); Ue=u1(X,Y)
        U_fa=fa_solver(n,f1,bz); U_cr=cr_solver(n,f1,bz); U_fc=facr_solver(n,f1,bz)
        print(f"{n:5d}|{np.max(np.abs(U_fa-Ue)):12.2e}|{np.max(np.abs(U_cr-Ue)):12.2e}|"
              f"{np.max(np.abs(U_fc-Ue)):12.2e}|{np.max(np.abs(U_fa-U_cr)):10.1e}")


def test_convergence():
    print("\n"+"="*75); print("CONVERGENCE RATE (expect ~2.0)"); print("="*75)
    def ue(x,y): return np.sin(2*np.pi*x)*np.sin(3*np.pi*y)
    def fe(x,y): return 13*np.pi**2*np.sin(2*np.pi*x)*np.sin(3*np.pi*y)
    def bz(x,y): return 0.0
    ns=[9,17,33,65,129,257,513]; errs=[]
    print(f"\n  {'n':>5s}{'h':>8s}{'error':>14s}{'rate':>6s}")
    print("  "+"-"*38)
    for n in ns:
        h=1/(n-1); U=fa_solver(n,fe,bz); X,Y=np.meshgrid(np.linspace(0,1,n),np.linspace(0,1,n),indexing='ij')
        e=np.max(np.abs(U-ue(X,Y))); errs.append(e)
        r=f"{np.log(errs[-2]/errs[-1])/np.log(ns[-1]/ns[-2]):.2f}" if len(errs)>1 else ""
        print(f"  {n:5d}{h:8.5f}{e:14.6e}{r:>6s}")


def test_performance():
    print("\n"+"="*75); print("PERFORMANCE BENCHMARK"); print("="*75)
    def frhs(x,y): return 2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bcz(x,y): return 0.0
    ns=[33,65,129,257,513]
    solvers=[("FA",lambda n:fa_solver(n,frhs,bcz)),("CR",lambda n:cr_solver(n,frhs,bcz)),
             ("FACR",lambda n:facr_solver(n,frhs,bcz))]
    print(f"\n  {'n':>5s}",end="")
    for nm,_ in solvers: print(f"|{nm:>14s}",end="")
    print(); print("  "+"-"*(6+16*len(solvers)))
    for n in ns:
        print(f"  {n:5d}",end="")
        for _,fn in solvers:
            t0=time.time();fn(n);print(f"|{time.time()-t0:14.4f}s",end="")
        print()


def test_facr_l():
    print("\n"+"="*75); print("FACR(l) VARIANTS"); print("="*75)
    def frhs(x,y): return 2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)
    def bcz(x,y): return 0.0
    n=513; N=n-2; lopt=max(0,int(np.round(np.log2(max(1,np.log2(N))))))
    print(f"\n  Grid n={n}, optimal l~{lopt}\n")
    print(f"  {'l':>3s}|{'Time(ms)':>10s}|{'Error':>12s}|Note")
    print("  "+"-"*48)
    X,Y=np.meshgrid(np.linspace(0,1,n),np.linspace(0,1,n),indexing='ij'); Uex=np.sin(np.pi*X)*np.sin(np.pi*Y)
    for lv in range(0,min(9,int(np.log2(N)))):
        t0=time.time();U=facr_solver(n,frhs,bcz,l=lv);ms=(time.time()-t0)*1000
        e=np.max(np.abs(U-Uex)); note="Pure FA" if lv==0 else ("***OPT***" if lv==lopt else "")
        print(f"  {lv:3d}|{ms:10.2f}|{e:12.2e}|{note}")


if __name__=="__main__":
    test_correctness()
    test_convergence()
    test_performance()
    test_facr_l()
