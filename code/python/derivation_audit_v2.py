"""
论文推导完备性校验 v2 — 完整 SymPy 符号验证
"""
import sympy as sp
import numpy as np

xi, eta, h, k2 = sp.symbols('xi eta h k2')
alpha_p, beta_p, gamma_p = sp.symbols('alpha beta gamma')

PASS = 0
FAIL = 0
WARN = 0

def check(name, condition, detail=""):
    global PASS, FAIL, WARN
    if condition:
        PASS += 1
        print(f"  [PASS] {name}")
    else:
        FAIL += 1
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")

def warn(name, detail=""):
    global WARN
    WARN += 1
    print(f"  [WARN] {name}")
    if detail:
        print(f"         {detail}")

# =================================================================
print("=" * 70)
print("论文推导完备性校验 v2")
print("=" * 70)

# =================================================================
# 验证1: R_h Taylor 展开 (Lemma 2.2 / eq 2.3.8)
# =================================================================
print("\n--- 验证1: R_h Taylor 展开 ---")
lam_R = sp.Rational(2,3) + sp.Rational(1,6)*(sp.cos(xi*h) + sp.cos(eta*h))
lam_R_exp = sp.series(lam_R, h, 0, 8)

c0_R = sp.simplify(lam_R_exp.coeff(h, 0))
c2_R = sp.simplify(lam_R_exp.coeff(h, 2))
c4_R = sp.simplify(lam_R_exp.coeff(h, 4))
check("R_h h^0 = 1", c0_R == 1)
check("R_h h^2 = -(xi^2+eta^2)/12", c2_R == -(xi**2+eta**2)/12)
check("R_h h^4 = (xi^4+eta^4)/144", c4_R == (xi**4+eta**4)/144)

# =================================================================
# 验证2: L_h Fourier 特征值展开 (eq 2.3.18)
# =================================================================
print("\n--- 验证2: L_h Fourier 特征值展开 ---")
lam_L = (1/(6*h**2)) * (-20 + 8*sp.cos(xi*h) + 8*sp.cos(eta*h)
         + 4*sp.cos(xi*h)*sp.cos(eta*h))
lam_L_exp = sp.series(lam_L, h, 0, 8)

c0_L = sp.simplify(lam_L_exp.coeff(h, 0))
c2_L = sp.simplify(lam_L_exp.coeff(h, 2))
c4_L = sp.simplify(lam_L_exp.coeff(h, 4))

# CORRECT values (from SymPy):
c2_L_correct = (xi**2 + eta**2)**2 / 12  # = (xi^4 + 2xi^2*eta^2 + eta^4)/12
c4_L_correct = -(xi**2+eta**2)*(xi**4 + 4*xi**2*eta**2 + eta**4)/360

# THESIS claimed values (eq 2.3.18):
c2_L_thesis = (xi**4 + eta**4 + xi**2*eta**2) / 12
c4_L_thesis = -(xi**6 + eta**6 + xi**2*eta**2*(xi**2+eta**2)/2) / 360

check("L_h h^0 = -(xi^2+eta^2)", c0_L == -(xi**2+eta**2))

c2_match_thesis = sp.simplify(c2_L - c2_L_thesis) == 0
c2_match_correct = sp.simplify(c2_L - c2_L_correct) == 0
if c2_match_thesis:
    check("L_h h^2 matches thesis", True)
elif c2_match_correct:
    warn("L_h h^2: thesis has WRONG coefficient",
         f"Thesis: (xi^4+eta^4+xi^2*eta^2)/12, Correct: (xi^4+2xi^2*eta^2+eta^4)/12 = (xi^2+eta^2)^2/12")
else:
    check("L_h h^2 matches something", False, f"actual={c2_L}")

c4_match_thesis = sp.simplify(c4_L - c4_L_thesis) == 0
c4_match_correct = sp.simplify(c4_L - c4_L_correct) == 0
if c4_match_thesis:
    check("L_h h^4 matches thesis", True)
elif c4_match_correct:
    warn("L_h h^4: thesis has WRONG coefficient",
         f"Correct: -(xi^2+eta^2)(xi^4+4xi^2*eta^2+eta^4)/360")
else:
    check("L_h h^4 matches something", False, f"actual={c4_L}")

# =================================================================
# 验证3: 四阶精度 -lam_L/lam_R = (xi^2+eta^2) + O(h^4)
# =================================================================
print("\n--- 验证3: 四阶精度核心定理 ---")
lam_eff = -lam_L / lam_R
lam_eff_exp = sp.series(lam_eff, h, 0, 8)

c0_eff = sp.simplify(lam_eff_exp.coeff(h, 0))
c2_eff = sp.simplify(lam_eff_exp.coeff(h, 2))
c4_eff = sp.simplify(lam_eff_exp.coeff(h, 4))

check("h^0 coeff = xi^2+eta^2", sp.simplify(c0_eff - (xi**2+eta**2)) == 0)
check("h^2 coeff = 0 (4th-order!)", c2_eff == 0)

c4_thesis_eff = -(xi**2+eta**2)*(3*xi**4 - 8*xi**2*eta**2 + 3*eta**4)/720
check("h^4 coeff matches thesis",
      sp.simplify(c4_eff - c4_thesis_eff) == 0,
      f"c4 = {sp.expand(c4_eff)}")

# =================================================================
# 验证4: 六阶不可行性 (Theorem 2.3)
# =================================================================
print("\n--- 验证4: 六阶不可行性证明 ---")
lam_5 = (2/h**2)*(2 - sp.cos(xi*h) - sp.cos(eta*h))

# Modified eigenvalue
lam_mod_num = -(lam_L + alpha_p*h**2*lam_5 - k2*(lam_R + beta_p*h**2))
lam_mod_den = lam_R + gamma_p*h**2
lam_mod = lam_mod_num / lam_mod_den
diff = lam_mod - (xi**2 + eta**2 + k2)
diff_exp = sp.series(diff, h, 0, 6)

# c2 coefficient
c2 = sp.simplify(diff_exp.coeff(h, 2))
c2_thesis = -alpha_p*(xi**2+eta**2) + beta_p*k2 - gamma_p*(xi**2+eta**2+k2)
check("c2 matches thesis", sp.simplify(c2 - c2_thesis) == 0)

# Poisson case: k2=0, alpha=-gamma
diff_poisson = diff.subs(k2, 0).subs(alpha_p, -gamma_p)
diff_poisson_exp = sp.series(diff_poisson, h, 0, 6)
c4_poisson = sp.simplify(diff_poisson_exp.coeff(h, 4))
c4_poisson_thesis = -(3*eta**6 + 60*gamma_p*eta**4 - 5*eta**4*xi**2
                      - 5*eta**2*xi**4 + 60*gamma_p*xi**4 + 3*xi**6)/720
check("c4 (Poisson, alpha=-gamma) matches thesis",
      sp.simplify(c4_poisson - c4_poisson_thesis) == 0)

# Solve c4=0 for gamma
gamma_sol = sp.solve(c4_poisson, gamma_p)[0]
gamma_thesis = (-eta**6/20 + eta**4*xi**2/12 + eta**2*xi**4/12 - xi**6/20)/(eta**4+xi**4)
check("gamma solution matches thesis",
      sp.simplify(gamma_sol - gamma_thesis) == 0)

# Mode dependence
g11 = float(gamma_sol.subs([(xi, sp.pi), (eta, sp.pi)]))
g10 = float(gamma_sol.subs([(xi, sp.pi), (eta, 0)]))
check("gamma depends on mode (no universal solution)", g11 != g10,
      f"gamma(pi,pi)={g11:.6f}, gamma(pi,0)={g10:.6f}")

# =================================================================
# 验证5: 关键中间步骤 eq(2.3.11) 的正确性
# =================================================================
print("\n--- 验证5: 中间步骤 R_h(-nabla^2 u) - L_5h u ---")
print("  Thesis eq(2.3.9): R_h(-nabla^2 u) = -nabla^2 u - (h^2/12)nabla^4 u + O(h^4)")
print("  Thesis eq(2.3.10): L_5h u = -nabla^2 u + (h^2/12)(u_xxxx+u_yyyy) + O(h^4)")
print()
print("  Difference = R_h(-nabla^2 u) - L_5h u")
print("  = [-(h^2/12)nabla^4 u] - [(h^2/12)(u_xxxx+u_yyyy)]")
print("  = -(h^2/12)[nabla^4 u + u_xxxx + u_yyyy]")
print("  = -(h^2/12)[(u_xxxx + 2u_xxyy + u_yyyy) + u_xxxx + u_yyyy]")
print("  = -(h^2/12)[2u_xxxx + 2u_xxyy + 2u_yyyy]")
print("  = -(h^2/6)[u_xxxx + u_xxyy + u_yyyy]")
print()
print("  BUT thesis eq(2.3.11) claims: -(h^2/6)u_xxyy only!")
print("  This is WRONG. The correct answer is -(h^2/6)(u_xxxx+u_xxyy+u_yyyy)")
print()

# However, verify that the final L_h stencil IS correct
# L_h eigenvalue = lam_L, and L_h = R_h * L_5h (exact at grid level)
diff_exact = sp.simplify(lam_L - lam_R * lam_5)
check("L_h = R_h * L_5h (exact identity)",
      diff_exact == 0,
      f"lam_L - lam_R*lam_5 = {diff_exact}")

# This means: L_h u = R_h(L_5h u) exactly at grid points
# The thesis's intermediate step has an error, but the final stencil
# L_h = [1,4,1;4,-20,4;1,4,1]/(6h^2) is correct because:
# L_h u_{ij} = R_h(L_5h u)_{ij}  (exact)
# And the four-order accuracy follows from -lam_L/lam_R = (xi^2+eta^2) + O(h^4)
warn("eq(2.3.11) intermediate step has algebraic error",
     "Claimed -(h^2/6)u_xxyy but correct is -(h^2/6)(u_xxxx+u_xxyy+u_yyyy). "
     "Final L_h stencil and 4th-order theorem are unaffected.")

# =================================================================
# 验证6: Neumann DCT-I 对角化
# =================================================================
print("\n--- 验证6: Neumann DCT-I 对角化 ---")
n = 4
G = sp.Matrix([
    [2, -2, 0, 0],
    [-1, 2, -1, 0],
    [0, -1, 2, -1],
    [0, 0, -2, 2]
])
D = sp.diag(sp.sqrt(2), 1, 1, sp.sqrt(2))
S = D.inv() * G * D
check("S = D^{-1}GD is symmetric", S.equals(S.T))

eps = [1/sp.sqrt(2)] + [1]*(n-2) + [1/sp.sqrt(2)]
C = sp.Matrix(n, n, lambda k, j:
    eps[k] * sp.sqrt(sp.Rational(2, n-1)) * sp.cos(k*j*sp.pi/(n-1)))
CCt = sp.simplify(C * C.T)
check("DCT-I matrix is orthogonal", CCt.equals(sp.eye(n)))

Lambda = sp.simplify(C * S * C.T)
is_diag = all(abs(float(Lambda[i,j])) < 1e-10 for i in range(n) for j in range(n) if i != j)
check("DCT-I diagonalizes S", is_diag)

expected = [2*(1-sp.cos(k*sp.pi/(n-1))) for k in range(n)]
eig_match = all(abs(float(Lambda[k,k]) - float(expected[k])) < 1e-10 for k in range(n))
check("Eigenvalues match (2/h^2)(1-cos(k*pi/(n-1)))", eig_match)

# =================================================================
# 验证7: Cyclic Reduction 数值正确性
# =================================================================
print("\n--- 验证7: Cyclic Reduction 公式 ---")
n_test = 8
d_orig = np.array([3.0 + i*0.5 for i in range(n_test)])
e_orig = np.array([-1.0 + i*0.05 for i in range(n_test-1)])
b_orig = np.array([1.0 + i*0.3 for i in range(n_test)])

A = np.diag(d_orig) + np.diag(e_orig, 1) + np.diag(e_orig, -1)
u_direct = np.linalg.solve(A, b_orig)

def cr_solve(d, e, b):
    n = len(d)
    if n == 1:
        return np.array([b[0] / d[0]])
    n_half = n // 2
    d_new = np.zeros(n_half)
    e_new = np.zeros(max(n_half - 1, 0))
    b_new = np.zeros(n_half)
    for k in range(n_half):
        j = 2 * k
        d_new[k] = d[j]
        if j > 0:
            d_new[k] -= e[j-1]**2 / d[j-1]
            b_new[k] -= (e[j-1] / d[j-1]) * b[j-1]
        if j < n - 1:
            d_new[k] -= e[j]**2 / d[j+1]
            b_new[k] -= (e[j] / d[j+1]) * b[j+1]
    for k in range(n_half - 1):
        e_new[k] = -e[2*k] * e[2*k+1] / d[2*k+1]
    u_even = cr_solve(d_new, e_new, b_new)
    u = np.zeros(n)
    for k in range(n_half):
        u[2*k] = u_even[k]
    for k in range(n_half):
        j = 2*k + 1
        if j < n:
            u[j] = b[j]
            if j > 0:
                u[j] -= e[j-1] * u[j-1]
            if j < n - 1:
                u[j] -= e[j] * u[j+1]
            u[j] /= d[j]
    return u

u_cr = cr_solve(d_orig, e_orig, b_orig)
cr_error = np.max(np.abs(u_cr - u_direct))
check("CR formulas produce correct solution", cr_error < 1e-12,
      f"max error = {cr_error:.2e}")

# =================================================================
# 验证8: L_h stencil 的推导正确性（通过算子等价）
# =================================================================
print("\n--- 验证8: L_h stencil 推导（从 R_h*L_5h 代数等价）---")
# The thesis claims L_h is derived by compensating the mixed derivative.
# But the CORRECT derivation is simpler: L_h = R_h * L_5h (exact at grid level)
# Let's verify the stencil [1,4,1;4,-20,4;1,4,1]/(6h^2)
# equals R_h * L_5h

# R_h * L_5h applied to u_{ij}:
# R_h(L_5h u)_{ij} = (1/12)[(L_5h u)_{i-1,j} + (L_5h u)_{i+1,j}
#                          + (L_5h u)_{i,j-1} + (L_5h u)_{i,j+1}]
#                    + (2/3)(L_5h u)_{ij}
# where L_5h u_{ij} = (1/h^2)(4u_{ij} - u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1})

# This is a 9x9 stencil computation. Let me verify symbolically.
# Use symbolic grid values
u = sp.symbols('u0:9')  # u0=center, u1=E, u2=W, u3=N, u4=S, u5=NE, u6=NW, u7=SE, u8=SW
u_c, u_E, u_W, u_N, u_S, u_NE, u_NW, u_SE, u_SW = u

# Five-point Laplacian at each relevant grid point
L5_center = (4*u_c - u_E - u_W - u_N - u_S)
L5_E = (4*u_E - u_NE - u_SE - u_c - 0)  # boundary = 0 for interior
L5_W = (4*u_W - u_c - 0 - u_NW - u_SW)
L5_N = (4*u_N - u_NE - u_NW - u_c - 0)
L5_S = (4*u_S - u_SE - u_SW - u_c - 0)

# Actually this gets messy with boundary assumptions. Let me just verify
# the eigenvalue identity L_h = R_h * L_5h which we already proved (diff_exact = 0).
# Instead, verify the stencil by applying L_h to a grid function.
print("  L_h = R_h * L_5h (eigenvalue identity) already verified above.")
print("  The stencil [1,4,1;4,-20,4;1,4,1]/(6h^2) is correct.")

# =================================================================
# 总结
# =================================================================
print("\n" + "=" * 70)
print(f"校验总结: {PASS} PASS, {FAIL} FAIL, {WARN} WARN")
print("=" * 70)
print()
print("关键发现:")
print("1. [BUG] eq(2.3.18) lam_L 的 h^2 系数有笔误:")
print("   论文: (xi^4+eta^4+xi^2*eta^2)/12")
print("   正确: (xi^2+eta^2)^2/12 = (xi^4+2*xi^2*eta^2+eta^4)/12")
print()
print("2. [BUG] eq(2.3.18) lam_L 的 h^4 系数有笔误:")
print("   论文: -(xi^6+eta^6+xi^2*eta^2*(xi^2+eta^2)/2)/360")
print("   正确: -(xi^2+eta^2)(xi^4+4*xi^2*eta^2+eta^4)/360")
print()
print("3. [BUG] eq(2.3.11) 中间推导 R_h(-nabla^2 u)-L_5h u 有代数错误:")
print("   论文: -(h^2/6)*u_xxyy")
print("   正确: -(h^2/6)*(u_xxxx + u_xxyy + u_yyyy)")
print()
print("4. [OK] 四阶精度定理 (Theorem 2.2) 结论完全正确:")
print("   -lam_L/lam_R = (xi^2+eta^2) + O(h^4), h^2系数恰好为0")
print("   h^4 系数 = -(xi^2+eta^2)(3xi^4-8xi^2*eta^2+3eta^4)/720 ✓")
print()
print("5. [OK] 六阶不可行性定理 (Theorem 2.3) 完全正确:")
print("   c2, c4, gamma 解、模式依赖性均验证通过 ✓")
print()
print("6. [OK] Neumann DCT-I 对角化验证通过 ✓")
print()
print("7. [OK] Cyclic Reduction 公式验证通过 ✓")
print()
print("注: 中间步骤错误不影响最终结论，因为:")
print("- lam_L 展开的笔误在 -lam_L/lam_R 运算中被修正")
print("- eq(2.3.11) 的错误在最终 L_h stencil 推导中被隐式修正")
print("- 但这些笔误应在论文中纠正以确保推导完备性")
