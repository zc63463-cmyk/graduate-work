#!/usr/bin/env python3
"""
Verify FFT9 Fourier symbol expansion using SymPy.
Unified notation: use continuous wavenumbers xi=p*pi, eta=q*pi,
 then alpha = xi*h, beta = eta*h.

Goal: verify that the h^2 term cancels in -lambda_L / lambda_R,
giving O(h^4) interior Fourier symbol consistency.
"""
import sympy as sp

# Define symbols
xi, eta, h = sp.symbols('xi eta h', real=True)
alpha = xi * h
beta = eta * h

# ============================================================
# L_h Fourier symbol (9-point stencil)
# lambda_L = [-20 + 8(cos alpha + cos beta) + 4 cos alpha cos beta] / (6*h^2)
# Note: L_h approximates Delta, so lambda_L's leading term is negative.
# ============================================================
cos_a = sp.cos(alpha)
cos_b = sp.cos(beta)
lambda_L = (-20 + 8*(cos_a + cos_b) + 4*cos_a*cos_b) / (6*h**2)

# Expand cos(alpha) and cos(beta) to O(h^6)
# cos(xi*h) = 1 - (xi*h)^2/2 + (xi*h)^4/24 - (xi*h)^6/720 + O(h^8)
cos_a_series = sp.series(cos_a, h, 0, 8).removeO()
cos_b_series = sp.series(cos_b, h, 0, 8).removeO()

lambda_L_sub = lambda_L.subs({cos_a: cos_a_series, cos_b: cos_b_series})
lambda_L_expanded = sp.simplify(sp.series(lambda_L_sub, h, 0, 7).removeO())

print("=============================================")
print("L_h Fourier symbol lambda_L (expanded to O(h^6)):")
print("lambda_L = ")
sp.pprint(lambda_L_expanded)
print()

# -lambda_L should have leading term +(xi^2 + eta^2).
neg_lambda_L = -lambda_L_sub
neg_lambda_L_expanded = sp.simplify(sp.series(neg_lambda_L, h, 0, 7).removeO())

print("-lambda_L = ")
sp.pprint(neg_lambda_L_expanded)
print()

# Extract coefficients
print("Coefficient of h^0 (should be xi^2 + eta^2):")
coeff_h0 = sp.simplify(neg_lambda_L_expanded.coeff(h, 0))
print("  ", coeff_h0)
print("  Match: ", sp.simplify(coeff_h0 - (xi**2 + eta**2)) == 0)
print()

print("Coefficient of h^2 (should be -(xi^2 + eta^2)^2 / 12):")
coeff_h2 = sp.simplify(neg_lambda_L_expanded.coeff(h, 2))
print("  ", coeff_h2)
# Verify
target_h2 = -(xi**2 + eta**2)**2 / 12
print("  Target: ", target_h2)
print("  Match: ", sp.simplify(coeff_h2 - target_h2) == 0)
print()

print("Coefficient of h^4 (should be (xi^2+eta^2)(xi^4+4xi^2*eta^2+eta^4)/360):")
coeff_h4 = sp.simplify(neg_lambda_L_expanded.coeff(h, 4))
print("  ", coeff_h4)
target_h4 = (xi**2 + eta**2)*(xi**4 + 4*xi**2*eta**2 + eta**4) / 360
print("  Target: ", target_h4)
print("  Match: ", sp.simplify(coeff_h4 - target_h4) == 0)
print()

# ============================================================
# R_h Fourier symbol
# lambda_R = 2/3 + 1/6*(cos alpha + cos beta)
# Actually: R_h weights [0,1,0; 1,8,1; 0,1,0]/12
# Fourier symbol: (1/12)*(4(cos alpha + cos beta) + 4*1)  <- Wait.
# Correct: R_h = (1/12) * [u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1}] + (8/12) u_{i,j}
# Fourier symbol: (1/12)*(2*cos(alpha) + 2*cos(beta)) + 8/12
#                 = (1/6)*(cos alpha + cos beta) + 2/3
# ============================================================
lambda_R = (sp.Rational(1,6))*(cos_a + cos_b) + sp.Rational(2,3)

lambda_R_sub = lambda_R.subs({cos_a: cos_a_series, cos_b: cos_b_series})
lambda_R_expanded = sp.simplify(sp.series(lambda_R_sub, h, 0, 7).removeO())

print("=============================================")
print("R_h Fourier symbol lambda_R (expanded to O(h^6)):")
print("lambda_R = ")
sp.pprint(lambda_R_expanded)
print()

# Verify expansion
print("lambda_R = 1 - (xi^2+eta^2)*h^2/12 + (xi^4+eta^4)*h^4/144 + O(h^6)")
print("Coefficient of h^0 (should be 1):")
coeff_r_h0 = sp.simplify(lambda_R_expanded.coeff(h, 0))
print("  ", coeff_r_h0)
print()

print("Coefficient of h^2 (should be -(xi^2+eta^2)/12):")
coeff_r_h2 = sp.simplify(lambda_R_expanded.coeff(h, 2))
print("  ", coeff_r_h2)
target_r_h2 = -(xi**2 + eta**2) / 12
print("  Target: ", target_r_h2)
print("  Match: ", sp.simplify(coeff_r_h2 - target_r_h2) == 0)
print()

# ============================================================
# Effective eigenvalue: -lambda_L / lambda_R
# ============================================================
eff = -lambda_L_sub / lambda_R_sub
eff_expanded = sp.simplify(sp.series(eff, h, 0, 7).removeO())

print("=============================================")
print("Effective eigenvalue -lambda_L / lambda_R (expanded):")
sp.pprint(eff_expanded)
print()

# Check if h^2 term cancels
coeff_eff_h2 = sp.simplify(eff_expanded.coeff(h, 2))
print("Coefficient of h^2 in effective eigenvalue (should be 0 for O(h^4) symbol consistency):")
print("  ", coeff_eff_h2)
print("  Is zero?", sp.simplify(coeff_eff_h2) == 0)
print()

# The h^0 term should be xi^2 + eta^2
coeff_eff_h0 = sp.simplify(eff_expanded.coeff(h, 0))
print("Coefficient of h^0 (should be xi^2 + eta^2):")
print("  ", coeff_eff_h0)
print("  Match: ", sp.simplify(coeff_eff_h0 - (xi**2 + eta**2)) == 0)
print()

# The h^4 term should be non-zero (this is the truncation error)
coeff_eff_h4 = sp.simplify(eff_expanded.coeff(h, 4))
print("Coefficient of h^4 (truncation error term, should be non-zero):")
print("  ", coeff_eff_h4)
print()

print("=============================================")
print("CONCLUSION:")
print("  -lambda_L = xi^2+eta^2 - ((xi^2+eta^2)^2/12) h^2 + ...")
print("  lambda_R = 1 - (xi^2+eta^2)h^2/12 + ...")
print("  -lambda_L / lambda_R = xi^2+eta^2 + O(h^4)")
print("  => FFT9 has fourth-order interior Fourier symbol consistency.")
print("=============================================")
