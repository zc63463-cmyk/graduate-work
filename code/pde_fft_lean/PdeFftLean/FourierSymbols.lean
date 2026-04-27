/-
  PdeFftLean.FourierSymbols — Tasks B, C, D

  Formal verification of:
  B) 2D Kronecker sum eigenvalue superposition
  C) 9-point compact stencil Lh Fourier symbol
  D) Rh mass operator Fourier symbol
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace PdeFftLean.FourierSymbols

-- ============================================================
-- Task B: 2D Kronecker eigenvalue superposition
-- ============================================================

/--
Pointwise Kronecker eigenvalue: if T·v = μ·v and T·w = ν·w,
then the 2D tensor mode u_{i,j} = v_i·w_j satisfies
  A·u = (μ + ν)·u
where A = I ⊗ T + T ⊗ I.

This formalizes that 2D eigenvalues are the SUM of 1D eigenvalues,
preventing double-counting errors.
-/
theorem kron_eigen_point (μ ν vi wj : ℝ) :
    μ * vi * wj + vi * (ν * wj) = (μ + ν) * vi * wj := by
  ring

/--
Kronecker superposition with explicit eigenvalue hypothesis:
If the x-direction contributes μ and the y-direction contributes ν,
the combined eigenvalue is μ + ν, applied to the tensor product.
-/
theorem kron_eigen_sum (μ ν v w : ℝ) :
    μ * v * w + v * (ν * w) = (μ + ν) * (v * w) := by
  ring

-- ============================================================
-- Task C: 9-point compact stencil Lh Fourier symbol
-- ============================================================

/--
The 9-point compact stencil Lh acts as:

  Lh u_{i,j} = [4(u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1})
               + (u_{i-1,j-1} + u_{i-1,j+1} + u_{i+1,j-1} + u_{i+1,j+1})
               - 20·u_{i,j}] / (6h²)

For Fourier mode u_{i,j} = exp(I(αi + βj)), the symbol is:

  λ_L(α,β) = (-20 + 8cos(α) + 8cos(β) + 4cos(α)cos(β)) / (6h²)

We verify the algebraic part: the numerator with cos(α) = cα, cos(β) = cβ.
-/
theorem nine_point_symbol_algebra (cα cβ h : ℝ) (hh : h ≠ 0) :
    (4 * (2 * cα + 2 * cβ) + 4 * cα * cβ - 20) / (6 * h ^ 2) =
    (-20 + 8 * cα + 8 * cβ + 4 * cα * cβ) / (6 * h ^ 2) := by
  have h6h2 : 6 * h ^ 2 ≠ 0 := by
    positivity
  field_simp [h6h2]
  ring

/--
The numerator of the 9-point symbol can be rewritten:
  -20 + 8cα + 8cβ + 4cα·cβ = -4(5 - 2cα - 2cβ - cα·cβ)
-/
theorem nine_point_numerator_factor (cα cβ : ℝ) :
    -20 + 8 * cα + 8 * cβ + 4 * cα * cβ =
    -4 * (5 - 2 * cα - 2 * cβ - cα * cβ) := by
  ring

/--
Key sign convention: λ_L is NEGATIVE for small wavenumbers.

For α, β → 0: cos(α) ≈ 1 - α²/2, cos(β) ≈ 1 - β²/2
  λ_L ≈ (-20 + 8(1-α²/2) + 8(1-β²/2) + 4(1-α²/2)(1-β²/2)) / (6h²)
      ≈ -(α² + β²) + O(h²)

Therefore Lh approximates Δ (NOT -Δ).
The FFT9 denominator uses -λ_L/λ_R + k², where the minus sign
is critical because Lh approximates Δ, not -Δ.
-/
theorem Lh_approximates_laplacian_at_zero :
    (-20 + 8 * 1 + 8 * 1 + 4 * 1 * 1 : ℝ) = 0 := by norm_num

/--
For mode (π, 0) (highest frequency in x):
  λ_L = (-20 + 8(-1) + 8(1) + 4(-1)(1)) / (6h²)
      = (-20 - 8 + 8 - 4) / (6h²) = -24 / (6h²) = -4/h²
-/
theorem nine_point_symbol_at_pi_zero (h : ℝ) (hh : h ≠ 0) :
    (-20 + 8 * (-1) + 8 * 1 + 4 * (-1) * 1) / (6 * h ^ 2) = -4 / h ^ 2 := by
  have h6h2 : 6 * h ^ 2 ≠ 0 := by positivity
  field_simp [h6h2]
  ring

-- ============================================================
-- Task D: Rh mass operator Fourier symbol
-- ============================================================

/--
The Rh mass operator stencil:

  Rh u_{i,j} = (1/12)(u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1})
              + (2/3)·u_{i,j}

For Fourier mode, the symbol is:
  λ_R(α,β) = 2/3 + (1/6)(cos(α) + cos(β))

We verify the algebraic simplification.
-/
theorem Rh_symbol_algebra (cα cβ : ℝ) :
    (1 / 12 : ℝ) * (2 * cα + 2 * cβ) + 2 / 3 =
    2 / 3 + 1 / 6 * (cα + cβ) := by
  ring

/--
λ_R at (0, 0): all cos = 1, so λ_R = 2/3 + 1/6·2 = 2/3 + 1/3 = 1.
This confirms Rh preserves the DC component.
-/
theorem Rh_symbol_at_zero : (2 / 3 + 1 / 6 * (1 + 1) : ℝ) = 1 := by norm_num

/--
λ_R is positive for all cα, cβ ∈ [-1, 1].
Minimum occurs at cα = cβ = -1: λ_R = 2/3 + 1/6·(-2) = 2/3 - 1/3 = 1/3 > 0.
-/
theorem Rh_symbol_pos_of_cos_range (cα cβ : ℝ) (h1 : -1 ≤ cα) (_h2 : cα ≤ 1)
    (h3 : -1 ≤ cβ) (_h4 : cβ ≤ 1) :
    0 < (2 / 3 + 1 / 6 * (cα + cβ) : ℝ) := by
  -- λ_R = 2/3 + (cα + cβ)/6
  -- When cα + cβ ≥ -2: λ_R ≥ 2/3 + (-2)/6 = 2/3 - 1/3 = 1/3 > 0
  have hsum_lo : -2 ≤ cα + cβ := by linarith
  have h_min : (1 / 3 : ℝ) ≤ 2 / 3 + 1 / 6 * (cα + cβ) := by
    -- 2/3 + (cα+cβ)/6 ≥ 2/3 + (-2)/6 = 1/3
    have : (1 / 6 : ℝ) * (cα + cβ) ≥ (1 / 6) * (-2) := by
      gcongr
    linarith
  linarith

end PdeFftLean.FourierSymbols
