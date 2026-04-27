/-
  PdeFftLean.FivePoint — Task A

  Formal verification of the 1D Dirichlet Laplacian eigenvalue formula.
  
  For the tridiagonal matrix T = tridiag(-1, 2, -1), the sine modes
  v_i = sin(i·θ) are eigenvectors with eigenvalue 2 - 2·cos(θ).

  Core identity:
    -sin(x - θ) + 2·sin(x) - sin(x + θ) = (2 - 2·cos(θ))·sin(x)

  This prevents double-counting the diagonal entry as 4 in 2D.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace PdeFftLean.FivePoint

-- ============================================================
-- Section 1: Core trigonometric eigenvalue identity
-- ============================================================

/--
Core identity: -sin(x-θ) + 2·sin(x) - sin(x+θ) = (2 - 2·cos(θ))·sin(x)

This is the fundamental lemma that the 1D discrete Laplacian T = tridiag(-1,2,-1)
has eigenvalues 2 - 2·cos(θ) with eigenvectors sin(i·θ).

Proof strategy:
  sin(x-θ) = sin(x)cos(θ) - cos(x)sin(θ)
  sin(x+θ) = sin(x)cos(θ) + cos(x)sin(θ)
  LHS = -sin(x)cos(θ) + cos(x)sin(θ) + 2sin(x) - sin(x)cos(θ) - cos(x)sin(θ)
      = 2sin(x) - 2sin(x)cos(θ)
      = (2 - 2cos(θ))sin(x)
-/
theorem sine_eigen_core_real (x θ : ℝ) :
    -Real.sin (x - θ) + 2 * Real.sin x - Real.sin (x + θ) =
    (2 - 2 * Real.cos θ) * Real.sin x := by
  -- Expand sin(x-θ) and sin(x+θ) using addition formulas
  rw [Real.sin_sub, Real.sin_add]
  -- Now we have:
  -- -(sin x cos θ - cos x sin θ) + 2 sin x - (sin x cos θ + cos x sin θ)
  -- = -sin x cos θ + cos x sin θ + 2 sin x - sin x cos θ - cos x sin θ
  -- = 2 sin x - 2 sin x cos θ
  -- = (2 - 2 cos θ) sin x
  ring_nf

-- ============================================================
-- Section 2: Eigenvalue as 2(1 - cos θ) = 4 sin²(θ/2)
-- ============================================================

/--
The eigenvalue 2 - 2·cos(θ) equals 2·(1 - cos(θ)).
This is a trivial algebraic restatement.
-/
theorem eigenval_two_one_minus_cos (θ : ℝ) :
    2 - 2 * Real.cos θ = 2 * (1 - Real.cos θ) := by ring

/--
The eigenvalue 2 - 2·cos(θ) is non-negative for all θ.
This follows from cos(θ) ≤ 1, so 2 - 2cos(θ) ≥ 0.
-/
theorem eigenval_nonneg (θ : ℝ) :
    0 ≤ 2 - 2 * Real.cos θ := by
  have : Real.cos θ ≤ 1 := Real.cos_le_one θ
  linarith

/--
The eigenvalue 2 - 2·cos(θ) is positive when cos(θ) < 1.
This covers all Fourier modes θ = pπ/(N+1) with p = 1,...,N,
where cos(θ) < 1 strictly.
-/
theorem eigenval_pos_of_cos_lt_one (θ : ℝ) (h : Real.cos θ < 1) :
    0 < 2 - 2 * Real.cos θ := by linarith

-- ============================================================
-- Section 3: Preventing double-counting error
-- ============================================================

/--
The 1D eigenvalue is 2 - 2·cos(θ), NOT 4.
For θ = π/2 (the midpoint mode), eigenvalue = 2, not 4.
This lemma prevents the common error of writing the 2D eigenvalue
as 4 + 4 = 8 instead of the correct (2-2cosθ) + (2-2cosη).
-/
theorem eigenval_not_four_at_pi_half :
    2 - 2 * Real.cos (Real.pi / 2) = 2 := by
  rw [Real.cos_pi_div_two]
  ring

end PdeFftLean.FivePoint
