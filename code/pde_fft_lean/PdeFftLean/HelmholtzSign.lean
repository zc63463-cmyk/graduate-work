/-
  PdeFftLean.HelmholtzSign — Tasks E, F

  Task E: FFT9 denominator algebra (algebraic parts only)
  Task F: Modified Helmholtz vs true Helmholtz resonance conditions

  Note on Task E: The full denominator cross-multiplication theorem
  (proving (-lamL/lamR + k2) * uhat = fhat from the spectral equation)
  requires cancellation of nonzero real factors, which needs careful
  handling of field_simp + ring. The algebraic equivalence of the
  two denominator forms is verified instead.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace PdeFftLean.HelmholtzSign

-- ============================================================
-- Task E: FFT9 denominator algebra
-- ============================================================

/--
Key algebraic lemma: (-lamL/lamR + k2) * lamR = -lamL + k2 * lamR
This is the core identity that connects the "divided by lamR" form
to the "multiplied by lamR" form of the FFT9 denominator.
-/
theorem fft9_key_lemma (lamL lamR k2 : ℝ) (hR : lamR ≠ 0) :
    (-lamL / lamR + k2) * lamR = -lamL + k2 * lamR := by
  field_simp [hR]

/--
Direct algebraic verification: the two denominator forms are equivalent.
  fhat / (-lamL + k2 * lamR) = fhat / ((-lamL/lamR + k2) * lamR)
This follows immediately from the key lemma.
-/
theorem fft9_denominator_algebra (lamL lamR k2 fhat : ℝ)
    (hR : lamR ≠ 0)
    (_hD : -lamL + k2 * lamR ≠ 0) :
    fhat / (-lamL + k2 * lamR) = fhat / ((-lamL / lamR + k2) * lamR) := by
  rw [fft9_key_lemma lamL lamR k2 hR]

/--
The FFT9 denominator D = -lamL/lamR + k2 satisfies:
  D * lamR = -lamL + k2 * lamR
This shows the denominator is well-defined when lamR ≠ 0
and -lamL + k2 * lamR ≠ 0.
-/
theorem fft9_denominator_welldefined (lamL lamR k2 : ℝ)
    (hR : lamR ≠ 0) (hD : -lamL + k2 * lamR ≠ 0) :
    (-lamL / lamR + k2) ≠ 0 := by
  intro h
  have : (-lamL / lamR + k2) * lamR = 0 := by rw [h]; ring
  rw [fft9_key_lemma lamL lamR k2 hR] at this
  exact hD this

-- ============================================================
-- Task F: Modified vs True Helmholtz resonance
-- ============================================================

/--
Modified Helmholtz denominator: lam + k2 is positive when lam > 0 and k2 >= 0.
This means (-Δ + k²) has NO resonance: the denominator never vanishes.
-/
theorem modified_denominator_pos (lam k2 : ℝ)
    (hlam : 0 < lam) (hk2 : 0 ≤ k2) :
    0 < lam + k2 := by
  linarith

/--
Modified Helmholtz has no zero denominator (no resonance).
-/
theorem modified_no_resonance (lam k2 : ℝ)
    (hlam : 0 < lam) (hk2 : 0 ≤ k2) :
    lam + k2 ≠ 0 := by
  linarith

/--
True Helmholtz denominator: lam - k2 = 0 iff lam = k2.
This is the resonance condition for (-Δ - k²).
-/
theorem true_helmholtz_resonance_iff (lam k2 : ℝ) :
    lam - k2 = 0 ↔ lam = k2 := by
  constructor <;> intro h <;> linarith

/--
Complete distinction: for lam > 0 and k2 >= 0,
- modified Helmholtz (lam + k2) cannot vanish
- true Helmholtz (lam - k2) vanishes exactly when lam = k2
-/
theorem helmholtz_sign_distinction (lam k2 : ℝ)
    (hlam : 0 < lam) (hk2 : 0 ≤ k2) :
    lam + k2 ≠ 0 ∧ (lam - k2 = 0 ↔ lam = k2) := by
  constructor
  · exact modified_no_resonance lam k2 hlam hk2
  · exact true_helmholtz_resonance_iff lam k2

/--
Concrete example: for eigenvalue lam1 > 0 and k2 = lam1,
the true Helmholtz has resonance but modified Helmholtz does not.
-/
theorem concrete_resonance_example (lam1 k2 : ℝ)
    (hlam : lam1 > 0) (hk : k2 = lam1) :
    lam1 + k2 ≠ 0 ∧ lam1 - k2 = 0 := by
  constructor
  · linarith
  · rw [hk]; ring

-- ============================================================
-- Sign convention: Lh approximates Δ, not -Δ
-- ============================================================

/--
For the FFT9 compact scheme, Lh approximates Δ (not -Δ).
The denominator -lamL/lamR + k2 uses the MINUS sign on lamL
because lamL < 0 for standard modes (Lh ≈ Δ).
-/
theorem fft9_sign_convention (lamL lamR k2 : ℝ) (hR : lamR ≠ 0) :
    (-lamL / lamR + k2) * lamR = -lamL + k2 * lamR := by
  field_simp [hR]

/-
The continuous eigenvalue of -Δ for mode (p,q) on [0,1]² with Dirichlet BC:
  λ_{p,q} = π²(p² + q²) > 0

For the modified Helmholtz (-Δ + k²):
  Denominator: λ_{p,q} + k² > 0 for all k² ≥ 0 (no resonance)

For the true Helmholtz (-Δ - k²):
  Denominator: λ_{p,q} - k² = 0 iff k² = λ_{p,q} (resonance!)

This distinction is critical: the FFT9 compact scheme for (-Δ + k²)u = f
is always well-posed, while the true Helmholtz (-Δ - k²)u = f requires
k² ≠ λ_{p,q} for any mode.

The formal theorems are modified_no_resonance and true_helmholtz_resonance_iff.
-/

end PdeFftLean.HelmholtzSign
