/-
  SixthOrderImpossibility — Lean 4 Formalization (Verified)

  Scope: this file verifies the algebraic core used in the thesis proof for a
  restricted three-parameter nine-point correction family with mode-independent
  constants. It does not formalize PDE convergence, boundary closure, Fourier
  symbol derivation from stencils, or numerical implementation correctness.

  Algebraic obstruction: after the necessary c₂ conditions are imposed, the h⁴
  coefficient contains
    c₄(ξ,η;γ) = -(3η⁶ + 60γη⁴ - 5η⁴ξ² - 5η²ξ⁴ + 60γξ⁴ + 3ξ⁶)/720.
  The conditions c₄(1,0;γ)=0 and c₄(1,1;γ)=0 require incompatible values of
  the same mode-independent constant γ.

  This file uses only Lean 4 core (Int + omega), no Mathlib needed.
  All theorems are checked by the Lean 4 kernel.
-/

-- ============================================================
-- Section 1: c₄ numerator polynomial
-- ============================================================

/--
The numerator of 720·c₄(ξ,η;γ), evaluated over integers.
  N(ξ,η,γ) = 3η⁶ + 60γη⁴ - 5η⁴ξ² - 5η²ξ⁴ + 60γξ⁴ + 3ξ⁶
c₄ = 0 ⟺ N = 0 (since denominator 720 ≠ 0)
-/
def c4_num (ξ η γ : Int) : Int :=
  3 * η ^ 6 + 60 * γ * η ^ 4 - 5 * η ^ 4 * ξ ^ 2
  - 5 * η ^ 2 * ξ ^ 4 + 60 * γ * ξ ^ 4 + 3 * ξ ^ 6

/--
Evaluating c₄ at mode (1,0): N(1,0,γ) = 60γ + 3.
Key fact: γ must satisfy 60γ+3=0 for c₄(1,0)=0, forcing γ = -1/20.
-/
theorem c4_num_mode_10 (γ : Int) : c4_num 1 0 γ = 60 * γ + 3 := by
  unfold c4_num
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul,
             Int.mul_zero, Int.zero_mul, Int.add_zero, Int.zero_add]
  omega

/--
Evaluating c₄ at mode (1,1): N(1,1,γ) = 120γ - 4.
Key fact: γ must satisfy 120γ-4=0 for c₄(1,1)=0, forcing γ = 1/30.
-/
theorem c4_num_mode_11 (γ : Int) : c4_num 1 1 γ = 120 * γ - 4 := by
  unfold c4_num
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul]
  omega

/--
Evaluating c₄ at mode (0,1): N(0,1,γ) = 60γ + 3 (symmetric to (1,0)).
-/
theorem c4_num_mode_01 (γ : Int) : c4_num 0 1 γ = 60 * γ + 3 := by
  unfold c4_num
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul,
             Int.mul_zero, Int.zero_mul, Int.add_zero, Int.zero_add]
  omega

/--
Evaluating c₄ at mode (2,1): N(2,1,γ) = 1020γ + 95.
Third independent constraint confirming no universal γ exists.
-/
theorem c4_num_mode_21 (γ : Int) : c4_num 2 1 γ = 1020 * γ + 95 := by
  unfold c4_num
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul]
  omega

-- ============================================================
-- Section 2: c₂ coefficient analysis (Poisson case)
-- ============================================================

/--
The h² coefficient of the modified format's truncation error.
  c₂(α,ξ,η,β,k²,γ) = -α(ξ²+η²) + β·k² - γ(ξ²+η²+k²)
-/
def c2_coeff (α ξ η β k2 γ : Int) : Int :=
  -α * (ξ^2 + η^2) + β * k2 - γ * (ξ^2 + η^2 + k2)

/--
c₂ at mode (1,0) with k²=0: c₂ = -(α+γ).
-/
theorem c2_at_10 (α γ : Int) : c2_coeff α 1 0 0 0 γ = -(α + γ) := by
  unfold c2_coeff
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul,
             Int.mul_zero, Int.zero_mul, Int.add_zero, Int.zero_add]
  omega

/--
c₂ at mode (0,1) with k²=0: c₂ = -(α+γ).
-/
theorem c2_at_01 (α γ : Int) : c2_coeff α 0 1 0 0 γ = -(α + γ) := by
  unfold c2_coeff
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul,
             Int.mul_zero, Int.zero_mul, Int.add_zero, Int.zero_add]
  omega

/--
c₂ at mode (1,1) with k²=0: c₂ = -2(α+γ).
-/
theorem c2_at_11 (α γ : Int) : c2_coeff α 1 1 0 0 γ = -2 * (α + γ) := by
  unfold c2_coeff
  simp only [Int.pow_succ, Int.pow_zero, Int.mul_one, Int.one_mul]
  omega

/--
If c₂ = 0 for mode (1,0) in the Poisson case, then α = -γ.
(A single mode suffices since c₂ = -(α+γ) is independent of (ξ,η).)
-/
theorem c2_poisson_forces_alpha_neg_gamma (α γ : Int)
    (h10 : c2_coeff α 1 0 0 0 γ = 0) :
    α + γ = 0 := by
  rw [c2_at_10] at h10
  omega

-- ============================================================
-- Section 3: Main impossibility theorems
-- ============================================================

/--
No integer γ satisfies both 60γ+3=0 and 120γ-4=0.
Direct: 2·(60γ+3) - (120γ-4) = 10 ≠ 0.
-/
theorem no_integer_gamma :
    ¬∃ (γ : Int), 60 * γ + 3 = 0 ∧ 120 * γ - 4 = 0 := by
  intro ⟨γ, h1, h2⟩; omega

/--
Core impossibility: the linear system
  60·a + 3·d = 0  and  120·a - 4·d = 0
has no solution with d ≠ 0 (clearing denominators for rational γ).
-/
theorem no_rational_gamma_core :
    ¬∃ (a d : Int), d ≠ 0 ∧ 60 * a + 3 * d = 0 ∧ 120 * a - 4 * d = 0 := by
  intro ⟨a, d, hd, h1, h2⟩
  have h10d : 10 * d = 0 := by omega
  -- d ≠ 0 and 10*d = 0 is impossible over Int
  have hd12 : d ≥ 1 ∨ d ≤ -1 := by omega
  cases hd12 with
  | inl h => omega
  | inr h => omega

/--
No rational γ = a/d with d > 0 makes c₄ zero for all modes.
-/
theorem no_rational_gamma :
    ¬∃ (a d : Int), d > 0 ∧
      60 * a + 3 * d = 0 ∧ 120 * a - 4 * d = 0 := by
  intro ⟨a, d, hd, h1, h2⟩
  have h10d : 10 * d = 0 := by omega
  omega

/--
★ Algebraic core theorem ★ No integer γ makes c₄(ξ,η,γ) = 0 for all test modes.

This supports the restricted local Fourier-symbol obstruction for the
three-parameter correction family; it is not a full PDE convergence theorem.

Proof:
  c₄(1,0,γ) = 0  requires  60γ + 3 = 0
  c₄(1,1,γ) = 0  requires  120γ - 4 = 0
  2·(60γ+3) - (120γ-4) = 10 ≠ 0, contradiction.
-/
theorem sixth_order_impossible :
    ¬∃ (γ : Int), ∀ (ξ η : Int), c4_num ξ η γ = 0 := by
  intro ⟨γ, hγ⟩
  have h1 := hγ 1 0
  have h2 := hγ 1 1
  rw [c4_num_mode_10] at h1
  rw [c4_num_mode_11] at h2
  omega

/--
Complete Poisson impossibility: c₂=0 forces α=-γ, then c₄=0 is impossible.
-/
theorem poisson_no_sixth_order :
    ¬∃ (α γ : Int), α + γ = 0 ∧
      ∀ (ξ η : Int), c4_num ξ η γ = 0 := by
  intro ⟨α, γ, hαγ, hc4⟩
  exact sixth_order_impossible ⟨γ, hc4⟩

/--
Even with c₂ correctly eliminated (α=-γ), c₄ cannot be made zero for all modes.
This is the algebraic obstruction for the restricted three-parameter family.
-/
theorem c4_is_the_obstruction (γ : Int)
    (h : ∀ (ξ η : Int), c4_num ξ η γ = 0) : False :=
  sixth_order_impossible ⟨γ, h⟩

-- ============================================================
-- Section 4: Full formalization outline (requires Mathlib)
-- ============================================================

/-
With Mathlib (import Mathlib.Data.Real.Basic, Mathlib.Tactic),
the full formalization over ℝ would establish:

  def c4_num (ξ η γ : ℝ) : ℝ :=
    3 * η ^ 6 + 60 * γ * η ^ 4 - 5 * η ^ 4 * ξ ^ 2
    - 5 * η ^ 2 * ξ ^ 4 + 60 * γ * ξ ^ 4 + 3 * ξ ^ 6

  theorem sixth_order_impossible_poisson :
      ¬∃ (γ : ℝ), ∀ (ξ η : ℝ), c4_num ξ η γ = 0 := by
    intro ⟨γ, hγ⟩
    have h1 := hγ 1 0; have h2 := hγ 1 1
    simp [c4_num] at h1 h2; ring_nf at h1 h2; linarith

  theorem sixth_order_impossible_helmholtz (k2 : ℝ) (hk2 : k2 ≠ 0) :
      ¬∃ (α β γ : ℝ),
        (∀ (ξ η : ℝ), -α*(ξ^2+η^2) + β*k2 - γ*(ξ^2+ena^2+k2) = 0) ∧
        (∀ (ξ η : ℝ), c4_num ξ η γ = 0) := by
    intro ⟨α, β, γ, hc2, hc4⟩
    have : α = -γ := by linarith [hc2 1 0, hc2 1 1]
    have : β = γ := by nlinarith [hc2 1 0, this]
    exact sixth_order_impossible_poisson ⟨γ, hc4⟩

Key mathematical insight: once the c₂=0 necessary relations are imposed, the
same c₄ polynomial obstruction applies to the restricted three-parameter family
in both Poisson and Helmholtz settings.
-/
