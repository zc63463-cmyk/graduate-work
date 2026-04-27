-- Mathlib availability test for SixthOrderImpossibility project
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

-- ✅ Test 1: Real number basic
example : (0 : ℝ) ≤ 1 := by linarith

-- ✅ Test 2: ring tactic
example (x : ℝ) : (x + 1) ^ 2 = x ^ 2 + 2 * x + 1 := by ring

-- ✅ Test 3: norm_num with rationals
example : (1 / 20 : ℝ) = 0.05 := by norm_num

-- ✅ Test 4: Simple linear arithmetic over ℝ
example (γ : ℝ) (h : 60 * γ + 3 = 0) : γ = -3 / 60 := by linarith

-- ✅ Test 5: The core impossibility — 60γ+3=0 and 120γ-4=0 cannot hold simultaneously
-- Need nlinarith (nonlinear) or manual elimination
example (γ : ℝ) (h1 : 60 * γ + 3 = 0) (h2 : 120 * γ - 4 = 0) : False := by
  have : 2 * (60 * γ + 3) = (120 * γ - 4) + 10 := by ring
  linarith

-- ✅ Test 6: c₄ numerator polynomial over ℝ
def c4_num (ξ η γ : ℝ) : ℝ :=
  3 * η ^ 6 + 60 * γ * η ^ 4 - 5 * η ^ 4 * ξ ^ 2
  - 5 * η ^ 2 * ξ ^ 4 + 60 * γ * ξ ^ 4 + 3 * ξ ^ 6

-- c₄(1,0,γ) = 60γ + 3
theorem c4_num_mode_10_real (γ : ℝ) : c4_num 1 0 γ = 60 * γ + 3 := by
  unfold c4_num
  ring

-- c₄(1,1,γ) = 120γ - 4
theorem c4_num_mode_11_real (γ : ℝ) : c4_num 1 1 γ = 120 * γ - 4 := by
  unfold c4_num
  ring

-- ★ Main theorem over ℝ: no real γ makes c₄ zero for all modes
theorem sixth_order_impossible_real :
    ¬∃ (γ : ℝ), ∀ (ξ η : ℝ), c4_num ξ η γ = 0 := by
  intro ⟨γ, hγ⟩
  have h1 := hγ 1 0
  have h2 := hγ 1 1
  rw [c4_num_mode_10_real] at h1
  rw [c4_num_mode_11_real] at h2
  -- 2*(60γ+3) - (120γ-4) = 10 ≠ 0
  have : 2 * (60 * γ + 3) - (120 * γ - 4) = 10 := by ring
  linarith

-- ★ Poisson impossibility over ℝ
def c2_coeff (α ξ η β k2 γ : ℝ) : ℝ :=
  -α * (ξ^2 + η^2) + β * k2 - γ * (ξ^2 + η^2 + k2)

-- c₂=0 at (1,0) with k²=0 forces α = -γ
theorem c2_poisson_forces_neg_gamma (α γ : ℝ)
    (h : c2_coeff α 1 0 0 0 γ = 0) : α = -γ := by
  unfold c2_coeff at h
  linarith

-- Complete Poisson impossibility
theorem poisson_no_sixth_order_real :
    ¬∃ (α γ : ℝ), α + γ = 0 ∧
      ∀ (ξ η : ℝ), c4_num ξ η γ = 0 := by
  intro ⟨α, γ, hαγ, hc4⟩
  exact sixth_order_impossible_real ⟨γ, hc4⟩

-- ============================================================
-- Section 2: Helmholtz impossibility (k² ≠ 0)
-- ============================================================

-- c₂ at mode (1,0) with general k²: c₂ = -α - γ(1 + k²) + β·k²
-- For c₂=0 at (1,0): α = -γ(1+k²) + β·k²
-- c₂ at mode (1,1) with general k²: c₂ = -2α - γ(2 + k²) + β·k²
-- For c₂=0 at (1,1): 2α = -γ(2+k²) + β·k²

-- Key lemma: c₂=0 at both (1,0) and (1,1) forces α = -γ
-- (same as Poisson case, independent of k²)
theorem c2_helmholtz_forces_neg_gamma (α β γ k2 : ℝ)
    (h10 : c2_coeff α 1 0 β k2 γ = 0)
    (h11 : c2_coeff α 1 1 β k2 γ = 0) :
    α + γ = 0 := by
  unfold c2_coeff at h10 h11
  -- h10: -α + β·k2 - γ(1 + k2) = 0  =>  α = -γ - γ·k2 + β·k2
  -- h11: -2α + β·k2 - γ(2 + k2) = 0 =>  2α = -2γ - γ·k2 + β·k2
  -- Subtracting: α = -γ
  have h_diff : 2 * (-α + β * k2 - γ * (1 + k2))
      - (-2 * α + β * k2 - γ * (2 + k2)) = 0 := by linarith
  ring_nf at h_diff
  linarith

-- ★ Helmholtz impossibility: no (α, β, γ) makes both c₂=0 and c₄=0 for all modes
-- Note: the impossibility holds for ALL k², including k²=0 (Poisson case)
theorem helmholtz_no_sixth_order_real (k2 : ℝ) :
    ¬∃ (α β γ : ℝ),
      (∀ (ξ η : ℝ), c2_coeff α ξ η β k2 γ = 0) ∧
      (∀ (ξ η : ℝ), c4_num ξ η γ = 0) := by
  intro ⟨α, β, γ, hc2, hc4⟩
  -- c₂=0 forces α = -γ (from both modes)
  have hαγ : α + γ = 0 := c2_helmholtz_forces_neg_gamma α β γ k2 (hc2 1 0) (hc2 1 1)
  -- Then c₄ is impossible (same as Poisson)
  exact sixth_order_impossible_real ⟨γ, hc4⟩

-- ============================================================
-- Section 3: c₂ mode-dependence proof (Helmholtz case)
-- ============================================================

-- Even for a single parameter γ, c₂=0 cannot hold for all modes when k²≠0
-- because c₂(1,0) = -α - γ(1+k²) + β·k² depends on k²,
-- while c₂(1,1) = -2α - γ(2+k²) + β·k² also depends on k²
-- But the key constraint is c₂ at modes (1,0) and (0,1) force α=-γ,
-- and then c₂ at mode (1,1) gives -2γ + β·k² - γ(2+k²) = β·k² - γ(4+k²)
-- For this to vanish for ALL modes, we'd need β and γ to be mode-dependent

-- c₂ at mode (1,0) requires α = -γ + (β-γ)·k² (when k²≠0, depends on k²)
-- This means α must change with k², confirming mode-dependence
theorem c2_mode_dependent_helmholtz (α β γ k2 : ℝ)
    (h10 : c2_coeff α 1 0 β k2 γ = 0) :
    α = -γ + (β - γ) * k2 := by
  unfold c2_coeff at h10
  linarith

-- c₂ at mode (0,1) gives same constraint (by symmetry)
theorem c2_mode_01_helmholtz (α β γ k2 : ℝ)
    (h01 : c2_coeff α 0 1 β k2 γ = 0) :
    α = -γ + (β - γ) * k2 := by
  unfold c2_coeff at h01
  linarith

-- c₂ at mode (2,0) requires -4α + β·k² - γ(4+k²) = 0
-- Combined with (1,0): -α + β·k² - γ(1+k²) = 0
-- Subtracting: -3α - 3γ = 0, i.e., α = -γ
-- This is consistent with the Poisson result
theorem c2_helmholtz_mode_20_forces_neg_gamma (α β γ k2 : ℝ)
    (h10 : c2_coeff α 1 0 β k2 γ = 0)
    (h20 : c2_coeff α 2 0 β k2 γ = 0) :
    α + γ = 0 := by
  unfold c2_coeff at h10 h20
  -- h10: -α + β·k2 - γ(1+k2) = 0
  -- h20: -4α + β·k2 - γ(4+k2) = 0
  -- Subtract: -3α - 3γ = 0
  have h_sub : (-4 * α + β * k2 - γ * (4 + k2))
      - (-α + β * k2 - γ * (1 + k2)) = 0 := by linarith
  ring_nf at h_sub
  linarith
