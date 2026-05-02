# Optional Patch Plan

## Patch Principles

- Keep changes mathematically minimal.
- Do not add experiments.
- Do not restore old figures.
- Do not change the thesis main line.
- Do not modify Python or Lean code unless a later paper-code audit finds implementation mismatch.
- Keep Lean 4 scoped to algebraic verification, not a full PDE proof.

## Proposed Edits

| Priority | Issue ID | File | Edit Type | Minimal Change | Rerun Needed |
|---|---|---|---|---|---|
| 1 | A-001 | `thesis/chapters/2_math_preliminary.tex` | Formula sign convention | In the FFT9 derivation, keep `L_h approx \Delta`; replace `L_h u = R_h(-\nabla^2u)+O(h^4)` with `-L_h u = R_h(-\nabla^2u)+O(h^4)` and replace `L_h = L_h^{(5)} + h^2 C_h/6` with `-L_h = L_h^{(5)} - h^2 C_h/6` when `L_h^{(5)}` denotes the positive five-point `-\Delta` operator. | No |
| 2 | A-002 | `thesis/chapters/2_math_preliminary.tex` | Formula sign | Replace the three-parameter effective eigenvalue numerator with `-\hat{\lambda}_L - \alpha h^2\hat{\lambda}_5 + \sigma(\hat{\lambda}_R+\beta h^2)`. | No |
| 3 | A-003 | `thesis/chapters/appendix_lean4.tex` | Appendix formula alignment | Replace the effective-eigenvalue formula with the corrected three-parameter fraction, including `\beta`, or mark the current expression as schematic and defer to Chapter 2. | No |
| 4 | A-004 | `thesis/chapters/4_helmholtz.tex` | Theorem condition | Change the indefiniteness theorem to state: SPD for `\kappa^2<\lambda_{\min}`, singular at `\kappa^2=\lambda_{p,q}`, indefinite for `\lambda_{\min}<\kappa^2<\lambda_{\max}`, and negative definite for `\kappa^2>\lambda_{\max}`. | No |
| 5 | A-004 | `thesis/main.tex`, `thesis/chapters/1_introduction.tex`, `thesis/chapters/2_math_preliminary.tex`, `thesis/chapters/7_conclusion.tex` | Scope wording | Replace broad "true Helmholtz 矩阵不定" wording with "可能不定/在离散谱内移位时不定/存在近共振风险" where the sentence is not explicitly restricted to the indefinite parameter range. | No |
| 6 | B-001 | `thesis/chapters/4_helmholtz.tex` | Local sign | If keeping `\kappa^2=\lambda_{p^*,q^*}+\delta`, change the target denominator to `d_{p^*,q^*}=-\delta`; otherwise define `\kappa^2=\lambda_{p^*,q^*}-\delta`. | No |
| 7 | C-001 | `thesis/chapters/2_math_preliminary.tex` | Wording precision | Replace "对任意非零模式 ... 不为零" with "该系数不恒为零；例如取 `\xi=\eta` 时非零". | No |
| 8 | C-002 | `thesis/chapters/3_fft_direct.tex` | Complexity wording | Replace "第三项趋于常数" with "第三项降为 `O(N^2)`，总复杂度由 `O(N^2 log log N)` 项主导" or equivalent wording. | No |

## Deferred / Review First

- No generated patch requires GPT-5.5 Pro review before drafting.
- If the advisor wants the thesis to keep the shorthand "true Helmholtz 不定", add a parenthetical qualification rather than removing the warning: "在本文关注的穿越离散谱/近共振参数区间可能不定".

## Code / Experiment Impact

- Code changes needed: No.
- Experiment reruns needed: No.
- Figures or CSV changes needed: No.
