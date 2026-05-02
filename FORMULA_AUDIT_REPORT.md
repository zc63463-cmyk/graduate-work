# Formula Audit Report

## Scope

- Files audited:
  - `thesis/main.tex`
  - `thesis/chapters/1_introduction.tex`
  - `thesis/chapters/2_math_preliminary.tex`
  - `thesis/chapters/3_fft_direct.tex`
  - `thesis/chapters/4_helmholtz.tex`
  - `thesis/chapters/5_gmres.tex`
  - `thesis/chapters/6_experiments.tex`
  - `thesis/chapters/7_conclusion.tex`
  - `thesis/chapters/appendix_lean4.tex`
  - `README.md`
  - `FINAL_POLISH_REPORT.md`
- Files not audited:
  - Python/Lean source files were not part of this requested formula-only input set.
- Assumptions:
  - `L_h` in FFT9 is intended to approximate `\Delta`, as explicitly stated in Chapter 2 and used by the raw denominator `-\hat{\lambda}_L+\sigma\hat{\lambda}_R`.
  - `\hat{\lambda}_5 = 2h^{-2}(2-\cos\alpha-\cos\beta)` is the positive five-point `-\Delta` symbol.

## Executive Summary

- A issues: 4
- B issues: 1
- C issues: 2
- D confirmations: 9

The unified PDE framework is mostly stable: `(-\nabla^2+\sigma)u=f`, modified Helmholtz as `\sigma=+\kappa^2`, true Helmholtz as `\sigma=-\kappa^2`, five-point denominators, nonhomogeneous Dirichlet signs, FFT9 raw denominator, Neumann DCT-I handling, GMRES scope, FACR-like complexity boundary, and Lean 4 scope boundaries are all broadly aligned.

The main remaining risks are local but important: Chapter 2 still has FFT9 derivation sign drift around `L_h` and the three-parameter effective eigenvalue; the Lean appendix repeats a mismatched effective-eigenvalue formula; Chapter 4 states true Helmholtz indefiniteness with an insufficient condition.

## Findings

### A-001: FFT9 `L_h` Sign Drift in the Derivation

- Severity: A
- File: `thesis/chapters/2_math_preliminary.tex:269`, `thesis/chapters/2_math_preliminary.tex:317`, `thesis/chapters/2_math_preliminary.tex:324`, `thesis/chapters/2_math_preliminary.tex:341`
- Current writing:
  - `L_h^{(5)}u=-\nabla^2u+...`
  - `L_h u = R_h(-\nabla^2u)+O(h^4)`
  - `L_h = L_h^{(5)} + h^2 C_h/6`
- Why risky:
  - The final FFT9 scheme correctly uses `L_h approx \Delta` and `(-L_h+\sigma R_h)u=R_h f`.
  - The intermediate derivation instead treats `L_h` as if it approximated `-\Delta`, which contradicts the stencil and the later Fourier symbol `\hat{\lambda}_L<0`.
- Correct version:
  - If `L_h^{(5)}` denotes the positive five-point negative Laplacian, then
    `-L_h = L_h^{(5)} - h^2 C_h/6`
    and
    `-L_h u = R_h(-\nabla^2u)+O(h^4)`.
  - Equivalently, `L_h = -L_h^{(5)} + h^2 C_h/6`.
- Minimal fix:
  - Rewrite the local derivation so `L_h approx \Delta` is maintained throughout.
  - Keep the final scheme `(-L_h+\sigma R_h)u=R_h f` and the raw denominator unchanged.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### A-002: Three-Parameter Effective Eigenvalue Has the Wrong Sign

- Severity: A
- File: `thesis/chapters/2_math_preliminary.tex:507`
- Current writing:
  - `\lambda_{\mathrm{mod}} = (\hat{\lambda}_L + \alpha h^2\hat{\lambda}_5 - \sigma(\hat{\lambda}_R+\beta h^2))/(\hat{\lambda}_R+\gamma h^2)`
- Why risky:
  - The operator displayed two lines above is `(-L_h - \alpha h^2 L_h^{(5)})u + \sigma(R_h+\beta h^2 I)u`.
  - Its Fourier numerator should start with `-\hat{\lambda}_L`, not `+\hat{\lambda}_L`.
  - The subsequent `c_2` formula matches the corrected sign, so the theorem currently contains a visible formula inconsistency.
- Correct version:
  - `\lambda_{\mathrm{mod}} = (-\hat{\lambda}_L - \alpha h^2\hat{\lambda}_5 + \sigma(\hat{\lambda}_R+\beta h^2))/(\hat{\lambda}_R+\gamma h^2)`.
- Minimal fix:
  - Replace the displayed effective eigenvalue formula and ensure `\hat{\lambda}_5` is described as the positive five-point `-\Delta` symbol.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### A-003: Lean Appendix Formula Does Not Match the Three-Parameter Family

- Severity: A
- File: `thesis/chapters/appendix_lean4.tex:17`
- Current writing:
  - `\hat{\lambda}_{\text{eff}} = ((-\hat{\lambda}_L+\sigma\hat{\lambda}_R)/\hat{\lambda}_R) * ...`
  - The formula omits `\beta`.
- Why risky:
  - The appendix says it is summarizing the three-parameter family, but the formula does not contain all three parameters and does not cleanly match the operator in Chapter 2.
  - This can make the Lean appendix look broader or less precise than the actual algebraic verification.
- Correct version:
  - Use the same corrected fraction as Chapter 2:
    `(-\hat{\lambda}_L - \alpha h^2\hat{\lambda}_5 + \sigma(\hat{\lambda}_R+\beta h^2))/(\hat{\lambda}_R+\gamma h^2)`.
- Minimal fix:
  - Replace the local formula or mark it as a schematic and point to the exact formula in Chapter 2.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### A-004: True Helmholtz Indefiniteness Condition Is Too Broad

- Severity: A
- File:
  - `thesis/chapters/4_helmholtz.tex:135`
  - Related broad wording: `thesis/main.tex:91`, `thesis/chapters/1_introduction.tex:20`, `thesis/chapters/2_math_preliminary.tex:19`, `thesis/chapters/7_conclusion.tex:23`
- Current writing:
  - `当 \kappa^2 > \lambda_{\min} 时，true Helmholtz 矩阵 A_h 是不定的`
  - Several summaries say true Helmholtz matrix is simply `不定`.
- Why risky:
  - For finite Dirichlet discretization, eigenvalues are `d_{p,q}=\lambda_{p,q}-\kappa^2`.
  - If `\kappa^2<\lambda_{\min}`, the matrix is still SPD.
  - If `\lambda_{\min}<\kappa^2<\lambda_{\max}`, it is indefinite.
  - If `\kappa^2>\lambda_{\max}`, all eigenvalues are negative and the matrix is negative definite.
- Correct version:
  - State the cases explicitly and use "可能不定/存在近共振风险" in global summaries.
- Minimal fix:
  - Replace the theorem statement with the interval condition `\lambda_{\min}<\kappa^2<\lambda_{\max}` for indefiniteness, include singular cases at eigenvalues, and soften cross-chapter summary wording.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### B-001: Near-Resonance `\delta` Sign Is Reversed

- Severity: B
- File: `thesis/chapters/4_helmholtz.tex:164`
- Current writing:
  - `\kappa^2 = \lambda_{p^*,q^*}+\delta`
  - `d_{p^*,q^*}=\delta`
- Why risky:
  - With `d=\lambda-\kappa^2`, the displayed definition gives `d=-\delta`.
  - The qualitative and norm-based conclusions using `|\delta|` remain correct, so this is not an experiment risk.
- Correct version:
  - Either write `d_{p^*,q^*}=-\delta`, or define `\kappa^2=\lambda_{p^*,q^*}-\delta` if the desired denominator is `+\delta`.
- Minimal fix:
  - Change the sign in the local modal formula and the following `\hat u` expression.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### C-001: FFT9 Fourth-Order Coefficient Wording Is Stronger Than Needed

- Severity: C
- File: `thesis/chapters/2_math_preliminary.tex:466`
- Current writing:
  - `对任意非零模式 ... h^4 项系数不为零`
- Why risky:
  - The coefficient is enough to show the `h^4` term is not identically zero.
  - As a continuous polynomial in wave-number ratios, it can vanish for special nonzero ratios, even though integer DST modes may avoid those ratios.
- Correct version:
  - Use "该系数不恒为零" or "一般不为零；例如取 `\xi=\eta` 时非零".
- Minimal fix:
  - Weaken the sentence; no formula or data change.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

### C-002: FACR Classical Complexity Wording Is Slightly Loose

- Severity: C
- File: `thesis/chapters/3_fft_direct.tex:166`
- Current writing:
  - `第三项趋于常数`
- Why risky:
  - In the displayed expression, the third term includes an `N^2` factor and is `O(N^2)`, not literally constant.
  - The important conclusion remains correct: the current FACR-like implementation is `O(N^2 log N)`, while classical optimized FACR is only background.
- Correct version:
  - Say the third term becomes `O(N^2)` and the total optimized classical FACR cost is governed by the `O(N^2 log log N)` part.
- Minimal fix:
  - Adjust one explanatory sentence.
- Code or experiment rerun needed:
  - No code change; no experiment rerun.
- Suggest asking GPT-5.5 Pro:
  - No.

## Confirmed Items

- Unified PDE sign convention is consistent across the audited text.
- Modified Helmholtz is consistently `\sigma=+\kappa^2`; true Helmholtz is consistently `\sigma=-\kappa^2`.
- No prohibited Kronecker block `B=tridiag(-1,4,-1)` was found.
- Five-point eigenvalues and condition-number asymptotics match the requested formulas.
- Five-point nonhomogeneous Dirichlet boundary terms use the plus sign on the RHS.
- FFT9 uses `L_h approx \Delta` in the final scheme and raw denominator `-\hat{\lambda}_L+\sigma\hat{\lambda}_R`.
- FFT9 nonhomogeneous Dirichlet correction uses `+L_IB g_B - sigma R_IB g_B` and the bottom-edge adjacent formula has the requested signs.
- FFT9 Fourier expansion uses continuous wave numbers `xi=p*pi`, `eta=q*pi`, `alpha=xi h`, `beta=eta h`.
- Neumann DCT-I discussion includes ghost-point symmetry, endpoint weights, compatibility, zero-mode handling, and weighted mean normalization.
- GMRES is limited to an unpreconditioned Krylov baseline; Chapter 6 states absolute residual tolerance.
- The positive-definite GMRES convergence estimate is not applied to true Helmholtz or raw Neumann Poisson zero-mode cases.
- FACR-like current implementation is stated as `O(N^2 log N)`, with `O(N^2 log log N)` only as classical background.
- Lean 4 is not described as a complete PDE proof; it is scoped to algebraic core checks.

## GPT-5.5 Pro Review Queue

No issue from this audit requires GPT-5.5 Pro adjudication. The A/B items above are local sign or spectral-condition corrections with standard derivations.
