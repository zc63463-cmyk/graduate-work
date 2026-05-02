# Paper-Code Consistency Report

This is a read-only audit. No thesis, code, test, experiment, CSV, or figure file was modified.

## Executive Summary

- Matrix rows reviewed: 14
- A-class issues: 1
- B-class issues: 2
- C-class issues: 1
- D-class confirmations: 10
- Solver code generally matches the final Chapter 3/Chapter 6 framework. The main A-class problem is a Chapter 2 FFT9 derivation/sign inconsistency relative to the final scheme and Python implementation.

## A-Class Issues

### A-01: Chapter 2 FFT9 derivation signs do not fully match the implemented final scheme

- Location: `thesis/chapters/2_math_preliminary.tex:315`, `thesis/chapters/2_math_preliminary.tex:322`, `thesis/chapters/2_math_preliminary.tex:507`
- Current paper writing: Chapter 2 contains an FFT9 derivation/effective-eigenvalue path that uses signs inconsistent with the final statement that `L_h` approximates `Delta` and the final raw scheme `(-L_h + sigma R_h)u = R_h f`.
- Code evidence: `code/python/helmholtz_solver.py:993-1050` solves the equivalent sign-flipped system `(L_h - sigma R_h)u = -R_h f` with denominator `lambda_L - sigma lambda_R`; `code/python/verify_fft9_expansion.py:107-113` verifies `-lambda_L/lambda_R`.
- Test evidence: `code/python/tests/test_04_fft9_vs_spsolve.py:9-17` explicitly documents the sign-flipped equivalent implementation and compares against sparse matrix rows.
- Risk: A reader can conclude a different compact operator or denominator from Chapter 2 than the one implemented and used in experiments.
- Correct writing: Chapter 2 should consistently state `L_h` approximates `Delta`, `(-L_h + sigma R_h)u = R_h f`, raw Fourier denominator `-lambda_L + sigma lambda_R`, and effective expansion `-lambda_L/lambda_R`.
- Minimal fix: Paper-only formula correction in Chapter 2. No code change appears necessary.
- Need code change: no
- Need experiment rerun: no

## B-Class Issues

### B-01: Exp04 GMRES figure can hide max-iteration failures

- Location: `code/python/experiments/exp04_modified_vs_true.py:331-349`
- Current code writing: the GMRES iteration plot includes all rows with `gmres_iters > 0`; it does not filter or mark `gmres_success == False`.
- Evidence: `code/python/experiments/results/exp04_modified_vs_true.csv` contains failed rows with capped iteration counts such as 501. The figure `code/python/experiments/figures/exp04_gmres_iters_vs_sigma.png` plots the largest grid size series without a visible failure marker.
- Paper status: `thesis/chapters/6_experiments.tex:229` uses the safer `n=33` table and does not overclaim full iteration-method performance, so this is an evidence-presentation risk rather than a central formula error.
- Risk: The PNG can be read as "501 successful iterations" instead of "hit the cap or failed to meet tolerance".
- Correct writing/behavior: Either filter to successful rows, use a different marker for failures, or state in the caption that capped failures are included.
- Minimal fix: Plot or caption update only.
- Need code change: yes, if regenerating the figure from script
- Need experiment rerun: yes, regenerate Exp04 figure after the plot change

### B-02: Exp05 baseline rows hard-code `gmres_success=True`

- Location: `code/python/experiments/exp05_true_helmholtz_resonance.py:176-190`
- Current code writing: baseline rows append `'gmres_success': True` regardless of the GMRES solver result.
- Evidence: `code/python/experiments/results/exp05_resonance.csv` can therefore label baseline rows as successful even when the iteration count reaches the cap or the residual is large.
- Paper status: the near-resonance rows and figure support the Chapter 6 near-resonance conclusion; this issue affects baseline metadata quality.
- Risk: CSV consumers may treat failed or capped baseline GMRES runs as successful.
- Correct writing/behavior: Set `gmres_success` from `info['success']`, or from the same residual/tolerance criterion used in the GMRES solver.
- Minimal fix: Code metadata fix in Exp05 script, then regenerate the CSV.
- Need code change: yes
- Need experiment rerun: yes, regenerate Exp05 CSV; figure regeneration is safe if the script rewrites it

## C-Class Issues

### C-01: `test_07_krylov_baselines.py` header overstates iteration-ordering coverage

- Location: `code/python/tests/test_07_krylov_baselines.py:21-23`
- Current code writing: the module docstring claims it compares GMRES iteration counts with ordering `modified < Poisson < true`.
- Contradicting implementation: `code/python/tests/test_07_krylov_baselines.py:202-205` states that single-mode RHS is not meaningful for iteration-count comparison and only verifies convergence sanity.
- Paper status: Chapter 6 uses Gaussian RHS for GMRES behavior and does not rely on this stale test header.
- Risk: Low. The stale comment can confuse future maintainers but does not affect reported thesis evidence.
- Correct writing/behavior: Update the docstring to say this test checks solver applicability/convergence sanity, not iteration-count ordering.
- Minimal fix: Comment-only update.
- Need code change: yes, comment only
- Need experiment rerun: no

## D-Class Confirmations

- Five-point nonhomogeneous Dirichlet correction is consistent: paper and code use `F += bc/h^2` in FA, GMRES, and the cyclic-reduction prototype.
- Exp02 formula is consistent across paper, `utils.py`, experiment script, tests, and CSV generation: the RHS uses the sine mode plus `sigma`, not `(lambda+sigma)u`.
- Chapter 3 FFT9 raw denominator and Python implementation are consistent once the implementation's global sign-flipped system is accounted for.
- FFT9 boundary correction signs are consistent between Chapter 3 and code under the same sign convention.
- FFT9 fourth-order verification is scoped to Dirichlet; Neumann and mixed cases warn or fall back in code.
- GMRES is consistently presented as an unpreconditioned baseline using absolute residual stopping.
- Exp04 uses Gaussian non-eigenfunction RHS, not a single eigenmode.
- Exp05 near-resonance CSV and figure exist and support the near-resonance conclusion.
- FACR-like complexity is consistently described as current `O(N^2 log N)`, with classical `O(N^2 log log N)` only as background.
- README and FINAL_POLISH_REPORT are consistent on 52 PDF pages, six experiments, six CSV outputs, and seven core PNG figures.

## Files Needing Patch If You Choose to Fix

- `thesis/chapters/2_math_preliminary.tex` for the A-class FFT9 sign/effective-eigenvalue correction.
- `code/python/experiments/exp04_modified_vs_true.py` for failure-aware GMRES plotting.
- `code/python/experiments/exp05_true_helmholtz_resonance.py` for baseline `gmres_success` metadata.
- `code/python/tests/test_07_krylov_baselines.py` for a stale comment/docstring update.

No new experiments are required by this audit. Two existing experiment artifacts would need regeneration only if their scripts are patched: Exp04 figure and Exp05 CSV/figure.
