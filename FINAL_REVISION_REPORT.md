# Final Revision Report

## 1. Final Chapter 6 Structure

Chapter 6 has been reorganized into one setup/validation section plus five core experiments:

1. `6.1 实验设置、可复现性与实现校验`
   - `exp00_fft_vs_sparse.py` is kept as implementation validation only.
   - It verifies FFT and sparse direct paths solve the same discrete system.
   - It is not counted as a core experiment and does not claim continuous PDE accuracy.

2. `6.2 实验一：Dirichlet 离散格式收敛性验证`
   - Combines `exp01_convergence.py` and `exp02_nonhom_bc.py`.
   - Keeps the Fourier/eigenfunction track as a frequency-formula sanity check.
   - Uses the polynomial manufactured solution as the main Dirichlet convergence evidence.
   - Keeps the nonhomogeneous Dirichlet figures in the body as boundary visual sanity checks.

3. `6.3 实验二：Neumann 与 mixed 边界处理`
   - Uses `exp03_neumann_mixed.py`.
   - Emphasizes ghost-point + DCT-I/DST-I symmetrization, flux residual, Dirichlet boundary residual, weighted mean, and zero-mode handling.

4. `6.4 实验三：精度--成本对比与复杂度实证`
   - Uses `exp06_accuracy_cost.py`.
   - Promotes the accuracy-cost comparison as the core evidence for the "FFT fast solvers vs iterative method" thesis theme.

5. `6.5 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为`
   - Integrates `exp07_spectral_denominator_maps.py`, `exp04_modified_vs_true.py`, the condition check, and GMRES residual history.

6. `6.6 实验五：True Helmholtz 近共振模态放大`
   - Uses `exp05_true_helmholtz_resonance.py`.
   - Integrates the original near-resonance scan, multimode modal projection check, dominant-mode projection, and near-resonance GMRES residual history.

## 2. Modified Files

- `thesis/chapters/6_experiments.tex`
- `thesis/main.tex`
- `thesis/chapters/1_introduction.tex`
- `thesis/chapters/7_conclusion.tex`
- `README.md`
- `FINAL_POLISH_REPORT.md`
- `VISUALIZATION_SUPPLEMENT_REPORT.md`
- `thesis/main.pdf`
- Regenerated experiment outputs from exp04, exp05, exp06, and exp07 scripts.

## 3. Core Figures And Tables

### Experiment 1: Dirichlet Convergence

- `fig:exp01_convergence`
- `tab:exp01_poly_convergence`
- `fig:exp02_nonhom_bc`
- `tab:exp02_nonhom_bc`
- `fig:exp02_temperature_field_comparison`
- `fig:exp06_complex_fields`
- `fig:exp06_complex_convergence`

### Experiment 2: Neumann/Mixed Boundary Processing

- `fig:exp03_neumann_mixed`
- `fig:exp03_neumann_mixed_fields`
- `tab:exp03_neumann_mixed`

### Experiment 3: Accuracy-Cost And Complexity

- `fig:exp06_accuracy_cost`
- `fig:exp06_time_scaling`

### Experiment 4: Helmholtz Spectrum And GMRES

- `fig:exp07_spectral_denominator_maps`
- `fig:exp04_spectral`
- `fig:exp04_condition_check`
- `fig:exp04_gmres`
- `tab:exp04_mod_true`
- `fig:exp04_gmres_history`

### Experiment 5: True Helmholtz Near-Resonance

- `fig:true_helmholtz_near_resonance_summary`
- `tab:exp05_resonance`
- `fig:exp05_multimode_resonance_summary`
- `fig:exp05_dominant_mode_projection`
- `fig:exp05_resonance_gmres_history`

## 4. Original Experiment 3 Figures

The original nonhomogeneous Dirichlet figures were preserved in Chapter 6, but folded into Experiment 1:

- `fig:exp02_nonhom_bc`: now described as a nonhomogeneous Dirichlet boundary subcase checking RHS/boundary correction signs and templates.
- `fig:exp02_temperature_field_comparison`: retained as a physical-space visual sanity check for moving boundary contributions to the RHS.

The multimode manufactured solution figures were also retained in Experiment 1:

- `fig:exp06_complex_fields`
- `fig:exp06_complex_convergence`

They are explicitly presented as Dirichlet visual/sanity-check supplements, not as a new core experiment.

## 5. Downgraded Content

- `exp00_fft_vs_sparse.py`: downgraded from a core experiment to implementation validation.
- Fourier eigenfunction convergence track: downgraded to frequency-formula sanity check.
- Nonhomogeneous Dirichlet temperature/field figures: downgraded to boundary visual sanity checks while remaining in the main text.
- Repeated convergence discussion: consolidated under Experiment 1.

## 6. Key Scope Limits

- GMRES results are for unpreconditioned restarted GMRES, primarily GMRES(30), and do not represent all Krylov or preconditioned methods.
- FFT9 is implemented and verified only for Dirichlet boundary conditions.
- FACR-like current implementation remains `O(N^2 log N)` and does not claim classical FACR `O(N^2 log log N)` complexity.
- The accuracy-cost timing scope is not a strict kernel-to-kernel benchmark; direct solvers include internal RHS/transform setup, while GMRES core timing excludes matrix/RHS setup.
- Near-resonance conclusions are discrete-spectrum experiments and are not direct continuous-spectrum limit claims.
- The condition-number equivalence between spectral denominator indicator and `cond_2(A)` is limited to the current Dirichlet five-point system with orthogonal DST-I diagonalization.

## 7. Validation Commands

Executed from repository root:

```powershell
PYTHONPATH=code/python python code/python/experiments/exp04_modified_vs_true.py
PYTHONPATH=code/python python code/python/experiments/exp05_true_helmholtz_resonance.py
PYTHONPATH=code/python python code/python/experiments/exp06_accuracy_cost.py
PYTHONPATH=code/python python code/python/experiments/exp07_spectral_denominator_maps.py
```

All four scripts completed successfully and regenerated their expected CSV/figure outputs.

Python tests:

```powershell
cd code/python
python -m pytest -q
```

Result:

```text
103 passed, 2 skipped, 2 warnings
```

The two warnings are the expected Neumann Poisson compatibility warnings in tests.

LaTeX:

```powershell
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Result:

```text
Output written on main.pdf (66 pages).
```

## 8. LaTeX Log Scan

Final scan status:

- no fatal error
- no undefined citation/reference
- no multiply-defined label
- no overfull hbox
- no missing character

Remaining non-blocking warnings:

- SimSun bold italic font substitution.
- hyperref PDF bookmark token warnings for math symbols.

## 9. Final PDF

- `thesis/main.pdf`
- Final page count: `66 pages`

## 10. Final Status

READY TO SUBMIT
