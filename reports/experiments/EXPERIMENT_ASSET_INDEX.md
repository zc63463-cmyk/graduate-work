# Experiment Asset Index

Phase 0 asset freeze for Chapter 6.

This index was built by scanning the current source files and generated assets only:

- `thesis/chapters/6_experiments.tex`
- `thesis/figures/`
- `code/python/experiments/`
- `code/python/experiments/results/`
- recent reports: `FINAL_REVISION_REPORT.md`, `RECENT_EXPERIMENT_UPDATE_REPORT.md`, `RISKMAP_PATCH_REPORT_V2.md`, `PATCH4_REPORT.md`, `PATCH5_REPORT.md`

No experiment, algorithm, thesis text, CSV, PNG, or PDF was modified while building this index.

## 1. Final Chapter 6 Structure

The current Chapter 6 structure is:

1. `6.1 实验设置、可复现性与实现校验`
2. `6.2 实验一：Dirichlet 离散格式收敛性验证`
3. `6.3 实验二：Neumann 与 mixed 边界处理`
4. `6.4 实验三：精度--成本对比与复杂度实证`
5. `6.5 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为`
6. `6.6 实验五：True Helmholtz 近共振模态放大`
7. `6.7 实验总结`

Core experiment count: 5.

`exp00` is retained in Section 6.1 as implementation validation, not as a numbered core experiment.

## 2. Experiment-Level Asset Map

### 6.1 Implementation Validation

Status: not a core experiment.

Scripts:

- `code/python/experiments/exp00_fft_vs_sparse.py`

CSV:

- `code/python/experiments/results/exp00_fft_vs_sparse.csv`

Thesis figures:

- None.

Thesis tables:

- `tab:exp00_sparse_consistency`: FFT solver vs sparse direct consistency, `n=65`.

Scientific question answered:

- Do FFT-based direct solver paths and sparse direct solve the same discrete linear system to numerical precision?

Questions not answered:

- Does the continuous PDE discretization converge?
- Is FFT9 more accurate than five-point methods?
- Does this validate GMRES behavior or runtime performance?

Notes:

- This is implementation validation only.
- It supports discrete solver-path consistency, not continuous numerical accuracy.

### 6.2 Experiment 1: Dirichlet Discretization Convergence

Scripts:

- `code/python/experiments/exp01_convergence.py`
- `code/python/experiments/exp02_nonhom_bc.py`
- `code/python/experiments/generate_visualization_supplements.py`
- `code/python/experiments/exp06_complex_manufactured_visualization.py`

CSV:

- `code/python/experiments/results/exp01_convergence.csv`
- `code/python/experiments/results/exp02_nonhom_bc.csv`
- `code/python/experiments/results/exp06_complex_manufactured_convergence.csv`

Thesis figures:

- `fig:exp01_convergence`: `thesis/figures/exp01_convergence.png`
- `fig:exp02_nonhom_bc`: `thesis/figures/exp02_nonhom_bc.png`
- `fig:exp02_temperature_field_comparison`: `thesis/figures/exp02_temperature_field_comparison.png`
- `fig:exp06_complex_fields`: `thesis/figures/exp06_complex_manufactured_fields.png`
- `fig:exp06_complex_convergence`: `thesis/figures/exp06_complex_manufactured_convergence.png`

Thesis tables:

- `tab:exp01_poly_convergence`: polynomial manufactured solution convergence order, `sigma=10`.
- `tab:exp02_nonhom_bc`: nonhomogeneous Dirichlet representative data, `n=65`.

Scientific questions answered:

- Do five-point FA/CR/FACR-like solvers show second-order convergence for Dirichlet problems?
- Does FFT9 show fourth-order convergence for smooth Dirichlet problems?
- Is the nonhomogeneous Dirichlet boundary correction sign/template handled correctly?
- Can the solvers reproduce smooth multimode Dirichlet manufactured solutions without the old max-error anomaly?

Questions not answered:

- Does FFT9 support Neumann or mixed boundaries?
- Does the Fourier eigenfunction track alone establish general convergence?
- Does the nonhomogeneous Dirichlet visualization prove a new independent core experiment?
- Does this section evaluate GMRES performance?

Positioning:

- Fourier eigenfunction track is a frequency-formula sanity check.
- Nonhomogeneous Dirichlet figures are boundary visual sanity checks.
- Multimode Dirichlet visualization is supplementary and not a new core experiment.

### 6.3 Experiment 2: Neumann And Mixed Boundary Handling

Scripts:

- `code/python/experiments/exp03_neumann_mixed.py`

CSV:

- `code/python/experiments/results/exp03_neumann_mixed.csv`

Thesis figures:

- `fig:exp03_neumann_mixed`: `thesis/figures/exp03_neumann_mixed_summary.png`
- `fig:exp03_neumann_mixed_fields`: `thesis/figures/exp03_neumann_mixed_fields.png`

Thesis tables:

- `tab:exp03_neumann_mixed`: Neumann and mixed representative FA results, `n=65`.

Scientific questions answered:

- Does the ghost-point + DCT-I five-point framework handle supported Neumann and mixed boundary cases consistently?
- Do representative supported cases show second-order behavior?
- Do boundary diagnostics close: Neumann flux residual, mixed Dirichlet residual, weighted mean, mean offset, and compatibility integral?

Questions not answered:

- Does FFT9 support Neumann/mixed fourth-order verification?
- Does this prove all possible boundary combinations?
- Does this compare solver performance?

Positioning:

- This experiment is about boundary constraints and zero-mode handling, not a new proof of FFT9 boundary capability.

### 6.4 Experiment 3: Accuracy--Cost Comparison And Complexity Evidence

Scripts:

- `code/python/experiments/exp06_accuracy_cost.py`

CSV:

- `code/python/experiments/results/exp06_accuracy_cost.csv`

Thesis figures:

- `fig:exp06_accuracy_cost`: `thesis/figures/exp06_accuracy_cost_error_time.png`
- `fig:exp06_time_scaling`: `thesis/figures/exp06_time_scaling.png`

Thesis tables:

- `tab:exp06_accuracy_cost_compact`: representative accuracy--cost values, `sigma=10`, `n=129`.

Scientific questions answered:

- Under the same smooth Dirichlet manufactured solution, how do FA, CR, FACR-like, FFT9, and unpreconditioned GMRES(30) compare in error vs time?
- Does FFT9 achieve significantly lower error at comparable FFT-direct cost in this smooth Dirichlet setting?
- How do current implementations scale in time against an `O(N^2 log N)` reference trend?

Questions not answered:

- Does this constitute a strict kernel-to-kernel benchmark?
- Does this compare all possible optimized implementations?
- Does this establish preconditioned Krylov performance?
- Does FACR-like achieve classical `O(N^2 log log N)` complexity?

Positioning:

- This is the core evidence for the thesis title's “FFT fast solver vs iterative method comparison”.
- Direct solver timing is solver-call timing including internal RHS/transform setup.
- GMRES timing is core solve timing with matrix/RHS setup excluded.
- The timing scope is deliberately stated and should not be read as a strict kernel benchmark.

### 6.5 Experiment 4: Modified/True Helmholtz Spectral Structure And GMRES Behavior

Scripts:

- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `code/python/experiments/exp04_modified_vs_true.py`

CSV:

- `code/python/experiments/results/exp07_spectral_denominator_summary.csv`
- `code/python/experiments/results/exp04_modified_vs_true.csv`
- `code/python/experiments/results/exp04_condition_check.csv`
- `code/python/experiments/results/exp04_gmres_history.csv`

Thesis figures:

- `fig:exp07_spectral_denominator_maps`: `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `fig:exp04_spectral`: `thesis/figures/exp04_min_denom_vs_sigma.png` and `thesis/figures/exp04_spectral_indicator_vs_sigma.png`
- `fig:exp04_condition_check`: `thesis/figures/exp04_condition_check.png`
- `fig:exp04_gmres`: `thesis/figures/exp04_gmres_iters_vs_sigma.png`
- `fig:exp04_gmres_history`: `thesis/figures/exp04_gmres_history.png`

Thesis tables:

- `tab:exp04_mod_true`: representative spectral indicator and GMRES iteration values, `n=33`; GMRES setting as in `fig:exp04_gmres`.

Scientific questions answered:

- How do frequency denominators differ between modified Helmholtz and true Helmholtz?
- How does true Helmholtz produce sign-changing denominators and small-denominator risk?
- In the Dirichlet five-point case, does the spectral denominator indicator agree with dense `cond_2(A)`?
- How does unpreconditioned GMRES(30) behave under modified/true Helmholtz Gaussian RHS tests?
- What does the residual history show beyond final iteration counts?

Questions not answered:

- Does the condition-number equivalence apply to Neumann ghost-point, mixed-boundary, non-normal, or non-orthogonally diagonalized systems?
- Does the risk map itself represent a condition number?
- Does this validate preconditioned GMRES, MINRES, or shifted-Laplacian methods?
- Does this establish all wave-number regimes?

Positioning:

- `exp07 + exp04` together form the Helmholtz spectral-structure and GMRES behavior experiment.
- The spectral denominator indicator equals `cond_2(A)` only for the current Dirichlet five-point orthogonal DST-I system.
- GMRES is unpreconditioned GMRES(30) with absolute residual tolerance.
- `501` in exp04 denotes capped count under `max_iter=500`.

### 6.6 Experiment 5: True Helmholtz Near-Resonance Modal Amplification

Scripts:

- `code/python/experiments/exp05_true_helmholtz_resonance.py`

CSV:

- `code/python/experiments/results/exp05_resonance.csv`
- `code/python/experiments/results/exp05_multimode_resonance.csv`
- `code/python/experiments/results/exp05_resonance_gmres_history.csv`

Thesis figures:

- `fig:true_helmholtz_near_resonance_summary`: `thesis/figures/exp05_near_resonance_summary.png`
- `fig:exp05_multimode_resonance_summary`: `thesis/figures/exp05_multimode_resonance_summary.png`
- `fig:exp05_dominant_mode_projection`: `thesis/figures/exp05_dominant_mode_projection.png`
- `fig:exp05_resonance_gmres_history`: `thesis/figures/exp05_resonance_gmres_history.png`

Thesis tables:

- `tab:exp05_resonance`: near-resonance representative data, `n=65`.

Scientific questions answered:

- Does true Helmholtz near a discrete eigenvalue exhibit small-denominator modal amplification?
- Does the target modal coefficient follow the `1/|delta|` law?
- Does the near-resonance solution become dominated by the target modal subspace in the tested Gaussian RHS setting?
- Does unpreconditioned GMRES(30) reach capped/not-converged behavior in these near-resonance cases?
- What does residual history show about stagnation?

Questions not answered:

- Does this describe the continuous spectrum limit?
- Does it apply to every RHS or every Helmholtz solver?
- Does it evaluate preconditioned Krylov methods?
- Does it show a general monotonic trend in energy/correlation for all RHS choices?

Positioning:

- `exp05` is the near-resonance modal amplification mechanism verification.
- The current Gaussian RHS already gives target subspace energy and shape correlation near 1 across the scan; panels (b)(c) diagnose dominance, not a large trend.
- `1001` denotes the capped count convention in this near-resonance GMRES scan.

### 6.7 Experiment Summary

Scripts / CSV / figures:

- No new assets.

Role:

- Summarizes the evidence from Sections 6.1--6.6.
- Keeps `exp00` downgraded as implementation validation.
- Preserves five core experiments.

## 3. Key Figure Index

| Figure label | Figure file(s) | Section / experiment | Axes, colors, markers | Reading conclusion | Notes / limitations |
|---|---|---|---|---|---|
| `fig:exp01_convergence` | `exp01_convergence.png` | 6.2 Experiment 1 | log-log; x=`h`; y=`||u-u_h||_infty`; method colors/markers; reference `O(h^2)` and `O(h^4)` lines | Five-point methods follow second order; FFT9 follows fourth order on Dirichlet smooth tests | Fourier eigenfunction track is only a frequency-formula sanity check; not GMRES evidence |
| `fig:exp02_nonhom_bc` | `exp02_nonhom_bc.png` | 6.2 nonhomogeneous Dirichlet subcase | log-log; x=`h`; y=`||u-u_h||_infty`; method curves; second/fourth-order references | Nonhomogeneous Dirichlet boundary correction gives expected convergence | Boundary visual/symbol check; not an independent core experiment |
| `fig:exp02_temperature_field_comparison` | `exp02_temperature_field_comparison.png` | 6.2 nonhomogeneous Dirichlet subcase | field panels for exact, FFT9 numerical, log-error, line cut | Physical-space boundary contribution handling is visually consistent | Visual sanity check only; convergence conclusion comes from tables/curves |
| `fig:exp06_complex_fields` | `exp06_complex_manufactured_fields.png` | 6.2 multimode Dirichlet visual check | 2x3 fields; exact/FA/FFT9 fields, log errors, line cut | FA and FFT9 reproduce multimode smooth structure; old max-error anomaly removed | Supplementary visualization; not a new core experiment |
| `fig:exp06_complex_convergence` | `exp06_complex_manufactured_convergence.png` | 6.2 multimode Dirichlet visual check | log-log; x=`h`; y=error; FA vs FFT9 | Multimode FA stays near second order; FFT9 stays near fourth order | Supports visual check; main convergence evidence remains `exp01` |
| `fig:exp03_neumann_mixed` | `exp03_neumann_mixed_summary.png` | 6.3 Experiment 2 | 1x3: (a) FA `L_infty` vs `h`; (b) boundary residual vs `h`; (c) zero-mode diagnostics vs `h`; log-log; `x` marker for mixed Dirichlet residual | Supported Neumann/mixed five-point solvers show convergence, boundary residual closure, and zero-mode diagnostics | Does not claim FFT9 Neumann/mixed support |
| `fig:exp03_neumann_mixed_fields` | `exp03_neumann_mixed_fields.png` | 6.3 Experiment 2 | 2x3 fields; exact, numerical, log error for Pure Neumann and Mixed `(N,D)` | Representative spatial fields and error distributions look consistent | Visual supplement; CSV/summary figure carry convergence/diagnostic claims |
| `fig:exp06_accuracy_cost` | `exp06_accuracy_cost_error_time.png` | 6.4 Experiment 3 | log-log; x=median solve time; y=`L_infty` error; marker/edge color=method; colorbar=`N^2`; black x=GMRES not converged | FFT9 reaches much lower error at comparable FFT-direct cost; GMRES30 can fail/not converge | Current implementation benchmark, not strict kernel-to-kernel timing |
| `fig:exp06_time_scaling` | `exp06_time_scaling.png` | 6.4 Experiment 3 | log-log; x=unknowns `N^2`; y=median solve-call time; method curves; `O(N^2 log N)` reference | Current FFT direct methods follow expected FFT-like scaling trend | FACR-like is not claimed to reach classical `O(N^2 log log N)` |
| `fig:exp07_spectral_denominator_maps` | `exp07_spectral_denominator_heatmaps.png` | 6.5 Experiment 4 | top: low-frequency `p,q<=12` risk map `R=-log10|d|`; green markers=(2,3)/(3,2); bottom: sorted smallest `|d|` rank plot; star markers for near targets | Modified: low risk; true-away: sign-changing but no small denominator hotspot; true-near: first two denominators are `(2,3)/(3,2)` at `10^-2` scale | Risk map is not a condition number; `d=0` contour is interpolated sign boundary, not a discrete zero denominator |
| `fig:exp04_spectral` | `exp04_min_denom_vs_sigma.png`, `exp04_spectral_indicator_vs_sigma.png` | 6.5 Experiment 4 | x=`sigma` symlog; y=min denominator or spectral denominator indicator; equation-type curves | Modified/true Helmholtz differ through denominator structure; true cases can approach small denominators | Figure uses max grid `n=129`; Table 9 uses `n=33`, so statuses should not be compared pointwise across grid sizes |
| `fig:exp04_condition_check` | `exp04_condition_check.png` | 6.5 Experiment 4 | log-log scatter; x=spectral denominator indicator; y=dense `cond_2(A)`; `y=x` reference; marker by `n` | For current Dirichlet five-point DST-I system, indicator equals `cond_2(A)` numerically | Equivalence is not generalized to other boundary/non-normal systems |
| `fig:exp04_gmres` | `exp04_gmres_iters_vs_sigma.png` | 6.5 Experiment 4 | x=`sigma` symlog; y=GMRES iterations; red x=failed/capped; line at 501=capped under `max_iter=500` | Unpreconditioned GMRES(30) gets harder near true-Helmholtz small denominators | Gaussian RHS baseline only; no preconditioner tested |
| `fig:exp04_gmres_history` | `exp04_gmres_history.png` | 6.5 Experiment 4 | 1x2 residual history; x=iteration; y=absolute/relative residual log scale; red x=capped end; horizontal absolute tolerance line | Residual history explains convergence/stagnation beyond final iteration counts | Stop criterion is absolute residual; relative residual is diagnostic only |
| `fig:true_helmholtz_near_resonance_summary` | `exp05_near_resonance_summary.png` | 6.6 Experiment 5 | 2x2: amplification vs `|delta|`, GMRES residual vs `|delta|`, raw fields for `delta=1e-1` and `1e-4`; red x=capped | Smaller `|delta|` amplifies the solution and leaves GMRES capped/not converged | Raw fields use same near-resonance branch; spatial shapes may remain similar |
| `fig:exp05_multimode_resonance_summary` | `exp05_multimode_resonance_summary.png` | 6.6 Experiment 5 | 2x2: (a) modal coefficients `|u_hat_23|`, `|u_hat_32|` vs `|delta|`; (b) target energy fraction; (c) shape correlation; (d) GMRES capped iterations; red x=capped | Modal amplitudes follow `1/|delta|`; target subspace dominates; GMRES hits capped count | Panels (b)(c) are near-1 dominance diagnostics, not large trend plots |
| `fig:exp05_dominant_mode_projection` | `exp05_dominant_mode_projection.png` | 6.6 Experiment 5 | fields: full solution, dominant projection, difference; common spatial axes | Near resonance, full solution is visually close to target modal projection | Specific to `(2,3)/(3,2)`, `delta=1e-4`, Gaussian RHS |
| `fig:exp05_resonance_gmres_history` | `exp05_resonance_gmres_history.png` | 6.6 Experiment 5 | 1x2 residual history; x=iteration; y=absolute/relative residual log scale; red x=capped end | Residual history shows stagnation under near-resonance detuning | Unpreconditioned GMRES(30); absolute tolerance `1e-10`; no preconditioner |

## 4. Table Index

| Table label | Section | CSV source | Purpose | Notes |
|---|---|---|---|---|
| `tab:exp00_sparse_consistency` | 6.1 | `exp00_fft_vs_sparse.csv` | FFT/direct solver path consistency | Implementation validation only |
| `tab:exp01_poly_convergence` | 6.2 | `exp01_convergence.csv` | Polynomial manufactured convergence orders | Main Dirichlet convergence table |
| `tab:exp02_nonhom_bc` | 6.2 | `exp02_nonhom_bc.csv` | Nonhomogeneous Dirichlet representative data | Boundary correction sanity check |
| `tab:exp03_neumann_mixed` | 6.3 | `exp03_neumann_mixed.csv` | Neumann/mixed representative FA data | Boundary/zero-mode diagnostics |
| `tab:exp06_accuracy_cost_compact` | 6.4 | `exp06_accuracy_cost.csv` | Representative time-error pairings | Timing scope differs between direct solvers and GMRES and is stated |
| `tab:exp04_mod_true` | 6.5 | `exp04_modified_vs_true.csv` | Spectral indicator and GMRES representative values | Table uses `n=33`; figures use larger grid where stated |
| `tab:exp05_resonance` | 6.6 | `exp05_resonance.csv` | Near-resonance representative data | `1001` is capped-count convention |

## 5. Artifact Positioning Tags

| Asset / experiment | Positioning |
|---|---|
| `exp00_fft_vs_sparse.py` | Implementation validation; not a core experiment |
| Fourier eigenfunction convergence track in `exp01` | Frequency formula sanity check |
| `exp02_nonhom_bc` and `exp02_temperature_field_comparison` | Nonhomogeneous Dirichlet boundary correction / visual sanity check |
| `exp06_complex_manufactured_visualization.py` | Multimode Dirichlet visual supplement; not a separate core experiment |
| `exp03_neumann_mixed.py` | Supported five-point Neumann/mixed boundary validation |
| `exp06_accuracy_cost.py` | Core evidence for “FFT fast solver vs iterative method comparison” |
| `exp07_spectral_denominator_maps.py` + `exp04_modified_vs_true.py` | Helmholtz spectral structure and GMRES behavior |
| `exp05_true_helmholtz_resonance.py` | Near-resonance modal amplification mechanism verification |

## 6. Global Limitations To Preserve

- FFT9 is implemented and verified only for Dirichlet boundary conditions.
- GMRES claims refer to unpreconditioned restarted GMRES, primarily GMRES(30), under the stated absolute residual tolerances.
- GMRES relative residual curves are diagnostic; stopping criteria are absolute residual tolerances where stated.
- FACR-like current implementation is treated as `O(N^2 log N)` because it still includes full-direction DST work; classical `O(N^2 log log N)` is not claimed for the current implementation.
- Accuracy--cost figures are current implementation benchmarks, not strict kernel-to-kernel timing comparisons.
- Spectral denominator indicator equals `cond_2(A)` only for the current Dirichlet five-point system with orthogonal DST-I diagonalization.
- Risk maps visualize small-denominator risk; they are not condition-number plots.
- Near-resonance experiments are discrete-spectrum experiments and should not be stated as continuous-spectrum limit results.

## 7. Current Asset Inventory Snapshot

Experiment scripts in `code/python/experiments/`:

- `exp00_fft_vs_sparse.py`
- `exp01_convergence.py`
- `exp02_nonhom_bc.py`
- `exp03_neumann_mixed.py`
- `exp04_modified_vs_true.py`
- `exp05_true_helmholtz_resonance.py`
- `exp06_accuracy_cost.py`
- `exp06_complex_manufactured_visualization.py`
- `exp07_spectral_denominator_maps.py`
- `generate_visualization_supplements.py`
- `utils.py`

CSV files in `code/python/experiments/results/`:

- `exp00_fft_vs_sparse.csv`
- `exp01_convergence.csv`
- `exp02_nonhom_bc.csv`
- `exp03_neumann_mixed.csv`
- `exp04_condition_check.csv`
- `exp04_gmres_history.csv`
- `exp04_modified_vs_true.csv`
- `exp05_multimode_resonance.csv`
- `exp05_resonance.csv`
- `exp05_resonance_gmres_history.csv`
- `exp06_accuracy_cost.csv`
- `exp06_complex_manufactured_convergence.csv`
- `exp07_spectral_denominator_summary.csv`

Current Chapter 6 thesis figures:

- `exp01_convergence.png`
- `exp02_nonhom_bc.png`
- `exp02_temperature_field_comparison.png`
- `exp03_neumann_mixed_summary.png`
- `exp03_neumann_mixed_fields.png`
- `exp04_min_denom_vs_sigma.png`
- `exp04_spectral_indicator_vs_sigma.png`
- `exp04_condition_check.png`
- `exp04_gmres_iters_vs_sigma.png`
- `exp04_gmres_history.png`
- `exp05_near_resonance_summary.png`
- `exp05_multimode_resonance_summary.png`
- `exp05_dominant_mode_projection.png`
- `exp05_resonance_gmres_history.png`
- `exp06_accuracy_cost_error_time.png`
- `exp06_time_scaling.png`
- `exp06_complex_manufactured_fields.png`
- `exp06_complex_manufactured_convergence.png`
- `exp07_spectral_denominator_heatmaps.png`

Supplemental or legacy-present figures not currently used as main indexed Chapter 6 evidence:

- `exp03_neumann_mixed.png`
- `exp05_resonance.png`
- `exp05_denominator_heatmap.png`

These should not be promoted back into core evidence without an explicit later decision.
