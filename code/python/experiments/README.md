# Chapter 6 Experiment Assets

This directory contains the Python entrypoints and generated evidence used by Chapter 6.

## Directory Layout

| Path | Purpose |
|---|---|
| `exp00_fft_vs_sparse.py` | Implementation validation: FFT solver path vs sparse direct solve on the same discrete system |
| `exp01_convergence.py` | Dirichlet five-point and FFT9 convergence study |
| `exp02_nonhom_bc.py` | Nonhomogeneous Dirichlet boundary visual sanity checks |
| `exp03_neumann_mixed.py` | Neumann and mixed boundary verification |
| `exp04_modified_vs_true.py` | Modified/true Helmholtz spectral indicator, condition check, and GMRES behavior |
| `exp05_true_helmholtz_resonance.py` | True Helmholtz near-resonance modal amplification |
| `exp06_accuracy_cost.py` | Accuracy-cost and timing comparison |
| `exp07_spectral_denominator_maps.py` | Low-frequency small-denominator risk maps |
| `results/` | CSV outputs used as evidence |
| `figures/` | Generated experiment figures copied into `thesis/figures/` when used in the thesis |

## Runtime Safety

The scripts assume the current `results/` and `figures/` layout. Do not move these subdirectories without updating scripts, thesis figure references, and rebuild instructions together.

## Standard Commands

From the repository root:

```powershell
PYTHONPATH=code/python python code/python/experiments/exp04_modified_vs_true.py
PYTHONPATH=code/python python code/python/experiments/exp05_true_helmholtz_resonance.py
PYTHONPATH=code/python python code/python/experiments/exp06_accuracy_cost.py
PYTHONPATH=code/python python code/python/experiments/exp07_spectral_denominator_maps.py
```

For routine project validation, prefer tests first:

```powershell
python -m pytest -q code/python/tests
```
