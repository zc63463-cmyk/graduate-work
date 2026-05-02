# FIG10 Final Visual Polish Report

## 1. Scope

本轮只微调第 6 章图 10 的可读性，不新增实验，不改变 exp07 的参数设置，不改变理论结论。

图 10 仍保持现有结构：

- 上排：低频窗口 `p,q <= 12` 的 small-denominator risk map；
- 下排：前 20 个最小 `|d_{p,q}|` 的排序图；
- 输出文件名保持 `exp07_spectral_denominator_heatmaps.png/pdf` 不变。

## 2. Files Modified

- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `thesis/chapters/6_experiments.tex`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.png`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `thesis/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/main.pdf`

## 3. Visual Changes

### Top Row

Panel titles were simplified to:

- `modified`
- `true-away`
- `true-near`

The case descriptions and `min |d|` values were moved into compact in-panel boxes, reducing title crowding while keeping the numerical summary visible.

The colorbar label was simplified to:

```text
risk = -log10 |d_{p,q}|
```

The plot does not describe risk as a condition number.

### d=0 Contour

The `d=0` contour was retained only as an auxiliary sign-boundary cue in the true-away panel.

The old label:

```text
black contour: d=0
```

was replaced with:

```text
interpolated sign boundary d=0
```

The thesis text now states that this contour is an interpolated sign boundary and does not mean a discrete modal denominator is exactly zero.

### Bottom Row

The first two near-resonance sorted-denominator points were emphasized with large green star markers with black edges.

A single annotation now labels both points together:

```text
(2,3)/(3,2), |d| ≈ 10^{-2}
```

A horizontal reference line `|d| = 10^{-2}` was added to make the near-resonance scale immediately visible.

## 4. Numerical Summary

From `code/python/experiments/results/exp07_spectral_denominator_summary.csv`:

| Case | min `|d|` | nearest mode |
|---|---:|---|
| modified Helmholtz, `sigma=+100` | `1.197382e+02` | `(1,1)` |
| true-away | `1.480069e+01` | `(1,2)` |
| true-near, `(2,3)/(3,2)` | `1.000000e-02` | `(2,3)` |

The sorted-denominator plot shows that the near-resonance case has its first two smallest denominators at the target degenerate pair `(2,3)/(3,2)`, both at the `10^{-2}` scale.

## 5. Thesis Text Update

The Figure 10 surrounding text now explicitly states:

- `R_{p,q}<0` only means `|d_{p,q}|>1`, i.e. low small-denominator risk, not “negative risk”.
- A `d=0` contour, when shown, is an interpolated sign boundary and not a discrete zero-denominator mode.
- The lower sorted plot identifies `(2,3)/(3,2)` as the first two smallest denominators in the near-resonance case.
- The risk map is only a visualization of small-denominator risk and is not a condition-number plot.

## 6. Commands Run

```powershell
PYTHONPATH=code/python python code/python/experiments/exp07_spectral_denominator_maps.py

cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

## 7. Validation

- exp07 script: success
- XeLaTeX/BibTeX/XeLaTeX/XeLaTeX: success
- `thesis/main.pdf`: generated successfully
- PDF pages: 62
- LaTeX log scan:
  - no fatal error
  - no undefined reference
  - no multiply-defined label
  - no overfull hbox

Non-blocking warning observed:

- SimSun bold italic font substitution warning.

## 8. Final Status

READY TO USE
