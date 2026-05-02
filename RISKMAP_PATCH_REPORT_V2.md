# RISKMAP Patch Report V2

## 1. Scope

本轮继续优化 exp07 谱分母可视化，只替换图 10 的表现形式，不新增实验、不改变 exp07 三个 case 的参数、不改变理论结论。

## 2. Files Modified

- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `code/python/experiments/results/exp07_spectral_denominator_summary.csv`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.png`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `thesis/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 3. Why the Full-Domain Risk Heatmap Was Replaced

全域 `127 x 127` risk map 中，高频模态的 `|d_{p,q}|` 普遍较大，导致大面积低风险背景占据视觉空间。真正有信息量的小分母区域集中在低频左下角，因此全域图虽然数学正确，但无法直接突出 near-resonance 热点。

新版图 10 改为：

- 第一行：低频窗口 `p,q <= 12` 的 risk map；
- 第二行：三种 case 的前 20 个最小 `|d_{p,q}|` 排序曲线。

这样同时展示“小分母出现在哪里”和“三种情形离零有多远”。

## 4. Quantity Plotted

低频 heatmap 仍使用小分母风险指标

```text
risk_{p,q} = -log10(max(|d_{p,q}|, 1e-16)).
```

该指标越大表示 `|d_{p,q}|` 越小，小分母风险越高。该图不是条件数图，也不声称 `risk` 等于 `cond_2(A)`。

## 5. Case Summary

| Case | Setting | min \|d\| |
| --- | --- | ---: |
| modified Helmholtz | `sigma=+100` | `1.197e+02` |
| true Helmholtz away from resonance | `sigma=-64.140291...` | `1.480e+01` |
| true Helmholtz near resonance | target `(2,3)/(3,2)`, `delta=1e-2` | `1.000e-02` |

The near-resonance sorted-denominator plot has its first two entries at:

- `(3,2): 1.000e-02`
- `(2,3): 1.000e-02`

Thus the first two minima are exactly the target degenerate modal pair `(2,3)/(3,2)`.

## 6. Thesis Integration

`thesis/chapters/6_experiments.tex` was updated around Figure 10:

- It now explains why full-domain `d_{p,q}` or `R_{p,q}` plots can hide the low-frequency small-denominator mechanism.
- The caption describes the two-row layout: low-frequency risk maps and sorted smallest-denominator comparison.
- The text states that modified Helmholtz remains far from zero, true-away is sign-changing but still away from resonance, and true-near pushes the target pair to the `10^{-2}` scale.

## 7. Validation

Commands run:

```powershell
$env:PYTHONPATH='code/python'; python code/python/experiments/exp07_spectral_denominator_maps.py
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Results:

- exp07 script: PASS
- XeLaTeX/BibTeX: PASS
- `main.pdf`: generated, 62 pages
- log scan: no fatal error, no undefined citation/reference, no multiply-defined label, no overfull hbox
- non-blocking warnings remain: SimSun font substitution and hyperref bookmark math-token warnings
- `git diff --check`: no whitespace errors; only existing CRLF normalization warnings

## 8. Final Status

READY TO USE
