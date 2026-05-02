# RISKMAP Patch Report

## 1. Scope

本轮只优化 exp07 的谱分母可视化表现，不新增实验、不改变三个 case 的参数设置、不改理论结论，也不重跑其他大规模实验。

## 2. Files Modified

- `code/python/experiments/exp07_spectral_denominator_maps.py`
- `code/python/experiments/results/exp07_spectral_denominator_summary.csv`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.png`
- `code/python/experiments/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/figures/exp07_spectral_denominator_heatmaps.png`
- `thesis/figures/exp07_spectral_denominator_heatmaps.pdf`
- `thesis/chapters/6_experiments.tex`
- `thesis/main.pdf`

## 3. Figure Variable

图 10 现在不直接绘制

```text
d_{p,q} = lambda_{p,q} + sigma
```

而绘制小分母风险指标

```text
risk_{p,q} = -log10(max(|d_{p,q}|, 1e-16)).
```

该值越大表示 `|d_{p,q}|` 越小，对应模态越接近共振，小分母风险越高。该图不是条件数图，只用于突出 small-denominator risk。

## 4. Case Summary

| Case | sigma / setting | min \|d\| | Visual interpretation |
| --- | --- | ---: | --- |
| modified Helmholtz | `sigma=+100` | `1.197e+02` | 全部分母为正，整体风险低，无局部亮热点。 |
| true Helmholtz away from resonance | `sigma=-64.140291...` | `1.480e+01` | 分母正负分裂，但仍远离零；保留黑色 `d=0` contour 作为符号分界。 |
| true Helmholtz near resonance | target `(2,3)/(3,2)`, `delta=1e-2` | `1.000e-02` | 在 `(2,3)/(3,2)` 附近出现明显热点，并用 marker 与低模态 inset 标出。 |

CSV 中保留原有列，并追加：

- `risk_min`
- `risk_max`
- `risk_at_nearest_mode`

## 5. Thesis Integration

第 6 章实验四中，图 10 的正文和 caption 已改为“符号结构与小分母风险图”的口径：

- 明确图中绘制 `R_{p,q}=-log10|d_{p,q}|`；
- 说明 risk 越大表示小分母风险越高；
- 说明 modified、true-away、true-near 三种情形的差异；
- 若出现黑色等值线，解释为 `d=0` 的符号分界；
- 不将 risk 图解释为条件数图。

## 6. Validation

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

## 7. Readability Improvement

旧版 symlog `d_{p,q}` 热图容易被大范围正值背景淹没，三个 panel 视觉上接近。新版 risk map 使用共享 `inferno` 色标，直接把小分母区域映射成亮热点：

- modified 图保持低风险、无热点；
- true-away 图显示风险有限，即使存在符号分裂；
- true-near 图在目标模态 `(2,3)/(3,2)` 附近给出清晰热点，并通过 inset 展示低模态区域。

Final status: READY TO USE
