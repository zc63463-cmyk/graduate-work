# FIG16 Presentation Patch Report

## 1. Scope

本轮只优化第 6 章实验五中图 16 的展示方式和正文解释。

未新增实验，未改变 near-resonance 参数，未重新扫描更多 `delta`，未更换 Gaussian RHS，未改写求解器。

## 2. Files Modified

- `code/python/experiments/exp05_true_helmholtz_resonance.py`
- `thesis/chapters/6_experiments.tex`
- `code/python/experiments/figures/exp05_multimode_resonance_summary.png`
- `code/python/experiments/figures/exp05_multimode_resonance_summary.pdf`
- `thesis/figures/exp05_multimode_resonance_summary.png`
- `thesis/figures/exp05_multimode_resonance_summary.pdf`
- `thesis/main.pdf`

## 3. Figure 16 Changes

图 16 仍保留 2x2 四个 panel：

1. **Panel (a): target pair amplification**
   - 分别标出目标简并模态 `(2,3)` 与 `(3,2)` 的模态系数；
   - 保留 log scale；
   - 保留参考线 `C/|delta|`；
   - 标注接近共振时约 `1000x` 的放大比例。

2. **Panel (b): target energy fraction**
   - 展示 target subspace energy fraction；
   - 标题改为 `target energy fraction (near 1)`；
   - 图中加入 `all values near 1` 提示；
   - y 轴聚焦在接近 1 的合理区间，用于说明目标子空间已经主导解场。

3. **Panel (c): dominant projection correlation**
   - 展示 dominant projection 与 full solution 的 shape correlation；
   - 标题改为 `dominant projection correlation (near 1)`；
   - 图中加入说明：full solution 已与 target subspace 高度一致。

4. **Panel (d): GMRES capped iterations**
   - 保留 GMRES(30) capped 汇总；
   - 标题改为 `GMRES capped iterations (all hit 1001)`；
   - 红叉表示 capped / not converged；
   - 图中提示详细停滞过程见 residual history figure。

## 4. Thesis Text Changes

第 6 章实验五中新增/强化了图 16 的读图说明：

- panel (a) 用于验证目标模态幅值随 `|delta|` 减小按 `1/|delta|` 放大；
- panel (b) 显示 target subspace energy 在当前 Gaussian RHS 设置下已接近 1；
- panel (c) 显示 dominant projection 与 full solution 的 shape correlation 已接近 1；
- panel (d) 只汇总 GMRES30 capped 结果，详细停滞过程见后续 residual history 图。

正文加入稳健表述：

> 在当前 Gaussian RHS 设置下，target subspace energy 与 shape correlation 在扫描区间内已接近 1，因此图 (b)(c) 的作用主要是说明 near-resonance 解场已被目标模态主导，而非展示显著的随 delta 变化趋势。

## 5. Data / Rerun Status

未重做 near-resonance 扫描。

图 16 仅基于已有文件重新绘制：

- `code/python/experiments/results/exp05_multimode_resonance.csv`

使用的重绘命令为：

```powershell
$env:PYTHONPATH='code/python'
python - <<'PY'
from pathlib import Path
import pandas as pd
from experiments import exp05_true_helmholtz_resonance as exp05

results_dir = Path(exp05.get_results_dir())
figures_dir = exp05.get_figures_dir()
df = pd.read_csv(results_dir / exp05.PATCH5_MULTIMODE_CSV)
exp05.plot_multimode_summary(df, figures_dir)
PY
```

没有运行完整 `exp05_true_helmholtz_resonance.py` 主实验，以避免重做 near-resonance 扫描。

## 6. Validation

LaTeX commands:

```powershell
cd thesis
xelatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

Results:

- XeLaTeX/BibTeX/XeLaTeX/XeLaTeX: success
- Output PDF: `thesis/main.pdf`
- PDF pages: 62
- Log scan:
  - no fatal error
  - no undefined reference
  - no multiply-defined label
  - no overfull hbox

Non-blocking warning:

- SimSun bold italic font substitution warning.

## 7. Final Status

READY TO USE
