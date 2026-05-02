from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "reports" / "ch6_experiment_report"
ASSET_DIR = OUT_DIR / "assets"
FIG_DIR = ROOT / "thesis" / "figures"
CSV_DIR = ROOT / "code" / "python" / "experiments" / "results"
HTML_PATH = OUT_DIR / "index.html"
BUILD_REPORT = OUT_DIR / "CH6_EXPERIMENT_REPORT_BUILD.md"
REFINEMENT_REPORT = OUT_DIR / "CH6_MATERIAL_REFINEMENT_REPORT.md"


FIGURES = [
    "exp01_convergence.png",
    "exp02_nonhom_bc.png",
    "exp02_temperature_field_comparison.png",
    "exp06_complex_manufactured_fields.png",
    "exp06_complex_manufactured_convergence.png",
    "exp03_neumann_mixed_summary.png",
    "exp03_neumann_mixed_fields.png",
    "exp06_accuracy_cost_error_time.png",
    "exp06_time_scaling.png",
    "exp07_spectral_denominator_heatmaps.png",
    "exp04_min_denom_vs_sigma.png",
    "exp04_spectral_indicator_vs_sigma.png",
    "exp04_condition_check.png",
    "exp04_gmres_iters_vs_sigma.png",
    "exp04_gmres_history.png",
    "exp05_near_resonance_summary.png",
    "exp05_multimode_resonance_summary.png",
    "exp05_dominant_mode_projection.png",
    "exp05_resonance_gmres_history.png",
]


CSVS = [
    "exp00_fft_vs_sparse.csv",
    "exp01_convergence.csv",
    "exp02_nonhom_bc.csv",
    "exp03_neumann_mixed.csv",
    "exp04_modified_vs_true.csv",
    "exp04_condition_check.csv",
    "exp04_gmres_history.csv",
    "exp05_resonance.csv",
    "exp05_multimode_resonance.csv",
    "exp05_resonance_gmres_history.csv",
    "exp06_accuracy_cost.csv",
    "exp06_complex_manufactured_convergence.csv",
    "exp07_spectral_denominator_summary.csv",
]


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value))


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def read_csv(name: str) -> list[dict[str, str]]:
    path = CSV_DIR / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def fmt(value: str | float | int | None) -> str:
    if value is None or value == "":
        return "-"
    try:
        x = float(value)
    except (TypeError, ValueError):
        return esc(value)
    if x == 0:
        return "0"
    ax = abs(x)
    if ax < 1e-3 or ax >= 1e4:
        return f"{x:.3e}"
    return f"{x:.4g}"


def table(headers: list[str], rows: list[list[object]], cls: str = "") -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = "\n".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def callout(title: str, body: str, kind: str = "info") -> str:
    return f"""
<div class="callout {esc(kind)}">
  <div class="callout-title">{esc(title)}</div>
  <div>{body}</div>
</div>
"""


def bullets(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{item}</li>" for item in items) + "</ul>"


def fig_card(filename: str | list[str], title: str, notes: dict[str, str], label: str = "") -> str:
    filenames = [filename] if isinstance(filename, str) else filename
    images = "\n".join(
        f'<img src="assets/{esc(name)}" alt="{esc(title)}">'
        for name in filenames
    )
    label_html = f'<div class="fig-label">{esc(label)}</div>' if label else ""
    note_html = "".join(
        f"<li><strong>{esc(k)}：</strong>{esc(v)}</li>"
        for k, v in notes.items()
    )
    return f"""
<figure class="figure-card">
  {label_html}
  <div class="figure-title">{esc(title)}</div>
  <div class="image-row">{images}</div>
  <ul class="note-list">{note_html}</ul>
</figure>
"""


def section(anchor: str, title: str, body: str) -> str:
    return f'<section id="{esc(anchor)}"><h2>{esc(title)}</h2>{body}</section>'


def find_learning_note() -> Path | None:
    for path in ROOT.glob("*.md"):
        try:
            text = path.read_text(encoding="utf-8-sig")
        except Exception:
            continue
        if "第六章数值实验学习笔记" in text[:800] and "频域分母" in text:
            return path
    downloads = Path.home() / "Downloads"
    if not downloads.exists():
        return None
    for path in downloads.glob("*.md"):
        try:
            text = path.read_text(encoding="utf-8-sig")
        except Exception:
            continue
        if "第六章数值实验学习笔记" in text[:500] and "频域分母" in text:
            return path
    return None


def copy_assets() -> tuple[list[str], list[str]]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    copied, missing = [], []
    for name in FIGURES:
        src = FIG_DIR / name
        if src.exists():
            shutil.copy2(src, ASSET_DIR / name)
            copied.append(name)
        else:
            missing.append(name)
    return copied, missing


def csv_status() -> tuple[list[str], list[str]]:
    present, missing = [], []
    for name in CSVS:
        if (CSV_DIR / name).exists():
            present.append(name)
        else:
            missing.append(name)
    return present, missing


def exp00_table() -> str:
    rows = read_csv("exp00_fft_vs_sparse.csv")
    compact = []
    for r in rows:
        if r.get("pass", "").lower() == "true":
            compact.append([
                esc(r.get("bc_type", "")),
                esc(r.get("sigma", "")),
                esc(r.get("equation_type", "")),
                esc(r.get("method", "")),
                fmt(r.get("error_vs_sparse") or r.get("linf_error")),
            ])
    return table(["边界", "sigma", "类型", "方法", "与 sparse direct 的误差"], compact[:10])


def exp01_table() -> str:
    rows = read_csv("exp01_convergence.csv")
    selected = []
    for r in rows:
        if r.get("rhs_type") == "polynomial" and r.get("sigma") in {"10", "10.0"}:
            if r.get("method") in {"FA", "FFT9"} and r.get("n") in {"33", "65", "129"}:
                selected.append([
                    esc(r.get("method", "")),
                    esc(r.get("n", "")),
                    fmt(r.get("error") or r.get("linf_error")),
                    fmt(r.get("rate") or r.get("observed_order_linf")),
                ])
    return table(["方法", "n", "L∞ 误差", "观测阶"], selected)


def exp03_table() -> str:
    rows = read_csv("exp03_neumann_mixed.csv")
    selected = []
    seen = set()
    for r in rows:
        key = (r.get("display"), r.get("n"), r.get("method"))
        if r.get("method") == "FA" and r.get("n") == "65" and key not in seen:
            seen.add(key)
            selected.append([
                esc(r.get("display", "")),
                esc(r.get("bc_type", "")),
                fmt(r.get("linf_error")),
                fmt(r.get("boundary_flux_linf")),
                fmt(r.get("weighted_mean_solution")),
            ])
    return table(["子情形", "边界", "L∞ 误差", "flux residual", "weighted mean"], selected)


def exp06_table() -> str:
    rows = read_csv("exp06_accuracy_cost.csv")
    selected = []
    for method in ["FA", "CR", "FACR", "FFT9", "GMRES30"]:
        candidates = [
            r for r in rows
            if r.get("method") == method and r.get("sigma") in {"10", "10.0"} and r.get("n") == "129"
        ]
        if candidates:
            r = candidates[0]
            selected.append([
                esc(method),
                esc(r.get("discretization", "")),
                fmt(r.get("time_s_median")),
                fmt(r.get("time_s_iqr")),
                fmt(r.get("linf_error")),
                esc(r.get("timing_scope", "")),
            ])
    return table(["方法", "离散", "median time", "IQR", "L∞ 误差", "timing scope"], selected, "small")


def exp07_table() -> str:
    rows = read_csv("exp07_spectral_denominator_summary.csv")
    selected = []
    for r in rows:
        selected.append([
            esc(r.get("case", "")),
            fmt(r.get("sigma")),
            fmt(r.get("d_min_abs")),
            esc(r.get("nearest_p", "")),
            esc(r.get("nearest_q", "")),
            esc(r.get("num_positive", "")),
            esc(r.get("num_negative", "")),
        ])
    return table(["case", "sigma", "min |d|", "nearest p", "nearest q", "# positive", "# negative"], selected)


def exp04_condition_table() -> str:
    rows = read_csv("exp04_condition_check.csv")
    if not rows:
        return ""
    max_rel = max(float(r.get("relative_difference", "0") or 0) for r in rows)
    return callout(
        "condition check 数值摘要",
        f"在当前 Dirichlet 五点系统中，dense cond₂(A) 与 spectral denominator indicator 的最大相对差约为 <strong>{max_rel:.3e}</strong>。这说明二者在该正交 DST-I 对角化系统下数值一致，但不能外推到 non-normal 或 mixed/Neumann 原始矩阵。",
        "ok",
    )


def exp05_table() -> str:
    rows = read_csv("exp05_multimode_resonance.csv")
    selected = []
    for r in rows:
        if r.get("delta") in {"0.0001", "1e-04", "1.0e-04"}:
            selected.append([
                esc(r.get("mode_label", "")),
                esc(r.get("target_modes", "")),
                fmt(r.get("target_uhat_norm")),
                fmt(r.get("modal_amplitude_relative_error")),
                fmt(r.get("shape_correlation")),
                fmt(r.get("gmres_final_abs_residual")),
            ])
    return table(["mode label", "target modes", "target |û|", "理论相对误差", "shape corr.", "GMRES final absres"], selected, "small")


def build_html() -> str:
    overview_table = table(
        ["实验", "核心问题", "主要证据", "一句话结论", "不回答的问题"],
        [
            ["exp00", "实现路径是否一致", "FFT vs sparse direct 表", "只是 implementation validation", "不回答连续 PDE 精度"],
            ["实验一", "Dirichlet 格式阶数是否达标", "收敛图、非齐次边界可视化", "五点二阶，FFT9 四阶", "不证明 FFT9 支持 Neumann/mixed"],
            ["实验二", "Neumann/mixed 边界是否闭合", "误差、flux、weighted mean", "边界残差和零模态诊断闭合", "不做性能比较"],
            ["实验三", "FFT 与 GMRES 如何对比", "error-time、time scaling", "FFT9 在光滑 Dirichlet 下更优", "不是严格 kernel benchmark"],
            ["实验四", "Helmholtz 为什么影响 GMRES", "risk map、condition check、residual history", "true Helmholtz 小分母导致谱风险", "不外推到所有系统"],
            ["实验五", "近共振具体放大哪个模态", "模态投影、dominant projection、GMRES history", "目标模态按 1/|delta| 放大", "不是连续谱极限结论"],
        ],
    )

    body = ""
    body += section(
        "proof-chain",
        "0. 第六章到底在证明什么",
        f"""
        {callout("主线", "第六章不是“跑几张图”，而是按 <strong>离散格式 → 边界处理 → 精度/成本 → 频域谱结构 → near-resonance 模态放大 → GMRES 停滞</strong> 的链条组织证据。")}
        <div class="chain">
          <span>实现校验</span><span>Dirichlet 收敛</span><span>边界处理</span><span>精度--成本</span><span>谱分母</span><span>模态放大</span><span>GMRES 停滞</span>
        </div>
        {overview_table}
        {callout("答辩时的总说法", "第六章先确认实现路径和离散精度，再比较 FFT 直接法与无预处理 GMRES30 的误差--时间表现，最后用频域分母解释 true Helmholtz near-resonance 下的模态放大和 GMRES 残差停滞。", "ok")}
        """,
    )

    body += section(
        "equations",
        "1. 统一方程与频域分母小抄",
        f"""
        <div class="formula-row">
          <div class="formula">(-∇² + σ)u = f</div>
          <div class="formula">d<sub>p,q</sub> = λ<sub>p,q</sub> + σ</div>
          <div class="formula">û<sub>p,q</sub> = f̂<sub>p,q</sub> / d<sub>p,q</sub></div>
        </div>
        {table(["类型", "sigma", "分母", "直觉"], [
            ["Poisson", "0", "d = λ", "标准椭圆问题"],
            ["modified Helmholtz", "+κ²", "d = λ + κ²", "分母全正且更远离 0"],
            ["true Helmholtz", "-κ²", "d = λ - κ²", "可能穿过离散谱，产生小分母"],
        ])}
        {callout("小分母因果链", "若某个 d<sub>p,q</sub> 很小，则 û<sub>p,q</sub>=f̂<sub>p,q</sub>/d<sub>p,q</sub> 被放大；这会增加解范数和系统病态程度，并使无预处理 GMRES30 更容易出现残差停滞。")}
        {callout("不能误读", "near-resonance 是当前有限维 Dirichlet 五点离散谱上的实验，不是连续谱极限结论。", "warn")}
        """,
    )

    body += section(
        "how-to-read",
        "2. 两个基本读图规则",
        f"""
        <div class="two-col">
          <div>
            <h3>收敛阶怎么看</h3>
            <p>若误差满足 ||u<sub>h</sub>-u||≈Ch<sup>p</sup>，则 log-log 图中的斜率就是收敛阶 p。网格加密一倍时，二阶误差约缩小 4 倍，四阶误差约缩小 16 倍。</p>
            <p>因此看收敛图时，不要只盯某个误差数字，而要看曲线是否贴近 O(h²) 或 O(h⁴) 参考线。</p>
          </div>
          <div>
            <h3>error-time 图怎么看</h3>
            <p>横轴是计算时间，纵轴是误差，越靠左下越好。FFT9 的优势来自四阶离散误差，而 GMRES30 曲线同时受到离散误差和代数残差影响。</p>
            <p>本论文的 timing scope 不是严格 kernel-to-kernel benchmark：direct solvers 计时包含 solver call 内部 RHS/transform setup，GMRES30 计时排除了 matrix/RHS setup。</p>
          </div>
        </div>
        """,
    )

    body += section(
        "exp00",
        "3. 实现校验：exp00 为什么不算核心实验",
        f"""
        <p>exp00 比较 FFT solver 与 sparse direct solver 的结果。如果误差在 10<sup>-14</sup> 量级，说明两条求解路径对齐到同一个离散线性系统。</p>
        {exp00_table()}
        {callout("它证明什么", "代码求解路径没有写错，FFT 路径与 sparse direct 在同一离散系统上对齐。", "ok")}
        {callout("它不证明什么", "它不证明连续 PDE 精度、不证明 FFT9 更准、不评价 GMRES 快慢，也不解释 true Helmholtz 难点。", "warn")}
        """,
    )

    body += section(
        "exp1",
        "4. 实验一：Dirichlet 离散格式收敛性验证",
        f"""
        <h3>为什么做</h3>
        <p>本实验回答两个基础问题：五点格式是否达到二阶，FFT9 是否达到四阶。主测试采用 polynomial manufactured solution，避免单一 Fourier 模态过于“顺着 FFT 方法”。</p>
        <h3>方法定位</h3>
        {table(["方法", "离散模板", "理论阶数", "边界范围"], [
            ["FA / CR / FACR-like", "五点格式", "O(h²)", "D/N/mixed 五点求解器"],
            ["FFT9", "九点紧致格式", "O(h⁴)", "本文只实现并验证 Dirichlet"],
        ])}
        {exp01_table()}
        {fig_card("exp01_convergence.png", "Dirichlet 收敛验证", {
            "看什么": "看 log-log 曲线斜率，而不是单个误差值。",
            "关键现象": "FA/CR/FACR-like 贴近二阶参考线，FFT9 贴近四阶参考线。",
            "理论对应": "五点截断误差为 O(h²)，FFT9 紧致修正达到 O(h⁴)。",
            "不要误读": "FFT9 四阶结论只限本文实现的 Dirichlet 情形。"
        }, "图 1")}
        <h3>非齐次 Dirichlet 子情形</h3>
        <p>非齐次 Dirichlet 主要检查边界值 g<sub>D</sub> 移到右端项时符号和模板处理是否正确。它保留在正文中作为 boundary visual sanity check，不再单独占用核心实验编号。</p>
        {fig_card(["exp02_nonhom_bc.png", "exp02_temperature_field_comparison.png"], "非齐次 Dirichlet 边界与温度场可视化", {
            "看什么": "解析场、数值场和误差分布是否一致，边界附近是否异常。",
            "关键现象": "非齐次边界贡献处理后，五点仍二阶，FFT9 仍四阶。",
            "理论对应": "边界已知值从离散模板移到右端项，符号必须正确。",
            "不要误读": "这两张图是边界 sanity check，收敛阶主证据仍是前面的收敛图和表。"
        }, "图 2--3")}
        {fig_card(["exp06_complex_manufactured_fields.png", "exp06_complex_manufactured_convergence.png"], "多模态 manufactured solution 可视化补充", {
            "看什么": "复杂空间结构下 exact、FA、FFT9 和误差分布是否合理。",
            "关键现象": "旧版 max error≈1 的异常已排除，FA 二阶、FFT9 四阶。",
            "理论对应": "多个正弦模态逐项构造 RHS，避免统一系数错误。",
            "不要误读": "这是可视化补充，不是新增核心实验。"
        }, "补充图")}
        {callout("答辩一句话", "实验一说明：Dirichlet 情形下五点格式稳定二阶，FFT9 稳定四阶；非齐次边界图只是确认边界修正项符号和模板处理没有错误。", "ok")}
        """,
    )

    body += section(
        "exp2",
        "5. 实验二：Neumann 与 mixed 边界处理",
        f"""
        <h3>为什么独立成实验</h3>
        <p>Neumann 给定的是法向导数，不是函数值。ghost-point 会改变边界行结构，Pure Neumann Poisson 还存在常数零模态，因此本实验重点不是重复五点二阶，而是验证边界约束、DCT-I 对称化和零模态处理。</p>
        {exp03_table()}
        {fig_card("exp03_neumann_mixed_summary.png", "Neumann/mixed boundary verification summary", {
            "看什么": "panel (a) 看 FA 在不同边界组合下仍二阶；panel (b) 看 flux 和 Dirichlet residual；panel (c) 看 zero-mode。",
            "关键现象": "误差收敛、边界残差和 weighted mean 诊断同时闭合。",
            "理论对应": "ghost-point + DCT-I 对称化给出可对角化求解路径。",
            "不要误读": "CR/FACR-like 的重合曲线保留在 CSV 中，图中只画 FA 以避免视觉冗余；不声称 FFT9 支持 Neumann/mixed 四阶。"
        }, "图 6")}
        {fig_card("exp03_neumann_mixed_fields.png", "Neumann/mixed 代表性解场与误差", {
            "看什么": "代表性 exact、numerical 和 log-error 空间结构。",
            "关键现象": "数值场重现了解析场，误差集中和边界行为合理。",
            "理论对应": "边界条件处理后的五点系统仍能正确求解。",
            "不要误读": "空间图只用于直观展示，收敛和边界闭合仍由 summary 图与 CSV 支撑。"
        }, "图 7")}
        {callout("答辩一句话", "实验二验证的是边界处理闭合：Neumann/mixed 下误差、边界 residual 和零模态诊断都合理，但不声明 FFT9 的 Neumann/mixed 四阶性。", "ok")}
        """,
    )

    body += section(
        "exp3",
        "6. 实验三：精度--成本对比与复杂度实证",
        f"""
        <h3>为什么这是核心实验</h3>
        <p>论文标题包含“FFT 快速求解器与迭代法对比研究”。如果没有 error-time 图，就只是分别验证了一些方法，而没有真正做对比。</p>
        {exp06_table()}
        {fig_card("exp06_accuracy_cost_error_time.png", "accuracy-cost：误差与计算时间", {
            "看什么": "横轴 median solve time，纵轴 L∞ error，二者均为 log scale；越靠左下越好。左 panel 是 sigma=0，右 panel 是 sigma=10。",
            "关键现象": "FFT9 在相近 FFT 复杂度下给出显著更低误差；五点 direct solvers 误差相近；GMRES30 更靠右且可能未达到容差。",
            "理论对应": "FFT9 的优势来自 Dirichlet 九点紧致四阶离散；FA/CR/FACR-like 与 GMRES30 都受五点二阶离散误差限制。",
            "不要误读": "FFT9 不是把同一个五点系统解得更准；timing scope 不是严格 kernel benchmark，GMRES setup 被排除，因此并未刻意偏向 FFT。"
        }, "图 8")}
        <div class="grid two">
          <div>
            <h3>坐标轴和视觉编码怎么读？</h3>
            <ul>
              <li><strong>横轴</strong>：median solve-call time，越左越快。</li>
              <li><strong>纵轴</strong>：L∞ error，越下越准。</li>
              <li><strong>形状/边框</strong>：算法类别；<strong>颜色条</strong>：未知量 $N^2$ 的对数规模。</li>
              <li><strong>连线</strong>：同一算法随网格加密的轨迹，通常往右下移动。</li>
            </ul>
          </div>
          <div>
            <h3>两张 panel 分别说明什么？</h3>
            {table(["panel", "参数", "含义"], [
                ["左图", "sigma=0", "Poisson 基准问题"],
                ["右图", "sigma=10", "正定 modified Helmholtz 基准问题"],
            ], "compact")}
            <p>这不是 true Helmholtz near-resonance 实验；near-resonance 的小分母和 GMRES 停滞由后续实验四、五解释。</p>
          </div>
        </div>
        {table(["方法", "离散/求解口径", "读图时的定位"], [
            ["FA / CR / FACR-like", "五点二阶 direct", "若线性系统解准，PDE 误差主要是 O(h²) 离散误差"],
            ["FFT9", "Dirichlet 九点紧致四阶 direct", "更低误差来自 O(h⁴) 离散格式，不是同一个五点系统更准"],
            ["GMRES30", "五点二阶 unpreconditioned restarted GMRES(30)", "误差同时受五点离散误差和代数残差影响"],
        ], "compact")}
        {callout("为什么“越靠左下越好”", "同样时间下谁更低谁更准；同样误差下谁更左谁更快。FFT9 的关键形态是横向没有付出数量级更高时间，纵向却降低了多个数量级误差。", "ok")}
        {callout("答辩时的边界", "这张图只说明当前规则 Dirichlet 光滑 manufactured solution、当前实现和当前计时口径下的精度--时间表现。不能说 FFT9 适用于所有边界，也不能说 GMRES 或预处理 Krylov 方法不行。", "warn")}
        <h3>常见追问怎么答</h3>
        {table(["问题", "稳妥回答"], [
            ["为什么 FFT9 更准？", "因为它是 Dirichlet 九点紧致四阶离散，误差为 O(h⁴)，而五点方法是 O(h²)。"],
            ["FFT9 是不是同一五点系统解得更准？", "不是。FFT9 求解九点紧致离散系统，优势来自更高阶离散格式。"],
            ["GMRES30 为什么不能达到 FFT9 误差？", "GMRES30 作用于五点离散矩阵，即使代数残差很小，PDE 误差仍受五点二阶离散限制。"],
            ["为什么 timing 不是严格 kernel benchmark？", "direct solvers 计时包含内部 RHS/transform setup；GMRES30 计时排除 matrix/RHS setup，因此它是当前实现趋势图。"],
        ], "compact")}
        {fig_card("exp06_time_scaling.png", "time scaling：当前实现的规模趋势", {
            "看什么": "不同方法随未知量规模增长的时间趋势。",
            "关键现象": "FFT 类直接法呈 O(N² log N) 量级趋势；GMRES30 作为无预处理基准成本增长更敏感。",
            "理论对应": "规则矩形上 DST/FFT 对角化带来快速直接求解。",
            "不要误读": "FACR-like 当前实现仍按 O(N² log N) 处理，不声称达到经典 FACR 的 O(N² log log N)。"
        }, "图 9")}
        {callout("答辩一句话", "实验三说明：在规则 Dirichlet 光滑问题上，FFT9 用同量级 FFT 成本获得更低误差；GMRES30 是无预处理基线，不能外推到预处理 Krylov 方法。", "ok")}
        """,
    )

    body += section(
        "exp4",
        "7. 实验四：Modified/True Helmholtz 的谱结构与 GMRES 行为",
        f"""
        <h3>核心问题</h3>
        <p>实验四解释 true Helmholtz 为什么比 modified Helmholtz 更难：根源是频域分母 d<sub>p,q</sub>=λ<sub>p,q</sub>+σ 的符号和大小。</p>
        {exp07_table()}
        {fig_card("exp07_spectral_denominator_heatmaps.png", "低频小分母 risk map 与最小分母排序", {
            "看什么": "上排看 p,q≤12 的 risk 热点在哪里；下排看前 20 个最小 |d| 的数量级。",
            "关键现象": "modified 的 min|d|≈1.197e2，true-away 的 min|d|≈1.480e1，true-near 的 (2,3)/(3,2) 降至 1.000e-2。",
            "理论对应": "R=-log10|d| 越大表示越接近小分母；小分母解释后续模态放大和 GMRES 困难。",
            "不要误读": "R<0 只表示 |d|>1、风险低；risk map 不是条件数图，sign-changing 也不等于共振。"
        }, "图 10")}
        <div class="grid two">
          <div>
            <h3>图 10 到底画什么？</h3>
            <p>它不是物理空间解场，也不是误差，而是五点 Dirichlet 离散系统的频域分母风险。DST 对角化后，每个模态满足 û<sub>p,q</sub>=f̂<sub>p,q</sub>/d<sub>p,q</sub>，所以 |d| 越小，该模态越容易被放大。</p>
            <p>之所以只画低频窗口，是因为高频模态的 |d| 普遍很大；若画完整 127×127 频域，真正危险的低频小分母会被大范围低风险背景淹没。</p>
          </div>
          <div>
            <h3>三种 case 的读图结论</h3>
            {table(["case", "min |d|", "读法"], [
                ["modified", "1.197e2", "分母全正且远离零，无热点"],
                ["true-away", "1.480e1", "已 sign-changing，但不是 near-resonance"],
                ["true-near", "1.000e-2", "(2,3)/(3,2) 是危险小分母"],
            ], "compact")}
            <p>下排排序图比颜色更直接：modified 和 true-away 的最低点仍远离 1，true-near 的前两个点突然降到 1e-2。</p>
          </div>
        </div>
        {callout("为什么 near-resonance 前两个点并列？", "因为正方形 Dirichlet 五点离散谱中 lambda_23^h=lambda_32^h。当 kappa² 接近这个特征值时，(2,3) 和 (3,2) 两个模态的分母同时变为约 delta=1e-2；第 3 个最小分母来自相邻模态，明显跳升。", "info")}
        {callout("答辩防误读", "risk map 只可视化 small-denominator risk，不等于条件数图；若有 d=0 contour，它只是插值符号分界，不表示离散零分母模态。true-away 的分母正负混合说明系统可能不定，但不等于近共振。", "warn")}
        {fig_card(["exp04_min_denom_vs_sigma.png", "exp04_spectral_indicator_vs_sigma.png"], "最小分母与 spectral denominator indicator", {
            "看什么": "看 min |d| 何时变小，以及 max|d|/min|d| 何时变大。",
            "关键现象": "true Helmholtz 在接近离散谱时指标迅速增大。",
            "理论对应": "小分母使谱分母指标变大。",
            "不要误读": "该指标用于解释 Dirichlet 五点系统的小分母风险，不能无条件外推。"
        }, "图 11")}
        {exp04_condition_table()}
        {fig_card("exp04_condition_check.png", "condition check：谱分母指标与 dense cond₂(A)", {
            "看什么": "横轴 spectral denominator indicator，纵轴 dense cond₂(A)，点是否贴近 y=x。",
            "关键现象": "点落在 y=x 附近。",
            "理论对应": "Dirichlet 五点矩阵由正交 DST-I 对角化，非奇异时 max|d|/min|d| 等于 2-范数条件数。",
            "不要误读": "这不是新条件数理论，只是当前正交对角化系统的数值校验。"
        }, "图 12")}
        {fig_card(["exp04_gmres_iters_vs_sigma.png", "exp04_gmres_history.png"], "GMRES iteration 与 residual history", {
            "看什么": "iteration 图看哪些点 capped；history 图看 residual 如何下降或停滞。",
            "关键现象": "modified 的某些参数残差下降更快，true/near-resonance 出现停滞或 capped。",
            "理论对应": "谱分母指标增大与无预处理 GMRES30 困难一致。",
            "不要误读": "501 是 max_iter=500 下 capped-count convention，不是成功收敛 501 步；停止准则是 absolute residual tolerance。"
        }, "图 13--14")}
        {callout("答辩一句话", "实验四说明：modified Helmholtz 分母全正且远离零；true Helmholtz 可能出现正负分裂和小分母，小分母会放大谱指标并使无预处理 GMRES30 更难把残差降到绝对容差。", "ok")}
        """,
    )

    body += section(
        "exp5",
        "8. 实验五：True Helmholtz near-resonance 模态放大",
        f"""
        <h3>实验四与实验五的区别</h3>
        <p>实验四说明“小分母存在，并影响 GMRES”；实验五进一步回答“小分母具体放大哪个模态、是否符合 û=f̂/(λ-κ²)、放大后的解场是否由目标模态主导”。</p>
        <div class="formula-row">
          <div class="formula">û<sub>p,q</sub> = f̂<sub>p,q</sub> / (λ<sub>p,q</sub> - κ²)</div>
          <div class="formula">κ² = λ<sub>target</sub> + δ</div>
          <div class="formula">|û<sub>target</sub>| ≈ |f̂<sub>target</sub>| / |δ|</div>
        </div>
        {fig_card("exp05_near_resonance_summary.png", "baseline near-resonance summary", {
            "看什么": "panel (a) 看 1/|delta| 放大；panel (b) 看 final absolute residual 是否低于容差；panel (c)(d) 看同一模态分支下 raw field 的幅值变化。",
            "关键现象": "delta 从 1e-1 到 1e-4 时，max |u| 约从 5.62e-1 增至 5.64e2，约放大 10^3 倍；GMRES 红叉表示 capped/not converged。",
            "理论对应": "目标分母 λ_target-(λ_target+delta)=-delta，因此代表模态系数按 1/|delta| 放大。",
            "不要误读": "下排形状相似不是矛盾：它们属于同一 near-resonance branch，主要差异是幅值；这是离散谱实验，不是连续谱极限。"
        }, "图 15")}
        <div class="grid two">
          <div>
            <h3>图 15 四个 panel 怎么读？</h3>
            {table(["Panel", "看什么", "回答的问题"], [
                ["(a)", "||u||₂(RMS)、|û₂,₃| 与 C/|delta|", "是否符合 1/|delta| 放大"],
                ["(b)", "final absolute residual 与 absolute tolerance", "无预处理 GMRES 是否达到容差"],
                ["(c)", "delta=1e-1 的 raw field", "离共振较远时解场幅值"],
                ["(d)", "delta=1e-4 的 raw field", "接近共振时解场幅值"],
            ], "compact")}
          </div>
          <div>
            <h3>下排 raw fields 的关键量级</h3>
            {table(["detuning", "min |d|", "max |u|", "读法"], [
                ["delta=1e-1", "1e-1", "5.62e-1", "幅值较小"],
                ["delta=1e-4", "1e-4", "5.64e2", "幅值约放大 10^3 倍"],
            ], "compact")}
            <p>两幅 raw fields 使用独立 colorbar，所以颜色强弱不能直接比较；真正要比较的是标题中的 max |u|。</p>
          </div>
        </div>
        {callout("横轴怎么读", "图 (a)(b) 的横轴是 |delta| 的对数坐标；从 1e-1 到 1e-4 表示越来越接近离散共振。", "info")}
        {callout("为什么下排形状相似？", "这是同一目标模态对 (2,3)/(3,2) 下的 detuning 对比。near-resonance 改变的主要是目标模态幅值，而不是模态形状本身，所以空间结构可以相似，但 max |u| 相差约 1000 倍。", "warn")}
        {callout("这张图不能证明什么", "它不证明所有 Helmholtz 求解器都会失败，也不说明预处理 Krylov 方法无效；结论限定为当前 Dirichlet 五点离散系统、当前 RHS 和无预处理 GMRES(30)。", "info")}
        {exp05_table()}
        {fig_card("exp05_multimode_resonance_summary.png", "multimode resonance summary", {
            "看什么": "panel (a) 当前专门看简并对 (2,3)/(3,2) 的两个模态系数；panel (b)(c)(d) 比较 (1,1)、(2,3)/(3,2)、(3,3) 三类目标集合。",
            "关键现象": "两个简并系数都沿 1/|delta| 放大；三组目标的 energy fraction 最小仍约 0.996，shape correlation 最小仍约 0.998；GMRES 全部 capped=1001。",
            "理论对应": "对目标模态有 û_target=f̂_target/(λ_target-κ²)=f̂_target/(-delta)，因此幅值按 1/|delta| 放大。",
            "不要误读": "panel (a) 不是三组模态同图比较，而是展开简并对内部两个系数；panel (b)(c) 接近平线不是缺少信号，而是说明 near-resonance 解已经被目标模态主导。"
        }, "图 16")}
        {table(["目标集合", "作用", "读图口径"], [
            ["(1,1)", "最低频单模态", "说明低频单模态也遵循小分母放大机制"],
            ["(2,3)/(3,2)", "简并双模态", "panel (a) 分别展示两个系数；二者高度可不同，但斜率同为 1/|delta|"],
            ["(3,3)", "较高频单模态", "说明机制不只出现在最低频或简并模态上"],
        ], "compact")}
        {callout("为什么 panel (b)(c) 接近 1 不是没信息", "panel (b) 是能量账本，说明完整解能量几乎都落在目标子空间；panel (c) 是图像相似度，说明 full solution 的空间形状几乎就是 dominant projection。它们的任务不是展示明显趋势，而是证明支配性已经建立。", "info")}
        {callout("和相邻图的关系", "图 10 说明小分母在哪里；图 16 验证不同目标模态被小分母放大并主导解场；dominant projection 图把 shape correlation≈1 画成物理空间对比；GMRES history 图再解释 capped 背后的残差停滞过程。", "ok")}
        {fig_card("exp05_dominant_mode_projection.png", "dominant mode projection", {
            "看什么": "左图是 full solution，中图是只保留 (2,3)/(3,2) 后重构的 dominant projection，右图是 full - dominant。",
            "关键现象": "full 与 dominant 的最大幅值均约 6.81e1，而差值最大约 6.89e-4，相对量级约 1e-5。",
            "理论对应": "near-resonance 下目标分母为 -delta，目标简并模态系数远大于其他模态，故物理空间解场由目标子空间主导。",
            "不要误读": "右图颜色明显不代表误差大；它使用独立色标，必须看标题中的 max |u| 数量级。"
        }, "图 17")}
        <div class="grid two">
          <div>
            <h3>dominant projection 怎么算？</h3>
            <ol>
              <li>先对完整解在 interior 网格上做正交 DST-I 展开，得到所有模态系数。</li>
              <li>只保留目标模态 $(2,3)$ 与 $(3,2)$，其余模态系数置零。</li>
              <li>反变换回物理空间，得到 u<sub>dominant</sub>。</li>
            </ol>
            <p>因此它不是另一个求解器的输出，而是一个诊断量：检查 near-resonance 的完整解能否被目标模态子空间解释。</p>
          </div>
          <div>
            <h3>这张图最关键的量级</h3>
            {table(["对象", "最大幅值", "含义"], [
                ["full solution", "$6.81\\times10^1$", "完整离散谱解"],
                ["dominant projection", "$6.81\\times10^1$", "仅保留 $(2,3)/(3,2)$ 后重构"],
                ["full - dominant", "$6.89\\times10^{-4}$", "未被目标子空间解释的剩余部分"],
            ], "compact")}
            <p>差值与主幅值之比约为 10<sup>-5</sup>，说明目标模态子空间解释了几乎全部空间结构。</p>
          </div>
        </div>
        {callout("右图为什么还有红蓝颜色？", "右图使用自己的 colorbar，色标量级约为 10<sup>-4</sup>；左、中图的色标量级约为 10<sup>1</sup>。所以右图有颜色不代表误差大。如果右图使用左、中图相同色标，几乎会是一片接近零的颜色。", "warn")}
        {callout("为什么差值不是严格 0？", "RHS 是偏心 Gaussian，不是纯目标模态，因此在非目标模态上也有投影；同时邻近模态的分母虽不为零，也可能被有限放大。差值约 10<sup>-5</sup> 的相对量级是可解释的非目标模态贡献，不是投影操作失败。", "info")}
        {fig_card("exp05_resonance_gmres_history.png", "near-resonance GMRES residual history", {
            "看什么": "横轴是 GMRES iteration，纵轴是 log scale residual；左图 absolute residual 对应实际停止准则，右图 relative residual 仅用于归一化比较。",
            "关键现象": "三条 near-resonance 曲线不是完全不动，而是先有限下降，随后进入平台；最终 residual 约 1.83e-1，仍远高于 1e-10 绝对残差容差。",
            "理论对应": "small denominator → modal amplification → near-singular/ill-conditioned discrete system → unpreconditioned GMRES(30) stagnation。",
            "不要误读": "capped=1001 是 max_iter=1000 下的 capped-count convention，不是成功收敛到第 1001 步；这也不代表所有 Helmholtz 求解器都不行。"
        }, "图 18")}
        {callout("这张图回答什么问题", "它把 summary 图里的 capped/not converged 拆成全过程：GMRES 是否完全不动、残差从哪里开始变平、capped 是终点标签还是残差平台化的结果。", "info")}
        {callout("为什么 residual 降了仍不算收敛", "本实验停止准则是 absolute residual tolerance = 1e-10；当前 history 的末端约为 1.83e-1，距离容差仍有多个数量级，因此只能解释为未满足当前停止准则。", "warn")}
        {callout("边界口径", "该图只说明当前 Dirichlet 五点离散、目标模态 (2,3)/(3,2)、当前 RHS、无预处理 restarted GMRES(30) 下的停滞过程；不是所有 Helmholtz 求解器或预处理 Krylov 方法的结论。", "warn")}
        {callout("答辩一句话", "实验五用离散模态投影公式直接验证：当 κ² 接近目标离散特征值时，目标模态按 1/|delta| 放大，并使当前无预处理 GMRES30 出现残差停滞。", "ok")}
        """,
    )

    body += section(
        "qa",
        "9. 最容易被问的问题",
        table(
            ["问题", "稳健回答"],
            [
                ["为什么 exp00 不算核心实验？", "它只验证同一离散系统的不同求解路径一致，不回答连续 PDE 精度。"],
                ["为什么 FFT9 比五点更准？", "FFT9 是 Dirichlet 九点紧致四阶格式；五点格式是二阶。"],
                ["FFT9 能用于 Neumann/mixed 吗？", "本文只实现并验证 Dirichlet FFT9；Neumann/mixed 四阶紧致边界模板是后续工作。"],
                ["GMRES 为什么慢？", "在 true Helmholtz near-resonance 下，小分母使系统病态，当前无预处理 GMRES30 残差下降困难。"],
                ["risk map 是条件数图吗？", "不是。risk map 是 small-denominator risk visualization；condition check 是另一张图。"],
                ["near-resonance 证明了什么？", "证明当前离散谱系统中目标模态按 1/|delta| 放大；不代表连续谱极限。"],
                ["为什么图 16 的 panel (b)(c) 是平的？", "因为 target subspace energy fraction 和 shape correlation 已接近 1，说明解已被目标模态主导；这不是缺少信号。"],
                ["GMRES 计时口径是否偏向 FFT？", "GMRES 的 matrix/RHS setup 被排除，计时条件更宽容；在这种口径下 FFT 仍占优，反而加强结论。"],
            ],
        ),
    )

    body += section(
        "dont-say",
        "10. 绝对不能乱说的话",
        f"""
        {table(["不要说", "应该说"], [
            ["GMRES 不适合 Helmholtz", "当前无预处理 GMRES30 在 true Helmholtz near-resonance 下容易停滞；预处理方法不在本文实验范围。"],
            ["FFT9 支持 Neumann/mixed 四阶", "本文 FFT9 只实现并验证 Dirichlet；Neumann/mixed 四阶紧致格式是后续方向。"],
            ["FACR-like 达到 O(N² log log N)", "经典 FACR 理论可达 O(N² log log N)，本文 FACR-like 当前实现仍按 O(N² log N) 处理。"],
            ["risk map 是条件数图", "risk map 只可视化 R=-log10|d| 的小分母风险；condition check 另行验证 Dirichlet 五点下的 cond₂ 等价。"],
            ["near-resonance 证明连续谱结论", "near-resonance 实验基于有限维 Dirichlet 五点离散谱。"],
            ["error-time 是严格公平 kernel benchmark", "该图展示当前实现下的整体基准趋势；GMRES setup 被排除，比较没有刻意偏向 FFT。"],
            ["图 16(b)(c) 展示随 delta 的显著趋势", "它们主要证明支配性：energy fraction 和 shape correlation 已接近 1，说明解由目标模态主导。"],
            ["exp04 和 exp05 的 capped 计数一样", "exp04 中 max_iter=500，capped 记为 501；exp05 中 max_iter=1000，capped 记为 1001。"],
        ])}
        """,
    )

    body += section(
        "review-plan",
        "11. 三天复习路线",
        f"""
        {table(["时间", "目标", "要掌握的句子"], [
            ["第一天", "只理解主线", "d=lambda+sigma；modified 分母全正，true 可能接近 0；小分母导致模态放大和 GMRES 停滞。"],
            ["第二天", "逐张读核心图", "图 1 看斜率，图 6 看三类诊断，error-time 看左下角，图 10 看热点和排序，图 16 看 1/|delta|。"],
            ["第三天", "练答辩话术", "每个实验都用“问题--方法--图--结论--限制”五步讲，不越界。"],
        ])}
        {callout("最终理解", "第六章的闭环是：FFT/DST/DCT 给出快速直接解；五点与 FFT9 分别提供二阶和四阶精度；Neumann/mixed 需要边界和零模态处理；accuracy-cost 兑现 FFT vs GMRES 对比；true Helmholtz 的小分母解释近共振模态放大和无预处理 GMRES30 停滞。", "ok")}
        """,
    )

    css = """
    :root {
      --bg: #f5f7fb;
      --paper: #ffffff;
      --ink: #1f2937;
      --muted: #64748b;
      --line: #dbe4ef;
      --blue: #1f5f9f;
      --teal: #02847d;
      --green: #258f61;
      --warn: #9a5b00;
      --warn-bg: #fff7e6;
      --ok-bg: #ecfdf3;
      --info-bg: #eef6ff;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: "Microsoft YaHei", "Noto Sans CJK SC", "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--bg);
      line-height: 1.68;
    }
    header {
      padding: 54px 7vw 38px;
      color: white;
      background: linear-gradient(120deg, #163250, #1f5f9f 58%, #02847d);
    }
    header h1 { margin: 0 0 12px; font-size: 36px; letter-spacing: 0; }
    header p { max-width: 980px; margin: 8px 0 0; color: #e9f4ff; font-size: 18px; }
    nav {
      position: sticky;
      top: 0;
      z-index: 5;
      padding: 10px 7vw;
      background: rgba(255,255,255,0.96);
      border-bottom: 1px solid var(--line);
      backdrop-filter: blur(8px);
    }
    nav a { margin-right: 16px; color: var(--blue); text-decoration: none; font-weight: 600; font-size: 14px; }
    main { max-width: 1180px; margin: 0 auto; padding: 30px 20px 70px; }
    section {
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 14px;
      margin: 22px 0;
      padding: 28px;
      box-shadow: 0 10px 28px rgba(30, 41, 59, 0.05);
    }
    h2 { margin: 0 0 18px; font-size: 27px; color: #173653; }
    h3 { margin: 22px 0 8px; color: #1f5f9f; font-size: 20px; }
    p { margin: 10px 0; }
    table { width: 100%; border-collapse: collapse; margin: 18px 0; font-size: 14px; }
    th, td { border: 1px solid var(--line); padding: 8px 10px; vertical-align: top; }
    th { background: #eaf1f8; color: #173653; text-align: left; }
    table.small { font-size: 12px; }
    .callout {
      margin: 16px 0;
      padding: 14px 16px;
      border-left: 5px solid var(--blue);
      border-radius: 10px;
      background: var(--info-bg);
    }
    .callout.ok { background: var(--ok-bg); border-left-color: var(--green); }
    .callout.warn { background: var(--warn-bg); border-left-color: var(--warn); }
    .callout-title { font-weight: 700; color: #173653; margin-bottom: 4px; }
    .chain {
      display: grid;
      grid-template-columns: repeat(7, minmax(0, 1fr));
      gap: 8px;
      margin: 18px 0;
    }
    .chain span {
      display: flex;
      min-height: 58px;
      align-items: center;
      justify-content: center;
      text-align: center;
      padding: 8px;
      border-radius: 999px;
      color: white;
      background: #1f5f9f;
      font-weight: 700;
      font-size: 13px;
    }
    .chain span:nth-child(n+5) { background: #02847d; }
    .formula-row {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: 12px;
      margin: 18px 0;
    }
    .formula {
      padding: 14px 16px;
      background: #f2f6fa;
      border: 1px solid var(--line);
      border-radius: 12px;
      text-align: center;
      font-size: 19px;
      font-weight: 700;
      color: #173653;
    }
    .two-col {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 18px;
    }
    .figure-card {
      margin: 22px 0;
      padding: 16px;
      border: 1px solid var(--line);
      border-radius: 14px;
      background: #fbfdff;
    }
    .fig-label { color: var(--teal); font-size: 13px; font-weight: 700; margin-bottom: 4px; }
    .figure-title { font-size: 18px; font-weight: 800; color: #173653; margin-bottom: 10px; }
    .image-row {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 12px;
      align-items: start;
    }
    .image-row img {
      width: 100%;
      background: white;
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 4px;
    }
    .note-list { margin: 12px 0 0; padding-left: 20px; }
    .note-list li { margin: 4px 0; }
    footer { max-width: 1180px; margin: 0 auto 40px; color: var(--muted); padding: 0 20px; }
    @media (max-width: 760px) {
      header { padding: 34px 20px; }
      header h1 { font-size: 28px; }
      section { padding: 20px; }
      .chain { grid-template-columns: 1fr; }
      nav { position: static; }
      nav a { display: inline-block; margin-bottom: 6px; }
    }
    """

    nav = """
    <nav>
      <a href="#proof-chain">主线</a>
      <a href="#equations">公式小抄</a>
      <a href="#exp1">实验一</a>
      <a href="#exp2">实验二</a>
      <a href="#exp3">实验三</a>
      <a href="#exp4">实验四</a>
      <a href="#exp5">实验五</a>
      <a href="#qa">答辩问答</a>
      <a href="#dont-say">禁忌表述</a>
    </nav>
    """

    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>第六章数值实验自学详解报告</title>
  <style>{css}</style>
</head>
<body>
  <header>
    <h1>第六章数值实验自学详解报告</h1>
    <p>FFT 快速求解器与无预处理 GMRES 对比：从收敛阶、边界处理、精度--成本，到 Helmholtz 小分母机制与 near-resonance 模态放大。</p>
    <p>定位：这是作者自学与答辩准备材料，不是论文正文；它解释“为什么做、图怎么看、结论是什么、不能乱说什么”。</p>
  </header>
  {nav}
  <main>{body}</main>
  <footer>Generated by reports/build_ch6_experiment_report.py. 所有图片来自 thesis/figures，所有数值来自 code/python/experiments/results。</footer>
</body>
</html>
"""


def write_build_reports(copied: list[str], missing_figs: list[str], present_csv: list[str], missing_csv: list[str], note_path: Path | None) -> None:
    status = "READY TO USE" if not missing_figs and not missing_csv else "READY WITH MISSING ASSETS"
    lines = [
        "# CH6 Experiment HTML Report Build",
        "",
        "## Outputs",
        f"- HTML: `{rel(HTML_PATH)}`",
        f"- Assets directory: `{rel(ASSET_DIR)}`",
        f"- Build script: `{rel(Path(__file__))}`",
        "",
        "## Source Note",
        f"- Learning note detected: `{note_path}`" if note_path else "- Learning note detected: not found; built from embedded refinement structure.",
        "",
        "## Figures",
        f"- Copied figures: {len(copied)}",
        "",
        "### Copied Figure Files",
    ]
    for name in copied:
        lines.append(f"  - `{name}`")
    lines.append("")
    lines.append(f"- Missing figures: {len(missing_figs)}")
    if missing_figs:
        lines.append("### Missing Figures")
        for name in missing_figs:
            lines.append(f"- `{name}`")
    lines += [
        "",
        "## CSV",
        f"- Present CSV files: {len(present_csv)}",
        f"- Missing CSV files: {len(missing_csv)}",
    ]
    for name in present_csv:
        lines.append(f"  - `{name}`")
    if missing_csv:
        lines.append("### Missing CSV")
        for name in missing_csv:
            lines.append(f"- `{name}`")
    lines += [
        "",
        "## Structure",
        "- Sections: 12",
        "- Core experiments explained: 5",
        "- exp00 is treated as implementation validation, not a core experiment.",
        "- The HTML rebuild itself did not modify experiment scripts, CSV data, or figures; thesis text updates were made separately in `thesis/chapters/6_experiments.tex`.",
        "",
        "## Final Status",
        status,
        "",
    ]
    BUILD_REPORT.write_text("\n".join(lines), encoding="utf-8")

    refine_lines = [
        "# CH6 Material Refinement Report",
        "",
        "## 1. Goal",
        "本轮根据《第六章实验理解笔记》重修第六章实验讲解材料。HTML 面向作者自学，PPT 面向答辩展示。",
        "",
        "## 2. HTML Refinement",
        "- 增加第六章证明链条、统一方程与频域分母小抄、收敛阶读图规则、error-time 读图规则。",
        "- 每个实验按“为什么做、数学对象、图怎么看、结论是什么、不要误读、答辩一句话”组织。",
        "- 增加常见答辩问题、绝对不能乱说的话和三天复习路线。",
        "- 同步新版笔记中关于图 16(b)(c) 接近平线、GMRES 计时口径和 exp04/exp05 capped 计数差异的说明。",
        "- 根据 dominant mode projection 图解读笔记，补充 full/projection/difference 的量级对比、独立色标解释和非目标模态剩余贡献说明。",
        "- 根据 true Helmholtz near-resonance summary 图解读笔记，补充图 15 四个 panel 的读图表、raw field 幅值放大和 capped/not-converged 说明。",
        "- 根据 exp07 spectral denominator heatmaps 图解读笔记，补充图 10 的 low-frequency risk map、sorted |d|、三类 case 数量级和 risk map 非条件数图说明。",
        "- 根据 exp05 resonance GMRES history 图解读笔记，补充图 18 的 residual history 读图逻辑；数量级采用当前 CSV 数据，即最终残差约 1.83e-1，而不是写入与当前数据不一致的平台数量级。",
        "- 根据 exp05 multimode resonance summary 图解读笔记，补充图 16 的三类目标集合、panel (a) 简并对系数、panel (b)(c) 接近 1 的支配性含义，以及 panel (d) 与 residual history 的分工。",
        "- 根据 exp06 accuracy-cost error-time 图解读笔记，补充图 8 的坐标轴/视觉编码、sigma=0/10 panel 含义、FFT9 四阶离散优势、GMRES30 五点无预处理基线和 timing scope 边界。",
        "",
        "## 3. Constraints Preserved",
        "- 论文正文仅同步增强关键图附近的读图解释文字，不改变实验数据、图像或理论结论。",
        "- 未新增实验。",
        "- 未重跑数据。",
        "- 未重绘实验图。",
        "- 所有图片仍来自 thesis/figures；所有数值仍来自 code/python/experiments/results。",
        "",
        "## 4. Source Note",
        f"- `{note_path}`" if note_path else "- 未自动定位到源笔记；已按用户计划和现有实验材料重修。",
        "",
        "## 5. Final Status",
        status,
        "",
    ]
    REFINEMENT_REPORT.write_text("\n".join(refine_lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    copied, missing_figs = copy_assets()
    present_csv, missing_csv = csv_status()
    note_path = find_learning_note()
    HTML_PATH.write_text(build_html(), encoding="utf-8")
    write_build_reports(copied, missing_figs, present_csv, missing_csv, note_path)
    print(f"Wrote {HTML_PATH}")
    print(f"Wrote {BUILD_REPORT}")
    print(f"Wrote {REFINEMENT_REPORT}")
    print(f"figures copied={len(copied)} missing={len(missing_figs)} csv_present={len(present_csv)} csv_missing={len(missing_csv)}")


if __name__ == "__main__":
    main()
