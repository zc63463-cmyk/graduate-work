#!/usr/bin/env python3
"""
综合验证脚本 — 论文项目最终检查
检查项：
1. LaTeX 编译无错误
2. Lean 4 形式化验证通过
3. 所有 Python 求解器可运行
4. 论文页数在合理范围
"""

import subprocess
import sys
import os
from pathlib import Path

# Force UTF-8 output
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
if sys.stderr.encoding != 'utf-8':
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

PROJECT_ROOT = Path(r"c:\Users\20564\Desktop\Graduate\论文收集")
THESIS_DIR = PROJECT_ROOT / "thesis"
LEAN_DIR = PROJECT_ROOT / "code" / "lean4_formalization"
PYTHON_DIR = PROJECT_ROOT / "code" / "python"

results = {}

def check(name, passed, detail=""):
    results[name] = (passed, detail)
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] {name}" + (f": {detail}" if detail else ""))

print("=" * 60)
print("Thesis Project Final Verification")
print("=" * 60)

# --- 1. LaTeX compilation ---
print("\n[1] LaTeX compilation check...")
try:
    result = subprocess.run(
        ["latexmk", "-xelatex", "-interaction=nonstopmode", "main.tex"],
        cwd=str(THESIS_DIR), capture_output=True, text=True, timeout=300,
        encoding='utf-8', errors='replace'
    )
    errors = [l for l in result.stdout.split("\n") if "!" in l and "error" in l.lower()]
    check("LaTeX 编译", len(errors) == 0, f"{len(errors)} 个错误")
    
    # Check PDF exists
    pdf_path = THESIS_DIR / "main.pdf"
    if pdf_path.exists():
        size_kb = pdf_path.stat().st_size / 1024
        check("PDF 生成", True, f"{size_kb:.0f} KB")
    else:
        check("PDF 生成", False, "未找到 PDF")
except Exception as e:
    check("LaTeX 编译", False, str(e))

# --- 2. Lean 4 verification ---
print("\n[2] Lean 4 formal verification...")
lean_path = os.path.expanduser("~\\.elan\\bin")
env = os.environ.copy()
env["PATH"] = lean_path + ";" + env.get("PATH", "")

try:
    lake_exe = os.path.join(lean_path, "lake.exe")
    # Full build
    result = subprocess.run(
        [lake_exe, "build"],
        cwd=str(LEAN_DIR), capture_output=True, text=True, timeout=600,
        env=env, encoding='utf-8', errors='replace'
    )
    check("Lean 4 全量编译", result.returncode == 0,
          "OK" if result.returncode == 0 else "FAILED")
    
    # MathlibTest specifically
    result_mt = subprocess.run(
        [lake_exe, "build", "SixthOrderImpossibility.MathlibTest"],
        cwd=str(LEAN_DIR), capture_output=True, text=True, timeout=600,
        env=env, encoding='utf-8', errors='replace'
    )
    check("Lean 4 Mathlib R-ver", result_mt.returncode == 0,
          "OK" if result_mt.returncode == 0 else "FAILED")
except Exception as e:
    check("Lean 4 验证", False, str(e))

# --- 3. Python solver smoke test ---
print("\n[3] Python solver smoke test...")
try:
    sys.path.insert(0, str(PYTHON_DIR))

    # Test FFT9 (import only)
    import fft9_complete
    check("FFT9 模块加载", True, "import OK")

    # Test CR (import only)
    import cyclic_reduction
    check("CR 模块加载", True, "import OK")

    # Test Helmholtz solver (import only)
    import helmholtz_solver
    check("Helmholtz 模块加载", True, "import OK")

    # Test GMRES (import only)
    import gmres_solver
    check("GMRES 模块加载", True, "import OK")
except Exception as e:
    check("Python 求解器", False, str(e))

# --- 4. Thesis content checks ---
print("\n[4] Thesis content check...")
ch_files = {
    "Ch1 引言": "1_introduction.tex",
    "Ch2 数学基础": "2_math_preliminary.tex",
    "Ch3 FFT直接法": "3_fft_direct.tex",
    "Ch4 Helmholtz": "4_helmholtz.tex",
    "Ch5 GMRES": "5_gmres.tex",
    "Ch6 实验": "6_experiments.tex",
    "Ch7 结论": "7_conclusion.tex",
}
total_lines = 0
for name, fname in ch_files.items():
    fpath = THESIS_DIR / "chapters" / fname
    if fpath.exists():
        lines = len(fpath.read_text(encoding='utf-8').splitlines())
        total_lines += lines
        check(name, lines > 50, f"{lines} 行")
    else:
        check(name, False, "文件缺失")

check("论文总行数", 500 < total_lines < 3000, f"{total_lines} 行")

# --- 5. BibTeX references ---
print("\n[5] BibTeX references check...")
bib_path = THESIS_DIR / "references.bib"
if bib_path.exists():
    bib_content = bib_path.read_text(encoding='utf-8')
    num_entries = bib_content.count("@")
    check("参考文献数量", num_entries >= 20, f"{num_entries} 篇")
    
    # Check key references
    key_refs = ["Sutmann2007", "SaadSchultz1986", "HoustisPapatheodorou1979a",
                "Swarztrauber1977", "Lynch1977"]
    for ref in key_refs:
        found = ref in bib_content
        check(f"引用 {ref}", found)
else:
    check("references.bib", False, "文件缺失")

# --- Summary ---
print("\n" + "=" * 60)
passed = sum(1 for v in results.values() if v[0])
total = len(results)
print(f"验证结果: {passed}/{total} 通过")
if passed == total:
    print("All checks passed!")
else:
    print("Some checks failed. See details above.")
    for name, (ok, detail) in results.items():
        if not ok:
            print(f"  FAIL: {name}: {detail}")
