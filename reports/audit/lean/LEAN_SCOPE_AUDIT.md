# Lean Scope Audit

## Verdict

总体评价：论文主文和附录的 Lean 范围总体保守，没有把 Lean 4 说成完整 PDE 收敛证明、程序正确性证明、DST/DCT 完整谱理论证明，也没有把贡献扩大成所有六阶格式不可行。

需要终稿前关注的点有三个：一是 Lean theorem 数量统计与当前源码不一致；二是第 2 章个别句子应从“九点紧致差分格式”收窄为“该三参数九点修正族”；三是 Helmholtz 情形的 `beta = gamma` 不宜写成已有命名 Lean 定理直接验证，除非补一个对应 theorem。

## Build And Escape Status

- `lake build`: pass.
- Warnings: unused `simp` arguments; proofwidgets dependency local changes notice.
- Forbidden proof escapes in active project source: none found for `sorry`, `admit`, `axiom`, `unsafe`.
- Current compiled theorem declaration inventory: 24 total, with 14 in `SixthOrderImpossibility.lean` and 10 in `SixthOrderImpossibility/MathlibTest.lean`.

## Acceptable Claims

| Location | Wording summary | Why acceptable |
|---|---|---|
| `appendix_lean4.tex:7-8` | Lean is used for auxiliary algebraic verification of the three-parameter nine-point correction family, not full PDE global error estimates. | Correctly limits Lean to algebraic core and restricted family. |
| `appendix_lean4.tex:72-77` | Appendix says full PDE global error estimates, Fourier-symbol convergence, and DST-I spectral theory are not covered. | Avoids major overclaim categories. |
| `appendix_lean4.tex:88-90` | Remark says Lean verifies algebraic coefficient contradiction, not full PDE global error estimate or spectral theory. | Thesis-safe and accurate. |
| `7_conclusion.tex:35-43` | Conclusion states the impossibility for a three-parameter nine-point correction family and calls Lean an auxiliary formal check. | Scope is mostly accurate; only theorem count should be updated. |
| `7_conclusion.tex:116-121` | Conclusion says Lean is only for algebraic core and does not cover PDE global error, DST/DCT spectral theory, or program correctness. | This is the clearest guardrail and should be kept. |
| `FINAL_POLISH_REPORT.md:18` | Report says the thesis no longer claims Lean proves full PDE global error. | Accurate. |

## Scope Risks

| Location | Current wording | Risk | Suggested conservative wording |
|---|---|---|---|
| `appendix_lean4.tex:57-59`; `7_conclusion.tex:42-43` | “整数域 12 定理 + 实数域 8 定理” | Current active source has 14 integer theorem declarations and 10 real theorem declarations. The count is stale or uses an unstated “core theorem” filter. | “Lean 源码包含若干整数域和实数域辅助定理；表中列出与论文命题直接相关的核心定理。” Or update counts to 14 and 10. |
| `2_math_preliminary.tex:568-569` | “四阶是九点紧致差分格式的自然精度上限。” | Can be read as all nine-point compact formats cannot reach sixth order. Lean supports the restricted three-parameter correction family only. | “对形如式~\\eqref{eq:modified_compact} 的三参数九点修正族而言，四阶是该族在统一参数约束下的自然精度上限。” |
| `2_math_preliminary.tex:533`; `2_math_preliminary.tex:588-593`; `7_conclusion.tex:39-43` | Helmholtz c2 analysis mentions `beta = gamma`, then says the conclusion is Lean-verified. | Current named theorem table directly proves `alpha = -gamma` and the final no-sixth-order obstruction; there is no named theorem explicitly stating `beta = gamma`. | “Lean 验证了 Helmholtz 情形下 c2 条件推出的 `alpha = -gamma` 约束及最终 c4 矛盾；`beta = gamma` 条件由正文代数推导给出。” |
| `appendix_lean4.tex:72-77` | 未覆盖列表 does not explicitly mention program implementation correctness. | Chapter 7 covers it, but the appendix itself could be clearer. | Add one bullet: “不验证 Python/实验程序实现正确性。” |
| `code/lean4_formalization/SixthOrderImpossibility.lean:14` | Source comment says Lean kernel `v4.20.0`. | `lean-toolchain` is `leanprover/lean4:v4.29.1`; comment is stale. | Update source comment to `v4.29.1` or remove exact version from the source comment. |
| `code/lean4_formalization/SixthOrderImpossibility.lean:154-155`; `:181-182` | Source docstrings say the theorem proves the 9-point compact format cannot achieve sixth order and is the definitive obstruction to sixth-order compact schemes. | This is broader than the safe thesis scope if quoted directly. | “This proves the c4 obstruction for the restricted three-parameter correction family encoded here.” |

## Required Guardrails

- Lean verifies algebraic coefficient identities and parameter contradictions only.
- Lean does not verify full PDE convergence, global error estimates, stability, boundary-condition analysis, DST/DCT spectral theory, or Python implementation correctness.
- The impossibility claim is restricted to the three-parameter nine-point correction family in `eq:modified_compact`.
- Do not state that Lean proves all sixth-order compact schemes are impossible.
- Do not require formalizing full PDE theory; the current appendix only needs accurate scope language.

## Final Scope Status

Status: mostly accurate with minor scope/count fixes needed.

No evidence found that the thesis currently claims Lean proves complete PDE convergence, all nine-point sixth-order impossibility, program correctness, or complete DST/DCT spectral theory. The strongest required edit is to update theorem counts and narrow the broad “九点紧致差分格式自然精度上限” sentence.
