# RUN_ALL_AUDITS

你是毕业论文终稿审查助手。请对我的项目进行一轮只读审查，不要直接修改论文或代码。

项目：
矩形区域上 Poisson/Helmholtz 方程的 FFT 快速求解器与迭代法对比研究。

统一框架：
(-Delta + sigma)u = f

其中：
sigma = 0          Poisson
sigma = +kappa^2   modified Helmholtz
sigma = -kappa^2   true Helmholtz

当前论文 V2.3 已将第 6 章收束为 6 个核心实验：
1. FFT 解与 sparse direct 解一致性；
2. 二阶与四阶收敛阶验证；
3. 非齐次 Dirichlet 边界验证；
4. Neumann 与 mixed 边界验证；
5. modified/true Helmholtz 谱指标与 Gaussian RHS 下 GMRES 行为；
6. true Helmholtz near-resonance 扫描。

请按以下顺序审查：

1. 使用 latex-formula-auditor 审查论文公式推导，输出：
   - FORMULA_AUDIT_MATRIX.md
   - FORMULA_AUDIT_REPORT.md

2. 使用 paper-code-consistency-checker 审查论文公式和 Python 实现一致性，输出：
   - PAPER_CODE_CONSISTENCY_MATRIX.md
   - PAPER_CODE_CONSISTENCY_REPORT.md

3. 使用 experiment-evidence-auditor 审查第 6 章实验是否支撑论文结论，输出：
   - EXPERIMENT_EVIDENCE_MATRIX.md
   - EXPERIMENT_AUDIT_REPORT.md
   - CLAIM_DOWNGRADE_SUGGESTIONS.md

4. 使用 lean-proof-mapper 审查 Lean 4 附录和论文命题映射，输出：
   - LEAN_PROOF_MAP.md
   - LEAN_SCOPE_AUDIT.md

5. 使用 gpt-pro-question-writer 只整理仍无法裁决的高价值问题，输出：
   - QUESTIONS_FOR_GPT55_PRO.md

所有问题按 A/B/C/D 分类：

A = 必须终稿前修；
B = 建议修，能降低风险；
C = 可答辩解释；
D = 不是问题或影响很小。

硬性要求：
1. 不要新增实验。
2. 不要恢复旧图。
3. 不要大改论文主线。
4. 不要把 Lean 4 说成完整 PDE 证明。
5. 不要把 FACR-like 说成 O(N^2 log log N)。
6. 不要把 modified Helmholtz 写成 lambda+kappa^2=0 共振。
7. 不要把 FFT9 写成完整支持 Neumann 四阶。
8. 不要把无预处理 GMRES 写成完整迭代法性能研究。
9. 只读审查，不直接修改文件。
10. 如果发现必须修复的问题，只写入 OPTIONAL_PATCH_PLAN.md。

最后输出：

# Audit Run Summary

## 1. Generated Files
列出生成的审查文件。

## 2. A-class Issues
列出所有必须修复问题。

## 3. Questions for GPT-5.5 Pro
说明 QUESTIONS_FOR_GPT55_PRO.md 中有几个问题。

## 4. Final Status
写：
READY TO SUBMIT
或
READY WITH MINOR FORMULA QUESTIONS
或
NOT READY
