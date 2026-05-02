# Lean Appendix Table Suggestion

This file suggests a thesis-safe replacement or supplement for the Lean appendix theorem table. It does not require formalizing full PDE theory and does not expand the Lean contribution.

## Suggested Public Appendix Table

| Lean theorem | Mathematical statement | Thesis location | Proof strategy | Scope | Is thesis wording accurate? yes/no/partial | Suggested wording |
|---|---|---|---|---|---|---|
| `c4_num_mode_10_real` | For the encoded c4 numerator, mode `(1,0)` gives `60 gamma + 3`. | `appendix_lean4.tex:45`; `2_math_preliminary.tex:577` | `ring` | Real algebraic identity. | yes | Keep. |
| `c4_num_mode_11_real` | For the encoded c4 numerator, mode `(1,1)` gives `120 gamma - 4`. | `appendix_lean4.tex:46`; `2_math_preliminary.tex:578` | `ring` | Real algebraic identity. | yes | Keep. |
| `sixth_order_impossible_real` | No real gamma makes the encoded c4 numerator vanish for all modes. | `appendix_lean4.tex:47`; `2_math_preliminary.tex:580-586` | mode substitution, `ring`, `linarith` | Algebraic c4 obstruction only. | yes | Keep but say “encoded c4 polynomial”. |
| `c2_poisson_forces_neg_gamma` | In the Poisson c2 expression, c2 zero at the tested mode forces `alpha = -gamma`. | `appendix_lean4.tex:48`; `2_math_preliminary.tex:544-545` | unfold, `linarith` | Algebraic c2 helper. | yes | Keep. |
| `poisson_no_sixth_order_real` | Under `alpha + gamma = 0`, no real gamma makes c4 vanish for all modes. | `appendix_lean4.tex:49`; `2_math_preliminary.tex:584-591` | reduce to c4 contradiction | Restricted three-parameter family after c2 elimination. | yes | Keep. |
| `c2_helmholtz_forces_neg_gamma` | Helmholtz c2 zero at two modes forces `alpha = -gamma`. | `appendix_lean4.tex:50`; `2_math_preliminary.tex:519-533` | subtract equations, `ring_nf`, `linarith` | Algebraic c2 helper; does not explicitly state `beta = gamma`. | partial | Wording should not claim this theorem alone proves `beta = gamma`. |
| `helmholtz_no_sixth_order_real` | No real `alpha,beta,gamma` make the encoded c2 and c4 conditions vanish for all modes. | `appendix_lean4.tex:51`; `2_math_preliminary.tex:568-593` | use c2 helper and c4 contradiction | Restricted three-parameter algebraic obstruction. | yes | Keep with “三参数修正族”. |
| `c2_mode_dependent_helmholtz` | From Helmholtz c2 zero at `(1,0)`, Lean derives `alpha = -gamma + (beta - gamma) k2`. | Add only if expanding the table. | unfold, `linarith` | Helper for Helmholtz parameter relation. | partial | Optional helper row. |
| `c2_mode_01_helmholtz` | Symmetric Helmholtz c2 helper at `(0,1)`. | Add only if expanding the table. | unfold, `linarith` | Helper theorem. | partial | Optional helper row. |
| `c2_helmholtz_mode_20_forces_neg_gamma` | Alternative Helmholtz c2 mode pair also forces `alpha = -gamma`. | Add only if expanding the table. | subtract equations, `ring_nf`, `linarith` | Helper theorem. | partial | Optional helper row. |

## Suggested Count Wording

Replace exact theorem-count wording unless you want to count every active declaration:

```latex
Lean 4 源码包含整数域与实数域两组辅助定理。
表~\ref{tab:lean4_theorems} 仅列出与论文命题直接对应的核心实数域定理；
整数域版本用于对同一系数矛盾进行独立的离散代数校验。
```

If exact counts are preferred, use:

```latex
当前源码中，整数域文件包含 14 个 theorem 声明，
实数域 Mathlib 文件包含 10 个 theorem 声明；
表中列出其中与论文主命题直接相关的核心定理。
```

## Suggested Scope Paragraph

```latex
Lean 4 形式化验证的是三参数九点修正族六阶不可达性证明中的
代数核心步骤：c_2 条件给出的参数约束，以及特定 Fourier 模式
下 c_4=0 条件之间的矛盾。该验证不覆盖完整 PDE 全局误差估计、
DST/DCT 谱理论证明、边界条件分析或 Python 程序实现正确性；
也不声称所有可能的六阶紧致格式均不可行。
```

## Suggested Fix For Chapter 2 Wording

For `thesis/chapters/2_math_preliminary.tex:568-569`, replace the broad sentence with:

```latex
该代数矛盾由 Lean 4 定理 \texttt{helmholtz\_no\_sixth\_order\_real}
辅助验证。因此，对形如~\eqref{eq:modified_compact} 的三参数九点修正族而言，
四阶是该族在统一常数参数约束下的自然精度上限。
```
