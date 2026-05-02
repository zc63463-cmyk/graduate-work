# Claim Downgrade Suggestions

以下建议只用于降低证据风险，不要求新增实验，不恢复旧图。

## `thesis/chapters/6_experiments.tex`

### Exp04 GMRES 图注

当前位置：`thesis/chapters/6_experiments.tex:217-218`

当前含义：Gaussian RHS 下无预处理 GMRES 迭代次数对比。

建议改为：

```latex
\caption{Gaussian RHS 下无预处理 GMRES 迭代次数对比；
达到迭代上限的点表示未在给定容差内收敛的 capped marker。}
```

### Exp04 结论边界

当前位置：`thesis/chapters/6_experiments.tex:243-247`

建议保留当前“不构成完整迭代法性能研究”的方向，可进一步收紧为：

```latex
该实验只在 Gaussian RHS 和所测试参数下观察无预处理 GMRES
对谱小分母的敏感性；结果不构成完整 Krylov 方法比较，
也不评价任何预处理策略的优劣。
```

### Exp05 预处理结论

当前位置：`thesis/chapters/6_experiments.tex:287-288`

建议把“说明 true Helmholtz 问题通常需要 shifted-Laplacian 等预处理技术”降级为：

```latex
这些结果提示，在 near-resonance 参数下，无预处理 GMRES
可能出现停滞或达到迭代上限；实际大规模 true Helmholtz 求解
通常需要考虑预处理技术，但本文未对具体预处理器进行实验比较。
```

### 实验总结中“本质差异”

当前位置：`thesis/chapters/6_experiments.tex:298`

建议可选降级为：

```latex
实验五和实验六分别从远离共振的谱指标和 near-resonance
扫描两方面说明 modified/true Helmholtz 在频域分母结构及
无预处理 GMRES 观察结果上的差异。
```

## `thesis/chapters/7_conclusion.tex`

### FFT 直接法“高效求解”

当前位置：`thesis/chapters/7_conclusion.tex:58-60`

建议改为：

```latex
在规则矩形区域和均匀网格上，FFT 直接法具有快速直接求解
对应离散系统的算法结构。第 6 章实验验证了离散一致性与
精度表现；本文不基于无统计重复计时得出广泛性能结论。
```

### GMRES 作为替代品

当前位置：`thesis/chapters/7_conclusion.tex:62-65`

建议改为：

```latex
无预处理 GMRES 对离散矩阵的谱性质较为敏感。
在本文规则区域测试中，它作为 Krylov 基线用于观察参数敏感性，
而不是作为 FFT 直接法的完整性能替代方案。
```

### Modified Helmholtz GMRES 下降趋势

当前位置：`thesis/chapters/7_conclusion.tex:69-72`

建议保持“本文 Gaussian RHS 扫描中”这一限定，也可写为：

```latex
在本文 Gaussian RHS 和采样参数下，modified Helmholtz 的
频域分母远离零，谱指标改善，并观察到无预处理 GMRES
迭代次数下降的趋势。
```

### True Helmholtz 与预处理器

当前位置：`thesis/chapters/7_conclusion.tex:73-78`

建议改为：

```latex
True Helmholtz（$\sigma=-\kappa^2$）：当 $\kappa^2$
接近离散 Laplacian 的某些特征值时，频域除数接近零，
近共振扫描实验显示解范数放大，且无预处理 GMRES
在给定迭代上限内未达到容差。实际应用通常需要考虑
MINRES/GMRES 与 shifted-Laplacian 等预处理思想结合，
但本文未对预处理器进行实验评价。
```

### Neumann 方法形容词

当前位置：`thesis/chapters/7_conclusion.tex:83-85`

建议改为：

```latex
Ghost-point 对称化 DCT-I 方法在本文测试的五点 Neumann
和 mixed 边界情形中支持统一的 FFT 框架求解，并呈现预期
二阶收敛。
```

## `README.md`

当前 README 对六个实验、六个 CSV、七张核心 PNG、FFT9 Dirichlet-only 范围、FACR-like 复杂度和旧图不恢复的说明一致。无必须降级项。

## `FINAL_POLISH_REPORT.md`

当前 FINAL_POLISH_REPORT 对第 6 章收束、旧实验不保留、GMRES 只作无预处理基线、FACR-like 复杂度边界的说明一致。无必须降级项。
