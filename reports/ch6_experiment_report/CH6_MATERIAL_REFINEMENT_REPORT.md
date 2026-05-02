# CH6 Material Refinement Report

## 1. Goal
本轮根据《第六章实验理解笔记》重修第六章实验讲解材料。HTML 面向作者自学，PPT 面向答辩展示。

## 2. HTML Refinement
- 增加第六章证明链条、统一方程与频域分母小抄、收敛阶读图规则、error-time 读图规则。
- 每个实验按“为什么做、数学对象、图怎么看、结论是什么、不要误读、答辩一句话”组织。
- 增加常见答辩问题、绝对不能乱说的话和三天复习路线。
- 同步新版笔记中关于图 16(b)(c) 接近平线、GMRES 计时口径和 exp04/exp05 capped 计数差异的说明。
- 根据 dominant mode projection 图解读笔记，补充 full/projection/difference 的量级对比、独立色标解释和非目标模态剩余贡献说明。
- 根据 true Helmholtz near-resonance summary 图解读笔记，补充图 15 四个 panel 的读图表、raw field 幅值放大和 capped/not-converged 说明。
- 根据 exp07 spectral denominator heatmaps 图解读笔记，补充图 10 的 low-frequency risk map、sorted |d|、三类 case 数量级和 risk map 非条件数图说明。
- 根据 exp05 resonance GMRES history 图解读笔记，补充图 18 的 residual history 读图逻辑；数量级采用当前 CSV 数据，即最终残差约 1.83e-1，而不是写入与当前数据不一致的平台数量级。
- 根据 exp05 multimode resonance summary 图解读笔记，补充图 16 的三类目标集合、panel (a) 简并对系数、panel (b)(c) 接近 1 的支配性含义，以及 panel (d) 与 residual history 的分工。
- 根据 exp06 accuracy-cost error-time 图解读笔记，补充图 8 的坐标轴/视觉编码、sigma=0/10 panel 含义、FFT9 四阶离散优势、GMRES30 五点无预处理基线和 timing scope 边界。

## 3. Constraints Preserved
- 论文正文仅同步增强关键图附近的读图解释文字，不改变实验数据、图像或理论结论。
- 未新增实验。
- 未重跑数据。
- 未重绘实验图。
- 所有图片仍来自 thesis/figures；所有数值仍来自 code/python/experiments/results。

## 4. Source Note
- `C:\Users\20564\Desktop\Graduate\论文收集\第六章实验理解笔记.md`

## 5. Final Status
READY TO USE
