# arXiv Literature Search Report for Thesis Topics

**Search Date:** 2026-04-29
**Focus Period:** 2020–2026
**Thesis Topics:** Fast FFT-based Poisson/Helmholtz solvers, compact 9-point 4th-order schemes, FACR/CR algorithms, GMRES for Helmholtz, Modified vs True Helmholtz, Lean 4 formal verification

---

## 1. Fast FFT-based Poisson Solvers on Rectangular Domains

### 1.1 A Scalable High-Order Multigrid-FFT Poisson Solver for Unbounded Domains on Adaptive Multiresolution Grids
- **Authors:** Gilles Poncelet, Jonathan Lambrechts, Thomas Gillis, Philippe Chatelain
- **Year:** 2025
- **arXiv ID:** 2512.08555
- **Summary:** Presents a flexible multigrid solver for the Poisson equation on collocated adaptive grids, using Fourier-based direct solvers with high-order compact stencils. Validated on up to 16,384 cores, demonstrating scalability and the advantage of compact stencils in reducing communication while improving accuracy.

### 1.2 FFT-based Free Space Poisson Solvers: Why Vico-Greengard-Ferrando Should Replace Hockney-Eastwood
- **Authors:** Junyi Zou, Eugenia Kim, Antoine J. Cerfon
- **Year:** 2021
- **arXiv ID:** 2103.08531
- **Summary:** Compares two FFT-based methods for the free-space Poisson equation, arguing that the Vico-Greengard-Ferrando algorithm achieves much higher accuracy with rapidly converging error compared to the standard Hockney-Eastwood method, while being equally efficient and simple to implement.

### 1.3 A FFT-based GMRES for Fast Solving of Poisson Equation in Concatenated Geometry
- **Authors:** Zichao Jiang, Jiacheng Lian, Zhuolin Wang
- **Year:** 2025
- **arXiv ID:** 2509.23180
- **Summary:** Extends FFT-based Poisson solvers beyond simple rectangular domains using domain decomposition. An FFT-based preconditioner accelerates GMRES for the interface system, achieving 40x speedup over sparse GMRES on cross-shaped domains with O(N log N) per-step complexity.

### 1.4 FFT-based Efficient Poisson Solver in Non-rectangular Domains
- **Authors:** (Wiley publication, not arXiv)
- **Year:** 2023
- **Source:** Comput. Anim. Virtual Worlds, DOI: 10.1002/cav.2185
- **Summary:** Proposes an FFT-based Poisson solver that extends beyond the traditional rectangular domain restriction, using a novel approach to handle complex geometries while maintaining FFT efficiency.

---

## 2. Compact Difference Schemes for Helmholtz Equations

### 2.1 Compact 9-Point Finite Difference Methods with High Accuracy Order and/or M-Matrix Property for Elliptic Cross-Interface Problems
- **Authors:** Qiwei Feng, Bin Han, Peter Minev
- **Year:** 2022 (revised 2023)
- **arXiv ID:** 2210.01290
- **Summary:** Develops 4th-order and 6th-order 9-point compact FD schemes on uniform Cartesian meshes for elliptic problems with piecewise constant coefficients across intersecting interfaces. The 6th-order special case yields M-matrix systems with proven convergence via the discrete maximum principle — directly relevant to 9-point stencil design.

### 2.2 Sixth Order Compact Finite Difference Method for 2D Helmholtz Equations with Singular Sources and Reduced Pollution Effect
- **Authors:** Qiwei Feng, Bin Han, Michelle Michelle
- **Year:** 2021 (revised 2023)
- **arXiv ID:** 2112.07154
- **Summary:** Presents a high-order compact FD method for 2D Helmholtz equations achieving 6th-order consistency for constant wavenumbers, with a pollution minimization strategy based on average truncation error of plane waves. Directly relevant to compact scheme design and pollution effect control in Helmholtz problems.

### 2.3 High Order Compact Schemes for Flux Type BCs
- **Authors:** Zhilin Li, Kejia Pan
- **Year:** 2021 (revised 2022)
- **arXiv ID:** 2109.05638
- **Summary:** Develops 4th-order compact (HOC) schemes for Robin and Neumann BCs on elliptic PDEs (Poisson, Helmholtz, diffusion-advection) in 2D/3D using undetermined coefficient methods. The resulting M-matrix property ensures the discrete maximum principle — highly relevant for compact scheme analysis on rectangular domains.

### 2.4 Fourth-Order Compact Finite Difference Schemes for Solving Biharmonic Equations with Dirichlet Boundary Conditions
- **Authors:** Kejia Pan, Jin Li, Zhilin Li, Kang Fu
- **Year:** 2024
- **arXiv ID:** 2409.01064
- **Summary:** Proposes genuine 4th-order compact FD schemes for biharmonic equations with Dirichlet BCs in 2D and 3D. Condition numbers grow as O(h⁻²), and the number of unknowns does not increase with dimensions — relevant methodology for compact scheme construction.

### 2.5 A Unified Approach to High-Order Compact Finite Difference Schemes for the Poisson and Helmholtz Equations
- **Authors:** (ScienceDirect publication)
- **Year:** 2026
- **Source:** J. Comput. Appl. Math., DOI: 10.1016/j.cam.2025.116451
- **Summary:** Presents a unified framework for deriving 4th-order compact schemes for both Poisson and Helmholtz equations, extending the Poisson scheme to general Helmholtz problems — directly aligned with the thesis topic.

---

## 3. GMRES Iterative Methods for Helmholtz Equations

### 3.1 Applying GMRES to the Helmholtz Equation with Strong Trapping
- **Authors:** Pierre Marchand, Jeffrey Galkowski, Alastair Spence, Euan A. Spence
- **Year:** 2021 (revised 2021)
- **arXiv ID:** 2102.05367
- **Summary:** First comprehensive study of frequency-dependence of GMRES iteration counts for Helmholtz boundary-integral equations under strong trapping. Proves upper bounds on how GMRES iterations grow with frequency when the problem is exponentially ill-conditioned through resonant frequencies.

### 3.2 Preconditioning of GMRES for Helmholtz Problems with Quasimodes
- **Authors:** Victorita Dolean, Pierre Marchand, Axel Modave, Timothée Raynaud
- **Year:** 2025
- **arXiv ID:** 2511.04512
- **Summary:** Derives new GMRES convergence bounds incorporating nonlinear residual behavior and relating convergence to harmonic Ritz values. Analyzes how quasimode eigenvalues hinder convergence and combines domain decomposition with deflation techniques using approximate eigenvectors for resonant regimes.

### 3.3 Iterative Solution of Helmholtz Problem with High-Order Isogeometric Analysis
- **Authors:** (ScienceDirect publication)
- **Year:** 2020
- **Source:** Comput. Methods Appl. Mech. Eng., DOI: 10.1016/j.cma.2020.112787
- **Summary:** Investigates GMRES for solving linear systems from isogeometric analysis of the Helmholtz equation, analyzing convergence behavior for high-frequency problems — relevant to GMRES performance characterization for indefinite Helmholtz systems.

---

## 4. FACR and Cyclic Reduction Algorithms

### 4.1 Theoretical Analysis of the Extended Cyclic Reduction Algorithm
- **Authors:** Xuhao Diao, Jun Hu, Suna Ma
- **Year:** 2022 (published 2023)
- **arXiv ID:** 2204.02068
- **Summary:** Fills a theoretical gap in the extended cyclic reduction (ECR) algorithm by Swarztrauber (1974), providing forward error analysis and proving that zeros of the matrix polynomial B_i^(r) are eigenvalues of a principal submatrix. Directly relevant to FACR algorithm foundations.

### 4.2 The Methods of Cyclic Reduction, Fourier Analysis and the FACR Algorithm (Classic Reference)
- **Authors:** P. N. Swarztrauber
- **Year:** 1977 (re-issued 2024)
- **Source:** SIAM Rev., DOI: 10.1137/1019071
- **Summary:** The foundational paper reviewing cyclic reduction and Fourier analysis methods along with the FACR algorithm combining both. While a classic, it remains the primary reference for FACR and was re-accessioned in 2024, indicating continued relevance.

### 4.3 Fast and Accurate Recursive Algorithms for Solving Cyclic Tridiagonal Systems
- **Authors:** (Springer publication)
- **Year:** 2025
- **Source:** J. Math. Chem., DOI: 10.1007/s10910-025-01716-x
- **Summary:** Presents recursive algorithms for cyclic tridiagonal systems using block LU factorization, extending the toolbox of methods related to cyclic reduction for structured linear systems.

> **Note:** The FACR/CR literature is dominated by classic works (Swarztrauber 1977, Buzbee-Golub-Nielson 1970, Hockney 1965). The 2022 Diao-Hu-Ma paper (2204.02068) is the most significant recent theoretical contribution, providing the first rigorous forward error analysis of the extended cyclic reduction algorithm.

---

## 5. Formal Verification of Numerical Methods (Lean / Proof Assistants)

### 5.1 Verifying Numerical Methods with Isabelle/HOL
- **Authors:** Dustin Bryant, Jonathan Julian Huerta y Munive, Simon Foster
- **Year:** 2025
- **arXiv ID:** 2511.20550
- **Summary:** Presents a framework for verifying numerical methods in Isabelle/HOL using Interaction Trees, with automated verification condition discharge and code generation from formal specs. Demonstrated on bisection and fixed-point iteration — the most directly relevant recent work on PDE solver verification methodology.

### 5.2 A Comprehensive Survey of the Lean 4 Theorem Prover
- **Authors:** Xichen Tang
- **Year:** 2025
- **arXiv ID:** 2501.18639
- **Summary:** Comprehensive survey of Lean 4 covering architecture, type system, metaprogramming, and applications in formal verification. Useful background for understanding Lean 4's capabilities for formalizing numerical methods and PDE solvers.

### 5.3 VeriNum: Formally Verified Numerical Methods (Project)
- **PI:** Andrew W. Appel (Princeton), David Bindel (Cornell), Jean-Baptiste Jeannin (UMich)
- **Year:** 2020–ongoing
- **URL:** https://verinum.org/
- **Summary:** NSF/DOE-funded project taking a layered approach to end-to-end verification of numerical software in Coq. Key results include: formally verified Jacobi iterative solver (2023), LAProof library for sparse linear algebra proofs (2023), formal proof of the Lax Equivalence Theorem for FD schemes (2021), and verified leapfrog ODE method (2022). FEM verification is in progress. **The Lax Equivalence Theorem proof and Jacobi solver verification are the closest existing work to verifying PDE solvers.**

### 5.4 Mathematical Formalized Problem Solving and Theorem Proving in Different Fields in Lean 4
- **Authors:** Xichen Tang
- **Year:** 2024
- **arXiv ID:** 2409.05977
- **Summary:** Explores LLM-assisted formalization of mathematical proofs in Lean 4, bridging natural language and computerized verification. While not PDE-specific, it demonstrates the growing feasibility of using Lean 4 for formalizing mathematical arguments relevant to numerical analysis.

---

## Gap Analysis & Thesis Relevance

| Thesis Topic | arXiv Coverage | Key Gap |
|---|---|---|
| FFT-based Poisson solvers | Well-covered (3 recent arXiv papers) | No recent paper combines FFT with 4th-order compact 9-point specifically for rectangular domains |
| Compact 9-point 4th-order schemes | Well-covered (4 papers) | Most focus on interface problems or biharmonic; limited work on standard Helmholtz with 9-point stencil |
| FACR / Cyclic Reduction | Sparse — only 1 recent theoretical paper | **Major gap**: No recent implementation or extension of FACR to 4th-order compact schemes |
| GMRES for Helmholtz | Moderate (2 key papers) | Limited work on GMRES convergence for *modified* vs *true* Helmholtz distinction |
| Lean 4 verification of numerical PDE | **Very sparse** — no direct work exists | **Major gap**: No one has formally verified FD schemes for PDEs in Lean 4; VeriNum uses Coq, Isabelle work is at early stages |
| Modified vs True Helmholtz | **No recent arXiv papers found** | **Major gap**: This distinction appears primarily in older textbooks and the thesis itself |

### Most Relevant Papers by Priority

1. **2210.01290** — 9-point compact FD with 4th/6th order (directly matches 9-point stencil thesis component)
2. **2112.07154** — 6th-order compact FD for Helmholtz with pollution control
3. **2109.05638** — HOC schemes for elliptic PDEs with flux BCs
4. **2204.02068** — Theoretical analysis of extended cyclic reduction
5. **2511.04512** — GMRES preconditioning for Helmholtz with quasimodes/resonance
6. **2511.20550** — Verifying numerical methods with Isabelle/HOL
7. **2512.08555** — Multigrid-FFT Poisson solver with high-order compact stencils
8. **2102.05367** — GMRES for Helmholtz with strong trapping
9. **2509.23180** — FFT-based GMRES for Poisson in complex geometry
10. **2501.18639** — Lean 4 comprehensive survey
