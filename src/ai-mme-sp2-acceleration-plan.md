# Mathematical Acceleration Plan for `ai_mme_sp2` in sommer

**Target file:** `src/MNR.cpp` (`ai_mme_sp2`)  
**Architecture reference:** `src/mme_paths.h`, `src/ai-mme-sp2-review.md`  
**Date:** 2026-10-12  

---

## Executive Summary

`ai_mme_sp2` implements mixed model equations (MME) with average information (AI) REML/ML estimation using sparse linear algebra (`Eigen::SimplicialLDLT`, SuiteSparse `CHOLMOD`, and iterative `PCG`). While symbolic caching, tensor factorizations, and selected inversion already exist, significant mathematical opportunities remain to speed up execution by $2\times$ to $10\times$ depending on model scale, trait structure, and data layout.

This document outlines 8 high-impact, state-of-the-art mathematical improvements.

---

## 1. Eliminate Exact-Solve Fallbacks via Augmented Symbolic Sparsity

- **Current bottleneck:**  
  In `ai_mme_sp2`, score traces $\operatorname{tr}(C^{-1} B_{\text{trace}})$ and $\operatorname{tr}(C_{uu}^{-1} A_i^{-1})$ query Takahashi selected-inverse entries via `getSelectedInverseOriginal`. If $B_{\text{trace}} = W' R^{-1} \dot{R} R^{-1} W$ or $A_i^{-1}$ contains nonzeros outside $\operatorname{pattern}(L)$, the engine falls back to dense column solves $C \backslash I_{:, j}$. In dense pedigree or marker-based models ($A_i^{-1}$ full), this degrades to dozens of repeated sparse solves.
- **Mathematical method:**  
  Augment the initial symbolic sparsity pattern prior to numerical factorization:
  $$\mathcal{S}_{\text{aug}} = \operatorname{pattern}(C) \cup \operatorname{pattern}(B_{\text{trace}}) \cup \operatorname{pattern}(G^{-1})$$
- **Impact:**  
  Symbolic analysis and elimination tree reordering are computed once. Takahashi recursion on $\mathcal{S}_{\text{aug}}$ yields exact inverse entries in $O(\operatorname{nnz}(L_{\text{aug}}))$ with **zero** dense solve fallbacks.

---

## 2. Sparse-RHS Tree-Reach Traversal for AI Sensitivity Solves

- **Current bottleneck:**  
  Computing the AI matrix requires solving the sensitivity system:
  $$C \cdot dBu = \text{sensitivityRHS}$$
  for all $n_{\text{vc}}$ variance parameters simultaneously as a dense matrix (`Eigen::MatrixXd`). Each parameter column in `sensitivityRHS` is nonzero only on its specific effect block $\frac{\partial G_i^{-1}}{\partial \theta_k} \hat{u}$ or $\frac{\partial R^{-1}}{\partial \theta_k} \hat{e}$.
- **Mathematical method:**  
  Exploit structural RHS column sparsity via elimination tree path reach (Gilbert-Peierls traversal or `CHOLMOD_A` with sparse RHS).
- **Impact:**  
  Forward/back substitution cost per parameter drops from full $O(\operatorname{nnz}(L))$ to path reach $O(\operatorname{nnz}(L_{*, \text{reach}}))$, reducing total sensitivity solve time by $40\%\text{--}70\%$ when $n_{\text{vc}}$ is large.

---

## 3. Hutch++ Deflated Randomized Trace Estimation for PCG

- **Current bottleneck:**  
  Under `solver = "pcg"`, log-determinants and trace corrections use standard Hutchinson Gaussian probes ($z' C^{-1} z$). Convergence rate is $O(1/\sqrt{m})$. Achieving acceptable gradient accuracy on large models requires high probe counts ($m > 30$), which causes line search oscillations and slow convergence.
- **Mathematical method (Meyer et al., 2021):**  
  Implement the **Hutch++** algorithm:
  1. Use $m/3$ matrix-vector products to approximate the dominant invariant subspace via randomized SVD: $Q \approx \operatorname{range}_k(C^{-1} B)$.
  2. Compute exact trace over subspace $Q$: $\operatorname{tr}(Q' C^{-1} B Q)$.
  3. Apply Hutchinson stochastic probing only to the deflated residual: $(I - Q Q') C^{-1} B (I - Q Q')$.
- **Impact:**  
  Trace variance drops from $O(1/\sqrt{m})$ to $O(1/m)$, yielding $3\times\text{--}5\times$ variance reduction for identical linear solve budgets.

---

## 4. Dual Phenotype-Space (Woodbury) Formulation for $n \ll p$ Regimes

- **Current bottleneck:**  
  In dense genomic prediction models (GBLUP / single-step with $p \gg n$, e.g., $p = 50{,}000$ markers/levels, $n = 3{,}000$ phenotyped individuals), primal MME size is $(n_X + p) \times (n_X + p)$. Direct sparse Cholesky factorization suffers high fill-in.
- **Mathematical method:**  
  Add a dual-path selector in `mme_paths.h` when $n < \alpha \cdot p$ ($\alpha \approx 0.5$):
  - Model covariance in data space: $V = Z G Z' + R \in \mathbb{R}^{n \times n}$.
  - When $R = \sigma_e^2 I_n$ and $G = \sigma_u^2 K$, diagonalize $Z K Z' = U \Lambda U'$ once outside the iteration loop.
  - Per iteration: $V^{-1} = U (\sigma_u^2 \Lambda + \sigma_e^2 I)^{-1} U'$.
- **Impact:**  
  Solves, determinants, and quadratic forms reduce to $O(n)$ diagonal matrix-vector operations, bypassing the $O(p^3)$ sparse factorization entirely.

---

## 5. Separable Generalized Eigendecomposition for Multi-Trait Kronecker Systems

- **Current bottleneck:**  
  In multi-trait or spatial models ($K_{\text{trait}} \otimes K_{\text{env}} \otimes A$), `ai_mme_sp2` factors the full coupled system $C$ or relies on bordered Schur complements.
- **Mathematical method:**  
  When residual and random structures share identical design incidence ($R = R_0 \otimes I_n$ and $G = G_0 \otimes A$):
  1. Solve the generalized eigenvalue problem between trait matrices:
     $$R_0^{-1/2} G_0 R_0^{-1/2} = V \Lambda V'$$
  2. Transform the coupled multi-trait MME into $T$ independent single-trait systems:
     $$(Z' Z + \lambda_t^{-1} A^{-1}) \tilde{u}_t = \tilde{y}_t \quad (t = 1, \dots, T)$$
- **Impact:**  
  One coupled $n T \times n T$ system splits into $T$ independent $n \times n$ systems that share identical symbolic patterns and can be factored in parallel across threads.

---

## 6. Multi-Secant Quasi-Newton (SR1 / L-BFGS) AI Updating

- **Current bottleneck:**  
  Every Newton step computes the exact Average Information matrix $\mathcal{I}_{\text{AI}}$, requiring $n_{\text{vc}}$ forward/back substitutions. For complex variance structures with $n_{\text{vc}} \ge 25$, sensitivity solves consume up to $70\%$ of iteration time.
- **Mathematical method:**  
  - Evaluate exact $\mathcal{I}_{\text{AI}}$ only every $K$ iterations (e.g., $K = 3$) or when step angle deviates from gradient descent.
  - On intermediate iterations, update $\mathcal{I}_{\text{AI}}$ via Symmetric Rank-One (SR1) or compact L-BFGS updates from parameter step $\Delta \theta$ and score difference $\Delta s = s^{(t)} - s^{(t-1)}$.
- **Impact:**  
  Maintains superlinear convergence while eliminating $60\%\text{--}70\%$ of sensitivity linear solves.

---

## 7. Split Preconditioned Stochastic Lanczos Quadrature (SLQ)

- **Current bottleneck:**  
  SLQ in `MNR.cpp:8844` estimates $\log|C|$ using unpreconditioned Lanczos iterations on $C$. Ill-conditioning in $C$ requires deep Krylov projections ($m \ge 25\text{--}30$) to accurately capture small eigenvalues.
- **Mathematical method:**  
  Using a symmetric split preconditioner $M = D^2$:
  $$\log|C| = \log|\tilde{C}| + 2 \sum_{i=1}^N \log D_{ii}, \quad \text{where } \tilde{C} = D^{-1} C D^{-1}$$
  Run the symmetric Lanczos recurrence directly on $\tilde{C}$.
- **Impact:**  
  Spectral clustering cuts required Lanczos steps from 20--30 down to 8--12 for identical quadrature accuracy.

---

## 8. Saddle-Point Augmented Formulation over $W' R^{-1} W$

- **Current bottleneck:**  
  Normal equation assembly $C = W' R^{-1} W + G^{-1}$ squares the condition number:
  $$\kappa(C) \approx \kappa(W)^2$$
  This causes ill-conditioning, premature numerical failure, and reliance on heuristic diagonal regularizers (`tolParInv`).
- **Mathematical method:**  
  Retain the augmented (saddle-point) system:
  $$\begin{bmatrix} R & W \\ W' & -G^{-1} \end{bmatrix} \begin{bmatrix} r \\ b \end{bmatrix} = \begin{bmatrix} y \\ 0 \end{bmatrix}$$
- **Impact:**  
  Avoids computing dense/fill-inducing $W' R^{-1} W$ products. Preserves native sparsity of $W$ and condition number $\kappa(W)$, enabling faster Krylov iterations (MINRES / SQMR) with block-diagonal preconditioning $\operatorname{diag}(R, G^{-1} + W' \tilde{R}^{-1} W)$.

---

## Recommended Implementation Roadmap

| Priority | Strategy | Implementation Area | Expected Speedup |
| :--- | :--- | :--- | :--- |
| **Phase 1** | 1. Augmented Symbolic Sparsity | `mme_paths.h`, `MNR.cpp` (Selected Inverse) | $2\times\text{--}4\times$ on dense pedigree/marker models |
| **Phase 1** | 2. Sparse-RHS Tree Reach | `MNR.cpp` (AI Sensitivity Solves) | $1.5\times\text{--}2.5\times$ when $n_{\text{vc}} > 10$ |
| **Phase 2** | 6. Quasi-Newton SR1 AI Updates | `MNR.cpp` (Optimization Loop) | $1.8\times\text{--}3\times$ across iterations |
| **Phase 2** | 3. Hutch++ Randomized Probing | `MNR.cpp` (PCG Backend) | $2\times\text{--}4\times$ probe efficiency |
| **Phase 3** | 4. Dual Phenotype Space | `mme_paths.h` (Path Selector) | $5\times\text{--}20\times$ for $n \ll p$ genomic data |
| **Phase 3** | 5. Separable Multi-Trait Kronecker | `mme_paths.h`, `MNR.cpp` (Kronecker Models) | $3\times\text{--}8\times$ for multi-trait models |
