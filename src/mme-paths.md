# MME Numerical Paths

`ai_mme_sp2()` selects algorithms from mathematical capabilities, not covariance
model names. `mme_paths.h` owns the pure planning rules. Numerical kernels,
factorization state, caches, and optimization safeguards remain in `MNR.cpp`.

## Decision Layers

1. **Factor capabilities:** descriptor evaluation and precision construction.
2. **Residual layout:** explicit or compact storage; complete or ragged blocks;
   section ownership and unique coordinates.
3. **Residual engine:** diagonal, grid, repeated Kronecker, or sparse.
4. **Data assembly:** sparse multiplication, block-local multiplication, or
   batched precision applications.
5. **Coefficient backend:** requested LDLT, CHOLMOD, or PCG; CHOLMOD may use exact
   block Schur, and eligible PCG may leave random priors matrix-free.
6. **Trace preparation:** selected inverse, inverse-block cache, or probes.
7. **Final inverse:** separate from score traces and controlled by `computeCi`.

Storage and engine selection are deliberately separate. Compact storage does
not imply grid eligibility, and a diagonal covariance does not imply diagonal
precision after a non-diagonal weight transformation.

## Residual Rules

| Property | Path |
| --- | --- |
| Sectioned or non-diagonal factors, except complete weighted blocks | Compact covariance storage |
| Compact, unweighted, ragged or sectioned, non-diagonal or sectioned | Grid candidate |
| Grid candidate with factors and acceptable rectangles in every block | Exact grid precision using missing-cell Schur complements |
| Structurally diagonal and not grid | Elementwise covariance inverse |
| Complete, non-diagonal, and not grid | Repeated Kronecker precision |
| Otherwise | Record-space sparse factorization |

Grid admission requires rectangle size <= 3 * observed + 1000 and missing count
<= 4000. Failure falls back to sparse residual handling, except sectioned
residuals retain their existing explicit unsupported-layout error. Sectioned
residuals with weights remain unsupported.

For assembly, diagonal effective precision or diagonal covariance with sparse
weight factors uses sparse products. Sparse weight factors require both stored
nonzero counts <= 4 * observations. Kronecker/grid precision with absent or
diagonal weights uses block-local designs. Other cases use batched applications.

## Solver Rules

| Request | Coefficient factorization | Accepted-state trace preparation |
| --- | --- | --- |
| LDLT | Cached Eigen symbolic/numeric LDLT | Takahashi selected inverse |
| CHOLMOD, REML | Eligible block Schur, otherwise supernodal CHOLMOD | Exact inverse blocks and selected-inverse/solve fallbacks |
| CHOLMOD, ML | Same coefficient dispatch; separate random-only D factor | D solves, without the REML C-inverse cache |
| PCG | Iterative solves and stochastic Lanczos log determinant | Hutchinson probes |

Matrix-free PCG additionally requires REML, `computeCi=0`, diagonal
preconditioning, random effects, and factor-wise precision for every random
term. A numerical precision fallback disables it for that trial. Nyström PCG
retains assembled C. PCG with ML or `computeCi=1` retains its existing error.
Even PCG fits may factor a generic sparse residual R: factorization-free refers
to the coefficient system C, not every matrix in the model.

ML uses D for determinant and score corrections, while joint sensitivity solves
and final BLUE/BLUP uncertainty continue to use C. `computeCi=0` skips final
inverse extraction, not the trace work needed during fitting.

## Block-Schur Planning

Random candidates are covariance partitions for structurally diagonal terms and
whole terms otherwise. Stored entries of C define the coupling graph. The pure
planner compares two strategies:

- Border the smaller group of each uncovered edge, with existing degree-based
  tie-breaking and deterministic edge order.
- Merge coupled groups into connected components and keep only fixed effects
  in the border.

The existing cost is retained:

    (5/3) b^3 + sum_g [(5/3) m_g^3 + 4 b m_g^2 + 2 m_g b^2]

Limits remain 2000 border effects, 8000 effects per group, and 1.25e8 doubles
for `2 b^2 + sum_g (2 m_g^2 + 3 m_g b)`. These limits are shared with the
CHOLMOD inverse-cache configuration. For a group <= 64, or one whose
effects * group size fits the inverse-cache budget, density is not tested.
Otherwise stored within-group density must be >= 0.1.

Cost ties still choose the border strategy. The selected strategy is checked
against the total memory and density limits; rejection uses CHOLMOD rather than
trying the other strategy. Numerical block-Cholesky failure remains an error.
Changing those fallback policies or tuning thresholds is a separate numerical
or performance change, not part of this behavior-preserving extraction.

## Precision And Trace Fallbacks

The descriptor precision ladder is unchanged: eligible FA/RR Woodbury,
native AR1/compound-symmetry/ANTE formulas, structurally diagonal reciprocals,
then generic SPD inversion. Eligibility includes numerical reconstruction and
conditioning checks. Random assembled/near-PD fallback and repeated-residual
block LDLT/jitter recovery remain in their owning numerical kernels.

Inverse lookup and exact solve fallback remain complementary: a missing
selected-inverse/cache entry must not be treated as zero. Partial trace sums
are discarded before exact batched solves. ML lookups use D indices, REML
lookups use C indices. Cached topology is reusable only when the stored sparse
pattern matches; numerical values must be refreshed.

## Phase Boundaries

Layout planning and constant design caches precede iteration. Each likelihood
trial evaluates covariance/precision, assembles or applies C, factors/solves,
and evaluates likelihood. Derivatives and coefficient inverse/probe preparation
follow acceptance. Final inverse extraction follows convergence. Planning must
not move accepted-state work onto rejected trials, or assume ragged block
positions are full covariance coordinates.

## Validation

`tests/native/test-mme-paths.cpp` checks dispatch truth tables, policy boundaries,
and both graph strategies; `tests/testthat/test-mme-paths.R` runs it when the
source tree and a C++ compiler are available.

For pre/post fit equivalence, run `benchmarks/mme_paths_regression.R` once with
the baseline installed package and a new reference RDS path, then with the
refactored package and the same path. It compares 130 fits across diagonal,
complete, ragged, spatial, sectioned, weighted, and dense-Gu models, including
LDLT/CHOLMOD/PCG, REML/ML, and supported final-inverse modes.