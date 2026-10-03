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

For CHOLMOD REML, a bordered block-chain engine is attempted before the dense
block-Schur planner. One coupled random term is partitioned into its covariance
coordinate blocks; all other effects form the border. The assembled C graph
must be a collection of paths, the local diagonal blocks must have density at
least 0.1, there must be at least three blocks, and the border must fit the
2000-effect limit. Admission uses exact stored connectivity, not model names or
thresholded numerical entries. Cycles, branching, or excessive storage use the
existing dense/sparse dispatch. Zero correlation may produce disconnected paths;
the subsequent numeric pattern is checked again if the correlation changes.

Block elimination factors each pivot block and the border Schur complement.
Backward selected inversion obtains diagonal and neighboring interior inverse
blocks, then adds the exact border correction. Nonadjacent inverse requests
within a path fall back to exact chain solves; different path components have
only the border correction. The storage estimate is
`3 b^2 + sum_g (4 m_g^2 + 5 m_g b) + 3 sum_edges m_g m_h` doubles and uses
the same configurable dense budget. `engineDiagnostics$blockChainActive`
identifies this path; `blockSchurActive` remains true for either bordered engine.
ML retains the existing dispatch. Final `computeCi=1` retains the one-time
LDLT subset extraction; `computeCi=2` uses the chain solve for the full inverse.

Before the block-chain attempt, a latent-factor Schur engine admits a single
random covariance factor with a validated Woodbury decomposition
`K = L L' + diag(psi)`. It requires positive, sufficiently well-conditioned
specifics, dense relationship precision, equal-width covariance-coordinate
blocks with disjoint incidence rows, and diagonal effective residual precision.
The original coefficient groups are independent conditional on a virtual
`rank * mainEffectSize` latent border. All other coefficients remain in the
real border. The conditional prior is
`[diag(1/psi), -diag(1/psi)L; -L'diag(1/psi), I+L'diag(1/psi)L] kron Ai / s`.
Eliminating the virtual border recovers the exact marginal C; log determinants
subtract the determinant of its virtual prior block. RHS entries for virtual
coefficients are zero and virtual solutions are not returned.

The existing marginal analytic score, AI updates and uncertainty reporting
are retained, including free loadings and specifics. Cross-environment random
score traces contract `Ai * F_g * S^-1` with `F_h` directly rather than forming
every dense inverse block. The public C and coefficient ordering remain
marginal; C is still explicitly assembled during optimization.
Admission uses `3 b^2 + mainSize^2 + sum_g (4 m_g^2 + 6 m_g b)` doubles,
where b includes virtual equations. This separate resource rule has no fixed
2000-effect latent-border cap. It does not change the limits of the ordinary
dense or chain planners. Correlated residuals, overlapping incidence blocks,
extra covariance factors, ML, or near-boundary specifics retain existing paths.
Numerical latent Cholesky failures fall back to the original marginal CHOLMOD
factorization for that evaluation. Diagnostics report `factorSchurActive` and
`factorSchurFallbacks`; `blockSchurActive` includes the latent engine.

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

The border limit remains 2000 effects. There is no fixed effect-count limit
per dense group. Each strategy must fit the dense-storage budget before its
cost is compared, using `2 b^2 + sum_g (2 m_g^2 + 3 m_g b)` doubles.
The default is 1.25e8 doubles (1000 decimal MB); set
`options(sommer.mme.denseMemoryMB=4096)` to allow a larger plan.
This estimates dense factor/cache storage, not peak process memory: sparse C,
covariance assembly, temporary matrices, and fitted output need additional memory.
Diagnostics report `denseMemoryMB` and `blockSchurEstimatedDenseMB`.
The CHOLMOD inverse cache separately requires `2 m^2 + effects * m`
for dense group storage and forward workspace to fit the same budget. For a group <= 64, or one whose
effects * group size fits the inverse-cache budget, density is not tested.
Otherwise stored within-group density must be >= 0.1.

Cost ties still choose the border strategy. The selected strategy is checked
against the density limit; rejection uses CHOLMOD rather than trying the other
strategy. Numerical block-Cholesky failure remains an error.
For coupled Kronecker terms, m is the full coefficient-group dimension, not
the main-effect dimension: the current dense engine factors that full matrix.
Main-effect-sized admission requires an actual banded or latent-factor
representation. The block-chain engine supplies this for path-structured
precisions; a Kronecker covariance alone does not split C into independent blocks.

Heuristics are component- and capability-specific. Random-effect coefficient
planning uses the precision coupling graph, group density, factor/cache storage,
and border cost. AR1 has banded precision even though its covariance is dense;
FA/RR generally have dense precision unless represented with latent factors;
the automatic latent Schur path supplies that conditional representation.
Neither should be classified using covariance density alone. Residual planning
uses observation layout, section ownership, weights, missing cells, and grid or
repeated-block workspaces through `ResidualPolicy`, independently of the dense
coefficient budget. The block-chain engine budgets its actual main-effect-sized
blocks; the latent-factor engine likewise uses its own estimate rather than
reuse the coupled-dense estimate.

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