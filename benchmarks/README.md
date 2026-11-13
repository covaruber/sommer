# Performance regression suites

The general MME runner measures solver/backend fits.

Install the revision being measured into an explicit library, then run:

```sh
R CMD INSTALL --preclean --library=/path/to/library .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  R_LIBS=/path/to/library Rscript benchmarks/run.R /path/to/new-results
```

The library must already exist. Choose a new output directory for each run.
Results are not committed.

The four seeded workloads cover serial covariance, a dense relationship,
incomplete observations, and factor-analytic covariance. Each compares LDLT
and CHOLMOD with both sensitivity and augmented-working AI assembly using
eight optimization iterations, one warmup, and three timed
repetitions and the benchmark API's default relative parity tolerance of `1e-5`.
This is a bounded performance/parity suite, **not a convergence
certificate or a comprehensive solver validation**. Convergence status and
errors remain in the output.

Outputs:

- `results.csv`: timing, convergence, parity, selected paths, and phase counters;
- `ai-parity.csv`: frozen-parameter AI, score and coefficient parity for both
  assemblies/backends, with median information-evaluation times (including
  preparation and factorization, not isolated warm-factor kernel timings);
- `benchmark.rds`: controls, workloads, and complete benchmark rows;
- `environment.rds`: installed package/library, dependencies, BLAS, threads,
  hardware/OS, and the current checkout's revision and dirty-file list.

The checkout revision does not prove which source produced an installed package.
Archive the source tarball or full dirty-tree snapshot alongside important runs.
Compare performance only on the same machine, BLAS, thread allocation, and
dependency versions. Investigate repeatable regressions rather than comparing
single timings across unrelated CI workers.

Frozen-parameter parity requires relative differences no greater than `1e-8`;
the complete-fit comparisons retain their `1e-5` threshold. A faster kernel
does not imply a faster fit if the number of likelihood/recovery evaluations
changes. Inspect `engineDiagnostics$aiAssembly` for absorption work and workspace,
and compare records and coefficient dimensions independently in larger opt-in
workloads. The default assembly remains sensitivity, regardless of benchmark
results on a particular model.

With `diagnostics=list(profile="force")`, the native profile also reports
`residualBlockDesignCacheBlocks` and `residualBlockDesignCacheVisits`. Residual
block design caches are built only for repeated or grid residual block assembly,
not for diagonal residuals, non-diagonal weights, or saddle-point equations.
Eligible caches are populated in one sparse-design traversal, preserving local
row coordinates and diagonal weights without rescanning shared design columns
for every residual block. The visit count is the number of sparse design
entries traversed; cache allocation and column lookup work are not included.
