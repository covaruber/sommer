# Essential performance regression suite

This directory intentionally contains only one runner. Numerical correctness,
native primitives, memory-policy boundaries, and specialized solver paths belong
in the automated tests, not in a second collection of validation scripts.

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
with CHOLMOD using eight optimization iterations, one warmup, and three timed
repetitions and the benchmark API's default relative parity tolerance of `1e-5`.
This is a bounded performance/parity suite, **not a convergence
certificate or a comprehensive solver validation**. Convergence status and
errors remain in the output.

Outputs:

- `results.csv`: timing, convergence, parity, selected paths, and phase counters;
- `benchmark.rds`: controls, workloads, and complete benchmark rows;
- `environment.rds`: installed package/library, dependencies, BLAS, threads,
  hardware/OS, and the current checkout's revision and dirty-file list.

The checkout revision does not prove which source produced an installed package.
Archive the source tarball or full dirty-tree snapshot alongside important runs.
Compare performance only on the same machine, BLAS, thread allocation, and
dependency versions. Investigate repeatable regressions rather than comparing
single timings across unrelated CI workers.
