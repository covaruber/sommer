# Contributing to sommer

## Architecture

The public `mmes()` and `mmes_control()` interface and model-independent
compatibility rules live in [R/FUN_mme_modules.R](R/FUN_mme_modules.R).
Private defaults/validation and the shared PCG transport mapping live in
[R/FUN_mme_control.R](R/FUN_mme_control.R). Model-dependent eligibility stays
with model preparation or the engine; the capability map describes it but is
not a replacement for runtime checks.

| Responsibility | Source |
|---|---|
| Formula evaluation, filtering, preparation, dispatch | [R/mmes.R](R/mmes.R), [model helpers](R/FUN_model_helpers.R), [weight preparation](R/FUN_model_weights.R) |
| Covariance term assembly and known precision | [R/FUN_vsm.R](R/FUN_vsm.R) |
| Covariance descriptor contract | [R/FUN_covfactor.R](R/FUN_covfactor.R) |
| Covariance families | `R/FUN_covariance_{basic,serial,factor,spatial}.R` |
| Public inference and prediction | [R/FUN_methods_mmes.R](R/FUN_methods_mmes.R), [R/predict_mmes.R](R/predict_mmes.R) |
| Gaussian Henderson optimizer and dispatch | [src/ai_mme_sp2.cpp](src/ai_mme_sp2.cpp) and specialized `mme_*.h` kernels |
| Direct observation-space fitting | [src/mme_direct.cpp](src/mme_direct.cpp) |
| Shared direct/solve covariance evaluation | [src/mme_direct_covariance.h](src/mme_direct_covariance.h) |
| Known-covariance PCG solving | [src/mme_solve.cpp](src/mme_solve.cpp) |
| Post-fit inverse/contrast covariance | [src/mme_inference.cpp](src/mme_inference.cpp) |
| Matrix utilities and relationships | [src/native_utils.cpp](src/native_utils.cpp), [src/native_relationships.cpp](src/native_relationships.cpp) |
| Legacy native engines and R compatibility | [src/MNR.cpp](src/MNR.cpp), [R/DEPRECATED_mmer.R](R/DEPRECATED_mmer.R) |
| Lifecycle hooks | [R/zzz.R](R/zzz.R) |

Covariance constructors supply the existing CovarianceFactor-v2 descriptor;
do not make the optimizer depend on constructor names. Preserve parameter order,
coordinates, level names, formula environments, and accepted-point convergence
semantics. Native prediction factors are operation-local: never save external
pointers in fitted objects.

The largest Henderson routine still needs incremental decomposition. Extract
cohesive layout, assembly, evaluation, and optimization operations separately;
do not combine a mechanical move with changed statistical or dispatch logic.

## Local validation

Install the declared dependencies before running tests. From the repository root:

```sh
Rscript --vanilla -e 'testthat::test_local(filter="mme-modules|control-contracts", stop_on_failure=TRUE)'
Rscript --vanilla -e 'testthat::test_local(stop_on_failure=TRUE)'
```

`test_local()` compiles the working-tree native code and loads the working-tree
R sources. Tests against an older installed package are not sufficient.
The native wrappers normally skip when source headers or a compiler are absent.
The Rcpp-based AI probes also require `RcppArmadillo` at test time; its
`LinkingTo` declaration only guarantees availability during package compilation.
It is also declared in `Suggests` because tests check its namespace at runtime.
These probes skip when that dependency is unavailable during installed-package
checks; the model-fitting tests in the same files still run.
The required source-tree CI job turns those missing prerequisites into failures:

```sh
SOMMER_REQUIRE_NATIVE=true Rscript --vanilla -e \
  'testthat::test_local(filter="mme-control|mme-paths|mme-solver|mme-saddle|pcg-numerics|ai-recovery|ai-variance-floor|reduced-ai", stop_on_failure=TRUE)'
```

| Native source | Automated test entry |
|---|---|
| `test-mme-control.cpp`, `test-mme-control-transport.cpp` | `test-mme-control-native.R` |
| `test-mme-paths.cpp` | `test-mme-paths.R` |
| `test-mme-sparse.cpp` | `test-mme-solver.R` |
| `test-mme-saddle.cpp` | `test-mme-saddle.R` |
| `test-mme-pcg.cpp`, `test-mme-pcg-estimators.cpp`, `test-mme-sufficient.cpp` | `test-pcg-numerics.R` |
| `test-mme-ai-recovery.cpp` | `test-ai-recovery.R`, `test-ai-variance-floor.R`, `test-reduced-ai.R` |

CUDA tests require suitable hardware, a CUDA-enabled build, and
`SOMMER_TEST_CUDA=true`. A CPU-only check does not certify CUDA behavior.
Scope global options, BLAS settings, temporary files, and working directories
so they are restored on both success and failure.
Tests unrelated to thread control should leave `execution$blasThreads=NULL`:
optional controller packages do not guarantee a controllable BLAS backend.
Set process-level thread limits before launching R when needed.
Child-process tests should set and restore the parent's environment so the child
inherits it; `system2(env=...)` is not portable for arbitrary Windows executables.
Solver-parity tests requiring tight parameter agreement should explicitly request
tight fitting tolerances rather than rely on the default stopping criteria.
Avoid near-machine-precision fitting tolerances and platform-specific iteration
counts: require certified convergence within the configured iteration budget,
and check parameter agreement, score, boundary values and likelihood monotonicity.

## Release checks and documentation

Build in a separate output directory, then check the resulting source tarball:

```sh
R CMD build /path/to/sommer
R CMD check --as-cran sommer_4.5.11.tar.gz
```

Use the tarball version from DESCRIPTION. Full release checks require
`pdflatex` and `makeindex` from a TeX installation, `checkbashisms` (provided by
Debian's `devscripts`), and `qpdf` on `PATH`. HTML Tidy (`tidy`) enables the optional
HTML manual check. Missing tools are check-environment problems, not reasons to
remove the manual, configure script, or vignettes.

On clusters, load the available TeX module before starting R. To test PDF
generation without recompiling the package:

```sh
R CMD Rd2pdf --no-preview --output=sommer-manual.pdf /path/to/sommer
```

Compiler flags inherited from `R CMD config CXXFLAGS` or user/site Makevars can
produce portability notes even when the package does not set those flags.
Use a portable R toolchain or a check-only `R_MAKEVARS_USER` file to remove
architecture-specific flags; do not override the user's compiler flags in the
package Makevars.

Incoming URL, DOI, ORCID, and current-time checks need working outbound network
access. Investigate TLS failures or unavailable services from another network
before changing identifiers or reporting them as broken. Installed-size notes
are informational; keep required documentation and libraries.

Do not suppress vignette building for the release gate. A quick check with
`--no-build-vignettes` does not certify generated articles.
CI checks Linux, Windows, macOS, previous-release R, and R-devel. That matrix
does not by itself establish that every dependency combination allowed by the
historical minimum R version in DESCRIPTION works.

Release validation also needs GCC and Clang undefined-behavior sanitizer builds:
ordinary test success does not detect null-pointer arithmetic in empty matrix
views. Run the tests against each freshly instrumented installation with
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`. Keep zero-column AI tangent,
scale-only, compressed-model, and fixed-only prediction cases in that coverage.
Test BLAS-dependent documentation on the available backends. Examples that
change runtime thread settings are opt-in, since a controller's presence does
not guarantee that the BLAS accepts those settings.

The project mixes manually maintained Rd/namespace entries and roxygen comments.
Do not regenerate all documentation blindly. Update the directly related help
pages and verify code/documentation consistency with `R CMD check`.
When Rcpp exports change, regenerate bindings with `Rcpp::compileAttributes(".")`.
Keep Unix and Windows header dependencies in sync.

To validate the solver decision tree:

```sh
Rscript --vanilla pkgdown/mmes-options/build-mmes-options.R
Rscript --vanilla pkgdown/mmes-options/check-mmes-options.R
```

When compatibility rules change, update the guide's rule contexts and relevant
decision nodes together; the guide checker validates rule and node references.

## Numerical and performance changes

Use independent small-matrix references where possible. Existing tests cover
direct/Henderson and LDLT/CHOLMOD parity, fixed-only models, structured/weighted
and incomplete residuals, dense relationships, FA/RR, PQL, PCG reproducibility,
prediction, and uncertainty. Extend the relevant tests rather than starting
another standalone validation-script collection.

Compare estimates, scores, information, predictions, output shapes, convergence,
and selected paths, not just final likelihood. Use fixed seeds and appropriate
accuracy tolerances for stochastic methods. Avoid platform-wide bitwise equality
requirements.

The [essential benchmark suite](benchmarks/README.md) records bounded workloads,
timing, parity, convergence, and environment metadata. Keep generated results and
installed libraries out of Git. Compare repeated runs on the same machine and
thread configuration; investigate repeatable regressions above baseline noise.
