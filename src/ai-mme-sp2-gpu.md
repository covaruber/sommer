# Hardware-Responsive MME Solver

`ai_mme_sp2.cpp` contains the single Henderson solver. Shared helpers live in
`mme_sp2_support.h`. CPU remains the default; eligible dense operations can
selectively use CUDA. The former `ai_mme_sp2_gpu()` copy has been removed.

## Startup and Controls

Package attachment checks accessible NVIDIA hardware first. Only if hardware
is detected does it check CUDA libraries. Startup messages distinguish no GPU,
missing libraries, a CPU-only package build, failed device initialization, and
a ready CUDA backend. Other GPU vendors do not enable this NVIDIA CUDA backend.

```r
options(sommer.mme.gpu = TRUE)
control <- mmes_control(solver = "cholmod", gpu = TRUE, gpuMinSize = 256L)
cpuControl <- mmes_control(gpu = FALSE)
```

Explicit per-fit controls override global defaults. `mmes(..., gpu=TRUE)` also
works through the existing control merger. Fits recheck readiness and stop if
GPU use was requested but unavailable. Direct inversion and `solveOnly` remain
CPU-only. No GPU-requested fit silently falls back entirely to CPU.

## Build

On Linux, install an NVIDIA CUDA toolkit containing cuBLAS, cuSOLVER, and the
CUDA runtime, plus a compatible NVIDIA driver. Then build from clean objects:

```sh
SOMMER_CUDA_HOME=/usr/local/cuda R CMD INSTALL --preclean .
```

For nonstandard installations, set `SOMMER_CUDA_CPPFLAGS` and
`SOMMER_CUDA_LIBS` instead. Configure checks compilation and linking, not GPU
availability. Normal builds do not require CUDA. CUDA libraries are loaded
dynamically rather than mandatory package dependencies. A CUDA-enabled build
can attach and fit on CPU even without GPU hardware or CUDA libraries.
Runtime lookup uses the system loader path, `SOMMER_CUDA_HOME`, and
`/usr/local/cuda/lib64`. Set the runtime loader path for custom installations.
Changing CUDA build flags requires `--preclean` to avoid stale object files.
The default Windows build remains CPU-only.

## Accelerated Operations

| Operation in the unified solver | Library | Potential benefit |
| --- | --- | --- |
| Dense group, latent, chain, banded, and border Cholesky | cuSOLVER `Dpotrf` | Large dense blocks |
| Block and border triangular solves with multiple RHS | cuBLAS `Dtrsm` | Wide solve panels and inverse work |
| Schur updates and chain/banded pivot products | cuBLAS `Dgemm` | Large dense matrix products |
| Grid-factor and missing-cell Cholesky | cuSOLVER `Dpotrf` | Large residual factors |

All CUDA kernels use double precision. `gpuMinSize=256` keeps small operations
on CPU; GEMM requires all three dimensions to meet this threshold. Set
`gpuMinSize=1` for correctness tests, not as a performance recommendation.
The existing `solver="auto"` policy remains unchanged. Eligible dense Schur
paths under `solver="cholmod"` offer the main targets. Sparse CHOLMOD/Eigen solves,
PCG products, covariance descriptors, assembly, and remaining trace operations
still run on CPU. No sparse coefficient matrix is densified just to use CUDA.

GPU allocation failures report the requested allocation size and stop. CUDA
runtime errors also stop. Numerical non-positive-definiteness still follows the
solver's existing model-optimization safeguards.

## Experiment

Use `mmes(..., engine=mmes_control(gpu=TRUE))`. Raw callers can pass `gpu` and
`gpuMinSize` to `sommer:::ai_mme_sp2()`. Use identical solver arguments for
CPU/GPU comparisons. `gpu=FALSE` exercises CPU kernels.

The returned `gpuDiagnostics` reports `enabled`, `minimumSize`, `choleskyCalls`,
`triangularSolveCalls`, and `productCalls`. Zero call counts mean the fit did not
exercise accelerated kernels, even if CUDA was enabled.

```r
workload <- generate_mme_workload(n_ind = 500L, n_levels = 4L,
	structure = "ar1", relationship = "dense", seed = 39L)
fit <- mmes(workload$fixed, random = workload$random, rcov = workload$rcov,
	data = workload$data, computeCi = 0L, getPEV = FALSE,
	engine = mmes_control(solver = "cholmod", gpu = TRUE,
		experimental = list(blockSchur = "force")))
fit$gpuDiagnostics
```

```r
Sys.setenv(SOMMER_TEST_CUDA = "true")
testthat::test_file("tests/testthat/test-mme-gpu.R")
```

This first prototype retains host matrices and transfers operands per operation.
Handles persist for a fit, but device factors and workspaces do not. Transfers
can outweigh kernel speedups. Profile complete fits before claiming improvement;
device-resident factor caches and fused solve/update operations are future work.
GPU correctness and speed require validation on real hardware.

## CRAN Checks

Normal source builds leave CUDA kernels disabled unless CUDA build variables
are explicitly set. No CUDA toolkit, GPU, or GPU-specific check flag is needed
for the default CPU package. Keep the unified solver and loader sources in the
tarball; `.Rbuildignore` excludes only developer Markdown notes under `src`.

Build vignettes and check the resulting source tarball:

```sh
R CMD build .
R CMD check --as-cran sommer_4.4.91.tar.gz
```

Install all suggested packages, including `flexiblas` and `RhpcBLASctl`, for a
complete check. `checkbashisms` and `qpdf` support shell and PDF checks. Disabling
forced suggested packages can help diagnose an incomplete local environment,
but does not establish that a complete CRAN check passes.

Hardware-independent GPU control/readiness tests run on CPU machines. Actual
GPU numerical tests run only when `SOMMER_TEST_CUDA=true` is set explicitly on
a CUDA-enabled build with accessible hardware. Do not exclude all tests or
native GPU source files to hide unrelated check failures.