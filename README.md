# sommer: Solving Mixed Model Equations in R

Structural multivariate-univariate linear mixed model solver for estimation of multiple random effects with unknown variance-covariance structures and known covariance among levels of random effects (e.g., pedigree and genomic relationship matrices) ([Covarrubias-Pazaran, 2016](https://doi.org/10.1371/journal.pone.0156744); [Maier et al., 2015](https://doi.org/10.1016/j.ajhg.2014.12.006); Jensen et al., 1997). REML estimates can be obtained using the Direct-Inversion Newton-Raphson and Direct-Inversion Average Information algorithms for the problems r x r (r being the number of records) or using the Henderson-based average information algorithm for the problem c x c (c being the number of coefficients to estimate). Spatial models can also be fitted using the two-dimensional spline functionality available.

## Installation

You can install the development version of `sommer` from GitHub:

``` r
remotes::install_github('covaruber/sommer')
```

### Optional GPU Kernels

The Henderson solver uses one implementation for CPU and optional CUDA dense
kernels. CPU is the default. Package startup checks NVIDIA hardware first, then
CUDA readiness, and reports whether GPU use can be enabled.

On Linux, install a compatible NVIDIA driver and CUDA toolkit, then build from
source:

```sh
SOMMER_CUDA_HOME=/path/to/cuda R CMD INSTALL --preclean .
```

Enable eligible kernels globally with `options(sommer.mme.gpu = TRUE)`, or pass
`engine = mmes_control(gpu = TRUE, gpuMinSize = 256L)` to `mmes()`. Set
`gpu = FALSE` in the control object to override the global option. The threshold
also defaults from `options(sommer.mme.gpuMinSize = 256L)`.

CUDA runtime, cuBLAS, and cuSOLVER libraries load dynamically. CPU fitting does
not require them, even in a CUDA-enabled build. Explicit GPU requests stop when
unavailable. Only eligible dense kernels are accelerated; sparse work remains
on CPU, and transfers can outweigh gains. Direct inversion and `solveOnly`
remain CPU-only. See `?mmes` and `?mmes_control` for requirements and the
`gpuDiagnostics` counters.

## Direct-inversion fitting

Use `mmes(..., henderson = FALSE)` for the dense observation-level covariance
engine. AI is the default, with on-demand EM/GEM recovery after unreliable
information solves, exhausted line searches, or repeated poor progress.
The `blend` and `blendWeight` controls have the same meaning as in Henderson
fitting: `blendWeight=0` requests AI only, `blendWeight=1` requests EM/GEM,
intermediate weights mix directions, and `blend="pxem"` selects scalar PX-GEM.
Scalar/vector schedules and whole-number exponential tapers are supported.
Random-effect variances can reach the original-response-scale numerical
floor, `max(1e-8, tolParInv)`, without freezing the remaining components.
The assembled observation covariance must remain positive definite.

Convergence requires a small feasible score as well as a small likelihood
change or AI proposal. A rejected update is not convergence: a stalled line
search warns and returns the last accepted fit with `convergence = FALSE`.
Reaching `nIters` also returns the last evaluated point, not an unevaluated
proposal. Inspect `fit$engineDiagnostics$direct` for the stop reason, score
certificate, and line-search halving count.

Direct EM uses exact observation-space posterior moments, not a diagonal AI
approximation or a Henderson solve. All direct covariance families support
EM/GEM and PX-GEM under ML and REML, including structured residuals, custom
covariance factors, and zero-specific reduced-rank terms. Eligible scalar,
diagonal, and unstructured random updates have closed-form M-steps; other
updates numerically increase the complete-data Q objective. ML profiles fixed
effects inside the residual Q objective; REML includes fixed-effect posterior
uncertainty. Final covariance standard errors still use clean AI.
Inspect `engineDiagnostics$blend` for Q gains and M-step work, and
`engineDiagnostics$recovery` for automatic recovery triggers and acceptance.
Exact EM workspace obeys `denseMemoryMB`; structured/PX updates can cost more
than AI despite needing fewer iterations. Henderson algorithms are unchanged.

## Home website

Want to know more about the sommer project or you need professional services? Visit https://covaruber.github.io/sommer/ 

## Using `sommer` code and the GPL license

Please be aware that the `sommer` project is open-source software released under the GNU General Public License (GPL). The GPL gives you the freedom to use, study, modify, and redistribute the software, subject to the conditions of the license.

In particular, if you copy, modify, or create a work based on GPL-covered `sommer` source code and **distribute that work**, the GPL imposes requirements on how the resulting work may be licensed and distributed. Among other things, recipients must retain the freedoms provided by the GPL, and when GPL-covered software is distributed in binary or other non-source form, the corresponding source code must also be made available in accordance with the GPL.

Therefore, the fact that `sommer` is free and open-source does **not** mean that its source code can be copied into another software package, modified, redistributed as a proprietary product, and kept secret. If your project contains or is based on GPL-covered `sommer` code, please carefully review the GPL requirements before distributing it.

The purpose of the `sommer` project is to encourage open collaboration, reproducibility, and continued improvement of statistical and quantitative-genetics software. Contributions and derivative work are welcome, but they should respect the freedoms and obligations established by the GPL.

Please read the full GPL license distributed with `sommer` before incorporating `sommer` source code into another project.
