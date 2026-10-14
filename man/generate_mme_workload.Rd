\name{generate_mme_workload}
\alias{generate_mme_workload}
\alias{mme_workload_grid}
\alias{run_mme_benchmark}
\alias{analyze_mme_crossovers}
\alias{profile_mme_workload}
\title{Experimental Benchmarking and Crossover Analysis Framework}
\description{
A suite of tools to instantiate parameterized mixed model workloads across factor levels,
dimensions, covariance families, and relationship structures; benchmark candidate solvers
with execution timings and memory profiling; verify numerical parity; and fit crossover
decision boundaries to register in solver paths.
}
\usage{
generate_mme_workload(name = NULL,
                      n_ind = 60L,
                      n_levels = 5L,
                      structure = c("ar1", "ar2", "us", "diag", "fa",
                                    "rr", "spatial", "sar", "car"),
                      relationship = c("none", "dense", "sparse"),
                      missing_rate = 0.0,
                      reml = TRUE,
                      computeCi = 0L,
                      seed = 42L)

mme_workload_grid(structures = c("ar1", "us", "diag"),
                  n_ind = c(40L, 100L),
                  n_levels = c(4L, 10L),
                  relationships = c("none", "dense"),
                  reml = TRUE,
                  computeCi = 0L)

run_mme_benchmark(workloads,
                  engines = NULL,
                  replications = 3L,
                  verifyParity = TRUE,
                  tolerance = 1.0e-5,
                  nIters = 5L,
                  verbose = TRUE)

analyze_mme_crossovers(benchmark,
                       baseline_engine = NULL,
                       candidate_engine = NULL)

profile_mme_workload(workload, engine = mmes_control(), nIters = 5L,
                     dominantShare = 0.3)
}
\arguments{
  \item{name}{Descriptive label for the workload.}
  \item{n_ind}{Number of unique individuals / genotypes. Default is 60.}
  \item{n_levels}{Number of factor levels (e.g. environments, timepoints, or sites). Default is 5.}
  \item{structure}{Covariance structure type. Supported: \code{"ar1"}, \code{"ar2"},
  \code{"us"}, \code{"diag"}, \code{"fa"}, \code{"rr"}, \code{"spatial"}, \code{"sar"}, \code{"car"}.}
  \item{relationship}{Individual relationship matrix: \code{"none"} (identity),
  \code{"dense"} (genomic relationship inverse), or \code{"sparse"} (pedigree inverse).}
  \item{missing_rate}{Fraction of missing response values (0.0 to 0.5). Default is 0.0.}
  \item{reml}{Logical; whether to fit by REML (TRUE) or ML (FALSE). Default is TRUE.}
  \item{computeCi}{Final inverse demand: 0 (none), 1 (Takahashi selected inverse), 2 (full inverse).}
  \item{seed}{Random seed for reproducible dataset simulation.}
  \item{structures}{Vector of covariance structure names for grid generation.}
  \item{relationships}{Vector of relationship types for grid generation.}
  \item{workloads}{A single \code{mme_workload} object or a named list of workloads.}
  \item{engines}{A named list of \code{\link{mmes_control}} objects to evaluate.}
  \item{replications}{Number of timed repetitions per engine (warmup excluded). Default is 3.}
  \item{verifyParity}{Logical; whether to assert that experimental candidates match
  the baseline engine in likelihood, variance components, and BLUPs within \code{tolerance}.}
  \item{tolerance}{Numerical tolerance for parity verification. Default is \code{1e-5}.}
  \item{nIters}{Maximum REML/ML iterations for benchmark evaluations. Default is 5.}
  \item{verbose}{Logical; whether to print progress messages.}
  \item{benchmark}{An object of class \code{"sommer_mme_benchmark"}.}
  \item{baseline_engine}{Name of the reference baseline engine in the benchmark results.}
  \item{candidate_engine}{Name of the candidate experimental engine in the benchmark results.}
  \item{workload}{A single \code{mme_workload} to profile.}
  \item{engine}{A \code{\link{mmes_control}} object; profiling is forced on.}
  \item{dominantShare}{Share of fit time above which a phase is flagged dominant. Default is 0.3.}
}
\value{
\code{generate_mme_workload} returns an object of class \code{"mme_workload"}.

\code{mme_workload_grid} returns a list of \code{"mme_workload"} objects.

\code{run_mme_benchmark} returns an object of class \code{"sommer_mme_benchmark"} containing
execution timings, iteration counts, active solver paths, and parity status.

\code{analyze_mme_crossovers} returns an object of class \code{"sommer_mme_crossover"} summarizing
speedup ratios and empirical decision boundaries.

\code{profile_mme_workload} returns an object of class \code{"sommer_mme_profile"} with a
per-phase seconds table, workload counters, and per-iteration seconds.
}
\seealso{
\code{\link{mmes}}, \code{\link{mmes_control}}
}
\examples{
# Generate a small 1D autoregressive workload
wl <- generate_mme_workload(
  name = "demo_ar1",
  n_ind = 20,
  n_levels = 4,
  structure = "ar1"
)

# Benchmark baseline vs experimental path
engines <- list(
  baseline = mmes_control(solver = "ldlt"),
  cholmod  = mmes_control(solver = "cholmod")
)
bm <- run_mme_benchmark(wl, engines = engines, replications = 1, verbose = FALSE)
print(bm)
}
