\name{generate_mme_workload}
\alias{generate_mme_workload}
\title{Generate Parameterized MME Workloads for Solver Evaluation}
\description{
Instantiates reproducible synthetic mixed model datasets and formulas spanning
factor dimensions, covariance structures, relationship densities, and tasks.
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
}
\value{
A list of class \code{"mme_workload"} containing data, formulas, REML and
computeCi settings, and workload metadata.
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
}
