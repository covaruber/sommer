\name{covparams_mmes}
\alias{covparams_mmes}
\title{Extract native-scale covariance parameters from an mmes fit}
\description{
Returns covariance parameters in the native, model-facing scale defined by
each CovarianceFactor descriptor, rather than the unconstrained working
coordinates or normalized variance ratios used internally by the optimizer.
}
\usage{
covparams_mmes(object, term = NULL)
}
\arguments{
  \item{object}{A fitted object of class \code{"mmes"}.}
  \item{term}{Optional covariance-term name or numeric index. By default,
all random and residual covariance terms are returned.}
}
\details{
Every CovarianceFactor descriptor carries a native reporting callback. Built-in
structures use it to report meaningful quantities such as level-specific
variances for \code{dsm()}, variances and covariances for \code{usm()}, and
variance-correlation parameters for correlation structures. Structural zeros
and internal coordinates such as log standard deviations, normalized Cholesky
coefficients, and variance ratios are not returned where the covariance model
provides a more direct native representation.

All built-in CovarianceFactor families provide native reporters, including
identity, diagonal, unstructured, AR(1--3), compound-symmetry, moving-average,
general-correlation, factor-analytic, antedependence, reduced-rank, Matern,
Toeplitz, SAR, and CAR structures. Heterogeneous AR and compound-symmetry
models report their dependence parameter together with level-specific
variances rather than an overall variance plus variance ratios.

Both fixed and estimated covariance parameters are included. For an arbitrary
Kronecker product, the product-level variance scale is absorbed by the first
covariance-shaping factor. Subsequent factors are necessarily reported on their
normalized relative scale because the product has only one identifiable overall
scale.

User-defined \code{ownm()} structures can supply a \code{native_report}
callback. If none is supplied, their overall variance and naturally transformed
parameter values are returned.

Use \code{covparams_mmes_se()} when native-scale delta-method standard errors
are also required.
}
\value{
A data frame with columns \code{term}, \code{factor}, \code{parameter}, and
\code{estimate}.
}
\examples{
\dontrun{
fit <- mmes(Yield ~ Env,
            random = ~ vsm(dsm(Env), ism(Name)),
            rcov = ~ units, data = DT_example)
covparams_mmes(fit)
fit$covParNative
}
}
\seealso{
\code{\link{covparams_mmes_se}}, \code{\link{mmes}}, \code{\link{vsm}},
\code{\link{ownm}}
}
