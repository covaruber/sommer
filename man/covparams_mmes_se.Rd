\name{covparams_mmes_se}
\alias{covparams_mmes_se}
\title{Native-scale covariance parameters and standard errors}
\description{
Returns descriptor-defined native covariance parameters with delta-method
standard errors computed from the fitted covariance-parameter uncertainty
matrix.
}
\usage{
covparams_mmes_se(object, term = NULL, rel_step = 1e-6)
}
\arguments{
  \item{object}{A fitted object of class \code{"mmes"}.}
  \item{term}{Optional covariance-term name or numeric index. By default,
all random and residual covariance terms are returned.}
  \item{rel_step}{Positive relative step used to numerically differentiate
each CovarianceFactor native reporting callback.}
}
\details{
The fitted \code{theta_se} matrix describes uncertainty in the reported
\code{covPar} coordinates. For each covariance structure,
\code{covparams_mmes_se} numerically evaluates the Jacobian of its descriptor's
native reporting callback and applies the multivariate delta method,
\deqn{V_{native} = J V_{covPar} J^{\mathsf T}.}
Consequently, coupled transformations use the complete covariance information.
For example, the standard error of a \code{dsm()} level variance accounts for
uncertainty and covariance in both the overall variance and its variance ratio.
Likewise, \code{usm()} standard errors account for all relevant normalized
Cholesky parameters.

Parameters marked fixed in the covariance descriptor contribute no independent
sampling uncertainty: their rows and columns in \code{theta_se} are excluded
before transformation. A native quantity involving both fixed and free
parameters can still have a nonzero standard error through its free parameters.
Both fixed and estimated native quantities are returned.

The same table is stored automatically in \code{object$covParNativeSE} by
\code{mmes()}.
}
\value{
A data frame containing \code{term}, \code{factor}, \code{parameter},
\code{estimate}, \code{StdError}, and \code{Zratio}. A zero standard error is
reported for a fully fixed native quantity, with \code{Zratio=NA}.
}
\examples{
\dontrun{
fit <- mmes(Yield ~ Env,
            random = ~ vsm(dsm(Env), ism(Name)),
            rcov = ~ units, data = DT_example)
covparams_mmes_se(fit)
fit$covParNativeSE
}
}
\seealso{
\code{\link{covparams_mmes}}, \code{\link{mmes}}, \code{\link{vsm}}
}
