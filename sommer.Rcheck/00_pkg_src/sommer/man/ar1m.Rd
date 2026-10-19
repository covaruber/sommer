\name{ar1m}
\alias{ar1m}
\title{First-order autoregressive covariance structure}
\description{
First-order autoregressive covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ar1m(x, rho = 0.30, fixed = FALSE,
     variance = c("homogeneous", "heterogeneous"), values = NULL)
}
\arguments{
  \item{x}{Ordered factor or design defining q ordered covariance levels.}
  \item{rho}{Starting AR(1) correlation in (-1,1).}
  \item{fixed}{For homogeneous variance, a logical indicating whether rho is fixed. For heterogeneous variance, a logical vector of length q: rho followed by q-1 variance-ratio parameters.}
  \item{variance}{Whether the AR correlation has homogeneous or heterogeneous marginal variances.}
  \item{values}{Optional q positive starting variances when \code{variance = "heterogeneous"}.}
}
\details{
The homogeneous AR(1) shape is \deqn{K_{ij}=\rho^{|i-j|}.} Heterogeneous mode uses \deqn{K=D R(\rho)D,} with the first variance fixed as the relative-scale reference. The working correlation coordinate is \eqn{\eta=\operatorname{atanh}(\rho)}, ensuring \eqn{|\rho|<1}.

Factor levels define the lag order, so numeric-looking levels must be in increasing numeric order; \code{ar1m()} warns otherwise. Inside \code{vsm()} one correlation is shared by every level of the other Kronecker factors. For a separate correlation per trial, wrap the \code{vsm()} term in \code{\link{dsumm}}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ar1m(time), ism(id))

# Separable row-by-range residual shared across trials:
# rcov = ~ vsm(dsm(trial), ar1m(range), ar1m(row), ism(units))
# Trial-specific correlations and variances:
# rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial)
}
}
\seealso{
\code{\link{ar2m}}, \code{\link{ar3m}}, \code{\link{toeplitzm}}, \code{\link{vsm}},
\code{\link{dsumm}}.
}
