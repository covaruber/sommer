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
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ar1m(time), ism(id))
}
}
\seealso{
\code{\link{ar2m}}, \code{\link{ar3m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}
