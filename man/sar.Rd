\name{sar}
\alias{sar}
\title{Simultaneous autoregressive spatial covariance structure}
\description{
Simultaneous autoregressive spatial covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
sar(x, W, rho = 0.10, fixed = FALSE)
}
\arguments{
  \item{x}{Factor or design defining q spatial levels.}
  \item{W}{Square spatial weights matrix aligned to the levels of x.  Row/column names are used when available.}
  \item{rho}{Starting SAR dependence parameter.}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
The simultaneous autoregressive covariance is constructed from \deqn{B=I-\rho W,\qquad M=B^{-1}B^{-\mathsf T},\qquad K=M/M_{11}.} For a general real weights matrix, \code{rho} is restricted to the conservative interval \eqn{(-1/r(W),1/r(W))}, where \eqn{r(W)} is the spectral radius.  A bounded-logit working coordinate enforces the interval.  The first covariance derivative is supplied analytically to the generic CovarianceFactor engine.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(sar(location, W), ism(genotype))
}
}
\seealso{
\code{\link{car}}, \code{\link{maternm}}, \code{\link{vsm}}.
}
