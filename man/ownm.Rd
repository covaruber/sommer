\name{ownm}
\alias{ownm}
\title{User-defined covariance structure}
\description{
User-defined covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ownm(x, K = NULL, fun = NULL, par = numeric(), fixed = NULL,
  dfun = NULL, par_names = NULL, native_report = NULL)
}
\arguments{
  \item{x}{Variable or design defining q covariance levels.}
  \item{K}{Optional known positive-definite q by q covariance matrix.  If supplied, the structure is fixed.}
  \item{fun}{Optional function fun(par) returning a finite q by q covariance matrix.}
  \item{par}{Starting parameter vector for fun.}
  \item{fixed}{Logical vector of length par indicating fixed parameters.}
  \item{dfun}{Optional derivative function dfun(par, k) returning the derivative of the raw covariance matrix with respect to parameter k.}
  \item{par_names}{Optional names for the parameters.}
  \item{native_report}{Optional function with arguments \code{scale},
\code{par}, \code{factor}, and \code{absorb_scale}. It must return a named
numeric vector giving the fitted parameters in the model's preferred native
scale. This callback is used by \code{covparams_mmes()}.}
}
\details{
There are two modes.  With \code{K}, \code{ownm} creates a fixed covariance factor.  With \code{fun}, the matrix is evaluated at every optimizer point.  The returned matrix is symmetrized and normalized by its first diagonal element so that \code{vsm} retains the unique overall variance scale.  If \code{dfun} is supplied, its derivative is normalized analytically using the quotient rule; otherwise a central numerical derivative of the normalized factor is used.  This is the general extension mechanism for covariance structures that do not require a native C++ evaluator. When \code{native_report} is omitted, the overall variance and naturally transformed \code{par} values are reported.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
K <- matrix(c(1,.3,.3,1),2,2)
vsm(ownm(group, K=K), ism(id))

cfun <- function(p) {
  r <- tanh(p[1]); matrix(c(1,r,r,1),2,2)
}
native <- function(scale, par, factor, absorb_scale=TRUE) {
  c(variance=scale, correlation=tanh(par[1]))
}
vsm(ownm(group, fun=cfun, par=0, native_report=native), ism(id))
}
}
\seealso{
\code{\link{vsm}}, \code{\link{covparams_mmes}}, \code{\link{rrcm}}, \code{\link{maternm}}, \code{\link{sar}}, \code{\link{car}}.
}
