\name{corhm}
\alias{corhm}
\title{Heterogeneous uniform-correlation covariance structure}
\description{
Heterogeneous uniform-correlation covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
corhm(x, rho = 0.10, values = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{rho}{Starting common correlation in the valid compound-symmetry interval.}
  \item{values}{Optional q positive starting variances.}
  \item{fixed}{Logical vector of length q: rho followed by q-1 variance-ratio parameters.}
}
\details{
The structure is \deqn{K=D C(\rho)D,} where \eqn{C(\rho)} has unit diagonal and common off-diagonal correlation and \eqn{D=\mathrm{diag}(1,\sqrt{r_2},\ldots,\sqrt{r_q})}.  The first variance is the scale reference; the remaining variances are estimated as positive ratios.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(corhm(environment), ism(genotype))
}
}
\seealso{
\code{\link{csm}}, \code{\link{corgm}}, \code{\link{vsm}}.
}
