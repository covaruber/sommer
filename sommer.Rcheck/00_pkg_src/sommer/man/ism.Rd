\name{ism}
\alias{ism}
\title{Identity covariance structure}
\description{
Identity covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ism(x)
}
\arguments{
  \item{x}{Factor, character vector, numeric design, or matrix defining covariance levels.}
}
\details{
\code{ism} returns \eqn{K=I}.  It has no covariance-shape parameters; the only unknown covariance scale is the \eqn{\sigma^2} supplied by \code{vsm}.  It is commonly used as the final term of \code{vsm} to indicate an independent main effect.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ism(id))
vsm(dsm(environment), ism(genotype))
}
}
\seealso{
\code{\link{vsm}}, \code{\link{mmes}}.
}
