\name{atm}
\alias{atm}
\title{Selected-level diagonal covariance structure}
\description{
Selected-level diagonal covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
atm(x, levs, values = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining all available levels.}
  \item{levs}{Character names or numeric positions of the levels retained in this covariance factor.}
  \item{values}{Optional positive starting variances for the selected levels.}
  \item{fixed}{Logical vector of length length(levs)-1 controlling the selected variance ratios.}
}
\details{
\code{atm} is a selected-level version of \code{dsm}.  The design is restricted to \code{levs}; observations outside those levels receive zeros for this factor.  The first selected level is the scale reference and the remaining selected levels are positive variance ratios.  This representation is most naturally used for random effects; residual factors must still assign exactly one covariance-product coordinate to every observation.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(atm(environment, c("E1","E3")), ism(genotype))
}
}
\seealso{
\code{\link{dsm}}, \code{\link{vsm}}, \code{\link{mmes}}.
}
