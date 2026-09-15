\name{fam}
\alias{fam}
\title{Factor-analytic covariance structure}
\description{
Factor-analytic covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
fam(x, k = 1L, loadings = NULL, specific = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{k}{Factor-analytic rank, satisfying 1 <= k < q.}
  \item{loadings}{Optional finite q by k starting loading matrix.}
  \item{specific}{Optional q positive starting specific variances.}
  \item{fixed}{Logical vector controlling all loading and specific-variance-ratio parameters.}
}
\details{
The factor-analytic shape is based on \deqn{M=\Lambda\Lambda^{\mathsf T}+\Psi,} where \eqn{\Psi} is diagonal positive, followed by \eqn{K=M/M_{11}}.  For rotational identification, the leading \eqn{k\times k} loading block is lower triangular and its diagonal loadings are positive.  The first specific variance is an internal reference; remaining specific variances are positive ratios.  This removes the otherwise redundant internal scale because \code{vsm} supplies the overall \eqn{\sigma^2}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(fam(environment, k=2), ism(genotype))
}
}
\seealso{
\code{\link{rrcm}}, \code{\link{usm}}, \code{\link{vsm}}.
}
