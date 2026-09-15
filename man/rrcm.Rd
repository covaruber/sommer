\name{rrcm}
\alias{rrcm}
\title{Reduced-rank covariance structure}
\description{
Reduced-rank covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
rrcm(x, k = 1L, loadings = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{k}{Reduced rank, satisfying 1 <= k < q.}
  \item{loadings}{Optional finite q by k starting loading matrix.}
  \item{fixed}{Logical vector controlling the free loading parameters.}
}
\details{
The reduced-rank structure uses \deqn{M=\Lambda\Lambda^{\mathsf T}+I_q,\qquad K=M/M_{11}.} The rank-\eqn{k} term captures the dominant covariance pattern and the identity term provides a common isotropic remainder, keeping the covariance strictly positive definite for precision-based Henderson calculations.  The leading \eqn{k\times k} loading block is lower triangular for rotational identification, and its diagonal loadings are positive.  \code{rrcm} is compiled through the generic CovarianceFactor callback interface rather than requiring model-specific solver code.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(rrcm(environment, k=2), ism(genotype))
}
}
\seealso{
\code{\link{fam}}, \code{\link{ownm}}, \code{\link{vsm}}.
}
