\name{fam}
\alias{fam}
\title{Factor-analytic covariance structure}
\description{
Factor-analytic covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
fam(x, k = 1L, loadings = NULL, specific = NULL, fixed = NULL,
    variance = c("heterogeneous", "homogeneous"), shared = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{k}{Factor-analytic rank, satisfying 1 <= k < q.}
  \item{loadings}{Optional finite q by k starting loading matrix.}
  \item{specific}{Optional q positive starting specific variances. With
\code{variance="homogeneous"}, supply either one positive value or q equal
values. The default is 0.5 for every level.}
  \item{fixed}{Logical vector controlling all loading and specific-variance-ratio parameters.}
  \item{variance}{Whether to estimate separate specific variances for each level
(\code{"heterogeneous"}, the default) or a common specific variance
(\code{"homogeneous"}).}
  \item{shared}{Optional non-empty loading-group name. Matching names on the
single random and residual FA terms request common trait loadings with
source-specific diagonal factor variances. Requires
\code{mmes_control(representation=list(equations="conditional"))} for AI,
or conditional EM through \code{update=list(method="em")} with that
representation. The earlier \code{update=list(method="conditional-fa")}
selector still selects EM. Ordinary equations reject shared tags rather than
fitting independent loadings.}
}
\details{
The factor-analytic shape is based on \deqn{M=\Lambda\Lambda^{\mathsf T}+\Psi,} where \eqn{\Psi} is diagonal positive, followed by \eqn{K=M/M_{11}}.  For rotational identification, the leading \eqn{k\times k} loading block is lower triangular and its diagonal loadings are positive.  The first specific variance is an internal reference; remaining specific variances are positive ratios.  This removes the otherwise redundant internal scale because \code{vsm} supplies the overall \eqn{\sigma^2}.

With \code{variance="heterogeneous"}, \code{fixed} has length
\eqn{qk-k(k-1)/2+q-1}: triangular loadings followed by q-1 specific-variance
ratios. With \code{variance="homogeneous"}, all specific ratios are one,
\eqn{\Psi=I_q} internally, and \code{fixed} has length \eqn{qk-k(k-1)/2}.
The common specific variance on the fitted scale is estimated, not fixed
at one, and is reported as \code{common_specific_variance}. This homogeneous
model reproduces the covariance family previously fitted by \code{rrm()}.
For a genuine reduced-rank model with all specific variances fixed at zero,
use \code{\link{rrm}}.

With matching \code{shared} groups, the conditional engine fits
\eqn{G=L D_g L^\mathsf{T}+\Psi_g} and
\eqn{R=L D_e L^\mathsf{T}+\Psi_e}, with diagonal positive factor variances
normalized to \eqn{D_g+D_e=I}. Both sources must be entirely free and have
the same trait levels and rank. Starting loading columns, after applying
each term's variance scale, must be positive scalar multiples across sources.
The common loading matrix and factor variance allocations are returned in
\code{fit$conditionalFA$state}. Marginal loading reports use independent
triangular rotations of each source, not the common basis.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(fam(environment, k=2), ism(genotype))
vsm(fam(environment, k=2, variance="homogeneous"), ism(genotype))
}
}
\seealso{
\code{\link{rrm}}, \code{\link{usm}}, \code{\link{vsm}}.
}
