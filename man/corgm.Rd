\name{corgm}
\alias{corgm}
\title{General positive-definite correlation structure}
\description{
General positive-definite correlation structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
corgm(x, theta = NULL, fixed = NULL,
      variance = c("homogeneous", "heterogeneous"), values = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{theta}{Optional positive-definite q by q starting covariance or correlation matrix.
  With \code{variance="heterogeneous"} its diagonal is used as starting
  variances when \code{values} is not given.}
  \item{fixed}{Logical vector flagging parameters held at their starting value:
  length q(q-1)/2 (correlations) for the homogeneous mode and
  q(q-1)/2 + q-1 (correlations, then variance ratios of levels 2..q) for the
  heterogeneous mode.}
  \item{variance}{\code{"homogeneous"} (default) gives a pure correlation
  matrix scaled by the single variance of \code{vsm}.
  \code{"heterogeneous"} adds one variance per level.}
  \item{values}{Optional positive starting variances (length q) for the
  heterogeneous mode. Only their ratios to the first level matter; the first
  level's variance is the \code{vsm} scale. When omitted, the starting ratios
  are taken from the data.}
}
\details{
An unrestricted SPD correlation matrix is represented using a unit-diagonal lower factor \eqn{A}.  The intermediate matrix \eqn{S=AA^{\mathsf T}} is standardized to correlation scale, \deqn{K_{ij}=S_{ij}/\sqrt{S_{ii}S_{jj}}.} The unrestricted working parameters are the strict lower-triangular entries of \eqn{A}.

With \code{variance="heterogeneous"} the shape is \eqn{K = DRD} with
\eqn{D=\mathrm{diag}(1, \sqrt{r_2}, \ldots, \sqrt{r_q})}, where
\eqn{r_l=\sigma^2_l/\sigma^2_1 = \exp(\phi_l)} are variance ratios relative to
the first level; the \code{vsm} scale is \eqn{\sigma^2_1}. The model is the
same as \code{\link{usm}} but parameterized as correlations and variances,
which allows individual variances to be fixed. This is what mixed-family
multi-trait models (\code{\link{familym}}) need: the residual variance of a
binomial or Poisson trait is fixed to one while the correlations and the
variances of Gaussian traits are estimated. Reported parameters
(\code{\link{covparams_mmes}}) are the correlations and the variances of each
level.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
data(DT_example, package = "enhancer")
DT <- DT_example
# genetic correlations between environments with a homogeneous variance
m1 <- mmes(Yield ~ Env, random = ~ vsm(corgm(Env), ism(Name)),
           rcov = ~ units, data = DT, verbose = FALSE)
covparams_mmes(m1, 1)

# heterogeneous variances (equivalent to usm, parameterized as corr/var)
m2 <- mmes(Yield ~ Env,
           random = ~ vsm(corgm(Env, variance = "heterogeneous"), ism(Name)),
           rcov = ~ units, data = DT, verbose = FALSE)
covparams_mmes(m2, 1)
covmatrix_mmes(m2, 1)$correlation
}
\seealso{
\code{\link{csm}}, \code{\link{usm}}, \code{\link{vsm}},
\code{\link{covmatrix_mmes}}, \code{\link{familym}}.
}
