\name{metricm}
\alias{metricm}
\title{Distance-based (metric) spatial correlation structures}
\description{
Exponential, power, gaussian, spherical and circular correlation structures
computed from distances between locations, isotropic or anisotropic. The
returned dimensionless covariance shape is intended for use inside
\code{\link{vsm}}, which supplies the single overall variance scale.
}
\usage{
metricm(x, model = c("exponential", "power", "gaussian", "spherical", "circular"),
        range = NULL, rho = NULL, anisotropy = c("none", "product", "geometric"),
        angle = pi/2, ratio = 1.5, metric = c("euclidean", "manhattan"),
        distance = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Numeric coordinate vector (1-D) or matrix/data frame whose rows give coordinates. Repeated coordinates share one covariance level.}
  \item{model}{Correlation function, see Details.}
  \item{range}{Starting range \eqn{\phi} (one per axis with product anisotropy). Defaults to the median positive distance.}
  \item{rho}{Starting correlation per unit distance for \code{model="power"}, in (0,1). Defaults to \eqn{e^{-1}} at the median distance.}
  \item{anisotropy}{\code{"none"}; \code{"product"}, one parameter per axis, \eqn{K=f_x(|dx|)f_y(|dy|)}; or \code{"geometric"}, coordinates rotated by \code{angle} and the second rotated axis divided by \code{ratio}. Anisotropy requires 2-D coordinates.}
  \item{angle}{Starting anisotropy angle in (0, pi).}
  \item{ratio}{Starting positive anisotropy ratio.}
  \item{metric}{Distance metric for isotropic models; \code{"manhattan"} is only valid (positive definite) for the exponential and power models.}
  \item{distance}{Optional symmetric q by q distance matrix with zero diagonal, aligned to the unique locations (isotropic models only).}
  \item{fixed}{Logical vector, one value per parameter (in the order range/rho, then angle and ratio), recycled if of length one.}
}
\details{
With \eqn{h=d/\phi}:
\tabular{llll}{
\code{model} \tab \eqn{K(d)} \tab validity \tab ASReml \cr
exponential \tab \eqn{e^{-h}} \tab any dimension \tab \code{exp}, \code{iexp} \cr
power \tab \eqn{\rho^{d}} \tab any dimension \tab \code{exp}, \code{iexp} (same parameterization) \cr
gaussian \tab \eqn{e^{-h^2}} \tab any dimension \tab \code{gau}, \code{igau} \cr
spherical \tab \eqn{1-\frac32h+\frac12h^3} for \eqn{h<1}, else 0 \tab up to 3-D \tab \code{isp} \cr
circular \tab \eqn{\frac{2}{\pi}(\arccos h-h\sqrt{1-h^2})} for \eqn{h<1}, else 0 \tab up to 2-D \tab \code{cir} \cr
}
Product anisotropy (exponential, power and gaussian) corresponds to ASReml
\code{aexp}/\code{agau}; geometric anisotropy corresponds to ASReml's
anisotropic Matern and is also available in \code{\link{maternm}}. The power
model is the exponential model with \eqn{\phi=-1/\log\rho}, and the
exponential model equals \code{maternm(nu=0.5)}.

Ranges and ratios are optimized on log scales, correlations and angles through
bounded logits, with analytic derivatives. The structure is dense over the q
locations, so it is intended for up to a few thousand locations. On a complete
row-by-column grid, product anisotropy is a Kronecker product of two 1-D
structures; \code{vsm(metricm(row, "power"), metricm(col, "power"), ism(units))}
fits the same model while keeping the compact grid engine. The gaussian model
is nearly singular for close locations; keep a residual (nugget) term. With
geometric anisotropy the angle is not identifiable when the ratio is close to
one, and the ratio can drift to extreme values when the data are isotropic.
}
\value{
A list with the incidence matrix \code{Z} and a CovarianceFactor descriptor \code{covFactor}, for use inside \code{vsm}.
}
\examples{
set.seed(1)
co <- cbind(x = runif(30, 0, 10), y = runif(30, 0, 10))
d <- data.frame(x = rep(co[, 1], each = 2), y = rep(co[, 2], each = 2), one = factor(1))
d$yield <- rnorm(60)
fit <- mmes(yield ~ 1, random = ~ vsm(metricm(cbind(x, y), model = "power"), ism(one)),
            data = d, verbose = FALSE)
covparams_mmes(fit)
}
\seealso{\code{\link{maternm}}, \code{\link{ar1m}}, \code{\link{sar}}, \code{\link{car}}, \code{\link{vsm}}.}
\keyword{models}
