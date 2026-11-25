\name{plot.mmes}
\alias{plot.mmes}
\alias{plot.mmer}
\title{plot form a LMM plot with mmes}
\description{\code{plot} method for class \code{"mmes"}.
}
\usage{
\method{plot}{mmes}(x,stnd=TRUE, \dots)
}
\arguments{
 \item{x}{an object of class \code{"mmes"}}
  \item{stnd}{if \code{TRUE}, residuals are studentized, i.e. divided by the square
  root of their conditional variance \code{diag(R - W Ci t(W))}, which holds for any
  residual structure \code{R} (weights, \code{W}, heterogeneous, structured or
  multi-trait). If \code{FALSE}, raw residuals are shown.}
 \item{\dots}{Further arguments to be passed}
}
\details{
  Four panels are drawn: residuals vs fitted values, normal Q-Q, scale-location and
  residuals vs leverage. Leverage is \code{diag(W Ci t(W) solve(R))} with \code{W=[X Z]}.
  \code{Ci} comes from \code{\link{postPEV}} when the fit did not store it (the
  default), so no direct inversion of the coefficient matrix is done. The
  Henderson and direct engines are supported, as well as Gaussian nonlinear curve
  fits with additive random effects (linearized at the solution). GLMM/Laplace,
  rotated, non-Gaussian nonlinear and joint or random-curve nonlinear fits are
  not supported and give an error. Non-scalar residual structures build a dense
  n x n residual covariance, which is limited by the memory setting.
}
\value{
  Invisibly, a list with the \code{fitted}, \code{residuals} (as plotted) and
  \code{leverage} values.
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\examples{
data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
head(DT)
m3 <- mmes(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)

}
\seealso{\code{\link{plot}}, \code{\link{mmes}}}
\keyword{models}

