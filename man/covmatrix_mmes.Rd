\name{covmatrix_mmes}
\alias{covmatrix_mmes}
\title{Covariance and correlation matrix of a fitted term with standard errors}
\description{
Returns the full covariance matrix of a random or residual term of a fitted
\code{\link{mmes}} model over its covariance coordinates (for example the
genetic trait-by-trait covariance matrix of a multi-trait model), together
with the implied correlation matrix and delta-method standard errors of every
element.
}
\usage{
covmatrix_mmes(object, term, se = TRUE, max.dim = 500L, rel_step = 1e-5)
}
\arguments{
  \item{object}{A fitted \code{\link{mmes}} model.}
  \item{term}{Name (as in \code{names(object$covStruct)}) or position of the
  term. The residual term is the last one.}
  \item{se}{If \code{TRUE}, add standard errors.}
  \item{max.dim}{Refuse terms whose covariance coordinate space is larger than
  this (the matrices are dense).}
  \item{rel_step}{Relative step of the central differences used for the
  Jacobian.}
}
\details{
The matrix is \eqn{\sigma^2 K_1 \otimes \cdots \otimes K_m} built from all the
covariance factors of the \code{vsm} term (the relationship matrix
\code{Gu} is not included), so for \code{vsm(usm(trait), ism(id))} it is the
trait covariance matrix, and for \code{vsm(dsm(env), usm(trait), ism(id))} the
environment-by-trait matrix. Row and column names are the coordinate levels.

Standard errors use the delta method,
\eqn{\mathrm{Var}(g(\hat\theta)) \approx J V J^\top}, where \eqn{V} is
\code{object$theta_se} (the covariance of the reported parameters, as used by
\code{\link{covparams_mmes_se}}) and \eqn{J} the Jacobian of each element with
respect to the reported parameters. Fixed parameters contribute no variance.
Diagonal standard errors of the correlation matrix are zero. Terms split by
\code{\link{dsumm}} sections are not supported.
}
\value{A list with \code{covariance} and \code{correlation} and, if
\code{se=TRUE}, \code{covariance.se} and \code{correlation.se}.}
\seealso{\code{\link{covparams_mmes}}, \code{\link{covparams_mmes_se}},
\code{\link{vpredict}}, \code{\link{stackTraits}}}
\examples{
set.seed(17)
n <- 60
DT <- data.frame(id = factor(rep(seq_len(n), each = 3)))
u <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, 0.35, 0.35, 1.4), 2))
DT$trait_a <- u[DT$id, 1] + rnorm(nrow(DT), sd = 0.7)
DT$trait_b <- u[DT$id, 2] + rnorm(nrow(DT), sd = 1)
DTL <- stackTraits(DT, traits = c("trait_a", "trait_b"), keep = "id")
fit <- mmes(value ~ trait,
            random = ~ vsm(usm(trait), ism(id)),
            rcov = ~ vsm(usm(trait), ism(record)),
            data = DTL, verbose = FALSE)
G <- covmatrix_mmes(fit, 1)       # genetic
G$covariance
G$correlation
G$correlation.se
R <- covmatrix_mmes(fit, 2)       # residual
R$correlation
}
\keyword{models}
