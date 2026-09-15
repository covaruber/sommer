\name{dsr}
\alias{dsr}
\title{Diagonal covariance structure for \code{mmer} and \code{vsr}}
\description{
Creates a diagonal variance-component structure for use by the \code{mmer} covariance-model interface, typically inside \code{\link{vsr}}. Each column or level represented by \code{x} receives its own diagonal variance component, while all off-diagonal covariance components are constrained to zero.
}
\usage{
dsr(x)
}
\arguments{
  \item{x}{
  A factor, character vector, numeric vector, or design/incidence matrix defining the covariance dimension. Factor and character inputs are expanded to incidence columns. A matrix is used directly.
  }
}
\details{
\code{dsr()} is part of the covariance-structure interface used by \code{mmer}/\code{vsr}; it is distinct from \code{\link{dsm}}, which belongs to the newer \code{mmes}/\code{vsm} CovarianceFactor interface.

The function returns a design matrix \eqn{Z} and a covariance-parameter constraint matrix \code{thetaC}. If \eqn{q} columns or levels are represented by \code{x}, then
\deqn{\mathrm{thetaC}=I_q.}
Consequently, the associated covariance matrix has the form
\deqn{\Sigma=\mathrm{diag}(\sigma_1^2,\ldots,\sigma_q^2),}
so that the \eqn{q} variances are estimated separately and all covariances are fixed to zero.

For a factor or character vector, missing observations are retained through the construction of the incidence matrix. If only one non-missing level is present, a one-column incidence matrix is constructed explicitly.

For a numeric non-factor input, \code{x} is treated as a single design column rather than expanded into factor levels.
}
\value{
A list with components:
\itemize{
  \item \code{Z}: the design/incidence matrix associated with the covariance dimension.
  \item \code{thetaC}: a diagonal matrix defining the diagonal covariance-component pattern used by \code{vsr}.
}
}
\seealso{
\code{\link{vsr}}, \code{\link{usr}}, \code{\link{csr}}, \code{\link{atr}},
\code{\link{dsm}}, \code{\link{mmer}}
}
\examples{
# Diagonal covariance structure across levels
D <- dsr(factor(c("A","B","A","C")))
D$Z
D$thetaC
}
