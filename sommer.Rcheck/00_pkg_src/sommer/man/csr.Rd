\name{csr}
\alias{csr}
\title{User-specified covariance-component structure for \code{mmer} and \code{vsr}}
\description{
Combines a design/incidence representation with a user-supplied covariance-component constraint matrix for use by the \code{mmer} covariance-model interface, typically inside \code{\link{vsr}}.
}
\usage{
csr(x, mm)
}
\arguments{
  \item{x}{
  A factor, character vector, numeric vector, or design/incidence matrix defining the covariance dimension. Factor and character inputs are expanded to incidence columns. A matrix is used directly.
  }
  \item{mm}{
  A square covariance-component constraint matrix. Its dimensions must correspond to the columns of the design matrix generated from \code{x}. The matrix is returned as \code{thetaC}, with its row and column names set to the column names of the resulting design matrix.
  }
}
\details{
\code{csr()} provides the general constraint-matrix interface underlying custom covariance patterns for \code{mmer}/\code{vsr}. Unlike \code{\link{dsr}} and \code{\link{usr}}, it does not construct a predefined diagonal or unstructured pattern. Instead, the supplied matrix \code{mm} determines the covariance-component structure through \code{thetaC}.

Conceptually, \code{thetaC} identifies which elements of the covariance matrix are associated with estimated variance-covariance components and which constraints are shared among elements. The interpretation of the entries therefore follows the \code{thetaC} convention used by \code{\link{vsr}}.

When \code{x} is a matrix it is used directly as \eqn{Z}. Otherwise, factors and character vectors are converted to incidence matrices, while numeric non-factor inputs are treated as a single design column.

The function assigns the column names of \eqn{Z} to both dimensions of \code{mm}; therefore \code{mm} should be conformable with the number of columns of the resulting design matrix.
}
\value{
A list with components:
\itemize{
  \item \code{Z}: the design/incidence matrix associated with the covariance dimension.
  \item \code{thetaC}: the user-supplied covariance-component constraint matrix, with dimensions named according to \code{Z}.
}
}
\seealso{
\code{\link{vsr}}, \code{\link{dsr}}, \code{\link{usr}}, \code{\link{atr}},
\code{\link{mmer}}
}
\examples{
# A custom two-level covariance-component pattern
x <- factor(c("A","B","A","B"))
mm <- matrix(c(1, 2,
               2, 1), 2, 2, byrow = TRUE)
C <- csr(x, mm)
C$Z
C$thetaC
}
