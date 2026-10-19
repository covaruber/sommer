\name{usr}
\alias{usr}
\title{Unstructured covariance structure for \code{mmer} and \code{vsr}}
\description{
Creates an unstructured variance-covariance component pattern for use by the \code{mmer} covariance-model interface, typically inside \code{\link{vsr}}.
}
\usage{
usr(x)
}
\arguments{
  \item{x}{
  A factor, character vector, numeric vector, or design/incidence matrix defining the covariance dimension. Factor and character inputs are expanded to incidence columns. A matrix is used directly.
  }
}
\details{
\code{usr()} is part of the covariance-structure interface used by \code{mmer}/\code{vsr}; it is distinct from \code{\link{usm}}, which belongs to the newer \code{mmes}/\code{vsm} CovarianceFactor interface.

After constructing the design matrix \eqn{Z}, the function calls \code{unsm(q)}, where \eqn{q} is the number of columns of \eqn{Z}, to create the \code{thetaC} pattern for an unstructured covariance matrix. Thus the associated covariance model permits a separate variance for every level and a separate covariance for every pair of levels:
\deqn{
\Sigma =
\begin{bmatrix}
\sigma_1^2 & \sigma_{12} & \cdots & \sigma_{1q}\\
\sigma_{12} & \sigma_2^2 & \cdots & \sigma_{2q}\\
\vdots & \vdots & \ddots & \vdots\\
\sigma_{1q} & \sigma_{2q} & \cdots & \sigma_q^2
\end{bmatrix}.
}
The exact component coding is stored in the matrix returned by \code{unsm(q)} and passed to the \code{vsr} machinery through \code{thetaC}.

For a factor or character vector, missing observations are retained through the construction of the incidence matrix. For a numeric non-factor input, \code{x} is treated as a single design column.
}
\value{
A list with components:
\itemize{
  \item \code{Z}: the design/incidence matrix associated with the covariance dimension.
  \item \code{thetaC}: the unstructured covariance-component pattern produced by \code{unsm()}.
}
}
\seealso{
\code{\link{vsr}}, \code{\link{dsr}}, \code{\link{csr}}, \code{\link{atr}},
\code{\link{usm}}, \code{\link{mmer}}
}
\examples{
# Unstructured covariance pattern across three levels
U <- usr(factor(c("A","B","A","C")))
U$Z
U$thetaC
}
