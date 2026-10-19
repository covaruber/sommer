\name{atr}
\alias{atr}
\title{Selected-level diagonal covariance structure for \code{mmer} and \code{vsr}}
\description{
Creates a diagonal covariance-component selector for specified columns or levels of a covariance dimension. It is intended for the \code{mmer} covariance-model interface, typically inside \code{\link{vsr}}.
}
\usage{
atr(x, levs)
}
\arguments{
  \item{x}{
  A factor, character vector, numeric vector, or design/incidence matrix defining the covariance dimension. Factor and character inputs are expanded to incidence columns. A matrix is used directly.
  }
  \item{levs}{
  Character vector identifying the columns or levels whose diagonal entries are selected. If omitted, all available column names are selected.
  }
}
\details{
\code{atr()} constructs a design matrix \eqn{Z} together with a diagonal \code{thetaC} selector. Let \eqn{q} be the number of columns of \eqn{Z}, and define
\deqn{
a_j =
\begin{cases}
1, & \text{if level }j\text{ is included in \code{levs}},\\
0, & \text{otherwise}.
\end{cases}
}
The returned constraint matrix is
\deqn{\mathrm{thetaC}=\mathrm{diag}(a_1,\ldots,a_q).}

Thus \code{atr()} is useful when a covariance component is to be associated only with a specified subset of levels or design columns. Levels not selected in \code{levs} receive zero on the corresponding diagonal of \code{thetaC}.

This function belongs to the \code{mmer}/\code{vsr} covariance-structure interface. The related \code{\link{atm}} function belongs to the newer \code{mmes}/\code{vsm} CovarianceFactor interface and uses a different parameterization.

For a matrix input, its column names define the available levels. For factor or character input, the incidence-matrix column names define the levels. If \code{levs} is missing, all available columns are selected.
}
\value{
A list with components:
\itemize{
  \item \code{Z}: the design/incidence representation of \code{x}.
  \item \code{thetaC}: a diagonal zero/one matrix selecting the levels specified by \code{levs}.
}
}
\seealso{
\code{\link{vsr}}, \code{\link{dsr}}, \code{\link{usr}}, \code{\link{csr}},
\code{\link{atm}}, \code{\link{mmer}}
}
\examples{
# Select levels B and C
A <- atr(factor(c("A","B","A","C")), levs = c("B","C"))
A$Z
A$thetaC
}
