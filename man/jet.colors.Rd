\name{jet.colors}
\alias{jet.colors}
\title{Generate a sequence of colors alog the jet colormap.}
\usage{
  jet.colors(n, alpha = 1)
}
\arguments{
  \item{n}{The number of colors to return.}

  \item{alpha}{The transparency value of the colors. See
  \code{?rgb} for details.}
}
\value{
  A vector of colors along the jet colorramp.
}
\description{
  \code{jet.colors(n)} generates a sequence of \eqn{n}
  colors from dark blue to cyan to yellow to dark red. It
  is similar to the default color schemes in Python's
  matplotlib or MATLAB.
}
\examples{
{
# Plot a colorbar with jet.colors
image(matrix(seq(100), 100), col=jet.colors(100))
}
}
\seealso{The core function of the package \code{\link{mmer}} }

