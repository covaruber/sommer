\name{bathy.colors}
\alias{bathy.colors}
\title{Generate a sequence of colors for plotting bathymetric data.}
\usage{
  bathy.colors(n, alpha = 1)
}
\arguments{
  \item{n}{The number of colors to return.}

  \item{alpha}{Alpha values to be passed to \code{rgb()}.}
}
\value{
  A vector of blue scale colors.
}
\description{
  \code{bathy.colors(n)} generates a sequence of \eqn{n}
  colors along a linear scale from light grey to pure blue.
}
\examples{
{
# Plot a colorbar using bathy.colors
image(matrix(seq(100), 100), col=bathy.colors(100))
}
}

