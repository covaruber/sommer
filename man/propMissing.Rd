\name{propMissing}
\alias{propMissing}
\title{Proportion of missing data}
\description{
  \code{propMissing} quick calculation of the proportion of missing data in a vector.
}
\usage{
  propMissing(x)
}

\arguments{
  \item{x}{vector of observations.}
}
\value{
  \describe{
    \item{$res}{ a numeric value with the proportion of missing data.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommec. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
x <- c(1:10, NA)
propMissing(x)
}

