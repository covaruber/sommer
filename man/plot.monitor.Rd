\name{pmonitor}
\alias{pmonitor}
\title{plot the change of VC across iterations}
\description{\code{plot} for monitoring.
}
\usage{
pmonitor(object, \dots)
}
\arguments{
 \item{object}{model object of class \code{"mmes"}}
 \item{\dots}{Further arguments to be passed to the plot function.}
}
\value{
  vector of plot
 }
\author{Giovanny Covarrubias}
\examples{
data(DT_yatesoats)
DT <- DT_yatesoats
head(DT)
m3 <- mmes(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
pmonitor(m3)
}
\seealso{\code{\link{plot}}, \code{\link{mmes}}}
\keyword{models}

