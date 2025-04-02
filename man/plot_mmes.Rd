\name{plot.mmes}
\alias{plot.mmes}
\title{plot form a LMM plot with mmes}
\description{\code{plot} method for class \code{"mmes"}.
}
\usage{
\method{plot}{mmes}(x,stnd=TRUE, \dots)
}
\arguments{
 \item{x}{an object of class \code{"mmes"}}
  \item{stnd}{argument for ploting the residuals to know if they should be standarized.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of plot
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\examples{
data(DT_yatesoats)
DT <- DT_yatesoats
head(DT)
m3 <- mmes(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
plot(m3)
}
\seealso{\code{\link{plot}}, \code{\link{mmes}}}
\keyword{models}

