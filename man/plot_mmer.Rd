\name{plot.mmer}
\alias{plot.mmer}
\title{plot form a LMM plot with mmer}
\description{\code{plot} method for class \code{"mmer"}.
}
\usage{
\method{plot}{mmer}(x,stnd=TRUE, \dots)
}
\arguments{
 \item{x}{an object of class \code{"mmer"}}
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
m3 <- mmer(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
plot(m3)
}
\seealso{\code{\link{plot}}, \code{\link{mmer}}}
\keyword{models}

