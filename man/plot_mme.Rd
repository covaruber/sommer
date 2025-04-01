\name{plot.mme}
\alias{plot.mme}
\title{plot form a LMM plot with mme}
\description{\code{plot} method for class \code{"mme"}.
}
\usage{
\method{plot}{mme}(x,stnd=TRUE, \dots)
}
\arguments{
 \item{x}{an object of class \code{"mme"}}
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
m3 <- mme(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
plot(m3)
}
\seealso{\code{\link{plot}}, \code{\link{mme}}}
\keyword{models}

