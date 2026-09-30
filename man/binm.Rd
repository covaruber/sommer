\name{binm}
\alias{binm}
\title{Binomial proportion response with its number of trials}
\description{
Builds a single proportion response for binomial \code{\link{mmes}} fits from
counts of successes and failures (or trials). The number of trials travels with
the proportion as an attribute and is used as the binomial prior weight, so the
single-response interface of \code{mmes} is preserved.
}
\usage{
binm(successes, failures = NULL, trials = NULL)
}
\arguments{
  \item{successes}{Numeric vector of successes.}
  \item{failures}{Numeric vector of failures. Provide either \code{failures} or \code{trials}.}
  \item{trials}{Numeric vector with the number of trials.}
}
\details{
The result is \code{successes/trials} with class \code{"binm"} and attribute
\code{"trials"}. Subsetting keeps the trials aligned, so the vector can be stored
as a data-frame column. Missing counts give a missing response, handled by the
\code{naMethodY} rules of \code{mmes}. Fitting a proportion without \code{binm}
treats every value as one trial, which is wrong unless all trials equal one.
}
\value{A numeric vector of proportions of class \code{"binm"}.}
\seealso{\code{\link{mmes}}}
\examples{
set.seed(1)
d <- data.frame(x = rnorm(60), g = factor(rep(1:12, each = 5)))
d$n <- sample(5:10, 60, TRUE)
d$s <- rbinom(60, d$n, plogis(0.3 * d$x))
fit <- mmes(binm(s, trials = n) ~ x, random = ~ g, data = d,
            family = binomial(), verbose = FALSE)
fit$b
}
\keyword{models}
