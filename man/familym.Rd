\name{familym}
\alias{familym}
\title{Different distribution families for the traits of a multi-trait model}
\description{
Assigns a \code{stats::family()} object to each trait of a long-format
multi-trait \code{\link{mmes}} model, e.g. a binary disease score together
with a Gaussian yield, so that their genetic and residual correlations can be
estimated jointly by penalized quasi-likelihood (PQL).
}
\usage{
familym(by, ...)
}
\arguments{
  \item{by}{Name (character string) of the column of \code{data} that
  identifies the traits, typically \code{"trait"} as created by
  \code{\link{stackTraits}}.}
  \item{\dots}{Family objects named by the levels of the \code{by} column,
  e.g. \code{Disease = binomial(), Yield = gaussian()}. Every level present in
  the data needs a family.}
}
\details{
The object is passed as the \code{family} argument of \code{mmes}. Each PQL
iteration builds the working response and working weights of every row with
the family of its trait, and fits a weighted Gaussian multi-trait model.

The residual variance of binomial, Poisson and negative binomial traits is not
a free parameter (their dispersion is one), while that of Gaussian and other
families is estimated. The residual structure must therefore allow individual
trait variances to be fixed: use
\code{rcov = ~ vsm(corgm(trait, variance = "heterogeneous"), ism(record))}
(trait correlations and heterogeneous variances, see \code{\link{corgm}}) or
\code{rcov = ~ vsm(dsm(trait), ism(units))} (uncorrelated residuals).
The first level of the trait factor must be a fixed-dispersion trait,
because its variance is the \code{vsm} scale; reorder the levels
otherwise. \code{mmes} then fixes the residual variance of every
fixed-dispersion trait to one and estimates the rest.

As for single-family PQL fits, the likelihood, \code{AIC} and \code{BIC} are
\code{NA}. The fitted object has a named \code{dispersion} vector (one per
trait: one for fixed-dispersion traits, the REML residual variance for the
others), and \code{fitted()} and \code{residuals()} apply the inverse link and
variance function of each trait to its rows. Random-effect variances are on
the link scale of each trait.
}
\value{An object of class \code{"sommer_familym"}.}
\seealso{\code{\link{mmes}}, \code{\link{corgm}}, \code{\link{stackTraits}},
\code{\link{covmatrix_mmes}}}
\examples{
set.seed(12)
ng <- 100; nr <- 4
G <- matrix(c(1, 0.4, 0.4, 0.8), 2)
u <- matrix(rnorm(ng * 2), ng) \%*\% chol(G)
d <- data.frame(id = factor(rep(paste0("g", 1:ng), each = nr)))
d$disease <- rbinom(nrow(d), 1, plogis(-0.5 + u[as.integer(d$id), 1]))
d$yield <- 10 + u[as.integer(d$id), 2] + rnorm(nrow(d))
L <- stackTraits(d, traits = c("disease", "yield"))

fit <- mmes(value ~ trait,
            random = ~ vsm(usm(trait), ism(id)),
            rcov = ~ vsm(corgm(trait, variance = "heterogeneous"), ism(record)),
            family = familym("trait", disease = binomial(), yield = gaussian()),
            data = L, verbose = FALSE)
fit$dispersion
covmatrix_mmes(fit, 1)$correlation   # genetic correlation (link scale)
head(fitted(fit))                    # probabilities for disease rows
}
\keyword{models}
