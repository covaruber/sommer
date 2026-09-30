\name{wald_mmes}
\alias{wald_mmes}
\alias{print.wald.mmes}
\title{Wald tests for fixed effects with optional small-sample denominator df}
\description{
Pseudo analysis of variance for the fixed terms of an \code{\link{mmes}} fit,
based on incremental or conditional Wald statistics, with chi-square,
residual, Satterthwaite or Kenward-Roger reference distributions.
}
\usage{
wald_mmes(object, ssType = c("incremental", "conditional"),
          denDF = c("none", "residual", "satterthwaite", "kr"), terms = NULL)
}
\arguments{
  \item{object}{a model of class \code{"mmes"}.}
  \item{ssType}{\code{"incremental"} tests each term adjusted for the terms
  preceding it in the fixed formula. \code{"conditional"} tests each term
  adjusted for every other term that does not contain it, so marginality is
  respected: a term is never adjusted for a term that contains it
  structurally (A:B contains A) or implicitly (its columns, together with the
  intercept, span the tested term, e.g. locations coded separately within
  regions).}
  \item{denDF}{\code{"none"} gives chi-square tests of the Wald statistic.
  \code{"residual"} uses F tests with n - rank(X) denominator df.
  \code{"satterthwaite"} and \code{"kr"} (Kenward-Roger) compute approximate
  denominator df; \code{"kr"} also uses the Kenward-Roger adjusted covariance
  of the fixed effects and scaled F statistic.}
  \item{terms}{optional subset of fixed terms (as named in \code{object$Dtable}).}
}
\details{
Let \eqn{\Phi} be the covariance of the fixed effects and \eqn{M=\Phi^{-1}}.
For a given order of the terms, the upper Cholesky factor \eqn{U} of \eqn{M}
defines \eqn{z=U b}; the incremental Wald statistic of a term is the sum of the
squared elements of \eqn{z} belonging to it, i.e. a test of \eqn{L b=0} with
\eqn{L} the corresponding rows of \eqn{U}. Conditional tests use the order
[terms not containing the tested term, tested term, terms containing it].

For Satterthwaite and Kenward-Roger df the model inputs are rebuilt from the
stored call and \eqn{X^\top V^{-1}X} is re-evaluated at perturbed covariance
parameters (on the reported scale of \code{covPar}); its numerical first and
second derivatives give the quantities \eqn{P_i} and the Kenward-Roger
adjustment. The covariance of the covariance parameters is \code{theta_se},
the inverse average information. The second-derivative term of \eqn{V}
enters through the Hessian of \eqn{X^\top V^{-1}X}; it vanishes for
variance-component models, where the method is the exact Kenward-Roger
approximation. Because the average (rather than expected) information is
used, df can differ slightly from software using the expected information.
The re-evaluation requires the data and objects used in the original call
to be available, REML estimation, and a Gaussian fit.

\code{anova(object)} with a single model returns \code{wald_mmes(object, ...)}.
}
\value{A data frame of class \code{"wald.mmes"} with columns \code{Df},
\code{denDF}, \code{Wald}, \code{F.value} and \code{p.value}, one row per term.}
\references{
Kenward MG and Roger JH (1997). Small sample inference for fixed effects from
restricted maximum likelihood. Biometrics 53: 983-997.
}
\seealso{\code{\link{mmes}}, \code{\link{predict.mmes}}, \code{\link{anova.mmes}}}
\examples{
data(DT_yatesoats, package = "enhancer")
m <- mmes(Y ~ V * N, random = ~ B + B:MP, rcov = ~ units,
          data = DT_yatesoats, verbose = FALSE)
wald_mmes(m)
wald_mmes(m, ssType = "conditional", denDF = "kr")
}
\keyword{models}
