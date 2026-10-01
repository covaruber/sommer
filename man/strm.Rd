\name{strm}
\alias{strm}
\title{Covariance structure spanning several random terms}
\description{
Combines two or more random terms built with \code{\link{vsm}} into one
random structure with an estimated covariance between the terms. Typical uses are
direct-maternal animal models, correlated
permanent-environment effects and indirect genetic effects.
}
\usage{
strm(..., cov = usm, Gu = NULL, labels = NULL, sigma2 = NULL,
     fixedSigma2 = FALSE)
}
\arguments{
  \item{\dots}{Two or more random \code{vsm()} terms, optionally named (names
  become the term labels), e.g. \code{dir = vsm(ism(id)), mat = vsm(ism(dam))}.}
  \item{cov}{A covariance constructor applied to the term index, e.g.
  \code{usm} (default), \code{dsm}, \code{csm}, \code{corgm}, or
  \code{function(x) fam(x, k = 1)}.}
  \item{Gu}{Precision matrix shared by all terms (\code{attr(Gu, "inverse") = TRUE}).
  If \code{NULL}, a \code{Gu} given inside the terms is used (all terms must
  then use the same one) or an identity over the union of the terms' levels.}
  \item{labels}{Term labels; default the argument names or \code{t1, t2, ...}.}
  \item{sigma2}{Starting overall variance scale; \code{NULL} lets \code{mmes} choose.}
  \item{fixedSigma2}{Fix the overall variance scale.}
}
\details{
For terms \eqn{u_1,\dots,u_T} with incidence matrices \eqn{Z_1,\dots,Z_T}
the structure is
\deqn{\mathrm{Var}(u_1,\dots,u_T)=\sigma^2\,K_{terms}\otimes K_{inner}\otimes A,}
where \eqn{K_{terms}} is given by \code{cov}, \eqn{K_{inner}} is the Kronecker
product of the covariance factors used inside the terms (e.g.
\code{dsm(env)}), and \eqn{A} is the relationship matrix. The requirements are
those of a direct product: all terms share the
coefficient levels (their levels are aligned on the union, or on the levels of
\code{Gu}), the same \code{Gu}, and the same inner covariance constructors with
the same levels; the parameters of the inner factors are shared. The
\code{sigma2} of individual terms is ignored. BLUPs are returned in one
\code{uList} element with a column per term (and inner level).
\code{covm(ran1, ran2)} is the two-term special case with \code{cov = usm}.

The terms do not need to be adjacent in the model
formula and their order is checked by level names rather than assumed.
}
\value{A random structure for the \code{random} formula of \code{\link{mmes}}.}
\seealso{\code{\link{covm}}, \code{\link{vsm}}, \code{\link{usm}}, \code{\link{mmes}}}
\examples{
set.seed(3)
ids <- paste0("i", 1:30)
A <- diag(30) + diag(0.2, 30); dimnames(A) <- list(ids, ids)
Ai <- solve(A); attr(Ai, "inverse") <- TRUE
d <- data.frame(id = factor(sample(ids, 200, TRUE), levels = ids),
                dam = factor(sample(ids[1:15], 200, TRUE), levels = ids))
d$y <- rnorm(200)
fit <- mmes(y ~ 1, random = ~ strm(dir = vsm(ism(id)), mat = vsm(ism(dam)), Gu = Ai),
            data = d, verbose = FALSE)
covparams_mmes(fit, 1)
}
\keyword{models}
