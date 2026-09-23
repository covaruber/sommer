\name{loadings_mmes}
\alias{loadings_mmes}
\title{Extract factor-analytic loadings from a fitted mmes model}
\description{
Reconstructs the normalized loadings matrix and specific variances of a
\code{fam()} or \code{rrcm()} covariance-shaping factor fitted inside
\code{vsm()}, on the same scale as the fitted covariance matrix returned in
\code{object$theta}.
}
\usage{
loadings_mmes(object, term = NULL, varianceScale=TRUE, rotation=TRUE)
}
\arguments{

  \item{object}{a fitted model of class \code{"mmes"}.}
  \item{term}{character name of the random term (matching a name in
\code{object$covStruct}) built with a single \code{fam()} or \code{rrcm()}
covariance-shaping factor. If \code{NULL}, the unique such term in the model
is used automatically; an error is raised if none or more than one exist.}
\item{varianceScale}{a logical argument to indicate if loadings should be returned 
in variance scale (multiplied by sqrt(sigma2)).}
\item{rotation}{a logical value to indicate if loadings should be rotated by its
singular vectors.}

}
\details{
\code{fam()}/\code{rrcm()} remove their otherwise redundant internal scale so
that \code{vsm()} owns the single overall variance; the loadings and specific
variances reported by \code{loadings_mmes} are re-normalized so that
\deqn{\Sigma = \sigma^2(\Lambda\Lambda^{\mathsf T}+\Psi)}
reproduces \code{object$theta[[term]]} exactly, where \eqn{\sigma^2} is
\code{object$covPar[[term]][1]}.
}
\value{
A list containing:
\describe{
\item{loadings}{a levels-by-factor matrix \eqn{\Lambda}.}
\item{specific}{a named vector of specific variances \eqn{\Psi} (identically
1 for every level in a \code{rrcm()} term).}
\item{sigma2}{the fitted overall variance scale.}
\item{model}{\code{"fa"} or \code{"rr"}.}
\item{term}{the resolved term name.}
}
}
\examples{
\dontrun{
mix <- mmes(BLUEs ~ trial,
            random = ~ vsm(fam(trial, 2), ism(genotype)),
            rcov = ~ units, data = dt)
loadings_mmes(mix)
}
}
\seealso{
\code{\link{scores_mmes}}, \code{\link{fam}}, \code{\link{rrcm}}, \code{\link{vsm}}, \code{\link{mmes}}.
}
