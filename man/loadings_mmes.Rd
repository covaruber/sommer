\name{loadings_mmes}
\alias{loadings_mmes}
\title{Extract factor-analytic loadings from a fitted mmes model}
\description{
Reconstructs the normalized loadings matrix and specific variances of a
\code{fam()} or \code{rrm()} covariance-shaping factor fitted inside
\code{vsm()}, on the same scale as the fitted covariance matrix returned in
\code{object$theta}.
}
\usage{
loadings_mmes(object, term = NULL, varianceScale=TRUE, rotation=TRUE)
}
\arguments{

  \item{object}{a fitted model of class \code{"mmes"}.}
  \item{term}{character name of the random term (matching a name in
\code{object$covStruct}) built with a single \code{fam()} or \code{rrm()}
covariance-shaping factor. If \code{NULL}, the unique such term in the model
is used automatically; an error is raised if none or more than one exist.}
\item{varianceScale}{Whether to return loadings multiplied by
\code{sqrt(sigma2)} and specific variances multiplied by \code{sigma2}.
When \code{FALSE}, both describe the dimensionless covariance shape.}
\item{rotation}{a logical value to indicate if loadings should be rotated by its
singular vectors.}

}
\details{
\code{fam()}/\code{rrm()} remove their otherwise redundant internal scale so
that \code{vsm()} owns the single overall variance; the loadings and specific
variances reported by \code{loadings_mmes} are re-normalized so that
\deqn{\Sigma = \sigma^2(\Lambda\Lambda^{\mathsf T}+\Psi)}
reproduces \code{object$theta[[term]]} exactly when
\code{varianceScale=FALSE}, where \eqn{\sigma^2} is
\code{object$covPar[[term]][1]}. With \code{varianceScale=TRUE}, the
returned loadings crossproduct plus the diagonal of the returned specific
variances directly reproduces that covariance.
}
\value{
A list containing:
\describe{
\item{loadings}{a levels-by-factor matrix \eqn{\Lambda}.}
\item{specific}{a named vector of specific variances \eqn{\Psi} (identically
zero for every level in a \code{rrm()} term).}
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
\code{\link{scores_mmes}}, \code{\link{fam}}, \code{\link{rrm}}, \code{\link{vsm}}, \code{\link{mmes}}.
}
