\name{meffects_mmes}
\alias{meffects_mmes}
\title{Marker effects recovered from a fitted relationship-based term}
\description{
Back-transforms the BLUPs of a random term fitted with a genomic relationship
matrix (GBLUP) into BLUPs of marker effects, optionally with standard errors
and back-solved GWAS-type tests. This is the post-fit equivalent of ASReml's
\code{mef} argument and does not require refitting the model.
}
\usage{
meffects_mmes(object, term, M, method = c("A.mat", "custom"),
              min.MAF = 0, Z = NULL, scale = NULL, blend = 0, Gu = NULL,
              se = FALSE, chunk = 1000)
}
\arguments{
  \item{object}{A fitted \code{\link{mmes}} model.}
  \item{term}{Name (as in \code{names(object$uList)}) or position of the random term fitted with \code{Gu}.}
  \item{M}{Marker matrix (individuals by markers, coded -1/0/1 as for
  \code{\link{A.mat}}) with rownames matching the levels of the term. Use the
  same matrix (same rows) that was used to build the relationship matrix,
  because allele frequencies are computed from it. Rows that are not levels of
  the term are ignored.}
  \item{method}{\code{"A.mat"} reproduces exactly how \code{\link{A.mat}}
  builds G from \code{M} (imputation, \code{min.MAF} and monomorphic filters,
  centring by \eqn{2p} and scaling by \eqn{2\sum p_jq_j}). \code{"custom"} uses
  \code{Z} and \code{scale}.}
  \item{min.MAF}{The \code{min.MAF} used in \code{A.mat}.}
  \item{Z}{For \code{method="custom"}: transformed marker matrix with rownames
  such that \eqn{G = scale \cdot ZZ^\top}.}
  \item{scale}{For \code{method="custom"}: the scalar \eqn{s} in \eqn{G = sZZ^\top}.}
  \item{blend}{\eqn{\lambda} if the fitted matrix was \eqn{(1-\lambda)G+\lambda I}.
  Small bending constants added only for invertibility can be left at 0.}
  \item{Gu}{Optional relationship matrix (or inverse, with \code{attr(Gu,"inverse")=TRUE}).
  By default the inverse used in the fit is rebuilt from the model inputs.}
  \item{se}{If \code{TRUE}, add standard errors, z statistics and p-values.}
  \item{chunk}{Number of markers processed per block when \code{se=TRUE}.}
}
\details{
If the fitted relationship matrix is \eqn{G=(1-\lambda)sZZ^\top+\lambda I}, the
BLUP of the marker effects is
\deqn{\hat a=(1-\lambda)s\,Z^\top G^{-1}\hat u,}
using the same \eqn{G^{-1}} as the fit, so \eqn{Z\hat a=\hat u-\lambda G^{-1}\hat u}
(exactly \eqn{\hat u} when \eqn{\lambda=0}). For structured terms such as
\code{vsm(dsm(env), ism(id), Gu=Gi)} effects are returned per coordinate
(environment, trait, ...).

With \code{se=TRUE}, following Gualdron-Duarte et al. (2014),
\deqn{\mathrm{Var}(\hat a)=H^\top(\sigma^2_c G - C^{uu})H,\qquad H=(1-\lambda)sG^{-1}Z,}
where \eqn{\sigma^2_c} is the variance of the coordinate and \eqn{C^{uu}} the
prediction error covariance of the term; \eqn{z=\hat a/\mathrm{SE}} gives
back-solved GWAS tests. Only diagonal elements are formed: the \eqn{C^{uu}}
part uses solves with the stored coefficient matrix in blocks of \code{chunk}
markers. The cost grows as \eqn{O(mn^2)} for m markers and n individuals.

The rebuild of \eqn{G^{-1}} requires the data used in the original call; supply
\code{Gu} otherwise. Markers removed by the filters are returned as \code{NA}.
}
\value{A data frame with columns \code{marker}, \code{coordinate} (structured
terms only), \code{effect} and, if \code{se=TRUE}, \code{se}, \code{z} and
\code{p.value}. Attributes \code{scale} and \code{markersKept}.}
\references{
Gualdron-Duarte JL, Cantet RJC, Bates RO, Ernst CW, Raney NE, Steibel JP (2014).
Rapid screening for phenotype-genotype associations by linear transformations
of genomic evaluations. BMC Bioinformatics 15: 246.
}
\seealso{\code{\link{A.mat}}, \code{\link{mmes}}, \code{\link{vsm}}}
\examples{
set.seed(7)
M <- matrix(sample(c(-1, 0, 1), 40 * 80, TRUE), 40)
rownames(M) <- paste0("g", 1:40); colnames(M) <- paste0("snp", 1:80)
G <- 0.95 * A.mat(M) + 0.05 * diag(40)
Gi <- solve(G); attr(Gi, "inverse") <- TRUE
d <- data.frame(id = factor(rep(rownames(M), each = 2), levels = rownames(M)))
d$y <- rnorm(80)
fit <- mmes(y ~ 1, random = ~ vsm(ism(id), Gu = Gi), data = d, verbose = FALSE)
head(meffects_mmes(fit, 1, M, blend = 0.05, se = TRUE))
}
\keyword{models}
