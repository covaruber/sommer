\name{rotate_back_mmes}
\alias{rotate_back_mmes}
\title{Back-transform rotated mmes random effects}
\usage{
rotate_back_mmes(object)
}
\arguments{
  \item{object}{An \code{mmes} model fitted with
  \code{vsm(..., rotation=TRUE)}.}
}
\description{
Back-transforms the selected relationship random effect from the solver
 eigenbasis to the original relationship-level basis.
}
\details{
If \eqn{Gu=U\Lambda U'} and \eqn{a} denotes the solver-basis effect, the
reported effect is \eqn{g=Ua}. Fitted \code{mmes} objects already expose this
original-basis effect in \code{uList}; this function provides the explicit
transformation from the retained engine output.
}
\value{
A matrix of original-basis random effects, with relationship levels in rows
and covariance-product coordinates in columns.
}
\seealso{\code{\link{mmes}}, \code{\link{vsm}}}
