\name{vsm}
\alias{vsm}
\title{Variance-structure model with arbitrary Kronecker covariance factors}
\description{
\code{vsm} constructs random-effect or residual covariance structures for
\code{mmes}.  Covariance-shaping terms are represented by CovarianceFactor v2
descriptors and may be combined in an arbitrary-depth Kronecker product.  A
single product-level variance, \eqn{\sigma^2}, is owned by \code{vsm}, which
avoids scale confounding among the individual covariance factors.
}
\usage{
vsm(..., Gu = NULL, sigma2 = 0.15, fixedSigma2 = FALSE,
    isFixed = FALSE, verbose = TRUE)
}
\arguments{
  \item{...}{One or more covariance-constructor terms.  The last term supplies
the main-effect incidence matrix.  All preceding terms are covariance-shaping
factors.  Examples include \code{dsm(environment)}, \code{ar1m(row)}, and
\code{ism(genotype)}.}
  \item{Gu}{Optional known precision matrix for the levels of the final
main-effect term.  For the Henderson solver, \code{Gu} must have row and column
names matching the main-effect levels and \code{attr(Gu, "inverse")} must be
defined.  If omitted, an identity precision matrix is used.}
  \item{sigma2}{Positive starting value for the single overall covariance
scale of this \code{vsm} term.}
  \item{fixedSigma2}{Logical indicating whether the product-level
\eqn{\sigma^2} is fixed at its starting value.}
  \item{isFixed}{Logical retained for compatibility.  If \code{TRUE}, return
the combined design matrix instead of the structured \code{vsm} object.}
  \item{verbose}{Logical controlling messages, including notification when
levels present in \code{Gu} are appended to the model matrix.}
}
\details{
Let \eqn{K_1,\ldots,K_m} be the dimensionless covariance shapes supplied by
the covariance constructors.  The covariance represented by one \code{vsm}
term is
\deqn{\Sigma = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes K_m).}
For a random effect with known relationship precision \code{Gu}, the factor
product determines the covariance among the crossed covariance coordinates and
\code{Gu} determines the relationship structure among the levels of the final
main effect.

The final argument in \code{...} is the main-effect incidence term.  Every term
before it is a covariance-shaping factor.  Thus
\code{vsm(dsm(location), ar1m(row), ism(genotype), Gu=Ainv)} represents a
location-by-row covariance shape crossed with genotype.  Kronecker ordering is
left to right: earlier factors are the slow index and later factors are the fast
index.

All optimizer coordinates are stored on unconstrained working scales.  Each
CovarianceFactor descriptor also carries its evaluator, derivative strategy,
reporting transformation, and trust-region caps.  Built-in high-use structures
may use native C++ covariance primitives, while newer or user-defined
structures can use generic R callbacks without requiring changes to
\code{ai_mme_sp2}.

For residual covariance models, the covariance-shaping factors must define
exactly one product coordinate for each observation.  \code{mmes} combines
that local coordinate with an independently constructed residual block index,
so a large observation-level identity design is not materialized.
}
\value{
A list containing the random-effect design blocks in \code{Z}, the precision
matrix in \code{Gu}, the flattened Kronecker descriptor in \code{covStruct},
the product design, and, for residual models, the local covariance-coordinate
index.  The \code{covStruct} object contains descriptor version 2
CovarianceFactor objects.
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)


str(with(DT, vsm(dsm(Env),ism(Name))))
str(with(DT, vsm(corm(Env),ism(Name))))
str(with(DT, vsm(toeplitzm(Env),ism(Name))))

}
\seealso{
\code{\link{mmes}}, \code{\link{ism}}, \code{\link{dsm}}, \code{\link{atm}}, \code{\link{usm}},
\code{\link{csm}}, \code{\link{ar1m}}, \code{\link{ar2m}}, \code{\link{ar3m}},
\code{\link{mam}}, \code{\link{corm}}, \code{\link{corvm}}, \code{\link{corhm}},
\code{\link{corgm}}, \code{\link{fam}}, \code{\link{antem}}, \code{\link{rrcm}},
\code{\link{maternm}}, \code{\link{toeplitzm}}, \code{\link{sar}}, \code{\link{car}},
and \code{\link{ownm}}.
}
