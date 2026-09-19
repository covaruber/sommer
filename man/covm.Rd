\name{covm}
\alias{covm}

\title{Covariance Between Two Random Effects}

\description{
  \code{covm} combines two random-effect structures created with
  \code{\link{vsm}} into a single correlated random-effect structure.
  It fits an unstructured \(2 \times 2\) covariance matrix between the
  two random effects while allowing them to have different incidence
  matrices. The two effects must act on the same coefficient space and
  use the same relationship or precision matrix. The resulting structure
  can be fitted with the \code{\link{mmes}} solver.
}

\usage{
covm(ran1, ran2, thetaC = NULL, theta = NULL,
     fixed = NULL, fixedSigma2 = FALSE,
     labels = c("ran1", "ran2"), tol = 1e-10)
}

\arguments{

  \item{ran1}{
    A random-effect structure returned by \code{\link{vsm}} for the
    first random effect. Currently, \code{covm} requires a simple
    \code{vsm} structure with one covariance-product coordinate, for
    example \code{vsm(ism(focal))} or
    \code{vsm(ism(focal), Gu = Ai)}.
  }

  \item{ran2}{
    A random-effect structure returned by \code{\link{vsm}} for the
    second random effect. It must have the same coefficient dimension
    and coefficient levels as \code{ran1}, although its incidence
    matrix may differ.
  }

  \item{thetaC}{
    Deprecated and not supported by the current CovarianceFactor
    parameterization. The former cell-wise constraint matrix cannot in
    general be translated exactly into the normalized-Cholesky
    parameterization used by \code{covm}. Use \code{fixedSigma2} and
    \code{fixed} instead.
  }

  \item{theta}{
    An optional symmetric \(2 \times 2\) positive-definite matrix
    containing initial values for the variance-covariance matrix of the
    two random effects. Values are supplied directly on the natural
    variance-covariance scale.

    The diagonal elements define the initial variances of the two random
    effects and the off-diagonal element defines their initial covariance.

    If \code{theta = NULL}, the default starting matrix is

    \preformatted{
    diag(2) * 0.15 + matrix(0.015, 2, 2)
    }

    giving initial variances of 0.165 and an initial covariance of 0.015.
  }

  \item{fixed}{
    An optional logical vector of length two indicating whether the two
    normalized-Cholesky coordinates describing the covariance structure
    should be fixed at their starting values. The coordinates correspond
    to the off-diagonal Cholesky element and the second diagonal Cholesky
    element, respectively. The default is \code{c(FALSE, FALSE)}, so both
    are estimated.
  }

  \item{fixedSigma2}{
    Logical value indicating whether the variance scale associated with
    the first random effect should be fixed at its starting value.
    The default is \code{FALSE}.
  }

  \item{labels}{
    Character vector of length two giving labels for the two random
    effects in the covariance descriptor. The default is
    \code{c("ran1", "ran2")}. The labels must be different and non-empty.
  }

  \item{tol}{
    Numerical tolerance used when checking positive definiteness and
    whether the relationship/precision matrices supplied by the two
    random effects are equal. The default is \code{1e-10}.
  }
}

\details{

  \code{covm} is intended for models in which two distinct random effects
  are expected to be correlated. A common example is the joint modeling
  of direct and indirect genetic effects.

  If the coefficient vector for the first random effect is \(u_1\) and
  that for the second is \(u_2\), \code{covm} represents their covariance
  structure as

  \deqn{
    Var([u_1', u_2']') = \sigma^2 K_{\mathrm{effect}} \otimes A,
  }

  where \(A\) represents the covariance structure associated with the
  coefficient levels, and \(K_{\mathrm{effect}}\) is a \(2 \times 2\)
  unstructured covariance shape describing the relationship between the
  two random effects.

  Internally, \code{covm} uses the same CovarianceFactor interface as
  \code{\link{vsm}}. The \(2 \times 2\) effect covariance is represented
  using a normalized-Cholesky parameterization. This guarantees a
  positive-definite covariance matrix during optimization.

  The variance of the first random effect is represented by the overall
  scale parameter \code{sigma2}. The remaining variance ratio and
  covariance are represented through the normalized-Cholesky factor.

  The two random effects may have different incidence matrices. However,
  they must refer to identical coefficient levels in the same order and
  must use the same relationship/precision matrix.

  When a relationship matrix is supplied through \code{Gu}, the
  Henderson solver requires it to be an inverse/precision matrix with
  \code{attr(Gu, "inverse") = TRUE}. For example:

  \preformatted{
  attr(Ai, "inverse") <- TRUE

  covm(
    vsm(ism(focal), Gu = Ai),
    vsm(ism(neighbour), Gu = Ai)
  )
  }

  The current implementation combines simple \code{vsm} random effects
  only. Each input must contain one covariance-product coordinate. More
  complex covariance-factor products should therefore not currently be
  supplied independently inside \code{ran1} and \code{ran2}.
}

\value{

  A list describing a correlated random-effect structure compatible with
  the CovarianceFactor-v2 interface used by \code{\link{vsm}} and
  \code{\link{mmes}}. The main elements are:

  \describe{

    \item{\code{Z}}{
      A list containing the two incidence matrices, one for each random
      effect.
    }

    \item{\code{Gu}}{
      The common sparse relationship precision matrix associated with
      the coefficient levels. It is stored with
      \code{attr(Gu, "inverse") = TRUE} for use by the Henderson solver.
    }

    \item{\code{covStruct}}{
      A CovarianceFactor-v2 descriptor containing the overall variance
      scale and the normalized-Cholesky representation of the
      \(2 \times 2\) unstructured covariance between the two random
      effects.
    }

    \item{\code{residualLocalIndex}}{
      \code{NULL}, since \code{covm} describes a random-effect structure
      rather than a residual covariance structure.
    }

    \item{\code{productDesign}}{
      \code{NULL} for the current implementation.
    }

    \item{\code{partitionsR}}{
      \code{NULL} for the current implementation.
    }

    \item{\code{covm}}{
      Logical value equal to \code{TRUE}, identifying the object as one
      constructed by \code{covm}.
    }

    \item{\code{covm_labels}}{
      The two labels used for the correlated random effects.
    }
  }
}

\references{

  Covarrubias-Pazaran G (2016).
  Genome assisted prediction of quantitative traits using the R package
  sommer. PLoS ONE 11(6).
  \doi{10.1371/journal.pone.0156744}

  Bijma, P. (2014).
  The quantitative genetics of indirect genetic effects: a selective
  review of modelling issues.
  Heredity, 112(1), 61--69.
}

\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_ige, package = "enhancer")
DT <- DT_ige

## Correlate two random effects with identity relationship matrices
covRes <- with(
  DT,
  covm(
    vsm(ism(focal)),
    vsm(ism(neighbour))
  )
)

str(covRes)


## Custom initial variance-covariance matrix
covRes2 <- with(
  DT,
  covm(
    vsm(ism(focal)),
    vsm(ism(neighbour)),
    theta = matrix(
      c(0.5, 0.1,
        0.1, 0.3),
      2, 2
    )
  )
)


## Relationship precision matrix
## Ai must have matching row and column names and be marked as inverse.
##
## attr(Ai, "inverse") <- TRUE
##
## covRes3 <- with(
##   DT,
##   covm(
##     vsm(ism(focal), Gu = Ai),
##     vsm(ism(neighbour), Gu = Ai)
##   )
## )

## See the DT_ige help page for a complete model-fitting example.
}

\seealso{
  \code{\link{vsm}} for constructing random-effect covariance structures
  and \code{\link{mmes}} for fitting mixed models using the Henderson
  solver.
}