\name{postPEV}
\alias{postPEV}

\title{
Post-fit prediction error variances and inverse coefficient matrix
}
\description{
Computes prediction error variances (PEVs) or the complete inverse of the mixed model coefficient matrix after a model has been fitted with \code{mmes}. For PEV-only calculations, the function uses the Takahashi sparse inverse subset approach, avoiding computation of the complete inverse of the coefficient matrix.
}
\usage{
postPEV(object, mode = 1L)
}

\arguments{
  \item{object}{
A fitted model object of class \code{"mmes"}. The object must contain the final sparse mixed model coefficient matrix and its associated scaling factor.
}
  \item{mode}{
Integer specifying the post-fit inverse calculation to perform. Use \code{0} to perform no inverse calculation and clear previously computed inverse-derived outputs, \code{1} to compute the prediction error variances using the Takahashi sparse inverse subset approach without forming the complete inverse, or \code{2} to compute the complete inverse of the coefficient matrix and the prediction error variances. The default is \code{1}.
}
}

\details{
The function is intended for post-processing models fitted with \code{mmes}, allowing prediction error variances or the complete inverse of the mixed model coefficient matrix to be calculated only when they are needed.

When \code{mode = 1}, the sparse mixed model coefficient matrix stored in the fitted object is factorized and the Takahashi equations are used to obtain the selected elements of its inverse required for the diagonal prediction error variances. The complete inverse is not formed.

When \code{mode = 2}, the complete inverse of the mixed model coefficient matrix is computed. The prediction error variances are then obtained from the corresponding diagonal elements of this inverse.

When \code{mode = 0}, no matrix factorization or inverse calculation is performed and previously stored inverse-derived results are cleared.

\code{\link{predict.mmes}} does not require any particular \code{mode}: it computes exact standard errors for arbitrary linear combinations of fixed and random effects directly from the stored coefficient matrix \code{C}, without using \code{Ci} at all (see \code{\link{predict.mmes}} for details). \code{postPEV()} is only needed when the diagonal prediction error variances (\code{uPevList}) or the complete inverse (\code{Ci}) are wanted for their own sake.

For this post-fit calculation to be available, the \code{mmes} model must have been fitted with a version of the Henderson mixed model solver that stores the final sparse coefficient matrix and its scaling factor in the fitted model object. Models can therefore be initially fitted with \code{computeCi = 0} and the PEVs or complete inverse calculated afterwards only if required.
}
\value{

An object of class \code{"mmes"} containing the original fitted model together with the requested post-fit inverse information.

With \code{mode = 1}, \code{uPevList} contains the prediction error variances for the random effects while \code{Ci} remains empty.

With \code{mode = 2}, \code{Ci} contains the complete inverse of the mixed model coefficient matrix and \code{uPevList} contains the corresponding prediction error variances for the random effects.

With \code{mode = 0}, inverse-derived outputs are cleared.

The returned object also contains \code{CiComputed} and \code{CiMode}, indicating whether the complete inverse was computed and which calculation mode was used, respectively.

}
\references{

Takahashi, K., Fagan, J., and Chin, M. S. (1973). Formation of a sparse bus impedance matrix and its application to short circuit study. Proceedings of the 8th PICA Conference, Minneapolis, Minnesota.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}
\examples{
####=========================================####
#### Fit a mixed model without computing PEVs
#### or the complete inverse during estimation
####=========================================####

# mod <- mmes(
#   fixed = Yield ~ 1,
#   random = ~ vsm(ism(Genotype)),
#   rcov = ~ units,
#   data = DT,
#   computeCi = 0
# )

####=========================================####
#### Compute PEVs afterwards using the
#### Takahashi sparse inverse subset
####=========================================####

# mod <- postPEV(mod, mode = 1)
# mod$uPevList

####=========================================####
#### Compute the complete inverse afterwards
####=========================================####

# mod <- postPEV(mod, mode = 2)
# mod$Ci
# mod$uPevList
}
