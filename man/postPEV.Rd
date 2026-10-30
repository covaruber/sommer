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
A fitted model object of class \code{"mmes"}. Gaussian fits must contain the final sparse mixed model coefficient matrix and its associated scaling factor. Native Laplace fits use the retained fixed/random designs and native conditional mode Hessian instead.
}
  \item{mode}{
Integer specifying the post-fit inverse calculation to perform. Use \code{0} to perform no inverse calculation and clear previously computed inverse-derived outputs, \code{1} to compute the prediction error variances without forming the complete inverse, or \code{2} to compute the complete inverse and prediction error variances. The default is \code{1}. Gaussian mode 1 uses a Takahashi sparse inverse subset; native Laplace mode 1 uses bounded sparse solve panels.
}
}

\details{
The function is intended for post-processing models fitted with \code{mmes}, allowing prediction error variances or the complete inverse to be calculated only when they are needed.

For Gaussian fits, when \code{mode = 1}, the sparse mixed model coefficient matrix stored in the fitted object is factorized and the Takahashi equations obtain the selected elements of its inverse required for diagonal prediction error variances. For native Laplace fits, the joint fixed/random observed-information matrix is reconstructed from the fitted design and conditional mode Hessian, then its needed inverse diagonal is obtained with bounded sparse solve panels. The complete inverse is not formed.

When \code{mode = 2}, the complete inverse is computed. For Laplace fits this is the joint fixed/random conditional observed-information inverse, with fitted covariance and family parameters treated as fixed. Its random-effect block supplies approximate PEVs; covariance-parameter uncertainty is not included. The calculation is guarded by \code{execution$denseMemoryMB}.

Single-family \code{\link{ordinal}} Laplace fits are an exception to this
conditional-on-parameters approximation: cutpoint locations cannot be ignored.
After fitting, the joint marginal-parameter covariance is propagated through
the conditional random modes and added to the conditional-mode covariance.
This includes estimated cutpoint, fixed-effect and covariance-parameter
uncertainty without double-counting fixed effects. Mode 1 obtains only random
diagonal variances with bounded solve panels; mode 2 returns this approximate
joint fixed/random covariance in \code{Ci}. The fixed block matches
\code{vcov(object)}. No fit-time numerical Hessian or PEV is introduced.
Nonconverged, near-boundary and invalid-curvature ordinal fits refuse PEV.
These are local approximations, not statistically calibrated intervals.

Single-family \code{\link{multinomial}} fits use the complete cross-category
observation curvature when reconstructing the joint fixed/random information.
Their PEVs condition on the fitted covariance parameters; no full inverse is
formed in mode 1.

When \code{mode = 0}, no matrix factorization or inverse calculation is performed and previously stored inverse-derived results are cleared.

\code{\link{predict.mmes}} does not require any particular \code{mode}: it computes exact standard errors for arbitrary linear combinations of fixed and random effects directly from the stored coefficient matrix \code{C}, without using \code{Ci} at all (see \code{\link{predict.mmes}} for details). \code{postPEV()} is only needed when the diagonal prediction error variances (\code{uPevList}) or the complete inverse (\code{Ci}) are wanted for their own sake.

For the sampling variance of the BLUP itself, use \code{\link{postVarU}}:
\eqn{\operatorname{Var}(\hat u)=G-\operatorname{Var}(u-\hat u)}.
PEV and BLUP sampling variance are different quantities. For Gaussian fits
neither includes uncertainty from estimating covariance parameters.
\code{postVarU()} leaves
existing PEV outputs unchanged.

For Gaussian fits, the final coefficient matrix must be retained with
\code{control = mmes_control(output = list(coefficientMatrix = TRUE))}. Native
Laplace fits retain their required designs and conditional Hessian; no
coefficient matrix option is required. For either engine,
\code{control$output$inverse = 0} defers PEV work until this function is called.
}
\value{

An object of class \code{"mmes"} containing the original fitted model together with the requested post-fit inverse information.

With \code{mode = 1}, \code{uPevList} contains random-effect prediction error variances. \code{Ci} may contain the selected inverse entries used to obtain those PEVs, but it does not contain the complete inverse.

With \code{mode = 2}, \code{Ci} contains the complete inverse of the mixed model coefficient matrix and \code{uPevList} contains the corresponding prediction error variances for the random effects.

With \code{mode = 0}, inverse-derived outputs are cleared.

The returned object also contains \code{CiComputed} and \code{CiMode}. \code{CiComputed} is TRUE only when the complete inverse was computed; \code{CiMode} records the selected mode.

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
#   control = mmes_control(output = list(inverse = 0))
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
