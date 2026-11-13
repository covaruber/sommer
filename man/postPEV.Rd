\name{postPEV}
\alias{postPEV}
\alias{vcov.mmes.glmm.laplace}

\title{
Post-fit prediction error variances and inverse coefficient matrix
}
\description{
Computes prediction error variances (PEVs) or joint fixed/random error covariance after a model has been fitted with \code{mmes}. PEV-only calculations avoid forming the complete coefficient covariance, using selected inversion for Henderson fits or bounded contrast panels for direct and native Laplace fits.
}
\usage{
postPEV(object, mode = 1L)
}

\arguments{
  \item{object}{
A fitted model object of class \code{"mmes"}. Henderson fits require the final sparse coefficient matrix and its scaling factor. Direct fits require retained post-fit designs and covariance inputs. Native Laplace fits use retained fixed/random designs and either a conditional mode Hessian or observation-backend covariance metadata.
}
  \item{mode}{
Integer specifying the post-fit inverse calculation to perform. Use \code{0} to perform no inverse calculation and clear previously computed inverse-derived outputs, \code{1} to compute the prediction error variances without forming the complete inverse, or \code{2} to compute the complete inverse and prediction error variances. The default is \code{1}. Henderson mode 1 uses a Takahashi sparse inverse subset; direct and native Laplace mode 1 use bounded solve panels.
}
}

\details{
The function is intended for post-processing models fitted with \code{mmes}, allowing prediction error variances or the complete inverse to be calculated only when they are needed.

Traditional nonlinear fits dispatch to \code{\link{post_nonlinear}} for
modes 1 and 2. They distinguish conditional-mode covariance from a
first-order empirical Bayes approximate PEV that includes estimated-parameter
mode sensitivity. Mode 1 retains diagonal approximate PEVs; mode 2 retains
the complete random covariance and parameter/mode cross-covariance in
\code{nonlinearInference}, not a linear-model \code{Ci}. Both currently
require budgeted full conditional-covariance workspace. Mode 0 clears all
saved nonlinear uncertainty without performing inference. The linear-model
inverse descriptions below do not apply to nonlinear fits.

For Henderson fits, when \code{mode = 1}, the sparse mixed model coefficient matrix stored in the fitted object is factorized and the Takahashi equations obtain the selected elements of its inverse required for diagonal prediction error variances. Native Laplace fits obtain the needed joint-information inverse diagonal using bounded solve panels. Their global random backend factors the joint matrix; their observation or explicitly enabled random block backend reconstructs operation-local mode solves and a fixed-effect Schur complement. Observation space does not construct the random-effect Hessian. The complete inverse is not formed.

For direct Gaussian and direct PQL fits, post-fit covariance is computed in
observation space from \eqn{V=ZGZ'+R}. The joint error covariance has blocks
\deqn{E_{\beta\beta}=(X'V^{-1}X)^{-1},\quad
E_{\beta u}=-E_{\beta\beta}X'V^{-1}ZG,\quad
E_{uu}=G-GZ'PZG,}
where \eqn{P=V^{-1}-V^{-1}XE_{\beta\beta}X'V^{-1}}.
Termwise prior products avoid constructing the full \eqn{G} or requiring its
inverse, including singular reduced-rank random covariance.
Mode 1 evaluates identity contrasts in panels controlled by
\code{control$output$panelSize}, storing only the diagonal in \code{Ci}.
Mode 2 returns the complete joint error covariance in \code{Ci}, including
fixed/random and cross-term entries, with the public coefficient ordering.
Neither mode constructs a Henderson coefficient matrix. Factors are local
to the call, never serialized pointers in the fit. Workspace is checked
against \code{control$execution$denseMemoryMB}; mode 2 requires quadratic
coefficient storage and can be refused even when mode 1 fits the budget.
PQL uses the final accepted working response covariance and weights; its PEVs
are link-scale working-model approximations. Gaussian and PQL results
condition on fitted covariance parameters.

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

\code{\link{predict.mmes}} does not require any particular \code{mode}: it computes standard errors for arbitrary linear combinations of fixed and random effects without using \code{Ci}, from stored Henderson \code{C} or the direct observation-space covariance operator. \code{postPEV()} is only needed when diagonal PEVs (\code{uPevList}) or the complete coefficient covariance (\code{Ci}) are wanted for their own sake.

For the sampling variance of the BLUP itself, use \code{\link{postVarU}}:
\eqn{\operatorname{Var}(\hat u)=G-\operatorname{Var}(u-\hat u)}.
PEV and BLUP sampling variance are different quantities. For Gaussian fits
neither includes uncertainty from estimating covariance parameters.
\code{postVarU()} leaves
existing PEV outputs unchanged.

For Henderson fits, the final coefficient matrix must be retained with
\code{control = mmes_control(output = list(coefficientMatrix = TRUE))}.
For direct fits, the same option retains sparse prepared designs,
relationship precisions, residual mappings and final residual weights instead;
no \code{C} is created. Older direct fits missing these inputs must be refitted.
Native
Laplace fits retain their required designs and conditional Hessian or
observation-backend inputs; no
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
