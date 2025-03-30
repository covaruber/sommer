\name{vsr}
\alias{vsr}
\title{variance structure specification}
\description{
  \code{vsr} is the main function to build the variance-covariance structure for the random effects to be fitted in the \code{\link{mmer}} solver. 
}
\usage{
  vsr(..., Gu=NULL, buildGu=TRUE, meN=1, meTheta=NULL, meThetaC=NULL, 
      sp=FALSE, isFixed=FALSE, verbose=TRUE)
}

\arguments{
  \item{...}{variance structure to be specified following the logic desired in the internal kronecker product. For example, if user wants to define a diagonal variance structure for the random effect 'genotypes'(g) with respect to a random effect 'environments'(e), this is:

  \code{var(g) = G.e @ I.g}

  being \code{G.e}  a matrix containing the variance covariance components for g (genotypes) in each level of e (environments), \code{I.g} is the covariance among levels of g (genotypes; i.e. relationship matrix), and @ is the kronecker product. This would be specified in the mmer solver as:

  \code{random=~vsr(dsr(e),g)}

  One strength of sommer is the ability to specify very complex structures with as many kronecker products as desired. For example:

  \code{var(g) = G.e @ G.f @ G.h @ I.g}

  is equivalent to

  \code{random=~vsr(e,f,h,g)}

  where different covariance structures can be applied to the levels of \code{e,f,h} (i.e. \code{\link{dsr}}, \code{\link{usc}}, \code{\link{csc}}, \code{\link{atr}} or a combination of these). For more examples please see the vignettes 'sommer.start' available in the package.}

  \item{Gu}{matrix with either the inverse (in case of mme=TRUE) or raw variance-covariance (mme=FALSE)
  values for the levels of the u.th random effect (e.g., the inverse of a 
  relationship matrix among individuals or any other known inverse covariance 
  matrix). If \code{NULL}, then an identity matrix is assumed.
  The Gu matrix can have more levels than the ones present in the random effect
  linked to it but not the other way around. Otherwise, an error message of missing
  level in Gu will be returned.}

    \item{buildGu}{a \code{TRUE/FALSE} statement to indicate if the Gu matrix should be built in R when the value for the argument \code{Gu=NULL}. Repeat, only when when the value for the argument \code{Gu} is equal to \code{NULL}. In some cases when the incidence matrix is wide (e.g. rrBLUP models) the covariance structure is a huge p x p matrix that can be avoided when performing matrix operations. By setting this argument to \code{FALSE} it allows to skip forming this covariance matrix.}

   \item{meN}{number of main effects in the variance structure. Is always counted from last to first.}

   \item{meTheta}{variance covariance matrix between the main effects desired. Just to be modified if the number of main effects is greater of 1 (e.g., indirect genetic effects).}

   \item{meThetaC}{constraints for the variance covariance matrix between the main effects desired. Just to be modified if the number of main effects is greater of 1 (e.g., indirect genetic effects).}

   \item{sp}{a TRUE/FALSE statement to indicate if the VC from this structure should be multiplied by the scale parameter added in the mmer function through the addScaleParam argument in the \code{mmer} function .}
   
   \item{isFixed}{a TRUE/FALSE statement to indicate if the vsr function is being used in the fixed part of the model. When TRUE, the function only returns the model matrix to avoid any error messages associated to returning all elements for a random effect. FALSE is the default since it is assumed to be used for a variance structure in a random effect.}
   
   \item{verbose}{a TRUE/FALSE statement to indicate if messages should be printed when special situations accur. For example, adding unphenotyped individuals to the incidence matrices when present in the relationship  matrices.}

}
\details{

...

}
\value{
  \describe{
    \item{$res}{a list with all neccesary elements (incidence matrices, known var-cov structures, unknown covariance structures to be estimated and constraints) to be used in the mmer solver.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

  Covarrubias-Pazaran G (2018) Software update: Moving the R package sommer to
  multivariate mixed models for genome-assisted prediction. doi:
  https://doi.org/10.1101/354639
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_example)
DT <- DT_example
DT=DT[with(DT, order(Env)), ]
A <- A_example

x <- as.character(unique(DT$Name))
DT <- droplevels(DT[which(!is.na(match(DT$Name, x[1:5]))),])
## ============================ ##
## example without structure
## ============================ ##
isr(DT$Name)
mix <- mmer(Yield~Env,
            random= ~ vsr(isr(Name)),
            rcov=~ units,
            nIters=10,
            data=DT)

## ============================ ##
## example to without structure but
## using covariance among levels in the
## random effect Name
## ============================ ##
Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
mix <- mmer(Yield~Env,
            random= ~ vsr(isr(Name), Gu=Ai),
            rcov=~ units,
            nIters=10,
            data=DT)
summary(mix)$varcomp
## ============================ ##
## example to use dsr() structure (DIAGONAL)
## ============================ ##
dsr(DT$Year)
mix <- mmer(Yield~Env,
            random= ~ vsr(dsr(Year),isr(Name)),
            rcov=~ vsr(dsr(Year),isr(units)),
            nIters=10,
            data=DT)
summary(mix)$varcomp
## ============================ ##
## example to use atc() structure (level-specific)
## ============================ ##
# unique(DT$Year)
# mix <- mmer(Yield~Env,
#             random= ~ vsr(atc(Year,c("2011","2012")),isr(Name)),
#             rcov=~ vsr(dsr(Year),isr(units)),
#             data=DT)
## ============================ ##
## example to use usc() structure (UNSTRUCTURED)
## ============================ ##
usc(DT$Year)
mix <- mmer(Yield~Env,
            random= ~ vsr(usc(Year),isr(Name)),
            rcov=~ vsr(dsr(Year),isr(units)),
            nIters = 10,
            data=DT)
## ============================ ##
## example using structure in fixed effect
## (notice the isFixed argument)
## ============================ ##
mix <- mmer(Yield~ vsr(atc(Env,"CA.2011"), isr(Name), isFixed = TRUE),
            rcov=~ units,
            nIters=10,
            data=DT)

}
\seealso{The core function of the package: \code{\link{mmer}}}
