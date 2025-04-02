\name{vs}
\alias{vs}
\title{variance structure specification}
\description{
  \code{vs} DEPRECATED NOW. was the main function to build the variance-covariance structure for the random effects to be fitted in the \code{\link{mmer}} solver. 
}
\usage{
  vs(..., Gu=NULL, Gti=NULL, Gtc=NULL, reorderGu=TRUE, buildGu=TRUE)
}

\arguments{
  \item{...}{variance structure to be specified following the logic desired in the internal kronecker product. For example, if user wants to define a diagonal variance structure for the random effect 'genotypes'(g) with respect to a random effect 'environments'(e), this is:
  
  \code{var(g) = G.e @ I.g}
  
  being \code{G.e}  a matrix containing the variance covariance components for g (genotypes) in each level of e (environments), \code{I.g} is the covariance among levels of g (genotypes; i.e. relationship matrix), and @ is the kronecker product. This would be specified in the mmer solver as:
  
  \code{random=~vs(dsr(e),g)}
  
  One strength of sommer is the ability to specify very complex structures with as many kronecker products as desired. For example:
  
  \code{var(g) = G.e @ G.f @ G.h @ I.g}
  
  is equivalent to
  
  \code{random=~vs(e,f,h,g)}
  
  where different covariance structures can be applied to the levels of \code{e,f,h} or a combination of these). For more examples please see the vignettes 'sommer.start' available in the package.}

  \item{Gu}{matrix with the known variance-covariance values for the levels of the u.th random effect (i.e. relationship matrix among individuals or any other known covariance matrix). If \code{NULL}, then an identity matrix is assumed. The Gu matrix can have more levels than the ones present in the random effect linked to it but not the other way around. Otherwise, an error message of missing level in Gu will be returned.}

  \item{Gti}{matrix with dimensions t x t (t equal to number of traits) with initial values of the variance-covariance components for the random effect specified in the .... argument. If \code{NULL} the program will provide the initial values. The values need to be scaled, see Details section.}

  \item{Gtc}{matrix with dimensions t x t (t equal to number of traits) of constraints for the variance-covariance components for the random effect specified in the ... argument according to the following rules:
  
  \code{0: not to be estimated}
  
  \code{1: estimated and constrained to be positive (i.e. variance component)}
  
  \code{2: estimated and unconstrained (can be negative or positive, i.e. covariance component)}
  
  \code{3: not to be estimated but fixed (value has to be provided in the Gti argument)}
  
  In the multi-response scenario if the user doesn't specify this argument the default is to build an unstructured matrix (using the \code{\link{unsm}}() function). This argument needs to be used wisely since some covariance among responses may not make sense. Useful functions to specify constraints are; \code{\link{diag}}(), \code{\link{unsm}}(), \code{\link{fixm}}().
  }
  
  \item{reorderGu}{a \code{TRUE/FALSE} statement if the Gu matrix should be reordered based on the names of the design matrix of the random effect or passed with the custom order of the user. This may be important when fitting covariance components in a customized fashion. Only for advanced users.}
  
    \item{buildGu}{a \code{TRUE/FALSE} statement to indicate if the Gu matrix should be built in R when the value for the argument \code{Gu=NULL}. Repeat, only when when the value for the argument \code{Gu} is equal to \code{NULL}. In some cases when the incidence matrix is wide (e.g. rrBLUP models) the covariance structure is a huge p x p matrix that can be avoided when performing matrix operations. By setting this argument to \code{FALSE} it allows to skip forming this covariance matrix.}
    
}
\details{
 When providing initial values in the \code{Gti} argument the user has to provide scaled variance component values. The user can provide values from a previous model by accessing the \code{sigma_scaled} output from an \code{mmer} model or if an specific value is desired the user can obtain the scaled value as:
 
 \code{m = x/var(y)}
 
 where \code{x} is the desired initial value and \code{y} is the response variable. You can find an example in the \code{\link{DT_cpdata}} dataset.
 
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
A <- A_example

## ============================ ##
## example to without structure
## ============================ ##

mix <- mmer(Yield~Env, 
            random= ~ vs(Name),
            rcov=~ vs(units),
            data=DT)

}
\seealso{The core function of the package: \code{\link{mmer}}}