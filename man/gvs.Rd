\name{gvs}
\alias{gvs}
\title{general variance structure specification}
\description{
  \code{gvs} function to build general variance-covariance structures for combination of random effects to be fitted in the \code{\link{mmer}} solver. 
}
\usage{
  gvs(..., Gu=NULL, Guc=NULL, Gti=NULL, Gtc=NULL, form=NULL)
}

\arguments{
  \item{...}{names of the random effects (variables in the dataset) to be used to create a general variance structure. For example, for 2 random effects (variables); mom and progeny, a model specified as:
  
  \code{gvs(mom, progeny)}
  
  will create a variance structure of the form:
  
\code{| sigma2.m   sigma.pm |}

\code{| sigma.pm   sigma2.p |}
  
  where not only variance components for each random effect will be estimated but also the covariance component between the 2 random effects is estimated. The user can also provide a numeric vector or matrix to be considered the design matrix for the ith random effect. More than two random effects can be provided.
  
  }

  \item{Gu}{list of matrices with the known variance-covariance values for the levels of the different random effects provided in the "..." argument (i.e. relationship matrix among individuals or any other known covariance matrix). If NULL, then an identity matrix is assumed. For example, a model with 2 random effects with covariance structure should be provided as:
  
  \code{gvs(mom, progeny, Gu=list(Am,Ap))}
  
  where \code{Am} and \code{Ap} are the relationship matrices for the random effects for mom and progeny respectively.
  
  }
  
  \item{Guc}{matrix with the constraints for the u random effects. This is used to specify which variance and covariance parameters between the 1 to 1 combinations of random effects should be estimated. For example, for 2 random effects the expected variance-covariance matrix expected to be estimated (when the default Guc=NULL) is and unstructured model:
  
  \code{| sigma2.m   sigma.pm |}

  \code{| sigma.pm   sigma2.p |}

but the user can constrain which parameters should be estimated. Providing:

\code{Guc=diag(2)}     would fit:

  \code{| sigma2.m  ...0... |}

  \code{| ...0...  sigma2.p |}

  }

  \item{Gti}{matrix with dimensions t x t (t equal to number of traits) with initial values of the variance-covariance components for the random effect specified in the .... argument. If the value is NULL the program will provide the initial values.}

  \item{Gtc}{matrix with dimensions t x t (t equal to number of traits) of constraints for the variance-covariance components for the random effect specified in the ... argument according to the following rules:
  
  \code{0: not to be estimated}
  
  \code{1: estimated and constrained to be positive (i.e. variance component)}
  
  \code{2: estimated and unconstrained (can be negative or positive, i.e. covariance component)}
  
  \code{3: not to be estimated but fixed (value has to be provided in the Gti argument)}
  
  In the multi-response scenario if the user doesn't specify this argument the default is to build an unstructured matrix (using the \code{\link{unsm}}() function). This argument needs to be used wisely since some covariance among responses may not make sense. Useful functions to specify constraints are; \code{\link{diag}}(), \code{\link{unsm}}(), \code{\link{uncm}}(), \code{\link{fixm}}().
  }
  
  \item{form}{an additional structure to specify a kronecker product on top of the general covariance structure defined in the ... argument.}
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

data(DT_ige)
DT <- DT_ige
Af <- A_ige
An <- A_ige
### Direct genetic effects model
# modDGE <- mmer(trait ~ block,
#                random = ~ focal,
#                rcov = ~ units, 
#                data = DT)
# summary(modDGE)$varcomp
# 
### Indirect genetic effects model without covariance between DGE and IGE
# modDGE <- mmer(trait ~ block,
#                random = ~focal + neighbour,
#                rcov = ~ units, 
#                data = DT)
# summary(modDGE)$varcomp
# 
### Indirect genetic effects model with covariance between DGE and IGE
# modIGE <- mmer(trait ~ block,
#               random = ~ gvs(focal, neighbour),
#               rcov = ~ units, iters=4, 
#               data = DT)
# summary(modIGE)$varcomp
#
### Indirect genetic effects model with covariance between DGE and IGE using relatioship matrices
# modIGEb <- mmer(trait ~ block,
#                random = ~ gvs(focal, neighbour, Gu=list(Af,An)),
#                rcov = ~ units, 
#                data = DT)
# summary(modIGEb)$varcomp
 

}
\seealso{The core function of the package: \code{\link{mmer}}}