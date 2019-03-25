\name{vs}
\alias{vs}
\title{variance structure specification}
\description{
  \code{vs} is the main function to build the variance-covariance structure for the random effects to be fitted in the \code{\link{mmer}} solver.
}
\usage{
  vs(..., Gu=NULL, Gt=NULL, Gtc=NULL)
}

\arguments{
  \item{...}{variance structure to be specified following the logic desired in the internal kronecker product. For example, if user wants to define a diagonal variance structure for the random effect 'genotypes'(g) with respect to a random effect 'environments'(e), this is:
  
  \code{var(g) = G.e @ I.g}
  
  being \code{G.e}  a matrix containing the variance covariance components for g (genotypes) in each level of e (environments), \code{I.g} is the covariance among levels of g (genotypes; i.e. relationship matrix), and @ is the kronecker product. This would be specified in the mmer solver as:
  
  \code{random=~vs(ds(e),g)}
  
  One strength of sommer is the ability to specify very complex structures with as many kronecker products as desired. For example:
  
  \code{var(g) = G.e @ G.f @ G.h @ I.g}
  
  is equivalent to
  
  \code{random=~vs(e,f,h,g)}
  
  where different covariance structures can be applied to the levels of \code{e,f,h} (i.e. \code{\link{ds}}, \code{\link{us}}, \code{\link{cs}}, \code{\link{at}} or a combination of these). For more examples please see the vignettes 'sommer.start' available in the package.}
  \item{Gu}{matrix with the known variance-covariance structure for the random effect levels (i.e. relationship matrix among individuals or any other known covariance structure). If NULL, then an identity matrix is assumed.}
  \item{Gt}{matrix with dimensions t x t (t equal to number of traits) with initial values of the variance-covariance components for the random effect specified in the .... argument. If NULL the program will provide the initial values.}
  \item{Gtc}{matrix with dimensions t x t (t equal to number of traits) of constraints for the variance-covariance components for the random effect specified in the ... argument according to the following rules:
  
  \code{0: not to be estimated}
  
  \code{1: estimated and constrained to be positive (i.e. variance component)}
  
  \code{2: estimated and unconstrained (can be negative or positive, i.e. covariance component)}
  
  \code{3: not to be estimated but fixed (value has to be provided in the Gt argument)}
  
  In the multi-response scenario if the user doesn't specify this argument the default is to build an unstructured matrix (using the \code{\link{unsm}}() function). This argument needs to be used wisely since some covariance among responses may not make sense.
  
  }
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
head(DT)

## ============================ ##
## example to use ds() structure (DIAGONAL)
## ============================ ##
ds(DT$Year)
mix <- mmer(Yield~Env, 
            random= ~ vs(ds(Year),Name),
            rcov=~ vs(ds(Year),units),
            data=DT)

## ============================ ##
## example to use at() structure (level-specific)
## ============================ ##
unique(DT$Year)
mix <- mmer(Yield~Env, 
            random= ~ vs(at(Year,c("2011","2012")),Name),
            rcov=~ vs(ds(Year),units),
            data=DT)

## ============================ ##
## example to use us() structure (UNSTRUCTURED)
## ============================ ##
us(DT$Year)
mix <- mmer(Yield~Env, 
            random= ~ vs(us(Year),Name),
            rcov=~ vs(ds(Year),units),
            data=DT)

## ============================ ##
## example to use cs() structure (CUSTOMIZED)
## ============================ ##
unique(DT$Year)
mm <- matrix(1,3,3); mm[1,3] <- mm[3,1] <- 0;mm #don't estimate cov 2011-2013
mix <- mmer(Yield~Env, 
            random= ~ vs(cs(Year,mm),Name),
            rcov=~ vs(ds(Year),units),
            data=DT)

## ============================ ##
## example to use overlay() + vs() structure
## ============================ ##
data("DT_halfdiallel")
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
A <- diag(7); colnames(A) <- rownames(A) <- 1:7;A # if you want to provide a covariance matrix
#### model using overlay
modh <- mmer(sugar~1, 
             random=~vs(overlay(femalef,malef), Gu=A) 
                     + genof,
              data=DT)

## ============================ ##
## example to use spl2D() + vs() structure
## ============================ ##            
# ### mimic two fields
# aa <- DT; bb <- DT
# aa$FIELD <- "A";bb$FIELD <- "B"
# set.seed(1234)
# aa$Yield <- aa$Yield + rnorm(length(aa$Yield),0,4)
# DT2 <- rbind(aa,bb)
# head(DT2)
# 
# mix <- mmer(Yield~1,
#             random=~vs(ds(FIELD),id, Gu=A) + 
#               vs(ds(FIELD),Rowf) +
#               vs(ds(FIELD),Colf) +
#                 vs(ds(FIELD),spl2D(Row,Col)),
#             rcov=~vs(ds(FIELD),units),
#             data=DT2)

}
\seealso{The core function of the package: \code{\link{mmer}}}