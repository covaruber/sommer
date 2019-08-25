\name{EM}
\alias{EM}
\title{Expectation Maximization Algorithm}
\description{
Univariate version of the expectation maximization (EM) algorithm.
}
\usage{
EM(y,X=NULL,ZETA=NULL,R=NULL,iters=30,draw=TRUE,silent=FALSE, 
   constraint=TRUE, init=NULL, forced=NULL, tolpar = 1e-04, 
   tolparinv = 1e-06)
}
\arguments{
    \item{y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects.}
  \item{ZETA}{an incidence matrix for random effects. This can be for one or more random effects. This NEEDS TO BE PROVIDED AS A LIST STRUCTURE. For example Z=list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3)) makes a 2 level list for 3 random effects. The general idea is that each random effect with or without its variance-covariance structure is a list, i.e. list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix. When moving to more than one random effect we need to make several lists that need to be inside another list. What we call a 2-level list, i.e. list(Z=Z1, K=K1) and list(Z=Z2, K=K2) would need to be put in the form; list(list(Z=Z1, K=K1),list(Z=Z1, K=K1)), which as can be seen, is a list of lists (2-level list). }
  \item{R}{a list of matrices for residuals, i.e. for longitudinal data. if not passed is assumed an identity matrix.}
  \item{draw}{a TRUE/FALSE value indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  \item{iters}{a scalar value indicating how many iterations have to be performed if the EM is performed. There is no rule of tumb for the number of iterations. The default value is 100 iterations or EM steps.}
    \item{constraint}{a TRUE/FALSE value indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way.}
      \item{init}{vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.}
      \item{forced}{a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.}
      \item{tolpar}{tolerance parameter for convergence in the models.}
      \item{tolparinv}{tolerance parameter for matrix inversion in the models.}
}
\details{
This algorithm is based on Searle (1993) and Bernanrdo (2010). This handles models of the form:



y = Xb + Zu + e




b ~ N[b.hat, 0]  ............zero variance because is a fixed term

u ~ N[0, K*sigma(u)]  .......where: K*sigma(u) = G

e ~ N[0, I*sigma(e)]  .......where: I*sigma(e) = R

y ~ N[Xb, var(Zu+e)]   ......where;

var(y) = var(Zu+e) = ZGZ+R = V which is the phenotypic variance

.

The function allows the user to specify the incidence matrices with their respective variance-covariance matrix in a 2 level list structure. For example imagine a mixed model with the following design:

.

fixed = only intercept....................................b ~ N[b.hat, 0]

random = GCA1 + GCA2 + SCA.................u ~ N[0, G]       

.

where G is:

.

|K*sigma(gca1).....................0..........................0.........|
    
|.............0.............S*sigma(gca2).....................0.........| = G
|.............0....................0......................W*sigma(sca)..|

.

The function is based on useing initial values for variance components, i.e.:

.

var(e) <- 100      
var(u1) <- 100      with incidence matrix Z1
var(u2) <- 100      with incidence matrix Z2
var(u3) <- 100      with incidence matrix Z3

.

and estimates the lambda(vx) values in the mixed model equations (MME) developed by Henderson (1975), i.e. consider the 3 random effects stated above, the MME are:

.

|...............X'*R*X...............X'*R*Z1.............X'*R*Z2...................X'*R*Z3 ..............|   |...X'Ry...|

|.............Z1'*R*X.........Z1'*R*Z1+K1*v1.....Z1'*R*Z2..................Z1'*R*Z3.............|   |...Z1'Ry...| 
|.............Z2'*R*X.............Z2'*R*Z1.............Z2'*R*Z2+K2*v2......Z2'*R*Z3.............|   |...Z2'Ry...|

|.............Z3'*R*X.............Z3'*R*Z1.............Z3'*R*Z2.............Z3'*R*Z3+K3*v3......|   |...Z3'Ry...|
.
..............................................................C.inv...................................................................RHS

.

where "*"" is a matrix product, R is the inverse of the var-cov matrix for the errors, Z1, Z2, Z3 are incidence matrices for random effects, X is the incidence matrix for fixed effects, K1,K2, K3 are the var-cov matrices for random effects and v1,v2,v3 are the estimates of variance components.
.
The algorithm can be summarized in the next steps:
.
1) provide initial values for the variance components
2) estimate the coefficient matrix from MME known as "C"
3) solve the mixed equations as theta = RHS * C.inv 
4) obtain new estimates of fixed (b's) and random effects (u's) called theta
5) update values for variance components according to formulas
6) steps are repeated for a number of iterations specified by the user, ideally is enough when no more variations in the estimates is found, in several problems that could take thousands of iterations, whereas in other 10 iterations could be enough.
}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{
\item{$var.com}{ a vector with the values of the variance components estimated}
\item{$V.inv}{ a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{ a vector with BLUPs for random effects}
\item{$Var.u.hat}{ a vector with variances for BLUPs}
\item{$PEV.u.hat}{ a vector with predicted error variance for BLUPs}
\item{$beta.hat}{ a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{ a vector with variances for BLUEs}
\item{$X}{ incidence matrix for fixed effects}
\item{$Z}{ incidence matrix for random effects, if not passed is assumed to be a diagonal matrix}
\item{$K}{ the var-cov matrix for the random effect fitted in Z}
}
}
\references{
Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.
Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

# ###################################################
# ###################################################
# # IMPORT DATA FOR  ESTIMATING 3 VARIANCE COMPONENTS
# ###################################################
# ###################################################
# ## Import phenotypic data on inbred performance
# ## Full data
# data(cornHybrid)
# hybrid2 <- cornHybrid$hybrid # extract cross data
# A <- cornHybrid$K # extract the var-cov K
# ############################################
# ############################################
# ## breeding values with 3 variance components
# ############################################
# ############################################
# y <- hybrid2$Yield
# X1 <- model.matrix(~ Location, data = hybrid2);dim(X1)
# Z1 <- model.matrix(~ GCA1 -1, data = hybrid2);dim(Z1)
# Z2 <- model.matrix(~ GCA2 -1, data = hybrid2);dim(Z2)
# Z3 <- model.matrix(~ SCA -1, data = hybrid2);dim(Z3)
# 
# K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)     
# ## Realized IBS relationships for set of parents 1
# K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)     
# ## Realized IBS relationships for set of parents 2
# S <- kronecker(K1, K2) ; dim(S)   
# ## Realized IBS relationships for cross (as the Kronecker product of K1 and K2)
# rownames(S) <- colnames(S) <- levels(hybrid2$SCA)
# 
# ETA <- list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=S))
# #ans <- EM(y=y, ZETA=ETA, iters=20)

}
\seealso{The core functions of the package \code{\link{mmer}} }