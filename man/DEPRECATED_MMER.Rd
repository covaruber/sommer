\name{mmer}
\alias{mmer}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations for \strong{r} records

}
\description{

The mmer function has been DEPRECATED but uses the Direct-Inversion Newton-Raphson or Average Information 

}

\usage{

mmer(fixed, random, rcov, data, weights, W, nIters=20, tolParConvLL = 1e-03, 
     tolParInv = 1e-06, init=NULL, constraints=NULL,method="NR", getPEV=TRUE,
     naMethodX="exclude", naMethodY="exclude",returnParam=FALSE, 
     dateWarning=TRUE,date.warning=TRUE,verbose=TRUE, reshapeOutput=TRUE, stepWeight=NULL,
     emWeight=NULL, contrasts=NULL)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, e.g.:
  
  \emph{response ~ covariate} for univariate models
  
  \emph{cbind(response.i,response.j) ~ covariate} for multivariate models
  
  }
  
  \item{random}{A formula specifying the name of the \strong{random effects}, e.g. \emph{random= ~ genotype + year}.
  
  }
  \item{rcov}{A formula specifying the name of the \strong{error term}, e.g., \emph{rcov= ~ units}.
  
  Special heterogeneous and special variance models and constraints for the residual part are the same used on the random term but the name of the random effect is always "units" which can be thought as a column with as many levels as rows in the data, e.g., \emph{rcov=~vsr(dsr(covariate),units)}
  
  }
  \item{data}{A data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{weights}{Name of the covariate for weights. To be used for the product R = Wsi*R*Wsi, where * is the matrix product, Wsi is the square root of the inverse of W and R is the residual matrix.}
  \item{W}{Alternatively, instead of providing a vector of weights the user can specify an entire W matrix (e.g., when covariances exist). To be used first to produce Wis = solve(chol(W)), and then calculate R = Wsi*R*Wsi.t(), where * is the matrix product, and R is the residual matrix. Only one of the arguments weights or W should be used. If both are indicated W will be given the preference.}
   \item{nIters}{Maximum number of iterations allowed.}
  \item{tolParConvLL}{Convergence criteria for the change in log-likelihood.}
  \item{tolParInv}{Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure.}
  \item{init}{Initial values for the variance components. By default this is NULL and initial values for the variance components are provided by the algorithm, but in case the user want to provide initial values for ALL var-cov components this argument is functional. It has to be provided as a list, where each list element corresponds to one random effect (1x1 matrix) and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits for such random effect. Initial values can also be provided in the Gti argument. Is highly encouraged to use the Gti and Gtc arguments of the vsr function instead of this argument, but these argument can be used to provide all initial values at once}
  \item{constraints}{When initial values are provided these have to be accompanied by their constraints. See the vsr function for more details on the constraints. Is highly encouraged to use the Gti and Gtc arguments of the vsr function instead of this argument but these argument can be used to provide all constraints at once.}
  \item{method}{This refers to the method or algorithm to be used for estimating variance components. Direct-inversion Newton-Raphson \strong{NR} and Average Information \strong{AI} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015).}
  \item{getPEV}{A TRUE/FALSE value indicating if the program should return the predicted error variance and variance for random effects. This option is provided since this can take a long time for certain models where p is > n by a big extent.}
   \item{naMethodX}{One of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the X matrices for fixed effects with the median value. If "exclude" is selected it will get rid of all rows with missing values for the X (fixed) covariates. The default is "exclude". The "include" option should be used carefully.}
   \item{naMethodY}{One of the three possible values; "include", "include2" or "exclude" (default) to treat the observations in response variable to be used in the estimation of variance components. The first option "include" will impute the response variables for all rows with the median value, whereas "include2" imputes the responses only for rows where there is observation(s) for at least one of the responses (only available in the multi-response models). If "exclude" is selected (default) it will get rid of rows in response(s) where missing values are present for at least one of the responses.}
   \item{returnParam}{A TRUE/FALSE value to indicate if the program should return the parameters to be used for fitting the model instead of fitting the model.}
   \item{dateWarning}{A TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package.}
   \item{date.warning}{A TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package. This argument will be removed soon, just left for backcompatibility.}
  \item{verbose}{A TRUE/FALSE value to indicate if the program should return the progress of the iterative algorithm.}
  \item{reshapeOutput}{A TRUE/FALSE value to indicate if the output should be reshaped to be easier to interpret for the user, some information is missing from the multivariate models for an easy interpretation.}
  \item{stepWeight}{A vector of values (of length equal to the number of iterations) indicating the weight used to multiply the update (delta) for variance components at each iteration. If NULL the 1st iteration will be multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument can help to avoid that variance components go outside the parameter space in the initial iterations which doesn't happen very often with the NR method but it can be detected by looking at the behavior of the likelihood. In that case you may want to give a smaller weight to the initial 8-10 iterations.}
  \item{emWeight}{A vector of values (of length equal to the number of iterations) indicating with values between 0 and 1 the weight assigned to the EM information matrix. And the values 1 - emWeight will be applied to the NR/AI information matrix to produce a joint information matrix.}
  \item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}
    
}

\details{

The use of this function requires a good understanding of mixed models. Please review the 'sommer.quick.start' vignette and pay attention to details like format of your random and fixed variables (e.g. character and factor variables have different properties when returning BLUEs or BLUPs, please see the 'sommer.changes.and.faqs' vignette). 

\strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

vignette("v1.sommer.quick.start")

vignette("v2.sommer.changes.and.faqs")

vignette("v3.sommer.qg")

vignette("v4.sommer.gxe")


\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications. 


\strong{Additional Functions}

Additional functions for genetic analysis have been included such as relationship 
matrix building (\code{\link{A.mat}}, \code{\link{D.mat}}, \code{\link{E.mat}}, 
\code{\link{H.mat}}), build a genotypic hybrid marker matrix (\code{\link{build.HMM}}), 
plot of genetic maps (\code{\link{map.plot}}), and manhattan plots (\code{\link{manhattan}}).
If you need to build a pedigree-based relationship matrix use the \code{getA} function
from the pedigreemm package. 

\strong{Bug report and contact}

If you have any technical questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com 

If you have any bug report please go to https://github.com/covaruber/sommer or send me an email to address it asap, just make sure you have read the vignettes carefully before sending your question.

\strong{Models Enabled}

For details about the models enabled and more information about the covariance structures please check the help page of the package (\code{\link{sommer}}).

}

\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{Vi}{the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1}
\item{P}{the projection matrix  Vi - [Vi*(X*Vi*X)^-*Vi] }
\item{sigma}{a list with the values of the variance-covariance components with one list element for each random effect.}
\item{sigma_scaled}{a list with the values of the scaled variance-covariance components with one list element for each random effect.}
\item{sigmaSE}{Hessian matrix containing the variance-covariance for the variance components. SE's can be obtained taking the square root of the diagonal values of the Hessian.}
\item{Beta}{a data frame for trait BLUEs (fixed effects).}
\item{VarBeta}{a variance-covariance matrix for trait BLUEs}
\item{U}{a list (one element for each random effect) with a data frame for trait BLUPs.}
\item{VarU}{a list (one element for each random effect) with the variance-covariance matrix for trait BLUPs.}
\item{PevU}{a list (one element for each random effect) with the predicted error variance matrix for trait BLUPs.}
\item{fitted}{ Fitted values y.hat=XB}
\item{residuals}{Residual values e = Y - XB}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{percChange}{The percent change of variance components across iterations. There should be one column less than the number of iterations. Calculated as percChange = ((x_i/x_i-1) - 1) * 100 where i is the ith iteration.}
\item{dL}{The vector of first derivatives of the likelihood with respect to the ith variance-covariance component.}
\item{dL2}{The matrix of second derivatives of the likelihood with respect to the i.j th variance-covariance component.}
\item{method}{The method for extimation of variance components specified by the user.}
\item{call}{Formula for fixed, random and rcov used.}
\item{constraints}{contraints used in the mixed models for the random effects.}
\item{constraintsF}{contraints used in the mixed models for the fixed effects.}
\item{data}{The dataset used in the model after removing missing records for the response variable.}
\item{dataOriginal}{The original dataset used in the model.}
\item{terms}{The name of terms for responses, fixed, random and residual effects in the model.}
\item{termsN}{The number of effects associated to fixed, random and residual effects in the model.}
\item{sigmaVector}{a vectorized version of the sigma element (variance-covariance components) to match easily the standard errors of the var-cov components stored in the element sigmaSE.}
\item{reshapeOutput}{The value provided to the mmer function for the argument with the same name.}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   

Covarrubias-Pazaran G. 2018. Software update: Moving the R package sommer to multivariate mixed models for genome-assisted prediction. doi: https://doi.org/10.1101/354639

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics 23 (2018): 52-71.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42:355-360.

}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_example)
DT <- DT_example
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmer(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)


}