\name{mmer}
\alias{mmer}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations in \strong{R}

}
\description{

Sommer is a structural multivariate-univariate linear mixed model solver for multiple random effects allowing the specification and/or estimation of variance covariance structures. REML estimates can be obtained using the Direct-Inversion Newton-Raphson, Average Information and Efficient Mixed Model Association algorithms coded in C++ using the Armadillo library to optimize dense matrix operations common in genomic selection models. Sommer was designed for genomic prediction and genome wide association studies (GWAS), but also functional as a regular mixed model program. These algorithms are \strong{intended to be used for problems of the type p > n and/or dense matrices}. For problems with sparse covariance structures, or problems of the type n > p, the MME-based algorithms are faster and we recommend to shift to the use of such software (e.g., lme4, breedR, or asreml-R).

\strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

\strong{vignette("v1.sommer.quick.start")} 

\strong{vignette("v2.sommer.changes.and.faqs")}

\strong{vignette("v3.sommer.qg")}

\strong{vignette("v4.sommer.gxe")}

}

\usage{

mmer(fixed, random, rcov, data, weights, iters=20, tolpar = 1e-03, 
     tolparinv = 1e-06, init=NULL, constraints=NULL,method="NR", getPEV=TRUE,
     na.method.X="exclude", na.method.Y="exclude",return.param=FALSE, 
     date.warning=TRUE,verbose=TRUE, reshape.output=TRUE, stepweight=NULL,
     emupdate=NULL)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:
  
  \emph{response ~ covariate} for univariate models
  
  \emph{cbind(response.i,response.j) ~ covariate} for multivariate models
  
  The \code{\link{fcm}()} function can be used to constrain fixed effects in multi-response models.
  
  }
  
  \item{random}{A formula specifying the name of the \strong{random effects}, i.e. \emph{random= ~ genotype + year}.
  
  Useful functions can be used to fit heterogeneous variances and other special models (\emph{see 'Special Functions' in the Details section for more information}):
  
  \code{\link{vs}(...,Gu,Gti,Gtc)} is the main function to specify variance models and special structures for random effects. On the ... argument you provide the unknown variance-covariance structures (i.e. us,ds,at,cs) and the random effect where such covariance structure will be used (the random effect of interest). Gu is used to provide known covariance matrices among the levels of the random effect, Gti initial values and Gtc for constraints. Auxiliar functions for building the variance models are:
  
  ** \code{\link{ds}(x)}, \code{\link{us}(x)}, \code{\link{cs}(x)} and \code{\link{at}(x,levs)} can be used to specify unknown diagonal, unstructured and customized unstructured and diagonal covariance structures to be estimated by REML.
  
  ** \code{\link{unsm}(x)}, \code{\link{uncm}(x)}, \code{\link{fixm}(x)} and \code{\link{diag}(x)} can be used to build easily matrices to specify constraints in the Gtc argument of the \code{\link{vs}()} function.
  
  ** \code{\link{overlay}()}, \code{\link{spl2Da}()}, \code{\link{spl2Db}()}, and \code{\link{leg}()} functions can be used to specify overlayed of design matrices of random effects, two dimensional spline and random regression models within the \code{\link{vs}()} function.
  
    \code{\link{gvs}(...,Gu,Guc,Gti,Gtc)} is an alternative function to specify general variance structures between different random effects. An special case in the indirect genetic effect models. Is similar to the vs function but in the ... argument the different random effects are provided.
  
  }
  \item{rcov}{A formula specifying the name of the \strong{error term}, i.e. \emph{rcov= ~ units}.
  
  Special heterogeneous and special variance models and constraints for the residual part are the same used on the random term but the name of the random effect is always "units" which can be thought as a column with as many levels as rows in the data, i.e. \emph{rcov=~vs(ds(covariate),units)}
  
  }
  \item{data}{A data frame containing the variables specified in the formulas for response, fixed, and random effects.}
  \item{weights}{Name of the covariate for weights. To be used for the product R = Wsi*R*Wsi, where * is the matrix product, Wsi is the square root of the inverse of W and R is the residual matrix.}
   \item{iters}{Maximum number of iterations allowed.}
  \item{tolpar}{Convergence criteria for the change in log-likelihood.}
  \item{tolparinv}{Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure.}
  \item{init}{Initial values for the variance components. By default this is NULL and initial values for the variance components are provided by the algorithm, but in case the user want to provide initial values for ALL var-cov components this argument is functional. It has to be provided as a list, where each list element corresponds to one random effect (1x1 matrix) and if multitrait model is pursued each element of the list is a matrix of variance covariance components among traits for such random effect. Initial values can also be provided in the Gti argument of the \link{vs} function. Is highly encouraged to use the Gti and Gtc arguments of the \link{vs} function instead of this argument, but these argument can be used to provide all initial values at once}
  \item{constraints}{When initial values are provided these have to be accompanied by their constraints. See the \link{vs} function for more details on the constraints. Is highly encouraged to use the Gti and Gtc arguments of the \link{vs} function instead of this argument but these argument can be used to provide all constraints at once.}
  \item{method}{This refers to the method or algorithm to be used for estimating variance components. Direct-inversion Newton-Raphson \strong{NR} and Average Information \strong{AI} (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2015).}
  \item{getPEV}{A TRUE/FALSE value indicating if the program should return the predicted error variance and variance for random effects. This option is provided since this can take a long time for certain models where p is > n by a big extent.}
   \item{na.method.X}{One of the two possible values; "include" or "exclude". If "include" is selected then the function will impute the X matrices for fixed effects with the median value. If "exclude" is selected it will get rid of all rows with missing values for the X (fixed) covariates. The default is "exclude". The "include" option should be used carefully.}
   \item{na.method.Y}{One of the three possible values; "include", "include2" or "exclude" (default) to treat the observations in response variable to be used in the estimation of variance components. The first option "include" will impute the response variables for all rows with the median value, whereas "include2" imputes the responses only for rows where there is observation(s) for at least one of the responses (only available in the multi-response models). If "exclude" is selected (default) it will get rid of rows in response(s) where missing values are present for at least one of the responses.}
   \item{return.param}{A TRUE/FALSE value to indicate if the program should return the parameters to be used for fitting the model instead of fitting the model.}
   \item{date.warning}{A TRUE/FALSE value to indicate if the program should warn you when is time to update the sommer package.}
  \item{verbose}{A TRUE/FALSE value to indicate if the program should return the progress of the iterative algorithm.}
  \item{reshape.output}{A TRUE/FALSE value to indicate if the output should be reshaped to be easier to interpret for the user, some information is missing from the multivariate models for an easy interpretation.}
  \item{stepweight}{A vector of values (of length equal to the number of iterations) indicating the weight used to multiply the update (delta) for variance components at each iteration. If NULL the 1st iteration will be multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument can help to avoid that variance components go outside the parameter space in the initial iterations which doesn't happen very often with the NR method but it can be detected by looking at the behavior of the likelihood. In that case you may want to give a smaller weight to the initial 8-10 iterations.}
  \item{emupdate}{A vector of values (of length equal to the number of iterations) indicating with 1's the iterations where the algorithm should calculate an EM update instead of an update based on first and second derivatives indicated with zeros (i.e. NR or AI methods).}
    
}

\details{

The use of this function requires a good understanding of mixed models. Please review the 'sommer.quick.start' vignette and pay attention to details like format of your random and fixed variables (i.e. character and factor variables have different properties when returning BLUEs or BLUPs, please see the 'sommer.changes.and.faqs' vignette). 

\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications. 

\strong{Special variance structures}

 \code{\link{vs}(\link{at}(x,levels),y)}   
 
 can be used to specify heterogeneous variance for the "y" covariate at specific levels of the covariate "x", i.e. \emph{random=~vs(at(Location,c("A","B")),ID)} fits a variance component for ID at levels A and B of the covariate Location.
  
 \code{\link{vs}(\link{ds}(x),y)}   
 
 can be used to specify a diagonal covariance structure for the "y" covariate for all levels of the covariate "x", i.e. \emph{random=~vs(ds(Location),ID)}  fits a variance component for ID at all levels of the covariate Location.
  
 \code{\link{vs}(\link{us}(x),y)}   
 
 can be used to specify an unstructured covariance structure for the "y" covariate for all levels of the  covariate "x", i.e. \emph{random=~vs(us(Location),ID)}  fits variance and covariance components for ID at all levels of the covariate Location.
 
  \code{\link{vs}(\link{overlay}(...,rlist=NULL,prefix=NULL))}   
 
 can be used to specify overlay of design matrices between consecutive random effects specified, i.e. \emph{random=~vs(overlay(male,female))} overlays (overlaps) the incidence matrices for the male and female random effects to obtain a single variance component for both effects. The `rlist` argument is a list with each element being a numeric value that multiplies the incidence matrix to be overlayed. See \code{\link{overlay}} for details.Can be combined with vs().
  
  \code{\link{vs}(\link{leg}(x,n),y)}
      
  can be used to fit a random regression model using a numerical variable \code{x} that marks the trayectory for the random effect \code{y}. The leg function can be combined with the special functions \code{ds}, \code{us} \code{at} and \code{cs}. For example \emph{random=~vs(leg(x,1),y)} or \emph{random=~vs(us(leg(x,1)),y)}.
  
  \code{\link{vs}(x,Gtc=\link{fcm}(v))}
      
  can be used to constrain fixed effects in the multi-response mixed models. This is a vector that specifies if the fixed effect is to be estimated for such trait. For example \emph{fixed=cbind(response.i, response.j)~vs(Rowf, Gtc=fcm(c(1,0)))} means that the fixed effect Rowf should only be estimated for the first response and the second should only have the intercept.
  
  \code{\link{gvs}(x,y)}  
  
  can be used to fit variance and covariance parameters between two or more random effects. For example, indirect genetic effect models.
  
  \code{\link{spl2Da}(x.coord, y.coord, at.var, at.levels))}
      
  can be used to fit a 2-dimensional spline (i.e. spatial modeling) using coordinates \code{x.coord} and \code{y.coord} (in numeric class) assuming a single variance component. The 2D spline can be fitted at specific levels using the \code{at.var} and \code{at.levels} arguments. For example \emph{random=~spl2Da(x.coord=Row,y.coord=Range,at.var=FIELD)}.
  
  \code{\link{spl2Db}(x.coord, y.coord, at.var, at.levels))}
      
  can be used to fit a 2-dimensional spline (i.e. spatial modeling) using coordinates \code{x.coord} and \code{y.coord} (in numeric class) assuming multiple variance components. The 2D spline can be fitted at specific levels using the \code{at.var} and \code{at.levels} arguments. For example \emph{random=~spl2Db(x.coord=Row,y.coord=Range,at.var=FIELD)}.

\strong{S3 methods}{

S3 methods are available for some parameter extraction such as \code{\link{fitted.mmer}}, \code{\link{residuals.mmer}}, \code{\link{summary.mmer}}, \code{\link{randef}}, \code{\link{coef.mmer}}, \code{\link{anova.mmer}}, \code{\link{plot.mmer}}, and \code{\link{predict.mmer}} to obtain adjusted means. In addition, the \code{\link{vpredict}} function (replacement of the pin function) can be used to estimate standard errors for linear combinations of variance components (i.e. ratios like h2). 

}

\strong{Additional Functions}

Additional functions for genetic analysis have been included such as relationship matrix building (\code{\link{A.mat}}, \code{\link{D.mat}}, \code{\link{E.mat}}, \code{\link{H.mat}}),  heritability (\code{\link{h2.fun}}), build a genotypic hybrid marker matrix (\code{\link{build.HMM}}), plot of genetic maps (\code{\link{map.plot}}), and manhattan plots (\code{\link{manhattan}}). If you need to build a pedigree-based relationship matrix use the \code{getA} function from the pedigreemm package. 

\strong{Bug report and contact}{

If you have any technical questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com 

If you have any bug report please go to https://github.com/covaruber/sommer or send me an email to address it asap, just make sure you have read the vignettes carefully before sending your question.

}
\strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the sommer package: 

* \code{\link{DT_halfdiallel}}, \code{\link{DT_fulldiallel}} and \code{\link{DT_mohring}} datasets have examples to fit half and full diallel designs. 

* \code{\link{DT_h2}} to calculate heritability

* \code{\link{DT_cornhybrids}} and \code{\link{DT_technow}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link{DT_wheat}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{DT_cpdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 

* \code{\link{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{DT_gryphon}} data contains an example of an animal model including pedigree information.

* \code{\link{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link{DT_legendre}} simulated dataset for random regression model.

* \code{\link{DT_sleepstudy}} dataset to know how to translate lme4 models to sommer models.

* \code{\link{DT_ige}} dataset to show how to fit indirect genetic effect models.

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
\item{reshape.output}{The value provided to the mmer function for the argument with the same name.}

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

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####

## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmer(Yield~Env,
             random= ~Name + vs(ds(Env),Name),
             rcov= ~ vs(ds(Env),units),
             data=DT)
summary(ans2)

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmer(Yield~Env,
             random=~ vs(us(Env),Name),
             rcov=~vs(ds(Env),units), 
             data=DT)
summary(ans3)

# ####==========================================####
# #### Multivariate homogeneous variance models ####
# ####==========================================####
# 
# ## Multivariate Compound simmetry (CS) model
# DT$EnvName <- paste(DT$Env,DT$Name)
# ans4 <- mmer(cbind(Yield, Weight) ~ Env,
#               random= ~ vs(Name, Gtc = unsm(2)) + vs(EnvName,Gtc = unsm(2)),
#               rcov= ~ vs(units, Gtc = unsm(2)),
#               data=DT)
# summary(ans4)
# 
# ####=============================================####
# #### Multivariate heterogeneous variance models  ####
# ####=============================================####
# 
# ## Multivariate Compound simmetry (CS) + Diagonal (DIAG) model
# ans5 <- mmer(cbind(Yield, Weight) ~ Env,
#               random= ~ vs(Name, Gtc = unsm(2)) + vs(ds(Env),Name, Gtc = unsm(2)),
#               rcov= ~ vs(ds(Env),units, Gtc = unsm(2)),
#               data=DT)
# summary(ans5)
# 
# ####===========================================####
# #### Multivariate unstructured variance models ####
# ####===========================================####
# 
# ans6 <- mmer(cbind(Yield, Weight) ~ Env,
#               random= ~ vs(us(Env),Name, Gtc = unsm(2)),
#               rcov= ~ vs(ds(Env),units, Gtc = unsm(2)),
#               data=DT)
# summary(ans6)
# 
# ####=========================================####
# ####=========================================####
# #### EXAMPLE SET 2
# #### 2 variance components
# #### one random effect with variance covariance structure
# ####=========================================####
# ####=========================================####
# 
# data("DT_cpdata")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# head(DT)
# GT[1:4,1:4]
# #### create the variance-covariance matrix
# A <- A.mat(GT)
# #### look at the data and fit the model
# mix1 <- mmer(Yield~1,
#              random=~vs(id, Gu=A) + Rowf,
#              rcov=~units,
#              data=DT)
# summary(mix1)$varcomp
# #### calculate heritability
# vpredict(mix1, h1 ~ V1/(V1+V3) )
# #### multi trait example
# mix2 <- mmer(cbind(Yield,color)~1,
#               random = ~ vs(id, Gu=A, Gtc = unsm(2)) + # unstructured at trait level
#                             vs(Rowf, Gtc=diag(2)) + # diagonal structure at trait level
#                                 vs(Colf, Gtc=diag(2)), # diagonal structure at trait level
#               rcov = ~ vs(units, Gtc = unsm(2)), # unstructured at trait level
#               data=DT)
# summary(mix2)
# 
# ####=========================================####
# #### EXAMPLE SET 3
# #### comparison with lmer, install 'lme4' 
# #### and run the code below
# ####=========================================####
# 
# #### lmer cannot use var-cov matrices so we will not 
# #### use them in this comparison example
# 
# library(lme4)
# library(sommer)
# data("DT_cornhybrids")
# DT <- DT_cornhybrids
# DTi <- DTi_cornhybrids
# GT <- GT_cornhybrids
# 
# fm1 <- lmer(Yield ~ Location + (1|GCA1) + (1|GCA2) + (1|SCA),
#             data=DT )
# out <- mmer(Yield ~ Location,
#              random = ~ GCA1 + GCA2 + SCA,
#              rcov = ~ units,
#              data=DT)
# summary(fm1)
# summary(out)
# ### same BLUPs for GCA1, GCA2, SCA than lme4
# plot(out$U$GCA1$Yield, ranef(fm1)$GCA1[,1])
# plot(out$U$GCA2$Yield, ranef(fm1)$GCA2[,1])
# vv=which(abs(out$U$SCA$Yield) > 0)
# plot(out$U$SCA$Yield[vv], ranef(fm1)$SCA[,1])
# 
# ### a more complex model specifying which locations
# head(DT)
# out2 <- mmer(Yield ~ Location,
#               random = ~ vs(at(Location,c("3","4")),GCA2) +
#                          vs(at(Location,c("3","4")),SCA),
#               rcov = ~ vs(ds(Location),units),
#               data=DT)
# summary(out2)


}