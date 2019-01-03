\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

\strong{So}lving \strong{M}ixed \strong{M}odel \strong{E}quations in \strong{R}
\cr

\if{html}{\figure{mai.png}{options: height="10\%" alt="Figure: mai.png"}}

}
\description{

Sommer is a structural multivariate-univariate linear mixed model solver for multiple random effects allowing the specification and/or estimation of variance covariance structures. ML/REML estimates can be obtained using the Direct-Inversion Newton-Raphson, Average Information and Efficient Mixed Model Association algorithms coded in C++ using the Armadillo library to optimize dense matrix operations common in genomic selection models. Sommer was designed for genomic prediction and genome wide association studies (GWAS) to include i.e. additive, dominance and epistatic relationship structures or other covariance structures, but also functional as a regular mixed model program.  

The sommer package has been developed to provide R users with open-source code to understand how most popular likelihood algorithms in mixed model analysis work, but at the same time allowing to perform their real analysis in diploid and polyploid organisms with small and medium-size data sets (< 10,000 observations for average computers given the computational burden carried by the direct-inversion algorithms). The package is focused in the p > n problem and dense covariance structures when the direct-inversion algorithm becomes faster than MME-based algorithms. \strong{ The core of the package is the \code{\link{mmer}} (formula-based) function} that fit the multivariate linear mixed models. This package returns variance-covariance components, BLUPs, BLUEs, residuals, fitted values, variances-covariances for fixed and random effects, etc. The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance (\code{\link{D.mat}}), and epistatic (\code{\link{E.mat}}) relationship matrices for diploid and polyploid organisms. It also provides flexibility to fit other genetic models such as full and half diallel models and random regression models. In addition the \code{\link{pin}} function can be used to estimate standard errors for linear combinations of variance components (i.e. ratios like h2). A good converter from letter code to numeric format is implemented in the function \code{\link{atcg1234}}, which supports higher ploidy levels than diploid. Recently, spatial modeling has been added added to sommer using the two-dimensional splines (\code{\link{spl2D}}) to face the lack of AR1 covariance structures in sommer.

Starting with version 3.0 the GWAS models that use the M argument have their own function \code{\link{GWAS}} and is not part of the \code{\link{mmer}} function.
}

\section{Keeping updated}{
The sommer package is updated on CRAN every 3-months due to CRAN policies but you can find the latest source at https://github.com/covaruber/sommer . This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/sommer")

This is recommended since bugs fixes will be immediately available in the GitHub source but not in CRAN until the next update.
}

\section{Tutorials}{

 \strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

\strong{vignette("sommer.start")} 

\strong{vignette("sommer")}

or visit \strong{https://covaruber.github.io}

}

\section{Getting started}{
The package has been equiped with several datasets to learn how to use the sommer package: 

* \code{\link{DT_halfdiallel}} and \code{\link{DT_fulldiallel}} datasets have examples to fit half and full diallel designs. 

* \code{\link{DT_h2}} to calculate heritability

* \code{\link{DT_cornhybrids}} and \code{\link{DT_technow}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link{DT_wheat}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link{DT_cpdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 

* \code{\link{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link{DT_gryphon}} data contains an example of an animal model including pedigree information.

* \code{\link{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link{DT_legendre}} simulated dataset for random regression model.

Other functions such as \code{\link{summary}}, \code{\link{fitted}}, \code{\link{randef}} (notice here is randef not ranef), \code{\link{anova}}, \code{\link{variogram}}, \code{\link{residuals}}, \code{\link{coef}} and \code{\link{plot}} applicable to typical linear models can also be applied to models fitted using the mmer-type of functions. Additional functions for genetic analysis have been included such as heritability (\code{\link{h2.fun}}), build a genotypic hybrid marker matrix (\code{\link{build.HMM}}), plot of genetic maps (\code{\link{map.plot}}), creation of manhattan plots (\code{\link{manhattan}}). If you need to use pedigree you need to convert your pedigree into a relationship matrix (use the `getA` function from the pedigreemm package).

Useful functions for analyzing field trials are included such as the \code{\link{spl2D}}, and \code{\link{fill.design}}.

}

\section{Differences of sommer >= 3.7 with previous versions}{

Since version 3.7 I have completly redefined the specification of the variance-covariance structures to provide more flexibility to the user. This has particularly helped the residual covariance structures and the easier combination of custom random effects and overlay models. I think that although this will bring some uncomfortable situations at the beggining, in the long term this will help users to fit better models. In esence, I have abandoned the asreml formulation (not the structures available) given it's limitations to combine some of the sommer structures but all covariance structures can now be fitted using the `vs` functions. 

}
\section{Differences of sommer >= 3.0 with previous versions}{

Since version 3.0 I have decided to focus in developing the multivariate solver and for doing this I have decided to remove the M argument (for GWAS analysis) from the mmer function and move it to it's own function \code{\link{GWAS}}.

Before the mmer solver had implemented the us(trait), diag(trait), at(trait) asreml formulation for multivariate models that allow to specify the structure of the trait in multivariate models. Therefore the MVM argument was no longer needed. After version 3.7 now the multi-trait structures can be specified in the \code{Gt} and \code{Gtc} arguments of the \code{\link{vs}} function.

The Average Information algorithm had been removed in the past from the package because of its instability to deal with very complex models without good initial values. Now after 3.7 I have brought it back after I noticed that starting with NR the first three iterations gives enough flexibility to the AI algorithm.

Keep in mind that sommer uses direct inversion (DI) algorithm which can be very slow for datasets with many observations (big 'n'). The package is focused in problems of the type p > n (more random effect(s) levels than observations) and models with dense covariance structures. For example, for experiment with dense covariance structures with low-replication (i.e. 2000 records from 1000 individuals replicated  twice with a covariance structure of 1000x1000) sommer will be faster than MME-based software. Also for genomic problems with large number of random effect levels, i.e. 300 individuals (n) with 100,000 genetic markers (p). On the other hand, for highly replicated trials with small covariance structures or n > p (i.e. 2000 records from 200 individuals replicated 10 times with covariance structure of 200x200) asreml or other MME-based algorithms will be much faster and I recommend you to use that software.

}

\section{Models Enabled}{


The core of the package is the \code{\link{mmer}} (formula-based) function which solve the mixed model equations. The functions are an interface to call the `NR` Direct-Inversion Newton-Raphson, `AI` Direct-Inversion Average Information (Tunnicliffe 1989; Gilmour et al. 1995; Lee et al. 2016). Since version 2.0 sommer can handle multivariate models. Following Maier et al. (2015), the multivariate (and by extension the univariate) mixed model implemented has the form:




\if{html}{\figure{form1.png}{options: height="10\%" alt="Figure: form1.png"}}



where y_i is a vector of trait phenotypes, \eqn{\beta_i} is a vector of fixed effects, u_i is a vector of random effects for individuals and e_i are residuals for trait i (i = 1,..., t). The random effects (u_1 ... u_i and e_i) are assumed to be normally distributed with mean zero. X and Z are incidence matrices for fixed and random effects respectively. The distribution of the multivariate response and the phenotypic variance covariance (V) are:




\if{html}{\figure{form2.png}{options: height="10\%" alt="Figure: form2.png"}}




where K is the relationship or covariance matrix for the kth random effect (u=1,...,k), and R=I is an identity matrix for the residual term. The terms \eqn{\sigma^2_{g_{i}}} and \eqn{\sigma^2_{\epsilon_{i}}} denote the genetic (or any of the kth random terms) and residual variance of trait i, respectively and \eqn{\sigma_{g_{_{ij}}}} and \eqn{\sigma_{\epsilon_{_{ij}}}} the genetic (or any of the kth random terms) and residual covariance between traits i and j (i=1,...,t, and j=1,...,t). The algorithm implemented optimizes the log likelihood:




\if{html}{\figure{form3.png}{options: height="10\%" alt="Figure: form3.png"}}




where || is the determinant of a matrix. And the REML estimates are updated using a Newton optimization algorithm of the form:




\if{html}{\figure{form4.png}{options: height="10\%" alt="Figure: form4.png"}}




Where, \eqn{\theta} is the vector of variance components for random effects and covariance components among traits, H^{-1} is the inverse of the Hessian matrix of second derivatives for the kth cycle, dL/d\eqn{\sigma^2_i} is the vector of first derivatives of the likelihood with respect to the variance-covariance components. The Eigen decomposition of the relationship matrix proposed by Lee and Van Der Werf (2016) was included in the Newton-Raphson algorithm to improve time efficiency. Additionally, the popular pin function to estimate standard errors for linear combinations of variance components (i.e. heritabilities and genetic correlations) was added to the package as well.


}

\section{GWAS Models}{

The GWAS models in the sommer package are enabled by using the M argument in the functions \code{\link{GWAS}}, which is expected to be a numeric marker matrix. Markers are treated as fixed effects according to the model proposed by Yu et al. (2006) for diploids, and Rosyara et al. (2016) (for polyploids). The matrices X and M are both fixed effects, but they are separated by 2 different arguments to distinguish factors such as environmental and design factors for the argument "X" and markers with "M".
  
  The genome-wide association analysis is based on the mixed model:

\deqn{y = X \beta + Z g + M \tau + e}

where \eqn{\beta} is a vector of fixed effects that can model both environmental factors and population structure.  
The variable \eqn{g} models the genetic background of each line as a random effect with \eqn{Var[g] = K \sigma^2}.  
The variable \eqn{\tau} models the additive SNP effect as a fixed effect.  The residual variance is \eqn{Var[\varepsilon] = I \sigma_e^2}
  
  When principal components are included (P+K model), the loadings are determined from an eigenvalue decomposition of the K matrix and are used in the fixed effect part.

The argument "P3D" introduced by Zhang et al. (2010) can be used with the P3D argument.  When P3D=FALSE, this function is equivalent to AI/NR with REML where the variance components are estimated for each SNP or marker tested (Kang et al. 2008).  When P3D=TRUE, it is equivalent to NR (Kang et al. 2010) where the assumption is that variance components for all SNP/markers are the same and therefore the variance components are estimated only once (and markers are tested in a WLS framework being the the weight matrix (M) the inverse of the phenotypic variance matrix (V)).  Therefore, P3D=TRUE option is faster but can underestimate significance compared to P3D=FALSE.

Multivariate GWAS are based in Covarrubias-Pazaran et al. (2018, In preparation), which adjusts betas for all response variables and then does the regular GWAS with such adjusted betas or marker effects.

For extra details about the methods please read the canonical papers listed in the References section.

}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com and send me an email with the link at cova_ruber@live.com.mx

I'll be glad to help or answer any question. I have spent a valuable amount of time developing this package. Please cite this package in your publication. Type 'citation("sommer")' to know how to cite it.


}
\references{

Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Covarrubias-Pazaran G. 2018. Software update: Moving the R package sommer to multivariate mixed models for genome-assisted prediction. doi: https://doi.org/10.1101/354639
  
Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

}
\keyword{ R package }
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
             rcov=~vs(us(Env),units), 
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


}