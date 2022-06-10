\name{spl2Db}
\alias{spl2Db}
\title{
Two-dimensional penalised tensor-product of marginal B-Spline basis.
}
\description{

Auxiliary function used for modelling the spatial or environmental effect as a two-dimensional penalised tensor-product (isotropic approach) based on Lee et al. (2013) and Rodriguez-Alvarez et al. (2018). \code{spl2Db} gets Tensor-Product P-Spline Mixed Model Incidence Matrices
    for use with \code{sommer} and its main function \code{mmer}. We thank Sue Welham for making the TPSbits package available to the community. If you're using this function for your research please cite her TPSbits package :) this is mostly a wrapper of her tpsmmb function to enable the use in sommer.

}
\usage{
spl2Db(x.coord,y.coord,at.var=NULL,at.levels=NULL,nsegments = c(10,10), 
       degree = c(3,3), penaltyord = c(2,2), nestorder = c(1,1),
       minbound=NULL, maxbound=NULL, method="Lee", what="bits")
}
\arguments{
  \item{x.coord}{vector of coordinates on the x-axis direction (i.e. row) to use in the 2 dimensional spline.}
  \item{y.coord}{vector of coordinates on the y-axis direction (i.e. range or column) to use in the 2 dimensional spline.}
  \item{at.var}{vector of indication variable where heterogeneous variance is required (e.g., a different spl2D for each field).}
  \item{at.levels}{character vector with the names of the leves for the at term that should be used, if missing all levels are used.}
  \item{nsegments}{numerical vector of length 2 containing the number of segments for each marginal (strictly \code{nsegments} - 1 is the number of internal knots in the domain of the covariate). Atomic values are also valid, being recycled. Default set to 10.}
  \item{degree}{numerical vector of length 2 containing the order of the polynomial of the B-spline basis for each marginal. Atomic values are also valid, being recycled. Default set to 3 (cubic B-splines).}
  \item{penaltyord}{numerical vector of length 2 containing the penalty order for each marginal. Atomic values are also valid, being recycled. Default set to 2 (second order). Currently, only second order penalties are allowed.}
  \item{nestorder}{numerical vector of length 2 containing the divisor of the number of segments (\code{nsegments}) to be used for the construction of the nested B-spline basis for the smooth-by-smooth interaction component. In this case, the nested B-spline basis will be constructed assuming a total of \code{nsegments}/\code{nestorder} segments. Default set to 1, which implies that nested basis are not used. See \code{SAP} for more details.}
  \item{minbound}{A list of length 2. The lower bound to be used for column and row dimensions respectively; default calculated as the minimum value for each dimension.}
  \item{maxbound}{A list of length 2. The upper bound to be used for column and row dimensions respectively; default calculated as the maximum value for each dimension.}
  \item{method}{A string. Method for forming the penalty; default="Lee" ie the penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The alternative method is "Wood" ie. the method from Wood et al (2012, Stat Comp 23, 341-360). This option is a research tool and requires further investigation.}
  \item{what}{one of two values; 'base' or 'bits' to return:
  
        base = matrix for columns cbind(TP.col,TP.row,TP.C.n,TP.R.n,TP.CR.n). To be used in the fixed part.
        
        bits = matrices for the tensor products. To be used in the random part.
  }
    
}
\details{

\strong{The following documentation is taken from the SpATS package. Please refer to this package and associated publications if you are interested in going deeper on this technique:}

Within the P-spline framework, anisotropic low-rank tensor-product smoothers have become the general approach for modelling multidimensional surfaces (Eilers and Marx 2003; Wood 2006). In the original SpATS package, was proposed to model the spatial or environmental effect by means of the tensor-product of B-splines basis functions. In other words, was proposed to model the spatial trend as a smooth bivariate surface jointly defined over the the spatial coordinates. Accordingly, the current function has been designed to allow the user to specify the spatial coordinates that the spatial trend is a function of. There is no restriction about how the spatial coordinates shall be specified: these can be the longitude and latitude of the position of the plot on the field or the column and row numbers. The only restriction is that the variables defining the spatial coordinates should be numeric (in contrast to factors).

As far as estimation is concerned, we have used in this package the equivalence between P-splines and linear mixed models (Currie and Durban, 2002). Under this approach, the smoothing parameters are expressed as the ratio between variance components. Moreover, the smooth components are decomposed in two parts: one which is not penalised (and treated as fixed) and one with is penalised (and treated as random). For the two-dimensional case, the mixed model representation leads also to a very interesting decomposition of the penalised part of the bivariate surface in three different components (Lee and Durban, 2011): (a) a component that contains the smooth main effect (smooth trend)	along one of the covariates that the surface is a function of (as, e.g, the x-spatial coordinate or column position of the plot in the field), (b) a component that contains the smooth main effect (smooth trend) along the other covariate (i.e., the y-spatial coordinate or row position); and (c) a smooth interaction component (sum of the linear-by-smooth interaction components and the smooth-by-smooth interaction component). 

The default implementation assumes two different smoothing parameters, i.e., one for each covariate in the smooth component. Accordingly, the same smoothing parameters are used for both, the main effects and the smooth interaction. However, this approach can be extended to deal with the ANOVA-type decomposition presented in Lee and Durban (2011). In their approach, four different smoothing parameters are considered for the smooth surface, that are in concordance with the aforementioned decomposition: (a) two smoothing parameter, one for each of the main effects; and (b) two smoothing parameter for the smooth interaction component. 

It should be noted that, the computational burden associated with the estimation of the two-dimensional tensor-product smoother might be prohibitive if the dimension of the marginal bases is large. In these cases, Lee et al. (2013) propose to reduce the computational cost by using nested bases. The idea is to reduce the dimension of the marginal bases (and therefore the associated number of parameters to be estimated), but only for the smooth-by-smooth interaction component. As pointed out by the authors, this simplification can be justified by the fact that the main effects would in fact explain most of the structure (or spatial trend) presented in the data, and so a less rich representation of the smooth-by-smooth interaction component could be needed. In order to ensure that the reduced bivariate surface is in fact nested to the model including only the main effects, Lee et al. (2013) show that the number of segments used for the nested basis should be a divisor of the number of segments used in the original basis (\code{nsegments} argument). In the present function, the divisor of the number of segments is specified through the argument \code{nestorder}. For a more detailed review on this topic, see Lee (2010) and Lee et al. (2013). The "PSANOVA" approach represents an alternative method. In this case, the smooth bivariate surface (or spatial trend) is decomposed in five different components each of them depending on a single smoothing parameter (see Lee et al., 2013).


}

\value{
List of length 7 elements:
\enumerate{
\item \code{data} = the input data frame augmented with structures required
    to fit tensor product splines in \code{asreml-R}. This data frame can be used
    to fit the TPS model.

    Added columns:
    \itemize{
    \item \code{TP.col}, \code{TP.row} = column and row coordinates
    \item \code{TP.CxR} = combined index for use with smooth x smooth term
    \item \code{TP.C.n} for n=1:{diff.c} = X parts of column spline for use
    in random model (where diff.c is the order of column differencing)
    \item \code{TP.R.n} for n=1:{diff.r} = X parts of row spline for use in
    random model (where diff.r is the order of row differencing)
    \item \code{TP.CR.n} for n=1:{(diff.c*diff.r)} = interaction between the
    two X parts for use in fixed model. The first variate is
    a constant term which should be omitted from the model when the constant
    (1) is present. If all elements are
    included in the model then the constant term should be omitted,
    eg. \code{y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...}
    \item when \code{asreml="grp"} or \code{"sepgrp"}, the spline basis
    functions are also added into the data frame. Column numbers for each
    term are given in the \code{grp} list structure.
    }
  
     
\item \code{fR} = Xr1:Zc
\item \code{fC} = Xr2:Zc
\item \code{fR.C} = Zr:Xc1
\item \code{R.fC} = Zr:Xc2 
\item \code{fR.fC} = Zc:Zr
\item \code{all} = Xr1:Zc | Xr2:Zc | Zr:Xc1 | Zr:Xc2 | Zc:Zr

}
}

\references{

Sue Welham (2021). TPSbits: Creates Structures to Enable Fitting and Examination of 2D Tensor-Product Splines using ASReml-R. R package version 1.0.0.
  
Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2018). SpATS: Spatial Analysis of Field Trials with Splines. R package version 1.0-9. https://CRAN.R-project.org/package=SpATS.

Rodriguez-Alvarez, M.X., et al. (2015) Fast smoothng parameter separaton n multdmensonal generalzed P-splnes: the SAP algorthm. Statistics and Computing 25.5: 941-957.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

}
\seealso{
\code{\link{mmer}}, \code{\link{spl2Da}}
}
\examples{
## ============================ ##
## example to use spl2Db() 
## ============================ ## 
data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# A <- A.mat(GT)
# mix <- mmer(Yield~1,
#             random=~vsr(id, Gu=A) +
#               vsr(Rowf) +
#               vsr(Colf) +
#               spl2Db(Row,Col),
#             rcov=~units,
#             data=DT)
# summary(mix)$varcomp
## ============================ ##
## mimic 2 fields
## ============================ ## 
# aa <- DT; bb <- DT
# aa$FIELD <- "A";bb$FIELD <- "B"
# set.seed(1234)
# aa$Yield <- aa$Yield + rnorm(length(aa$Yield),0,4)
# DT2 <- rbind(aa,bb)
# head(DT2)
# A <- A.mat(GT)
# mix <- mmer(Yield~1,
#             random=~vsr(dsr(FIELD),id, Gu=A) +
#               vsr(dsr(FIELD),Rowf) +
#               vsr(dsr(FIELD),Colf) +
#                 spl2Db(Row,Col,at.var=FIELD),
#             rcov=~vsr(dsr(FIELD),units),
#             data=DT2)
}
