\name{spl2Dc}
\alias{spl2Dc}
\title{
Two-dimensional penalised tensor-product of marginal B-Spline basis.
}
\description{

Auxiliary function used for modelling the spatial or environmental effect as a 
two-dimensional penalised tensor-product (isotropic approach) based on Lee et al. 
(2013) and Rodriguez-Alvarez et al. (2018). This is a modified wrapper of some 
portions of the SpATS package to build a single incidence matrix containing all
the columns from tensor products of the x and y coordinates and it fits such matrix 
as a single random effect. Then the heterogeneous covariances structure capabilities
of sommer can be used to enhance the model fit. You may be interested in reading and 
citing not only sommec but also Wageningen publications if using this 2D spline methodology.

}
\usage{
spl2Dc(x.coord,y.coord,at.var=NULL,at.levels=NULL, type="PSANOVA", 
      nsegments = c(10,10), penaltyord = c(2,2), degree = c(3,3), 
      nestorder = c(1,1), thetaC=NULL, theta=NULL, sp=FALSE)
}
\arguments{
  \item{x.coord}{vector of coordinates on the x-axis direction (i.e. row) to use in the 2 dimensional spline.}
  \item{y.coord}{vector of coordinates on the y-axis direction (i.e. range or column) to use in the 2 dimensional spline.}
  \item{at.var}{vector of indication variable where heterogeneous variance is required (e.g., a different spl2D for each field).}
  \item{at.levels}{character vector with the names of the leves for the at term that should be used, if missing all levels are used.}
  \item{type}{one of the two methods "PSANOVA" or "SAP". See details below.}
  \item{nsegments}{numerical vector of length 2 containing the number of segments for each marginal (strictly \code{nsegments} - 1 is the number of internal knots in the domain of the covariate). Atomic values are also valid, being recycled. Default set to 10.}
  \item{penaltyord}{numerical vector of length 2 containing the penalty order for each marginal. Atomic values are also valid, being recycled. Default set to 2 (second order). Currently, only second order penalties are allowed.}
  \item{degree}{numerical vector of length 2 containing the order of the polynomial of the B-spline basis for each marginal. Atomic values are also valid, being recycled. Default set to 3 (cubic B-splines).}
  \item{nestorder}{numerical vector of length 2 containing the divisor of the number of segments (\code{nsegments}) to be used for the construction of the nested B-spline basis for the smooth-by-smooth interaction component. In this case, the nested B-spline basis will be constructed assuming a total of \code{nsegments}/\code{nestorder} segments. Default set to 1, which implies that nested basis are not used. See \code{SAP} for more details.}
  \item{thetaC}{an optional matrix for constraints in the variance components.}
  \item{theta}{an optional matrix for initial values of the variance components.}
  \item{sp}{a TRUE/FALSE statement to indicate if the VC from this structure should be multiplied by the scale parameter added in the mmes function through the addScaleParam argument in the \code{mmes} function .}
     
}
\details{

\strong{The following documentation is taken from the SpATS package. Please refer to this package and associated publications if you are interested in going deeper on this technique:}

Within the P-spline framework, anisotropic low-rank tensor-product smoothers have become the general approach for modelling multidimensional surfaces (Eilers and Marx 2003; Wood 2006). In the original SpATS package, was proposed to model the spatial or environmental effect by means of the tensor-product of B-splines basis functions. In other words, was proposed to model the spatial trend as a smooth bivariate surface jointly defined over the the spatial coordinates. Accordingly, the current function has been designed to allow the user to specify the spatial coordinates that the spatial trend is a function of. There is no restriction about how the spatial coordinates shall be specified: these can be the longitude and latitude of the position of the plot on the field or the column and row numbers. The only restriction is that the variables defining the spatial coordinates should be numeric (in contrast to factors).

As far as estimation is concerned, we have used in this package the equivalence between P-splines and linear mixed models (Currie and Durban, 2002). Under this approach, the smoothing parameters are expressed as the ratio between variance components. Moreover, the smooth components are decomposed in two parts: one which is not penalised (and treated as fixed) and one with is penalised (and treated as random). For the two-dimensional case, the mixed model representation leads also to a very interesting decomposition of the penalised part of the bivariate surface in three different components (Lee and Durban, 2011): (a) a component that contains the smooth main effect (smooth trend)	along one of the covariates that the surface is a function of (as, e.g, the x-spatial coordinate or column position of the plot in the field), (b) a component that contains the smooth main effect (smooth trend) along the other covariate (i.e., the y-spatial coordinate or row position); and (c) a smooth interaction component (sum of the linear-by-smooth interaction components and the smooth-by-smooth interaction component). 

The original implementation of SpATS assumes two different smoothing parameters, i.e., one for each covariate in the smooth component. Accordingly, the same smoothing parameters are used for both, the main effects and the smooth interaction. However, this approach can be extended to deal with the ANOVA-type decomposition presented in Lee and Durban (2011). In their approach, four different smoothing parameters are considered for the smooth surface, that are in concordance with the aforementioned decomposition: (a) two smoothing parameter, one for each of the main effects; and (b) two smoothing parameter for the smooth interaction component.

It should be noted that, the computational burden associated with the estimation of the two-dimensional tensor-product smoother might be prohibitive if the dimension of the marginal bases is large. In these cases, Lee et al. (2013) propose to reduce the computational cost by using nested bases. The idea is to reduce the dimension of the marginal bases (and therefore the associated number of parameters to be estimated), but only for the smooth-by-smooth interaction component. As pointed out by the authors, this simplification can be justified by the fact that the main effects would in fact explain most of the structure (or spatial trend) presented in the data, and so a less rich representation of the smooth-by-smooth interaction component could be needed. In order to ensure that the reduced bivariate surface is in fact nested to the model including only the main effects, Lee et al. (2013) show that the number of segments used for the nested basis should be a divisor of the number of segments used in the original basis (\code{nsegments} argument). In the present function, the divisor of the number of segments is specified through the argument \code{nestorder}. For a more detailed review on this topic, see Lee (2010) and Lee et al. (2013). The "PSANOVA" approach represents an alternative method. In this case, the smooth bivariate surface (or spatial trend) is decomposed in five different components each of them depending on a single smoothing parameter (see Lee et al., 2013).

-----------------

As mentioned at the beginning, the piece of documentation stated above was taken completely from the SpATS package in order to provide a deeper explanation. In practice, sommec uses some pieces of code from SpATS to build the design matrix containing all the columns from tensor products of the x and y coordinates and it fits such matrix as a single random effect. As a result the same variance component is assumed for the linear, linear by linear, linear by spline, and spline by spline interactions. This results in a less flexible approach than the one proposed by Rodriguez-Alvarez et al. (2018) but still makes a pretty good job to model the spatial variation. Use under your own risk.

}
\references{

Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2018). SpATS: Spatial Analysis of Field Trials with Splines. R package version 1.0-9. https://CRAN.R-project.org/package=SpATS.

Rodriguez-Alvarez, M.X., et al. (2015) Fast smoothng parameter separaton n multdmensonal generalzed P-splnes: the SAP algorthm. Statistics and Computing 25.5: 941-957.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for Natural and Extraneous Variation in the Analysis of Field Experiments. Journal of Agricultural, Biological, and Environmental Statistics, 2, 269 - 293.

}
\seealso{
\code{\link{mmes}}
}
\examples{
## ============================ ##
## example to use spl2Dc() 
## ============================ ## 
data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# A <- A.mat(GT)
## ============================ ##
## mimic 2 fields
## ============================ ## 
# aa <- DT; bb <- DT
# aa$FIELD <- "A";bb$FIELD <- "B"
# set.seed(1234)
# aa$Yield <- aa$Yield + rnorm(length(aa$Yield),0,4)
# DT2 <- rbind(aa,bb)
# head(DT2)
# mix <- mmes(Yield~1,
#             random=~vsm(dsm(FIELD),ism(id), Gu=Ai) +
#               vsm(dsm(FIELD),ism(Rowf)) +
#               vsm(dsm(FIELD),ism(Colf)) +
#                 spl2Dc(Row,Col,at.var=FIELD),
#             rcov=~vsm(dsm(FIELD),ism(units)),
#             data=DT2)

}
