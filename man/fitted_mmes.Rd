\name{fitted.mmes}
\alias{fitted.mmes}
\alias{fitted.mmer}
\title{fitted form a LMM fitted with mmes}
\description{\code{fitted} method for class \code{"mmes"}.
}
\usage{
\method{fitted}{mmes}(object, type=c("response", "link"), \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{type}{For a non-Gaussian PQL fit, return fitted means on the response
 scale (default) or linear predictors on the link scale.}
 \item{\dots}{Further arguments to be passed to the mmes function}
}
\value{
  For Gaussian fits, fitted values of the form y.hat = Xb + Zu. For non-Gaussian
  PQL fits, response-scale conditional means by default, or Xb + Zu with
  \code{type="link"}.
 }
\author{Giovanny Covarrubias}
\examples{
# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmes(Yield~1,
#               random=~vsm(ism(id),Gu=A)
#                       + Rowf + Colf + spl2Dc(Row,Col),
#               rcov=~units,
#               data=DT)
# 
# ff=fitted(mix1)
# 
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# lattice::wireframe(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted,  
#           aspect=c(61/87,0.4), drape=TRUE,# col.regions = colfunc,
#           light.source=c(10,0,10))
# lattice::levelplot(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted, col.regions = colfunc)

}
\seealso{\code{\link{fitted}}, \code{\link{mmes}}}
\keyword{models}

