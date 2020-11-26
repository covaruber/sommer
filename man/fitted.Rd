\name{fitted.mmer}
\alias{fitted.mmer}
\title{fitted form a LMM fitted with mmer}
\description{\code{fitted} method for class \code{"mmer"}.
}
\usage{
\method{fitted}{mmer}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmer"}}
 \item{\dots}{Further arguments to be passed to the mmer function}
}
\value{
  vector of fitted
 }
\author{Giovanny Covarrubias \email{cova_ruber@live.com.mx}}
\examples{
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#               random=~vs(id,Gu=A)
#                       + Rowf + Colf + vs(spl2D(Row,Col)),
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
\seealso{\code{\link{fitted}}, \code{\link{mmer}}}
\keyword{models}

