\name{predict.mmer}
\alias{predict.mmer}
\title{Predict form a LMM fitted with mmer}
\description{\code{predict} method for class \code{"mmer"}.
}
\usage{
\method{predict}{mmer}(object, classify=NULL,
        \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmer"}}
 \item{classify}{an optional character string with the variables that define the margins of the multiway table to be predicted (see Details).}
 \item{\dots}{Further arguments to be passed to the model fit}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument classify). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 
  
  For predicted values the pertinent design matrices X and Z together with BLUEs (b) and BLUPs (u) are multiplied and added together.
  
  predicted.value equal Xb + Zu.1 + ... + Zu.n
  
  For computing standard errors for predictions the parts of the coefficient matrix:
  
  C11 equal X'ViX
  
  C12 equal 0 - [(PX)'Vi(GZ')'] 
  
  C22 equal PEV equal G - [Z'G[Vi - (VX*tXVXVX)]GZ']
  
  and are summarized as:
  
  standard.errors equal sqrt( rowSums(XC11X') + rowSums(2*(XC12Z')) + rowSums(ZC22Z') )
  
 }
\value{
  The data frame used for obtaining the predictions, jointly with the predicted values and the corresponding standard errors. 
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias \email{cova_ruber@live.com.mx}}
\seealso{\code{\link{predict}}, \code{\link{mmer}}}
\keyword{models}
\examples{
data(DT_yatesoats)
DT <- DT_yatesoats

m3 <- mmer(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
summary(m3)$varcomp

p0 <- predict.mmer(object=m3, classify = "N")
p0$pvals

p0 <- predict.mmer(object=m3, classify = "V")
p0$pvals

p0 <- predict.mmer(object=m3, classify = "B")
p0$pvals
}
