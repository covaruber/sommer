\name{predict.mmes}
\alias{predict.mmes}
\title{Predict form of a LMM fitted with mmes}
\description{\code{predict} method for class \code{"mmes"}.
}
\usage{
\method{predict}{mmes}(object, Dtable=NULL, D, \dots)
}
\arguments{
 \item{object}{a mixed model of class \code{"mmes"}}
 \item{Dtable}{a table specifying the terms to be included or averaged.
 
  An "include" term means that the model matrices for that fixed or random effect is filled with 1's for the positions where column names and row names match.
  
  An "include and average" term means that the model matrices for that fixed or random effect is filled with 1/#1's in that row.
  
  An "average" term alone means that all rows for such fixed or random effect will be filled with 1/#levels in the effect.
  
  If a term is not considered "include" or "average" is then totally ignored in the BLUP and SE calculation. 
  
  The default rule to invoke when the user doesn't provide the Dtable is to include and average all terms that match the argument D.
 
 }
 \item{D}{a character string specifying the variable used to extract levels for the rows of the D matrix and its construction. Alternatively, the D matrix (of class dgCMatrix) specifying the matrix to be used for the predictions directly.}
 \item{\dots}{Further arguments to be passed.}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument D). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 
  
  For predicted values the pertinent design matrices X and Z together with BLUEs (b) and BLUPs (u) are multiplied and added together.
  
  predicted.value equal Xb + Zu.1 + ... + Zu.n
  
  For computing standard errors for predictions the parts of the coefficient matrix:
  
  C11 equal (X.t() V.inv() X).inv()
  
  C12 equal 0 - [(X.t() V.inv() X).inv() X.t() V.inv() G Z] 
  
  C22 equal PEV equal G - [Z.t() G[V.inv() - (V.inv() X X.t() V.inv() X V.inv() X)]G Z.t()]
  
  In practive C equals ( W.t() V.inv() W ).inv()
  
  when both fixed and random effects are present in the inclusion set. If only fixed and random effects are included, only the respective terms from the SE for fixed or random effects are calculated.
  
 }
\value{

  \item{pvals}{the table of predictions according to the specified arguments.}
  
  \item{vcov}{the variance covariance for the predictions.}
  
  \item{D}{the model matrix for predictions as defined in Welham et al.(2004).}
  
  \item{Dtable}{the table specifying the terms to include and terms to be averaged.}
  
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias-Pazaran}
\seealso{\code{\link{predict}}, \code{\link{mmes}}}
\keyword{models}
\examples{

data(DT_yatesoats)
DT <- DT_yatesoats
m3 <- mmes(fixed=Y ~ V + N + V:N ,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)

#############################
## predict means for nitrogen
#############################
Dt <- m3$Dtable; Dt
# first fixed effect just average
Dt[1,"average"] = TRUE
# second fixed effect include
Dt[2,"include"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N")
pp$pvals

#############################
## predict means for variety
#############################

Dt <- m3$Dtable; Dt
# first fixed effect include
Dt[1,"include"] = TRUE
# second fixed effect just average
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="V")
pp$pvals

#############################
## predict means for nitrogen:variety
#############################
# prediction matrix D based on (equivalent to classify in asreml)
Dt <- m3$Dtable; Dt
# first fixed effect include and average
Dt[1,"include"] = TRUE
Dt[1,"average"] = TRUE
# second fixed effect include and average
Dt[2,"include"] = TRUE
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N:V")
pp$pvals

}
