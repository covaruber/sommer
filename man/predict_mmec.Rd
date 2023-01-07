\name{predict.mmec}
\alias{predict.mmec}
\title{Predict form of a LMM fitted with mmec}
\description{\code{predict} method for class \code{"mmec"}.
}
\usage{
\method{predict}{mmec}(object, Dtable=NULL, D, \dots)
}
\arguments{
 \item{object}{a mixed model of class \code{"mmec"}}
 \item{Dtable}{a table specifying the terms to be included or averaged.
 
  An "ignored" term means that the model matrix part for that fixed or random effect are not used in the BLUP and SE calculation.
 
 An "averaged" term means that the model matrices for that fixed or random effect is filled purely with the value 1/#levels for that term.
 
  An "include" term means that the model matrices for that fixed or random effect is filled with 1's for the positions where column names and row names match.
  
  An "include and average" term means that the model matrices for that fixed or random effect is filled with 1/#levels for the positions where column names and row names match.
 
 }
 \item{D}{a character string specifying the variable used to extract levels for the rows of the D matrix and its construction. Alternatively, the D matrix (of class dgCMatrix) specifying the matrix to be used for the predictions directly.}
 \item{\dots}{Further arguments to be passed.}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument classify). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 
  
  For predicted values the pertinent design matrices X and Z together with BLUEs (b) and BLUPs (u) are multiplied and added together.
  
  predicted.value equal Xb + Zu.1 + ... + Zu.n
  
  For computing standard errors for predictions the parts of the coefficient matrix:
  
  C11 equal (X'ViX)-
  
  C12 equal 0 - [(X'V-X)-X' V- GZ] 
  
  C22 equal PEV equal G - [Z'G[V- - (VX*tXVXVX)]GZ']
  
  and are summarized as:
  
  standard.errors equal sqrt( rowSums(XC11*X') + rowSums(2*(XC12*Z')) + rowSums(ZC22*Z') )
  
  when both fixed and random effects are present in the inclusion set. If only fixed and random effects are included, only the respective terms from the SE for fixed or random effects are calculated.
  
 }
\value{

  \item{pvals}{the table of predictions according to the specified arguments.}
  
  \item{hypertable}{the summary table specifying the ignored, include and average sets.}
  
  \item{model}{the mixed model used within predict.}
  
  \item{C11}{the inverse of the coefficient matrix corresponding to the fixed effects.}
  
  \item{C12}{the inverse of the coefficient matrix corresponding to the covariance between fixed and random effects.}
  
  \item{C22}{the inverse of the coefficient matrix corresponding to the random effects.}
  
  \item{Xextended}{the model matrix for fixed effects used to compute BLUPs and SEs.}
  
  \item{Zextended}{the model matrix for random effects used to compute BLUPs and SEs.}
  
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{predict}}, \code{\link{mmec}}}
\keyword{models}
\examples{

data(DT_yatesoats)
DT <- DT_yatesoats
m3 <- mmec(fixed=Y ~ V + N + V:N ,
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
