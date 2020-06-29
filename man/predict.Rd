\name{predict.mmer}
\alias{predict.mmer}
\title{Predict form a LMM fitted with mmer}
\description{\code{predict} method for class \code{"mmer"}.
}
\usage{
\method{predict}{mmer}(object, classify=NULL,
        RtermsToForce=NULL,
        FtermsToForce=NULL,
        \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmer"}}
 \item{classify}{an optional character string with the variables that define the margins of the multiway table to be predicted (see Details).}
 \item{RtermsToForce}{numeric vector to force the random effects to be used instead of the ones defined in the classify argument.}
\item{FtermsToForce}{numeric vector to force the fixed effects to be used instead of the ones defined in the classify argument.}
 \item{\dots}{Further arguments to be passed}
}
\details{
  This function allows to produce predictions, either specifying: (1) the data frame on which to obtain the predictions (argument newdata), or (2) those variables that define the margins of the multiway table to be predicted (argument classify). In the first case, all fixed and random components need to be present in the data frame. In the second case, predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model.
 }
\value{
  The data frame used for obtaining the predictions, jointly with the predicted values and the corresponding standard errors. 
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
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
p0$predictions

p0 <- predict.mmer(object=m3, classify = "V")
p0$predictions

p0 <- predict.mmer(object=m3, classify = "B")
p0$predictions
}
