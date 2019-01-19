\name{fcm}
\alias{fcm}
\title{fixed effect constraint indication matrix}
\description{
  \code{fcm} creates a matrix with the correct number of columns to specify a constraint in the fixed effects using the Gtc argument of the \code{\link{vs}} function.
}
\usage{
  fcm(x, reps=NULL)
}

\arguments{
  \item{x}{vector of 1's and 0's corresponding to the traits for which this fixed effect should be fitted. For example, for a trivariate model if the fixed effect "x" wants to be fitted only for trait 1 and 2 but not for the 3rd trait then you would use fcm(c(1,1,0)) in the Gtc argument of the vs() function.}
  \item{reps}{integer specifying the number of times the matrix should be repeated in a list format to provide easily the constraints in complex models that use the ds(), us() or cs() structures.}
}
\value{
  \describe{
    
    \item{$res}{ a matrix or a list of matrices with the constraints to be provided in the Gtc argument of the \code{\link{vs}} function.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
fcm(c(1,1,0))
fcm(c(0,1,1))
fcm(c(1,1,1))

fcm(c(1,1,1),2)

# ## model with Env estimated for both traits
# ans4 <- mmer(cbind(Yield, Weight) ~ Env,
#               random= ~ vs(Name) + vs(Env:Name),
#               rcov= ~ vs(units),
#               data=DT)
# summary(ans4)$betas
# ## model with Env only estimated for Yield
# ans4b <- mmer(cbind(Yield, Weight) ~ vs(Env, Gtc=fcm(c(1,0))),
#              random= ~ vs(Name) + vs(Env:Name),
#              rcov= ~ vs(units),
#              data=DT)
# summary(ans4b)$betas

}
\seealso{The function \code{\link{vs}} to know how to use \code{fcm} in the \code{\link{mmer}} solver.}