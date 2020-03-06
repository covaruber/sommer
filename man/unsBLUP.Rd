\name{unsBLUP}
\alias{unsBLUP}
\title{unsBLUP}
\description{
  \code{unsBLUP} takes a list of BLUPs from the U component of an mmer model fitted as unstructured to add the main variance BLUP and the covariance BLUPs to provide the right BLUPs for the main variances.
}
\usage{
  unsBLUP(blups)
}

\arguments{
  \item{blups}{a list of BLUPs from the U component of the model.}
}
\value{
  \describe{
    
    \item{$res}{a list of BLUPs from the U component of the model adjusting the main variance BLUPs by adding the covariance BLUPs .}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_example)
DT <- DT_example

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmer(fixed=Yield~Env-1,
             random=~ vs(us(Env),Name),
             rcov=~vs(us(Env),units), 
             data=DT)#, return.param = T)

u <- unsBLUP(ans3$U[1:6])
# adjusted BLUP vs non-adjusted BLUP
plot(u$`CA.2011:Name`$Yield, ans3$U$`CA.2011:Name`$Yield)

}
\seealso{The function \code{\link{vs}} to know how to use \code{unsBLUP} in the \code{\link{mmer}} solver. }