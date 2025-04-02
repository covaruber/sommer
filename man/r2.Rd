\name{r2}
\alias{r2}

\title{
Reliability
}
\description{
Calculates the reliability of BLUPs in a sommer model.
}
\usage{
r2(object, object2=NULL)
}

\arguments{
  \item{object}{
Model fitted with the mmes function.
}
\item{object2}{
An optional model identical to object in the first argument but fitted with the argument returnParam set to TRUE to access the relationship matrices from the fitted model.
}
}
\details{

The reliability method calculated is the classical animal model: R2=(G-PEV)/G

}
\value{
\describe{
\item{result}{a list with as many elements as random effects fitted containing reliabilities for individual BLUPs.}
}
}

\references{

Mrode, R. A. (2014). Linear models for the prediction of animal breeding values. Cabi.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}

\examples{
####=========================================####
#### Example population
####=========================================####
data(DT_example)
DT <- DT_example
head(DT)
ans1 <- mmes(Yield~Env, 
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
rel=r2(ans1)
}
\seealso{\code{\link{mmes}} -- the core function of the package }
