\name{simGECorMat}
\alias{simGECorMat}

\title{
Create a GE correlation matrix for simulation purposes. 
}
\description{
Makes a simple correlation matrix based on the number of environments and megaenvironments desired. 
}
\usage{
simGECorMat(nEnv,nMegaEnv,mu=0.7, v=0.2, mu2=0, v2=0.3)
}

\arguments{
  \item{nEnv}{
Number of environments to simulate. Needs to be divisible by the nMegaEnv argument.
}
\item{nMegaEnv}{
Number of megaenvironments to simulate.
}
\item{mu}{
Mean value of the genetic correlation within megaenvironments.
}
\item{v}{
variance in the genetic correlation within megaenvironments.
}
\item{mu2}{
Mean value of the genetic correlation between megaenvironments.
}
\item{v2}{
variance in the genetic correlation between megaenvironments.
}
}
\details{
Simple simulation of a correlation matrix for environments and megaenvironments.
}
\value{
G the correlation matrix

\describe{
\item{$G}{the correlation matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}

\examples{
simGECorMat(9,3)
}
\seealso{\code{\link{mmer}} -- the core function of the package }
