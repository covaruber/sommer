\name{AR1}
\alias{AR1}

\title{
Autocorrelation matrix of order 1.
}
\description{
Creates an autocorrelation matrix of order one with parameters specified. 
}
\usage{
AR1(x,rho=0.25)
}

\arguments{
  \item{x}{vector of the variable to define the factor levels for the AR1 covariance structure.}
  \item{rho}{rho value for the matrix.}
}
\details{
Specially useful for constructing covariance structures for rows and ranges to capture better the spatial variation trends in the field. The rho value is assumed fixed and values of the variance component will be optimized through REML.

}
\value{
If everything is defined correctly the function returns:
\describe{
\item{$nn}{the correlation matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}

\examples{
x <- 1:4
R1 <- AR1(x,rho=.25)
image(R1)

data(DT_sleepstudy)
DT <- DT_sleepstudy
head(DT)

# define the correlation between Days
D = with(DT, AR1(Days, rho=0.5))
subs = unique(DT$Subject)
# define the correlation between Subjects
S = diag(length(subs))
rownames(S) <- colnames(S) <- subs
# make the kronecker product
DS = kronecker(D,S, make.dimnames = TRUE)
# form the covariance matrix between units 
# this is assumes correlation between timepoints
DT$ds <- paste(DT$Days, DT$Subject, sep=":")
DS <- DS[DT$ds, DT$ds]
colnames(DS) <- rownames(DS) <- paste0("u",1:nrow(DS))
# fit the residual model
head(DT)

fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject, 
            rcov = ~vsm(ism(units), Gu=DS), # equivalent to Subject:ar1(Days)
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

# the matrix D can take any form: AR1, ARMA, or a custom correlation matrix

}
