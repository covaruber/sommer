\name{build.HMM}
\alias{build.HMM}

\title{
Build a hybrid marker matrix using parental genotypes
}
\description{
Uses the 2 marker matrices from both sets of parents and creates all possible combinations unless the user specifies which hybrid genotypes to build (custom.hyb argument). It returns the marker matrix for hybrids coded as additive (-1,0,1; homo,het,homo) and dominance (0,1,0; homo,het,homo).
}
\usage{
build.HMM(M1,M2, custom.hyb=NULL, return.combos.only=FALSE)
}

\arguments{
  \item{M1}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are not allowed.
}

  \item{M2}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are not allowed.
}

\item{custom.hyb}{
A data frame with columns 'Var1' 'Var2', 'hybrid' which specifies which hybrids should be built using the M1 and M2 matrices provided.
}
\item{return.combos.only}{
A TRUE/FALSE statement inicating if the function should skip building the geotype matrix for hybrids and only return the data frame with all possible combinations to be build. In case the user wants to subset the hybrids before building the marker matrix.
}

}
\details{

It returns the marker matrix for hybrids coded as additive (-1,0,1; homo,het,homo) and dominance (0,1,0; homo,het,homo).

}
\value{
It returns the marker matrix for hybrids coded as additive (-1,0,1; homo,het,homo) and dominance (0,1,0; homo,het,homo).

\describe{
\item{$HMM.add}{marker matrix for hybrids coded as additive (-1,0,1; homo,het,homo)}
\item{$HMM.dom}{marker matrix for hybrids coded as dominance (0,1,0; homo,het,homo)}
\item{$data.used}{the data frame used to build the hybrid genotypes}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. Plos One 9(1), doi:10.1371/journal.pone.0085792

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

}

\examples{
####=========================================####
#### use Technow data as example
####=========================================####
data(DT_technow)
M.flint <- Mf # Marker matrix Flint
M.dent <- Md # Marker matrix Dent

## first get all possible hybrids
res1 <- build.HMM(M.dent, M.flint, 
                  return.combos.only = TRUE)
head(res1$data.used)

## build the marker matrix for the first 50 hybrids
res2 <- build.HMM(M.dent, M.flint,
                 custom.hyb = res1$data.used[1:50,]
                 )
res2$HMM.add[1:5,1:5]
res2$HMM.dom[1:5,1:5]

## now you can use the A.mat(), D.mat() and E.mat() functions
# M <- res2$HMM.add
# A <- A.mat(M)
# D <- D.mat(M)

}
\seealso{\code{\link{mmer}}-- the core function of the package }
