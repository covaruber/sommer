\name{build.HMM}
\alias{build.HMM}

\title{
Build a hybrid marker matrix using parental genotypes from inbred individuals
}
\description{
Uses the 2 marker matrices from both sets of inbred or partially inbred parents and creates all possible combinations unless the user specifies which hybrid genotypes to build (custom.hyb argument). It returns the additive and dominance marker matrices (-1,0,1; homo,het,homo in additive and 0,1,0; homo,het,homo for dominance).
}
\usage{
build.HMM(M1,M2, custom.hyb=NULL, return.combos.only=FALSE,separator=":")
}

\arguments{
  \item{M1}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} inbreds and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are not allowed.
}

  \item{M2}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} inbreds and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are not allowed.
}

\item{custom.hyb}{
A data frame with columns 'Var1' 'Var2', 'hybrid' which specifies which hybrids should be built using the M1 and M2 matrices provided.
}

\item{return.combos.only}{
A TRUE/FALSE statement inicating if the function should skip building the geotype matrix for hybrids and only return the data frame with all possible combinations to be build. In case the user wants to subset the hybrids before building the marker matrix.
}

\item{separator}{
Any desired character to be used when pasting the male and female columns to assign the name to the hybrids.
}

}
\details{

It returns the marker matrix for hybrids coded as additive (-1,0,1; homo,het,homo) and dominance (0,1,0; homo,het,homo). This function is deviced for building marker matrices for hybrids coming from inbreds. If the parents are close to inbred >F5 you can try deleting the heterozygote calls (0's) and imputing those cells with the most common genotype (1 or -1). The expectation is that for mostly inbred individuals this may not change drastically the result but will make the results more interpretable. For non-inbred parents (F1 to F3) the cross of an F1 x F1 has many possibilities and is not the intention of this function to build genotypes for heterzygote x heterozygote crosses.

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
DT <- DT_technow
Md <- (Md_technow * 2) - 1
Mf <- (Mf_technow * 2) - 1


## first get all possible hybrids
res1 <- build.HMM(Md, Mf, 
                  return.combos.only = TRUE)
head(res1$data.used)

## build the marker matrix for the first 50 hybrids
res2 <- build.HMM(Md, Mf,
                 custom.hyb = res1$data.used[1:50,]
                 )
res2$HMM.add[1:5,1:5]
res2$HMM.dom[1:5,1:5]

## now you can use the A.mat(), D.mat() and E.mat() functions
# M <- res2$HMM.add
# A <- A.mat(M)
# D <- D.mat(M)

}
