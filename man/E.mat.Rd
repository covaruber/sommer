\name{E.mat}
\alias{E.mat}

\title{
Epistatic relationship matrix
}
\description{
Calculates the realized epistatic relationship matrix of second order (additive x additive, additive x dominance, or dominance x dominance) using hadamard products with the C++ Armadillo library.
}
\usage{
E.mat(X,nishio=TRUE,type="A#A",min.MAF=0.02)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{nishio}{
If TRUE Nishio ans Satoh. (2014), otherwise Su et al. (2012) (see Details in the D.mat help page).
}
\item{type}{
An argument specifying the type of epistatic relationship matrix desired. The default is the second order epistasis (additive x additive) type="A#A". Other options are additive x dominant (type="A#D"), or dominant by dominant (type="D#D").  
}
\item{min.MAF}{
Minimum minor allele frequency. The A matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
}

\details{
it is computed as the Hadamard product of the epistatic relationship matrix; E=A#A, E=A#D, E=D#D. 
}
\value{
The epistatic relationship matrix is returned.

}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. Plos One 9(1), doi:10.1371/journal.pone.0085792

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

}

\examples{
####=========================================####
####random population of 200 lines with 1000 markers
####=========================================####
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}

E <- E.mat(X, type="A#A") 
# if heterozygote markers are present can be used "A#D" or "D#D"
}
\seealso{The core functions of the package \code{\link{mmer}}}
