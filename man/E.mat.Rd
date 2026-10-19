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
If TRUE, Nishio and Satoh (2014) (orthogonal parameterization of Vitezica et al. 2013), otherwise Su et al. (2012) (see Details in the D.mat help page).
}
\item{type}{
An argument specifying the type of epistatic relationship matrix desired. The default is the second order epistasis (additive x additive) type="A#A". Other options are additive x dominant (type="A#D"), or dominant by dominant (type="D#D").  
}
\item{min.MAF}{
Minimum minor allele frequency, a finite value between zero and 0.5. Markers at or below this threshold and all constant markers are removed. The default is 0.02.
}
}

\details{
The result is the element-wise (Hadamard) product \eqn{E=A\circ A}, \eqn{E=A\circ D}, or \eqn{E=D\circ D}. Mean imputation, allele-frequency estimation and marker filtering are shared between the component matrices. Their definitions and input requirements are the same as in \code{\link{A.mat}} and \code{\link{D.mat}}. As in \code{A.mat}, the correction from inadvertent median imputation to mean imputation can change results for missing genotypes.
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
