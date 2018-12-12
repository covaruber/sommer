\name{E.mat}
\alias{E.mat}

\title{
Epistatic relationship matrix
}
\description{
Calculates the realized epistatic relationship matrix of second order (additive x additive, additive x dominance, or dominance x dominance).
}
\usage{
E.mat(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE, type="A#A", ploidy=2)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency. The A matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
\item{max.missing}{
Maximum proportion of missing data; default removes completely missing markers.
}
\item{impute.method}{
There are two options. The default is "mean", which imputes with the mean for each marker.  The "EM" option imputes with an EM algorithm (see details).
}
\item{tol}{
Specifies the convergence criterion for the EM algorithm (see details). 
}
\item{n.core}{
Specifies the number of cores to use for parallel execution of the EM algorithm (use only at UNIX command line).
}
\item{shrink}{
Set shrink=TRUE to use the shrinkage estimation procedure (see Details).  
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}
\item{type}{
An argument specifying the type of epistatic relationship matrix desired. The default is the second order epistasis (additive x additive) type="A#A". Other options are additive x dominant (type="A#D"), or dominant by dominant (type="D#D").  
}
\item{ploidy}{
The ploidy of the organism. The default is 2 which means diploid but higher ploidy levels are supported.} 
}
\details{
it is computed as the Hadamard product of the epistatic relationship matrix (A); E=A#A.
}
\value{
If return.imputed = FALSE, the \eqn{n \times n} epistatic relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$E}{the E matrix}
\item{$imputed}{the imputed marker matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Poland, J., J. Endelman et al. 2012. Genomic selection in wheat breeding using genotyping-by-sequencing. Plant Genome 5:103-113. doi: 10.3835/plantgenome2012.06.0006
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
