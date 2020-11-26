\name{D.mat}
\alias{D.mat}

\title{
Dominance relationship matrix
}
\description{
C++ implementation of the dominance matrix. Calculates the realized dominance relationship matrix. Can help to increase the prediction accuracy when 2 conditions are met; 1) The trait has intermediate to high heritability, 2) The population contains a big number of individuals that are half or full sibs (HS & FS).
}
\usage{
D.mat(X,nishio=TRUE,min.MAF=0,return.imputed=FALSE)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{nishio}{
If TRUE Nishio ans Satoh. (2014), otherwise Su et al. (2012). See references.
}
\item{min.MAF}{
Minimum minor allele frequency. The D matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}


}
\details{

The additive marker coefficients will be used to compute dominance coefficients as: Xd = 1-abs(X) for diploids.

For nishio method: the marker matrix is centered by subtracting column means \eqn{M= Xd - ms} where ms is the column means. Then \eqn{A=M M'/c}, where \eqn{c = 2 \sum_k {p_k (1-p_k)}}.

For su method: the marker matrix is normalized by subtracting row means \eqn{M= Xd - 2pq} where 2pq is the product of allele frequencies times 2. Then \eqn{A=M M'/c}, where \eqn{c = 2 \sum_k {2pq_k (1-2pq_k)}}.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$D}{the D matrix}
\item{$imputed}{the imputed marker matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. Plos One 9(1), doi:10.1371/journal.pone.0085792

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

}

\examples{
####=========================================####
#### EXAMPLE 1
####=========================================####
####random population of 200 lines with 1000 markers
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}

D <- D.mat(X)

}
\seealso{The core functions of the package \code{\link{mmer}}}
