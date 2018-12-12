\name{D.mat}
\alias{D.mat}

\title{
Dominance relationship matrix
}
\description{
Calculates the realized dominance relationship matrix. Can help to increase the prediction accuracy when 2 conditions are met; 1) The trait has intermediate to high heritability, 2) The population contains a big number of individuals that are half or full sibs (HS & FS).
}
\usage{
D.mat(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE, ploidy=2, return.Xd=FALSE, 
      method=3)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency. The D matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
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
\item{ploidy}{
The ploidy of the organism. The default is 2 which means diploid but higher ploidy levels are supported.}
\item{return.Xd}{
If TRUE, the function returns the marker matrix converted to a dominant matrix instead of the dominance relationship matrix. The -1 and 1 are converted to 0, and 0's are converted to 1's which is a common incidence matrix for dominance. By default the markers that did not have heterozygote calls are not returned which would mean returning columns of only zeros putting fresh users in a mathematical problem when trying to estimate the variance component. The default for this argument is FALSE, so the function returns a dominant relationship matrix instead of a marker matrix.}
\item{method}{
Method 1 is Endelman and Jannink (2012) used for A.mat, method 2 is  Nishio ans Satoh. (2014), method 3 is Su et al. (2012). See references.}

}
\details{

The additive marker coefficients will be used to compute dominance coefficients as: 1-abs(X) for diploids.

At high marker density, the relationship matrix is estimated as \eqn{D=W W'/c}, where \eqn{W_{ik} = 1 - \| X_{ik} \|}. By using a normalization constant of \eqn{c = 2 \sum_k {p_k q_k (1-p_k q_k)}}.

The EM imputation algorithm is based on the multivariate normal distribution and was designed for use with GBS (genotyping-by-sequencing) markers, which tend to be high density but with lots of missing data.  Details are given in Poland et al. (2012). The EM algorithm stops at iteration \eqn{t} when the RMS error = \eqn{n^{-1} \|A_{t} - A_{t-1}\|_2} < tol.  

At low marker density (m < n), shrinkage estimation can improve the estimate of the relationship matrix and the accuracy of GEBVs for lines with low accuracy phenotypes (Endelman and Jannink 2012).  The shrinkage intensity ranges from 0 (no shrinkage, same estimator as high density formula) to 1 (completely shrunk to \eqn{(1+f)I}). The shrinkage intensity is chosen to minimize the expected mean-squared error and printed to the screen as output. 

The shrinkage and EM options are designed for opposite scenarios (low vs. high density) and cannot be used simultaneously. 

When the EM algorithm is used, the imputed alleles can lie outside the interval [-1,1].  Polymorphic markers that do not meet the min.MAF and max.missing criteria are not imputed.
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

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

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
