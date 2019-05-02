\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix. Is a wrapper of the A.mat function from the rrBLUP package published by Endelman (2011). If using this function please cite Endelman (2011). 
}
\usage{
A.mat(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE, ploidy=2)
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
\item{ploidy}{
The ploidy of the organism. The default is 2 which means diploid but higher ploidy levels are supported. 
}
}
\details{
At high marker density, the relationship matrix is estimated as \eqn{A=W W'/c}, where \eqn{W_{ik} = X_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  By using a normalization constant of \eqn{c = 2 \sum_k {p_k (1-p_k)}}, the mean of the diagonal elements is \eqn{1 + f} (Endelman and Jannink 2012).

The EM imputation algorithm is based on the multivariate normal distribution and was designed for use with GBS (genotyping-by-sequencing) markers, which tend to be high density but with lots of missing data.  Details are given in Poland et al. (2012). The EM algorithm stops at iteration \eqn{t} when the RMS error = \eqn{n^{-1} \|A_{t} - A_{t-1}\|_2} < tol.  

At low marker density (m < n), shrinkage estimation can improve the estimate of the relationship matrix and the accuracy of GEBVs for lines with low accuracy phenotypes (Endelman and Jannink 2012).  The shrinkage intensity ranges from 0 (no shrinkage, same estimator as high density formula) to 1 (completely shrunk to \eqn{(1+f)I}). The shrinkage intensity is chosen to minimize the expected mean-squared error and printed to the screen as output. 

The shrinkage and EM options are designed for opposite scenarios (low vs. high density) and cannot be used simultaneously. 

When the EM algorithm is used, the imputed alleles can lie outside the interval [-1,1].  Polymorphic markers that do not meet the min.MAF and max.missing criteria are not imputed.
}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$A}{the A matrix}
\item{$imputed}{the imputed marker matrix}
}
}

\references{

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Poland, J., J. Endelman et al. 2012. Genomic selection in wheat breeding using genotyping-by-sequencing. Plant Genome 5:103-113. doi: 10.3835/plantgenome2012.06.0006
}

\examples{
####=========================================####
#### random population of 200 lines with 1000 markers
####=========================================####
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

A <- A.mat(X)

####=========================================####
#### take a look at the Genomic relationship matrix 
#### (just a small part)
####=========================================####
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(A[1:15,1:15], col = colfunc(100),Colv = "Rowv")
# str(hv)
}
\seealso{\code{\link{mmer}} -- the core function of the package }
