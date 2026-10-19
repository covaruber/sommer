\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix. Currently is the C++ implementation of van Raden (2008). 
}
\usage{
A.mat(X,min.MAF=0,return.imputed=FALSE)
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\}. Fractional (imputed) and missing values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency, a finite value between zero and 0.5. Markers at or below this threshold and all constant markers are removed.
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}
}
\details{

Missing values are replaced by the column (marker) mean. At least two individuals and one variable marker surviving the frequency filter are required. Genotypes must lie in \eqn{[-1,1]}; a marker with no observed genotypes cannot be imputed.

For the VanRaden method, let \eqn{p_k} be the column mean of \eqn{(X+1)/2} for each retained marker and \eqn{W_{ik}=X_{ik}+1-2p_k}. Then \eqn{A=WW'/c}, where \eqn{c=2\sum_k p_k(1-p_k)}.

Mean imputation is explicit. Earlier implementations inadvertently used the imputation helper's default median despite reporting mean imputation; results with missing genotypes can therefore differ.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$A}{the A matrix}
\item{$X}{the imputed marker matrix}
}
}

\references{

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

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
