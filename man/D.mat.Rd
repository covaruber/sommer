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
If TRUE, Nishio and Satoh (2014), which uses the orthogonal parameterization of Vitezica et al. (2013). Otherwise Su et al. (2012). See Details and references.
}
\item{min.MAF}{
Minimum minor allele frequency. The D matrix is not sensitive to rare alleles, so by default only monomorphic markers are removed.
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}


}
\details{

Missing values are replaced by the column (marker) mean before the matrix is computed. Markers with minor allele frequency less than or equal to \code{min.MAF} (and therefore all monomorphic markers) are removed. For each marker \eqn{k}, \eqn{p_k} is the frequency of the allele counted by the code 1, estimated as the column mean of \eqn{(X+1)/2}, and \eqn{q_k = 1-p_k}.

For the nishio method (Nishio and Satoh, 2014; Vitezica et al., 2013): genotypes are coded as \eqn{w_{ik} = -2q_k^2} for \eqn{x_{ik}=1}, \eqn{2p_kq_k} for \eqn{x_{ik}=0} and \eqn{-2p_k^2} for \eqn{x_{ik}=-1}, computed as \eqn{w_{ik} = -x_{ik}^2 + (p_k-q_k)x_{ik} + 2p_kq_k}. Fractional values are passed through the same formula; note that a mean-imputed value \eqn{x_{ik} = p_k-q_k} receives the heterozygote code \eqn{2p_kq_k}. Then \eqn{D = W W'/c}, where \eqn{c = \sum_k (2p_kq_k)^2}. Under Hardy-Weinberg equilibrium this coding is orthogonal to the additive coding.

For the Su method (Su et al., 2012): heterozygosity is coded as \eqn{h_{ik} = 1-|x_{ik}|} and centered as \eqn{m_{ik} = h_{ik} - 2p_kq_k}. Then \eqn{D = M M'/c}, where \eqn{c = \sum_k 2p_kq_k(1-2p_kq_k)}.

Versions of sommer earlier than 4.4.9 used a different formula for \code{nishio=TRUE}: \eqn{h_{ik}} centered by its column mean and divided by \eqn{\sum_k (2p_kq_k)^2}. That was not the Nishio and Satoh (2014) method, so results from those versions will differ.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} dominance relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$D}{the D matrix}
\item{$X}{the imputed marker matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. Plos One 9(1), doi:10.1371/journal.pone.0085792

Vitezica ZG, Varona L, Legarra A. 2013. On the Additive and Dominant Variance and Covariance of Individuals Within the Genomic Selection Scope. Genetics 195(4): 1223-1230. doi:10.1534/genetics.113.155176

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

