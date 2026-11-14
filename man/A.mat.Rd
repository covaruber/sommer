\name{A.mat}
\alias{A.mat}

\title{
Additive relationship matrix
}
\description{
Calculates the realized additive relationship matrix for diploid or polyploid
genotypes using the VanRaden method, with optional Endelman-Jannink
covariance shrinkage, or the average semivariance (ASV) method of
Feldmann et al. (2023). Shrinkage and ASV are implemented in C++.
}
\usage{
A.mat(X,min.MAF=0,return.imputed=FALSE,ploidy=NULL,shrink=FALSE,
      method=c("VanRaden","ASV"))
}

\arguments{
  \item{X}{
Matrix (\eqn{n \times m}) of unphased genotypes for \eqn{n} lines and \eqn{m} biallelic markers, 
coded as \{-1,0,1\} when \code{ploidy=NULL} (the default).
When \code{ploidy} is supplied, use allele dosages from zero to \code{ploidy},
as returned by \code{enhancer::atcg1234}. Fractional (imputed) and missing
values (NA) are allowed.
}
\item{min.MAF}{
Minimum minor allele frequency, a finite value between zero and 0.5. Markers at or below this threshold and all constant markers are removed.
}
\item{return.imputed}{
When TRUE, the imputed marker matrix is returned.  
}
\item{ploidy}{
Number of chromosome copies, a finite integer of at least two.
The default \code{NULL} preserves the legacy diploid \{-1,0,1\} coding.
An explicit value, including \code{ploidy=2}, selects dosage coding
from zero to \code{ploidy}. Ploidy is not inferred from observed dosage
levels because missing classes can underestimate it.
}
\item{shrink}{
Logical. When \code{TRUE}, estimate and apply the Endelman-Jannink (2012)
shrinkage intensity using a native C++ implementation of the EJ method
used by \code{rrBLUP::A.mat(shrink=TRUE)}. The default \code{FALSE}
preserves the unshrunken matrix. Regression-based shrinkage
(\code{method="REG"} in rrBLUP) is not implemented.
Available only with \code{method="VanRaden"}; combining
\code{shrink=TRUE} with \code{method="ASV"} is an error.
}
\item{method}{
\code{"VanRaden"} (the default) uses allele-frequency normalization.
\code{"ASV"} uses the average semivariance normalization in equation (2)
of Feldmann et al. (2023). Both methods use the same marker filters,
mean imputation and genotype coding.
}
}
\details{

Missing values are replaced by the column (marker) mean. At least two individuals and one variable marker surviving the frequency filter are required. Genotypes must lie in \eqn{[-1,1]} by default, or \eqn{[0,P]} for explicit ploidy \eqn{P}; a marker with no observed genotypes cannot be imputed.

For the VanRaden method, let \eqn{p_k} be the column mean of \eqn{(X+1)/2} for each retained marker and \eqn{W_{ik}=X_{ik}+1-2p_k}. Then \eqn{A=WW'/c}, where \eqn{c=2\sum_k p_k(1-p_k)}.

For explicit ploidy \eqn{P}, let \eqn{p_k} be the column mean of
\eqn{X/P} and \eqn{W_{ik}=X_{ik}-Pp_k}. The additive relationship matrix
is \eqn{A=WW'/c}, where \eqn{c=P\sum_k p_k(1-p_k)}. This uses the
binomial allele-dosage variance and does not model double reduction.
The returned imputed matrix retains the input dosage scale.

\code{method="ASV"} uses the column-centered retained marker matrix
\eqn{W=X-1_n\bar x'} and computes
\deqn{A_{\mathrm{ASV}} =
       \frac{WW'}{\mathrm{tr}(WW')/(n-1)}.}
The denominator is the sum of sample variances of the retained markers,
not the expected allele-dosage variance used by VanRaden. For
\eqn{P_n=I_n-1_n1_n'/n}, this gives
\eqn{\mathrm{tr}(P_n A_{\mathrm{ASV}})/(n-1)=1}; equivalently, the
average over unordered pairs of
\eqn{(A_{ii}+A_{jj}-2A_{ij})/2} equals one.
The mean diagonal is \eqn{(n-1)/n}, not one.
ASV changes the scale, not the relative relationships, of the unshrunken
VanRaden matrix computed from the same markers. It does not add
regularization or make a rank-deficient matrix invertible.
The normalization is invariant to a common shift or nonzero scale
change of all marker dosages, including equivalent diploid and
polyploid codings. Explicit ploidy is still needed to interpret allele
frequencies correctly for filtering. The paper describes diploid coding;
applying its same centered-crossproduct formula to polyploid dosage
data is a natural extension.

Mean imputation is explicit. Earlier implementations inadvertently used the imputation helper's default median despite reporting mean imputation; results with missing genotypes can therefore differ.

\code{shrink=TRUE} is intended for lower-density marker data. For the
\eqn{m} retained markers, let \eqn{\mu_i} be the row mean of \eqn{W},
\eqn{Z=W-\mu 1_m'}, \eqn{S=ZZ'/m}, and
\eqn{T=\mathrm{tr}(S)I_n/n}. Estimate
\deqn{\hat b^2 = \frac{1}{m^2}\sum_{k=1}^{m}
       \|z_k z_k' - S\|_F^2,\qquad
       \hat d^2 = \|S-T\|_F^2,}
where \eqn{z_k} is column \eqn{k} of \eqn{Z}. The intensity
\eqn{\delta} is \eqn{\hat b^2/\hat d^2} clipped to \eqn{[0,1]}.
If \eqn{\hat d^2=0}, the covariance already equals its target and the
intensity is reported as zero. The result is
\deqn{A_{\mathrm{EJ}} =
       \frac{(1-\delta)S+\delta T+\mu\mu'}{c/m}.}
This preserves the mean diagonal of the unshrunken relationship matrix.
It is not generally equivalent to blending the unshrunken matrix with
the identity: the row-mean contribution is retained without shrinkage.
The same estimator applies to explicit ploidy using its dosage-based
\eqn{W} and \eqn{c}.

Shrinkage is deterministic and uses the same mean imputation, strict
MAF cutoff and constant-marker removal as the unshrunken path.
The EJ calculation matches rrBLUP for the same retained and imputed
markers. rrBLUP has different default MAF/missingness filters and also
offers EM imputation; these options are not added here.

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

With \code{shrink=TRUE}, the returned matrix has a numeric
\code{"shrinkage"} attribute containing the estimated intensity in
\eqn{[0,1]}. No attribute is added when \code{shrink=FALSE}.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$A}{the A matrix}
\item{$X}{the imputed marker matrix}
}
}

\references{

Endelman, J.B., and J.-L. Jannink. 2012. Shrinkage estimation of the realized relationship matrix. G3:Genes, Genomes, Genetics. 2:1405-1413. doi: 10.1534/g3.112.004259

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Feldmann MJ, Covarrubias-Pazaran G, Piepho HP, Lipka AE (2023).
Complex traits and candidate genes: estimation of genetic variance
components across multiple genetic architectures. G3: Genes, Genomes,
Genetics 13(9):jkad148. doi:10.1093/g3journal/jkad148.
\url{https://pmc.ncbi.nlm.nih.gov/articles/PMC10468314/}

}

\examples{
####=========================================####
#### random population of 20 lines with 100 markers
####=========================================####
X <- matrix(rep(0,20*100),20,100)
for (i in 1:20) {
  X[i,] <- ifelse(runif(100)<0.5,-1,1)
}

A <- A.mat(X)

#### Average semivariance: mean pairwise semivariance equals one
A.ASV <- A.mat(X, method="ASV")
sum(diag(A.ASV)) / (nrow(A.ASV)-1)

#### Endelman-Jannink shrinkage for a lower-density marker subset
A.EJ <- A.mat(X[,1:100], shrink=TRUE)
attr(A.EJ, "shrinkage")

#### Tetraploid allele dosages from the bundled example data
data(DT_polyploid, package="enhancer")
numo <- enhancer::atcg1234(data=GT_polyploid[,1:100], ploidy=4,
                         silent=TRUE)
common <- intersect(DT_polyploid$Name, rownames(numo$M))
marks <- numo$M[common,,drop=FALSE]
A4 <- A.mat(marks, ploidy=4)
A4.EJ <- A.mat(marks, ploidy=4, shrink=TRUE)
A4.ASV <- A.mat(marks, ploidy=4, method="ASV")

####=========================================####
#### take a look at the Genomic relationship matrix 
#### (just a small part)
####=========================================####
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(A[1:15,1:15], col = colfunc(100),Colv = "Rowv")
# str(hv)
}
