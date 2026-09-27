# Dominance relationship matrix

C++ implementation of the dominance matrix. Calculates the realized
dominance relationship matrix. Can help to increase the prediction
accuracy when 2 conditions are met; 1) The trait has intermediate to
high heritability, 2) The population contains a big number of
individuals that are half or full sibs (HS & FS).

## Usage

``` r
D.mat(X,nishio=TRUE,min.MAF=0,return.imputed=FALSE)
```

## Arguments

- X:

  Matrix (\\n \times m\\) of unphased genotypes for \\n\\ lines and
  \\m\\ biallelic markers, coded as {-1,0,1}. Fractional (imputed) and
  missing values (NA) are allowed.

- nishio:

  If TRUE Nishio ans Satoh. (2014), otherwise Su et al. (2012). See
  references.

- min.MAF:

  Minimum minor allele frequency. The D matrix is not sensitive to rare
  alleles, so by default only monomorphic markers are removed.

- return.imputed:

  When TRUE, the imputed marker matrix is returned.

## Details

The additive marker coefficients will be used to compute dominance
coefficients as: Xd = 1-abs(X) for diploids.

For nishio method: the marker matrix is centered by subtracting column
means \\M= Xd - ms\\ where ms is the column means. Then \\A=M M'/c\\,
where \\c = 2 \sum_k {p_k (1-p_k)}\\.

For su method: the marker matrix is normalized by subtracting row means
\\M= Xd - 2pq\\ where 2pq is the product of allele frequencies times 2.
Then \\A=M M'/c\\, where \\c = 2 \sum_k {2pq_k (1-2pq_k)}\\.

## Value

If return.imputed = FALSE, the \\n \times n\\ additive relationship
matrix is returned.

If return.imputed = TRUE, the function returns a list containing

- \$D:

  the D matrix

- \$imputed:

  the imputed marker matrix

## References

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6):
doi:10.1371/journal.pone.0156744

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic
BLUP Method for Genomic Evaluation. Plos One 9(1),
doi:10.1371/journal.pone.0085792

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating
Additive and Non-Additive Genetic Variances and Predicting Genetic
Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers.
PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

## Examples

``` r
####=========================================####
#### EXAMPLE 1
####=========================================####
####random population of 200 lines with 1000 markers
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}

D <- D.mat(X)
```
