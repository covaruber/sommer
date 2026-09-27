# Covariance Between Two Random Effects

`covm` combines two random-effect structures created with
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) into a
single correlated random-effect structure. It fits an unstructured \\2
\times 2\\ covariance matrix between the two random effects while
allowing them to have different incidence matrices. The two effects must
act on the same coefficient space and use the same relationship or
precision matrix. The resulting structure can be fitted with the
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md) solver.

## Usage

``` r
covm(ran1, ran2, thetaC = NULL, theta = NULL,
fixed = NULL, fixedSigma2 = FALSE,
labels = c("ran1", "ran2"), tol = 1e-10)
```

## Arguments

- ran1:

  A random-effect structure returned by
  [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) for the
  first random effect. Currently, `covm` requires a simple `vsm`
  structure with one covariance-product coordinate, for example
  `vsm(ism(focal))` or `vsm(ism(focal), Gu = Ai)`.

- ran2:

  A random-effect structure returned by
  [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) for the
  second random effect. It must have the same coefficient dimension and
  coefficient levels as `ran1`, although its incidence matrix may
  differ.

- thetaC:

  Deprecated and not supported by the current CovarianceFactor
  parameterization. The former cell-wise constraint matrix cannot in
  general be translated exactly into the normalized-Cholesky
  parameterization used by `covm`. Use `fixedSigma2` and `fixed`
  instead.

- theta:

  An optional symmetric \\2 \times 2\\ positive-definite matrix
  containing initial values for the variance-covariance matrix of the
  two random effects. Values are supplied directly on the natural
  variance-covariance scale.

  “\` The diagonal elements define the initial variances of the two
  random effects and the off-diagonal element defines their initial
  covariance.

  If `theta = NULL`, the default starting matrix is

      ```

      diag(2) * 0.15 + matrix(0.015, 2, 2)

  “\` giving initial variances of 0.165 and an initial covariance of
  0.015. “\`

- fixed:

  An optional logical vector of length two indicating whether the two
  normalized-Cholesky coordinates describing the covariance structure
  should be fixed at their starting values. The coordinates correspond
  to the off-diagonal Cholesky element and the second diagonal Cholesky
  element, respectively. The default is `c(FALSE, FALSE)`, so both are
  estimated.

- fixedSigma2:

  Logical value indicating whether the variance scale associated with
  the first random effect should be fixed at its starting value. The
  default is `FALSE`.

- labels:

  Character vector of length two giving labels for the two random
  effects in the covariance descriptor. The default is
  `c("ran1", "ran2")`. The labels must be different and non-empty.

- tol:

  Numerical tolerance used when checking positive definiteness and
  whether the relationship/precision matrices supplied by the two random
  effects are equal. The default is `1e-10`.

## Details

`covm` is intended for models in which two distinct random effects are
expected to be correlated. A common example is the joint modeling of
direct and indirect genetic effects.

If the coefficient vector for the first random effect is \\u_1\\ and
that for the second is \\u_2\\, `covm` represents their covariance
structure as

\$\$ \mathrm{Var}(\[u_1^\prime, u_2^\prime\]^\prime) = \sigma^2
K\_{\mathrm{effect}} \otimes A, \$\$

where \\A\\ represents the covariance structure associated with the
coefficient levels, and \\K\_{\mathrm{effect}}\\ is a \\2 \times 2\\
unstructured covariance shape describing the relationship between the
two random effects.

Internally, `covm` uses the same CovarianceFactor interface as
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md). The \\2
\times 2\\ effect covariance is represented using a normalized-Cholesky
parameterization. This guarantees a positive-definite covariance matrix
during optimization.

The variance of the first random effect is represented by the overall
scale parameter `sigma2`. The remaining variance ratio and covariance
are represented through the normalized-Cholesky factor.

The two random effects may have different incidence matrices. However,
they must refer to identical coefficient levels in the same order and
must use the same relationship/precision matrix.

When a relationship matrix is supplied through `Gu`, the Henderson
solver requires it to be an inverse/precision matrix with
`attr(Gu, "inverse") = TRUE`. For example:

    attr(Ai, "inverse") <- TRUE

    covm(
    vsm(ism(focal), Gu = Ai),
    vsm(ism(neighbour), Gu = Ai)
    )

The current implementation combines simple `vsm` random effects only.
Each input must contain one covariance-product coordinate. More complex
covariance-factor products should therefore not currently be supplied
independently inside `ran1` and `ran2`.

## Value

A list describing a correlated random-effect structure compatible with
the CovarianceFactor-v2 interface used by
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) and
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md). The main
elements are:

- `Z`:

  A list containing the two incidence matrices, one for each random
  effect.

- `Gu`:

  The common sparse relationship precision matrix associated with the
  coefficient levels. It is stored with `attr(Gu, "inverse") = TRUE` for
  use by the Henderson solver.

- `covStruct`:

  A CovarianceFactor-v2 descriptor containing the overall variance scale
  and the normalized-Cholesky representation of the \\2 \times 2\\
  unstructured covariance between the two random effects.

- `residualLocalIndex`:

  `NULL`, since `covm` describes a random-effect structure rather than a
  residual covariance structure.

- `productDesign`:

  `NULL` for the current implementation.

- `partitionsR`:

  `NULL` for the current implementation.

- `covm`:

  Logical value equal to `TRUE`, identifying the object as one
  constructed by `covm`.

- `covm_labels`:

  The two labels used for the correlated random effects.

## References

Covarrubias-Pazaran G (2016). Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6).
[doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744)

Bijma, P. (2014). The quantitative genetics of indirect genetic effects:
a selective review of modelling issues. Heredity, 112(1), 61–69.

## Author

Giovanny Covarrubias-Pazaran

## See also

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) for
constructing random-effect covariance structures and
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md) for
fitting mixed models using the Henderson solver.

## Examples

``` r
data(DT_ige, package = "enhancer")
DT <- DT_ige

## Correlate two random effects with identity relationship matrices

covRes <- with(
DT,
covm(
vsm(ism(focal)),
vsm(ism(neighbour))
)
)

str(covRes)
#> List of 8
#>  $ Z                 :List of 2
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:1000] 305 313 542 567 605 619 626 717 953 131 ...
#>   .. .. ..@ p       : int [1:99] 0 9 14 20 22 25 28 32 86 91 ...
#>   .. .. ..@ Dim     : int [1:2] 1000 98
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:98] "id_1019" "id_1101" "id_1137" "id_1268" ...
#>   .. .. ..@ x       : num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:1000] 350 567 591 593 595 717 755 859 105 200 ...
#>   .. .. ..@ p       : int [1:99] 0 8 14 21 24 28 30 32 84 88 ...
#>   .. .. ..@ Dim     : int [1:2] 1000 98
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:98] "id_1019" "id_1101" "id_1137" "id_1268" ...
#>   .. .. ..@ x       : num [1:1000] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>  $ Gu                :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:98] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ p       : int [1:99] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 98 98
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:98] "id_1019" "id_1101" "id_1137" "id_1268" ...
#>   .. .. ..$ : chr [1:98] "id_1019" "id_1101" "id_1137" "id_1268" ...
#>   .. ..@ x       : num [1:98] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>   .. ..$ inverse: logi TRUE
#>  $ covStruct         :List of 12
#>   ..$ type              : chr "kron"
#>   ..$ par               : Named num [1:3] -1.80181 0.09091 -0.00415
#>   .. ..- attr(*, "names")= chr [1:3] "sigma2" "chol[ran2,ran1]" "chol_diag[ran2]"
#>   ..$ free              : logi [1:3] TRUE TRUE TRUE
#>   ..$ par_names         : chr [1:3] "sigma2" "chol[ran2,ran1]" "chol_diag[ran2]"
#>   ..$ factors           :List of 1
#>   .. ..$ :List of 18
#>   .. .. ..$ dim                  : int 2
#>   .. .. ..$ levels               : chr [1:2] "ran1" "ran2"
#>   .. .. ..$ par                  : num [1:2] 0.09091 -0.00415
#>   .. .. ..$ free                 : logi [1:2] TRUE TRUE
#>   .. .. ..$ par_names            : chr [1:2] "chol[ran2,ran1]" "chol_diag[ran2]"
#>   .. .. ..$ us_row               : int [1:2] 2 2
#>   .. .. ..$ us_col               : int [1:2] 1 2
#>   .. .. ..$ us_diag              : logi [1:2] FALSE TRUE
#>   .. .. ..$ model                : chr "us"
#>   .. .. ..$ evaluator            :List of 2
#>   .. .. .. ..$ backend: chr "native"
#>   .. .. .. ..$ op     : chr "us"
#>   .. .. ..$ derivative           :List of 2
#>   .. .. .. ..$ backend: chr "native"
#>   .. .. .. ..$ op     : chr "us"
#>   .. .. ..$ report               :List of 4
#>   .. .. .. ..$ backend  : chr "builtin"
#>   .. .. .. ..$ transform: chr [1:2] "identity" "exp"
#>   .. .. .. ..$ lower    : num [1:2] NA NA
#>   .. .. .. ..$ upper    : num [1:2] NA NA
#>   .. .. ..$ native_report        :List of 2
#>   .. .. .. ..$ backend: chr "R"
#>   .. .. .. ..$ fun    :function (scale, par, factor, absorb_scale = TRUE)  
#>   .. .. ..$ trust_cap            : num [1:2] 2 1
#>   .. .. ..$ structurally_diagonal: logi FALSE
#>   .. .. ..$ descriptor_version   : int 2
#>   .. .. ..$ par_start            : int 2
#>   .. .. ..$ par_end              : int 3
#>   .. .. ..- attr(*, "class")= chr [1:2] "sommer_covfactor" "list"
#>   ..$ dim               : int 2
#>   ..$ levels            : chr [1:2] "ran1" "ran2"
#>   ..$ scale_index       : int 1
#>   ..$ descriptor_version: int 2
#>   ..$ factor_interface  : chr "CovarianceFactor"
#>   ..$ parameterization  : chr "working"
#>   ..$ main_levels       : chr [1:98] "id_1019" "id_1101" "id_1137" "id_1268" ...
#>  $ residualLocalIndex: NULL
#>  $ productDesign     : NULL
#>  $ partitionsR       : NULL
#>  $ covm              : logi TRUE
#>  $ covm_labels       : chr [1:2] "ran1" "ran2"

## Custom initial variance-covariance matrix

covRes2 <- with(
DT,
covm(
vsm(ism(focal)),
vsm(ism(neighbour)),
theta = matrix(
c(0.5, 0.1,
0.1, 0.3),
2, 2
)
)
)

## Relationship precision matrix

## Ai must have matching row and column names and be marked as inverse.

##

## attr(Ai, "inverse") <- TRUE

##

## covRes3 <- with(

## DT,

## covm(

## vsm(ism(focal), Gu = Ai),

## vsm(ism(neighbour), Gu = Ai)

## )

## )

## See the DT_ige help page for a complete model-fitting example.
```
