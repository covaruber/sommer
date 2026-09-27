# Variance-structure model with arbitrary Kronecker covariance factors

`vsm` constructs random-effect or residual covariance structures for
`mmes`. Covariance-shaping terms are represented by CovarianceFactor v2
descriptors and may be combined in an arbitrary-depth Kronecker product.
A single product-level variance, \\\sigma^2\\, is owned by `vsm`, which
avoids scale confounding among the individual covariance factors.

## Usage

``` r
vsm(..., Gu = NULL, sigma2 = NULL, fixedSigma2 = FALSE,
  rotation = FALSE, isFixed = FALSE, verbose = TRUE)
```

## Arguments

- ...:

  One or more covariance-constructor terms. The last term supplies the
  main-effect incidence matrix. All preceding terms are
  covariance-shaping factors. Examples include `dsm(environment)`,
  `ar1m(row)`, and `ism(genotype)`.

- Gu:

  Optional known precision matrix for the levels of the final
  main-effect term. For `mmes`, `Gu` must have row and column names
  matching the main-effect levels and `attr(Gu, "inverse")` must be
  `TRUE`. If omitted, an identity precision matrix is used.

- sigma2:

  Positive starting value for the single overall covariance scale of
  this `vsm` term. If `NULL` (the default), `mmes` replaces it with a
  data-driven starting value computed from the response and fixed
  effects; supplying a value here always takes precedence and is never
  overridden.

- fixedSigma2:

  Logical indicating whether the product-level \\\sigma^2\\ is fixed at
  its starting value.

- rotation:

  Logical indicating whether the supplied `Gu` precision matrix should
  be eigen-decomposed for Lee–van der Werf rotation in `mmes`. Rotation
  is opt-in and currently requires a balanced Gaussian random-effect
  term.

- isFixed:

  Logical retained for compatibility. If `TRUE`, return the combined
  design matrix instead of the structured `vsm` object.

- verbose:

  Logical controlling messages, including notification when levels
  present in `Gu` are appended to the model matrix.

## Details

Let \\K_1,\ldots,K_m\\ be the dimensionless covariance shapes supplied
by the covariance constructors. The covariance represented by one `vsm`
term is \$\$\Sigma = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes
K_m).\$\$ For a random effect with known relationship precision `Gu`,
the factor product determines the covariance among the crossed
covariance coordinates and `Gu` determines the relationship structure
among the levels of the final main effect.

The final argument in `...` is the main-effect incidence term. Every
term before it is a covariance-shaping factor. Thus
`vsm(dsm(location), ar1m(row), ism(genotype), Gu=Ainv)` represents a
location-by-row covariance shape crossed with genotype. Kronecker
ordering is left to right: earlier factors are the slow index and later
factors are the fast index.

All optimizer coordinates are stored on unconstrained working scales.
Each CovarianceFactor descriptor also carries its evaluator, derivative
strategy, reporting transformation, and trust-region caps. Built-in
high-use structures may use native C++ covariance primitives, while
newer or user-defined structures can use generic R callbacks without
requiring changes to `ai_mme_sp2`.

For residual covariance models, the covariance-shaping factors must
define exactly one product coordinate for each observation. `mmes`
combines that local coordinate with an independently constructed
residual block index, so a large observation-level identity design is
not materialized.

When `rotation=TRUE`, `Gu` is decomposed as \$\$Gu=U\Lambda U',\$\$
where \\\Lambda\\ is diagonal. The aligned original precision, diagonal
precision, and eigenvectors are retained in the returned descriptor. The
`mmes` front end validates the observation layout and applies the
engine-specific rotation after observation filtering.

## Value

A list containing the random-effect design blocks in `Z`, the precision
matrix in `Gu`, the flattened Kronecker descriptor in `covStruct`, the
product design, and, for residual models, the local
covariance-coordinate index. When rotation is requested, `GuRot`
contains the diagonal precision and `rotation` contains its eigen
metadata. The `covStruct` object contains descriptor version 2
CovarianceFactor objects.

## See also

[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md),
[`ism`](https://covaruber.github.io/sommer/reference/ism.md),
[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md),
[`atm`](https://covaruber.github.io/sommer/reference/atm.md),
[`usm`](https://covaruber.github.io/sommer/reference/usm.md),
[`csm`](https://covaruber.github.io/sommer/reference/csm.md),
[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`ar2m`](https://covaruber.github.io/sommer/reference/ar2m.md),
[`ar3m`](https://covaruber.github.io/sommer/reference/ar3m.md),
[`mam`](https://covaruber.github.io/sommer/reference/mam.md),
[`corgm`](https://covaruber.github.io/sommer/reference/corgm.md),
[`fam`](https://covaruber.github.io/sommer/reference/fam.md),
[`antem`](https://covaruber.github.io/sommer/reference/antem.md),
[`rrcm`](https://covaruber.github.io/sommer/reference/rrcm.md),
[`maternm`](https://covaruber.github.io/sommer/reference/maternm.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`sar`](https://covaruber.github.io/sommer/reference/sar.md),
[`car`](https://covaruber.github.io/sommer/reference/car.md), and
[`ownm`](https://covaruber.github.io/sommer/reference/ownm.md).

## Examples

``` r
####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)
#>                   Name     Env Loc Year     Block Yield    Weight
#> 33  Manistee(MSL292-A) CA.2013  CA 2013 CA.2013.1     4 -1.904711
#> 65          CO02024-9W CA.2013  CA 2013 CA.2013.1     5 -1.446958
#> 66  Manistee(MSL292-A) CA.2013  CA 2013 CA.2013.2     5 -1.516271
#> 67            MSL007-B CA.2011  CA 2011 CA.2011.2     5 -1.435510
#> 68           MSR169-8Y CA.2013  CA 2013 CA.2013.1     5 -1.469051
#> 103         AC05153-1W CA.2013  CA 2013 CA.2013.1     6 -1.307167


str(with(DT, vsm(dsm(Env),ism(Name))))
#> List of 8
#>  $ Z                 :List of 3
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:43] 177 179 162 181 118 132 178 184 119 157 ...
#>   .. .. ..@ p       : int [1:42] 0 2 2 4 6 6 6 8 10 10 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:43] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:77] 30 83 68 84 48 69 70 86 133 134 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 6 8 10 12 14 16 18 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:77] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:65] 156 180 31 85 32 49 87 135 5 12 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 4 6 8 10 12 14 16 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:65] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>  $ Gu                :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:41] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ p       : int [1:42] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 41 41
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. ..@ x       : num [1:41] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>   .. ..$ inverse: logi TRUE
#>  $ GuRot             : NULL
#>  $ rotation          : NULL
#>  $ covStruct         :List of 13
#>   ..$ type              : chr "kron"
#>   ..$ par               : Named num [1:3] -1.9 0 0
#>   .. ..- attr(*, "names")= chr [1:3] "sigma2" "dsm(Env):variance_ratio[CA.2012]" "dsm(Env):variance_ratio[CA.2013]"
#>   ..$ free              : logi [1:3] TRUE TRUE TRUE
#>   ..$ par_names         : chr [1:3] "sigma2" "dsm(Env):variance_ratio[CA.2012]" "dsm(Env):variance_ratio[CA.2013]"
#>   ..$ factors           :List of 1
#>   .. ..$ :List of 15
#>   .. .. ..$ dim                  : int 3
#>   .. .. ..$ levels               : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. .. ..$ par                  : num [1:2] 0 0
#>   .. .. ..$ free                 : logi [1:2] TRUE TRUE
#>   .. .. ..$ par_names            : chr [1:2] "variance_ratio[CA.2012]" "variance_ratio[CA.2013]"
#>   .. .. ..$ model                : chr "diag"
#>   .. .. ..$ evaluator            :List of 2
#>   .. .. .. ..$ backend: chr "native"
#>   .. .. .. ..$ op     : chr "diag"
#>   .. .. ..$ derivative           :List of 2
#>   .. .. .. ..$ backend: chr "native"
#>   .. .. .. ..$ op     : chr "diag"
#>   .. .. ..$ report               :List of 4
#>   .. .. .. ..$ backend  : chr "builtin"
#>   .. .. .. ..$ transform: chr [1:2] "exp" "exp"
#>   .. .. .. ..$ lower    : num [1:2] NA NA
#>   .. .. .. ..$ upper    : num [1:2] NA NA
#>   .. .. ..$ native_report        :List of 2
#>   .. .. .. ..$ backend: chr "R"
#>   .. .. .. ..$ fun    :function (scale, par, factor, absorb_scale = TRUE)  
#>   .. .. ..$ trust_cap            : num [1:2] 1 1
#>   .. .. ..$ structurally_diagonal: logi TRUE
#>   .. .. ..$ descriptor_version   : int 2
#>   .. .. ..$ par_start            : int 2
#>   .. .. ..$ par_end              : int 3
#>   .. .. ..- attr(*, "class")= chr [1:2] "sommer_covfactor" "list"
#>   ..$ dim               : int 3
#>   ..$ levels            : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   ..$ scale_index       : int 1
#>   ..$ descriptor_version: int 2
#>   ..$ factor_interface  : chr "CovarianceFactor"
#>   ..$ parameterization  : chr "working"
#>   ..$ main_levels       : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   ..$ sigma2_is_default : logi TRUE
#>  $ residualLocalIndex: NULL
#>  $ productDesign     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:185] 3 7 52 76 97 99 111 118 119 122 ...
#>   .. ..@ p       : int [1:4] 0 43 120 185
#>   .. ..@ Dim     : int [1:2] 185 3
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. ..@ x       : num [1:185] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>  $ partitionsR       : NULL
str(with(DT, vsm(csm(Env),ism(Name))))
#> List of 8
#>  $ Z                 :List of 3
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:43] 177 179 162 181 118 132 178 184 119 157 ...
#>   .. .. ..@ p       : int [1:42] 0 2 2 4 6 6 6 8 10 10 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:43] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:77] 30 83 68 84 48 69 70 86 133 134 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 6 8 10 12 14 16 18 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:77] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:65] 156 180 31 85 32 49 87 135 5 12 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 4 6 8 10 12 14 16 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:65] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>  $ Gu                :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:41] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ p       : int [1:42] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 41 41
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. ..@ x       : num [1:41] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>   .. ..$ inverse: logi TRUE
#>  $ GuRot             : NULL
#>  $ rotation          : NULL
#>  $ covStruct         :List of 13
#>   ..$ type              : chr "kron"
#>   ..$ par               : Named num [1:2] -1.897 -0.405
#>   .. ..- attr(*, "names")= chr [1:2] "sigma2" "csm(Env):rho"
#>   ..$ free              : logi [1:2] TRUE TRUE
#>   ..$ par_names         : chr [1:2] "sigma2" "csm(Env):rho"
#>   ..$ factors           :List of 1
#>   .. ..$ :List of 16
#>   .. .. ..$ variance             : chr "homogeneous"
#>   .. .. ..$ dim                  : int 3
#>   .. .. ..$ levels               : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. .. ..$ par                  : Named num -0.405
#>   .. .. .. ..- attr(*, "names")= chr "eta_rho"
#>   .. .. ..$ free                 : logi TRUE
#>   .. .. ..$ par_names            : chr "rho"
#>   .. .. ..$ model                : chr "csm"
#>   .. .. ..$ evaluator            :List of 2
#>   .. .. .. ..$ backend: chr "native"
#>   .. .. .. ..$ op     : chr "cor_uniform"
#>   .. .. ..$ derivative           :List of 2
#>   .. .. .. ..$ backend : chr "numeric"
#>   .. .. .. ..$ rel_step: num 1e-06
#>   .. .. ..$ report               :List of 4
#>   .. .. .. ..$ backend  : chr "builtin"
#>   .. .. .. ..$ transform: chr "bounded_logit"
#>   .. .. .. ..$ lower    : num -0.5
#>   .. .. .. ..$ upper    : num 1
#>   .. .. ..$ native_report        :List of 2
#>   .. .. .. ..$ backend: chr "R"
#>   .. .. .. ..$ fun    :function (scale, par, factor, absorb_scale = TRUE)  
#>   .. .. ..$ trust_cap            : num 1
#>   .. .. ..$ structurally_diagonal: logi FALSE
#>   .. .. ..$ descriptor_version   : int 2
#>   .. .. ..$ par_start            : int 2
#>   .. .. ..$ par_end              : int 2
#>   .. .. ..- attr(*, "class")= chr [1:2] "sommer_covfactor" "list"
#>   ..$ dim               : int 3
#>   ..$ levels            : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   ..$ scale_index       : int 1
#>   ..$ descriptor_version: int 2
#>   ..$ factor_interface  : chr "CovarianceFactor"
#>   ..$ parameterization  : chr "working"
#>   ..$ main_levels       : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   ..$ sigma2_is_default : logi TRUE
#>  $ residualLocalIndex: NULL
#>  $ productDesign     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:185] 3 7 52 76 97 99 111 118 119 122 ...
#>   .. ..@ p       : int [1:4] 0 43 120 185
#>   .. ..@ Dim     : int [1:2] 185 3
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. ..@ x       : num [1:185] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>  $ partitionsR       : NULL
str(with(DT, vsm(toeplitzm(Env),ism(Name))))
#> List of 8
#>  $ Z                 :List of 3
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:43] 177 179 162 181 118 132 178 184 119 157 ...
#>   .. .. ..@ p       : int [1:42] 0 2 2 4 6 6 6 8 10 10 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:43] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:77] 30 83 68 84 48 69 70 86 133 134 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 6 8 10 12 14 16 18 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:77] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>   ..$ :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. .. ..@ i       : int [1:65] 156 180 31 85 32 49 87 135 5 12 ...
#>   .. .. ..@ p       : int [1:42] 0 2 4 4 6 8 10 12 14 16 ...
#>   .. .. ..@ Dim     : int [1:2] 185 41
#>   .. .. ..@ Dimnames:List of 2
#>   .. .. .. ..$ : NULL
#>   .. .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..@ x       : num [1:65] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. ..@ factors : list()
#>  $ Gu                :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:41] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ p       : int [1:42] 0 1 2 3 4 5 6 7 8 9 ...
#>   .. ..@ Dim     : int [1:2] 41 41
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. .. ..$ : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   .. ..@ x       : num [1:41] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>   .. ..$ inverse: logi TRUE
#>  $ GuRot             : NULL
#>  $ rotation          : NULL
#>  $ covStruct         :List of 13
#>   ..$ type              : chr "kron"
#>   ..$ par               : Named num [1:3] -1.9 0.1 0.1
#>   .. ..- attr(*, "names")= chr [1:3] "sigma2" "toeplitzm(Env):pacf[1]" "toeplitzm(Env):pacf[2]"
#>   ..$ free              : logi [1:3] TRUE TRUE TRUE
#>   ..$ par_names         : chr [1:3] "sigma2" "toeplitzm(Env):pacf[1]" "toeplitzm(Env):pacf[2]"
#>   ..$ factors           :List of 1
#>   .. ..$ :List of 15
#>   .. .. ..$ dim                  : int 3
#>   .. .. ..$ levels               : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. .. ..$ par                  : num [1:2] 0.1 0.1
#>   .. .. ..$ free                 : logi [1:2] TRUE TRUE
#>   .. .. ..$ par_names            : chr [1:2] "pacf[1]" "pacf[2]"
#>   .. .. ..$ evaluator            :List of 2
#>   .. .. .. ..$ backend: chr "R"
#>   .. .. .. ..$ fun    :function (par)  
#>   .. .. ..$ derivative           :List of 2
#>   .. .. .. ..$ backend : chr "numeric"
#>   .. .. .. ..$ rel_step: num 1e-06
#>   .. .. ..$ report               :List of 4
#>   .. .. .. ..$ backend  : chr "builtin"
#>   .. .. .. ..$ transform: chr [1:2] "tanh" "tanh"
#>   .. .. .. ..$ lower    : num [1:2] NA NA
#>   .. .. .. ..$ upper    : num [1:2] NA NA
#>   .. .. ..$ native_report        :List of 2
#>   .. .. .. ..$ backend: chr "R"
#>   .. .. .. ..$ fun    :function (scale, par, factor, absorb_scale = TRUE)  
#>   .. .. ..$ trust_cap            : num [1:2] 1 1
#>   .. .. ..$ structurally_diagonal: logi FALSE
#>   .. .. ..$ descriptor_version   : int 2
#>   .. .. ..$ model                : chr "toeplitz"
#>   .. .. ..$ par_start            : int 2
#>   .. .. ..$ par_end              : int 3
#>   .. .. ..- attr(*, "class")= chr [1:2] "sommer_covfactor" "list"
#>   ..$ dim               : int 3
#>   ..$ levels            : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   ..$ scale_index       : int 1
#>   ..$ descriptor_version: int 2
#>   ..$ factor_interface  : chr "CovarianceFactor"
#>   ..$ parameterization  : chr "working"
#>   ..$ main_levels       : chr [1:41] "A01143-3C" "AC00206-2W" "AC01151-5W" "AC03433-1W" ...
#>   ..$ sigma2_is_default : logi TRUE
#>  $ residualLocalIndex: NULL
#>  $ productDesign     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:185] 3 7 52 76 97 99 111 118 119 122 ...
#>   .. ..@ p       : int [1:4] 0 43 120 185
#>   .. ..@ Dim     : int [1:2] 185 3
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:3] "CA.2011" "CA.2012" "CA.2013"
#>   .. ..@ x       : num [1:185] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..@ factors : list()
#>  $ partitionsR       : NULL
```
