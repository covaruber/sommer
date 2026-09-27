# Simultaneous autoregressive spatial covariance structure

Simultaneous autoregressive spatial covariance structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
sar(x, W, rho = 0.10, fixed = FALSE)
```

## Arguments

- x:

  Factor or design defining q spatial levels.

- W:

  Square spatial weights matrix aligned to the levels of x. Row/column
  names are used when available.

- rho:

  Starting SAR dependence parameter.

- fixed:

  Logical indicating whether rho is fixed.

## Details

The simultaneous autoregressive covariance is constructed from
\$\$B=I-\rho W,\qquad M=B^{-1}B^{-\mathsf T},\qquad K=M/M\_{11}.\$\$ For
a general real weights matrix, `rho` is restricted to the conservative
interval \\(-1/r(W),1/r(W))\\, where \\r(W)\\ is the spectral radius. A
bounded-logit working coordinate enforces the interval. The first
covariance derivative is supplied analytically to the generic
CovarianceFactor engine.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`car`](https://covaruber.github.io/sommer/reference/car.md),
[`maternm`](https://covaruber.github.io/sommer/reference/maternm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(sar(location, W), ism(genotype))
} # }
```
