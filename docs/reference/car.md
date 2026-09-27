# Proper conditional autoregressive spatial covariance structure

Proper conditional autoregressive spatial covariance structure. The
returned dimensionless covariance shape is intended for use inside
`vsm`, which supplies the single overall variance scale.

## Usage

``` r
car(x, W, rho = 0.10, fixed = FALSE)
```

## Arguments

- x:

  Factor or design defining q spatial levels.

- W:

  Symmetric nonnegative adjacency/weights matrix with zero diagonal,
  aligned to the levels of x.

- rho:

  Starting proper-CAR dependence parameter.

- fixed:

  Logical indicating whether rho is fixed.

## Details

Let \\D=\mathrm{diag}(W\mathbf{1})\\. The proper CAR precision and
covariance are \$\$Q=D-\rho W,\qquad M=Q^{-1},\qquad K=M/M\_{11}.\$\$
The constructor requires positive row sums, so isolated levels are not
allowed. The admissible open interval for `rho` is obtained from the
eigenvalues of \\D^{-1/2}WD^{-1/2}\\ and is enforced through a
bounded-logit working coordinate. This is a proper, nonsingular CAR
model; an intrinsic singular CAR is not used by the current
precision-based Henderson implementation. The covariance derivative is
analytic.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`sar`](https://covaruber.github.io/sommer/reference/sar.md),
[`maternm`](https://covaruber.github.io/sommer/reference/maternm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(car(location, W), ism(genotype))
} # }
```
