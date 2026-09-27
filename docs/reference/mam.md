# Moving-average covariance structures of order one or two

Moving-average covariance structures of order one or two. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
mam(x, order = 1L, theta = NULL, fixed = NULL)
ma1m(x, theta = 0.15, fixed = FALSE)
ma2m(x, theta = c(0.15, 0.05), fixed = NULL)
```

## Arguments

- x:

  Ordered factor defining covariance levels.

- order:

  Moving-average order; currently 1 or 2.

- theta:

  Starting MA polynomial coefficients. Defaults to 0.15 for each order.

- fixed:

  Logical vector with one value per MA coefficient.

## Details

Using the conventional polynomial
\\e_t+\theta_1e\_{t-1}+\theta_2e\_{t-2}\\, the autocovariance at lag
\\h\\ is proportional to \\\sum\_{j=0}^{p-h}\theta_j\theta\_{j+h}\\ with
\\\theta_0=1\\. The covariance is normalized to correlation scale and is
zero beyond lag \\p\\. Invertibility of the MA polynomial is not
required merely to define this covariance matrix. `ma1m` and `ma2m` are
convenience wrappers.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(ma1m(time), ism(id))
vsm(ma2m(time), ism(id))
} # }
```
