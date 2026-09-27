# Predict factor-analytic scores from a fitted mmes model

Predicts latent factor scores for each level of the main random-effect
term crossed with a
[`fam()`](https://covaruber.github.io/sommer/reference/fam.md) or
[`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
covariance-shaping factor, combining the fitted loadings (see
[`loadings_mmes`](https://covaruber.github.io/sommer/reference/loadings_mmes.md)),
the fitted covariance matrix, and the BLUPs already stored in
`object$uList`.

## Usage

``` r
scores_mmes(object, term = NULL, method = c("regression", "bartlett"), 
            varianceScale=TRUE, rotation=TRUE)
```

## Arguments

- object:

  a fitted model of class `"mmes"`.

- term:

  character name of the random term built with a single
  [`fam()`](https://covaruber.github.io/sommer/reference/fam.md) or
  [`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
  covariance-shaping factor. If `NULL`, the unique such term in the
  model is used automatically.

- method:

  `"regression"` (default) uses the full fitted covariance \\\Sigma\\
  (Thomson's method); `"bartlett"` uses only the specific (residual)
  variances \\\Psi\\ (Bartlett's classic unbiased estimator).

- varianceScale:

  a logical argument to indicate if loadings should be returned in
  variance scale (multiplied by sqrt(sigma2)).

- rotation:

  a logical value to indicate if loadings should be rotated by its
  singular vectors.

## Details

Let \\U\\ be the levels-by-levels BLUP matrix in `object$uList[[term]]`
and \\\Lambda\\ the loadings from
[`loadings_mmes`](https://covaruber.github.io/sommer/reference/loadings_mmes.md).
The regression (Thomson) scores are \$\$F = U\Sigma^{-1}\Lambda,\$\$ and
the Bartlett scores are \$\$F = U\Psi^{-1}\Lambda(\Lambda^{\mathsf
T}\Psi^{-1}\Lambda)^{-1}.\$\$ Because \\U\\ already contains shrunken
(BLUP) effects, the resulting scores are similarly regularized.

## Value

A matrix with one row per level of the main random-effect term and one
column per latent factor.

## See also

[`loadings_mmes`](https://covaruber.github.io/sommer/reference/loadings_mmes.md),
[`fam`](https://covaruber.github.io/sommer/reference/fam.md),
[`rrcm`](https://covaruber.github.io/sommer/reference/rrcm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
mix <- mmes(BLUEs ~ trial,
            random = ~ vsm(fam(trial, 2), ism(genotype)),
            rcov = ~ units, data = dt)
scores_mmes(mix)
} # }
```
