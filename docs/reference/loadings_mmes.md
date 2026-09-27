# Extract factor-analytic loadings from a fitted mmes model

Reconstructs the normalized loadings matrix and specific variances of a
[`fam()`](https://covaruber.github.io/sommer/reference/fam.md) or
[`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
covariance-shaping factor fitted inside
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md), on the
same scale as the fitted covariance matrix returned in `object$theta`.

## Usage

``` r
loadings_mmes(object, term = NULL, varianceScale=TRUE, rotation=TRUE)
```

## Arguments

- object:

  a fitted model of class `"mmes"`.

- term:

  character name of the random term (matching a name in
  `object$covStruct`) built with a single
  [`fam()`](https://covaruber.github.io/sommer/reference/fam.md) or
  [`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
  covariance-shaping factor. If `NULL`, the unique such term in the
  model is used automatically; an error is raised if none or more than
  one exist.

- varianceScale:

  a logical argument to indicate if loadings should be returned in
  variance scale (multiplied by sqrt(sigma2)).

- rotation:

  a logical value to indicate if loadings should be rotated by its
  singular vectors.

## Details

[`fam()`](https://covaruber.github.io/sommer/reference/fam.md)/[`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
remove their otherwise redundant internal scale so that
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) owns the
single overall variance; the loadings and specific variances reported by
`loadings_mmes` are re-normalized so that \$\$\Sigma =
\sigma^2(\Lambda\Lambda^{\mathsf T}+\Psi)\$\$ reproduces
`object$theta[[term]]` exactly, where \\\sigma^2\\ is
`object$covPar[[term]][1]`.

## Value

A list containing:

- loadings:

  a levels-by-factor matrix \\\Lambda\\.

- specific:

  a named vector of specific variances \\\Psi\\ (identically 1 for every
  level in a
  [`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
  term).

- sigma2:

  the fitted overall variance scale.

- model:

  `"fa"` or `"rr"`.

- term:

  the resolved term name.

## See also

[`scores_mmes`](https://covaruber.github.io/sommer/reference/scores_mmes.md),
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
loadings_mmes(mix)
} # }
```
