# Back-transform rotated mmes random effects

Back-transforms the selected relationship random effect from the solver
eigenbasis to the original relationship-level basis.

## Usage

``` r
rotate_back_mmes(object)
```

## Arguments

- object:

  An `mmes` model fitted with `vsm(..., rotation=TRUE)`.

## Details

If \\Gu=U\Lambda U'\\ and \\a\\ denotes the solver-basis effect, the
reported effect is \\g=Ua\\. Fitted `mmes` objects already expose this
original-basis effect in `uList`; this function provides the explicit
transformation from the retained engine output.

## Value

A matrix of original-basis random effects, with relationship levels in
rows and covariance-product coordinates in columns.

## See also

[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)
