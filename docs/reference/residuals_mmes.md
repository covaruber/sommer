# Residuals form a GLMM fitted with mmes

`residuals` method for class `"mmes"`.

## Usage

``` r
# S3 method for class 'mmes'
residuals(object,
                          type=c("response", "deviance", "working"), ...)
```

## Arguments

- object:

  an object of class `"mmes"`

- type:

  For a non-Gaussian PQL fit, return response residuals (default),
  signed deviance residuals, or final IRLS working residuals.

- ...:

  Further arguments to be passed

## Value

For Gaussian fits, residuals of the form e = y - Xb - Zu. For
non-Gaussian PQL fits, residuals on the requested scale.

## Author

Giovanny Covarrubias

## See also

[`residuals`](https://rdrr.io/r/stats/residuals.html),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md)
