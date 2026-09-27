# fitted form a LMM fitted with mmes

`fitted` method for class `"mmes"`.

## Usage

``` r
# S3 method for class 'mmes'
fitted(object, type=c("response", "link"), ...)
```

## Arguments

- object:

  an object of class `"mmes"`

- type:

  For a non-Gaussian PQL fit, return fitted means on the response scale
  (default) or linear predictors on the link scale.

- ...:

  Further arguments to be passed to the mmes function

## Value

For Gaussian fits, fitted values of the form y.hat = Xb + Zu. For
non-Gaussian PQL fits, response-scale conditional means by default, or
Xb + Zu with `type="link"`.

## Author

Giovanny Covarrubias

## See also

[`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md)

## Examples

``` r
# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmes(Yield~1,
#               random=~vsm(ism(id),Gu=A)
#                       + Rowf + Colf + spl2Dc(Row,Col),
#               rcov=~units,
#               data=DT)
# 
# ff=fitted(mix1)
# 
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# lattice::wireframe(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted,  
#           aspect=c(61/87,0.4), drape=TRUE,# col.regions = colfunc,
#           light.source=c(10,0,10))
# lattice::levelplot(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted, col.regions = colfunc)
```
