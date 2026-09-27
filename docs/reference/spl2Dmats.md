# Get Tensor Product Spline Mixed Model Incidence Matrices

`spl2Dmats` gets Tensor-Product P-Spline Mixed Model Incidence Matrices
for use with `sommer` and its main function `mmes`. We thank Sue Welham
for making the TPSbits package available to the community. If you're
using this function for your research please cite her TPSbits package :)
this is mostly a wrapper of her tpsmmb function to enable the use in
sommer.

## Usage

``` r
spl2Dmats(
  x.coord.name,
  y.coord.name,
  data,
  at.name,
  at.levels, 
  nsegments=NULL,
  minbound=NULL,
  maxbound=NULL,
  degree = c(3, 3),
  penaltyord = c(2,2), 
  nestorder = c(1,1),
  method = "Lee"
)
```

## Arguments

- x.coord.name:

  A string. Gives the name of `data` element holding column locations.

- y.coord.name:

  A string. Gives the name of `data` element holding row locations.

- data:

  A dataframe. Holds the dataset to be used for fitting.

- at.name:

  name of a variable defining if the 2D spline matrices should be
  created at different units (e.g., at different environments).

- at.levels:

  a vector of names indicating which levels of the at.name variable
  should be used for fitting the 2D spline function.

- nsegments:

  A list of length 2. Number of segments to split column and row ranges
  into, respectively (= number of internal knots + 1). If only one
  number is specified, that value is used in both dimensions. If not
  specified, (number of unique values - 1) is used in each dimension;
  for a grid layout (equal spacing) this gives a knot at each data
  value.

- minbound:

  A list of length 2. The lower bound to be used for column and row
  dimensions respectively; default calculated as the minimum value for
  each dimension.

- maxbound:

  A list of length 2. The upper bound to be used for column and row
  dimensions respectively; default calculated as the maximum value for
  each dimension.

- degree:

  A list of length 2. The degree of polynomial spline to be used for
  column and row dimensions respectively; default=3.

- penaltyord:

  A list of length 2. The order of differencing for column and row
  dimensions, respectively; default=2.

- nestorder:

  A list of length 2. The order of nesting for column and row
  dimensions, respectively; default=1 (no nesting). A value of 2
  generates a spline with half the number of segments in that dimension,
  etc. The number of segments in each direction must be a multiple of
  the order of nesting.

- method:

  A string. Method for forming the penalty; default=`"Lee"` ie the
  penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The
  alternative method is `"Wood"` ie. the method from Wood et al (2012,
  Stat Comp 23, 341-360). This option is a research tool and requires
  further investigation.

## Value

List of length 7 elements:

1.  `data` = the input data frame augmented with structures required to
    fit tensor product splines in `asreml-R`. This data frame can be
    used to fit the TPS model.

    Added columns:

    - `TP.col`, `TP.row` = column and row coordinates

    - `TP.CxR` = combined index for use with smooth x smooth term

    - `TP.C.n` for n=1:(diff.c) = X parts of column spline for use in
      random model (where diff.c is the order of column differencing)

    - `TP.R.n` for n=1:(diff.r) = X parts of row spline for use in
      random model (where diff.r is the order of row differencing)

    - `TP.CR.n` for n=1:((diff.c\*diff.r)) = interaction between the two
      X parts for use in fixed model. The first variate is a constant
      term which should be omitted from the model when the constant (1)
      is present. If all elements are included in the model then the
      constant term should be omitted, eg.
      `y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...`

    - when `asreml="grp"` or `"sepgrp"`, the spline basis functions are
      also added into the data frame. Column numbers for each term are
      given in the `grp` list structure.

2.  `fR` = Xr1:Zc

3.  `fC` = Xr2:Zc

4.  `fR.C` = Zr:Xc1

5.  `R.fC` = Zr:Xc2

6.  `fR.fC` = Zc:Zr

7.  `all` = Xr1:Zc \| Xr2:Zc \| Zr:Xc1 \| Zr:Xc2 \| Zc:Zr

## Examples

``` r
data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix

M <- spl2Dmats(x.coord.name = "Col", y.coord.name = "Row", data=DT, nseg =c(14,21))
head(M$data)
#>               id Row Col Year      color  Yield FruitAver Firmness Rowf Colf
#> FIELD1.P003 P003   3   1 2014 0.10075269 154.67     41.93  588.917    3    1
#> FIELD1.P004 P004   4   1 2014 0.13891940 186.77     58.79  640.031    4    1
#> FIELD1.P005 P005   5   1 2014 0.08681502  80.21     48.16  671.523    5    1
#> FIELD1.P006 P006   6   1 2014 0.13408561 202.96     48.24  687.172    6    1
#> FIELD1.P007 P007   7   1 2014 0.13519278 174.74     45.83  601.322    7    1
#> FIELD1.P008 P008   8   1 2014 0.17406685 194.16     44.63  656.379    8    1
#>             FIELDINST TP.col TP.row TP.CxR    TP.C.1     TP.C.2    TP.R.1
#> FIELD1.P003    FIELD1      1      3    103 0.2425356 -0.3465516 0.2041241
#> FIELD1.P004    FIELD1      1      4    104 0.2425356 -0.3465516 0.2041241
#> FIELD1.P005    FIELD1      1      5    105 0.2425356 -0.3465516 0.2041241
#> FIELD1.P006    FIELD1      1      6    106 0.2425356 -0.3465516 0.2041241
#> FIELD1.P007    FIELD1      1      7    107 0.2425356 -0.3465516 0.2041241
#> FIELD1.P008    FIELD1      1      8    108 0.2425356 -0.3465516 0.2041241
#>                    TP.R.2    TP.CR.1       TP.CR.2     TP.CR.3       TP.CR.4
#> FIELD1.P003 -2.064187e-01 0.04950738 -5.006390e-02 -0.07073956  7.153475e-02
#> FIELD1.P004 -1.548141e-01 0.04950738 -3.754792e-02 -0.07073956  5.365106e-02
#> FIELD1.P005 -1.032094e-01 0.04950738 -2.503195e-02 -0.07073956  3.576738e-02
#> FIELD1.P006 -5.160468e-02 0.04950738 -1.251597e-02 -0.07073956  1.788369e-02
#> FIELD1.P007 -1.743562e-13 0.04950738 -4.228759e-14 -0.07073956  6.042343e-14
#> FIELD1.P008  5.160468e-02 0.04950738  1.251597e-02 -0.07073956 -1.788369e-02
# m1g <- mmes(Yield~1+TP.CR.2+TP.CR.3+TP.CR.4,
#             random=~Rowf+Colf+vsm(ism(M$fC))+vsm(ism(M$fR))+
#               vsm(ism(M$fC.R))+vsm(ism(M$C.fR))+vsm(ism(M$fC.fR))+
#               vsm(ism(id),Gu=A),
#             data=M$data, tolpar = 1e-6,
#             iters=30)
# 
# summary(m1g)$varcomp
```
