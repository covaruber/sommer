# Get Tensor Product Spline Mixed Model Incidence Matrices

`tpsmmbwrapper` is a wrapper of tpsmmb function from the TPSbits package
to avoid version dependencies but if you're using this function for your
research please cite the TPSbits package. This function is internally
used by the `spl2Dmatrices` function to get Tensor-Product P-Spline
Mixed Model Bits (design matrices) for use with `sommer`.

## Usage

``` r
tpsmmbwrapper(
  columncoordinates,
  rowcoordinates,
  data,
  nsegments=NULL,
  minbound=NULL,
  maxbound=NULL,
  degree = c(3, 3),
  penaltyord = c(2, 2),
  nestorder = c(1, 1),
  asreml = "mbf",
  eigenvalues = "include",
  method = "Lee",
  stub = NULL
)
```

## Arguments

- columncoordinates:

  A string. Gives the name of `data` element holding column locations.

- rowcoordinates:

  A string. Gives the name of `data` element holding row locations.

- data:

  A dataframe. Holds the dataset to be used for fitting.

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

- asreml:

  A string. Indicates the types of structures to be generated for use in
  asreml models; default `"mbf"`. The appropriate eigenvalue scaling is
  included within the Z matrices unless setting `scaling="none"` is
  used, and then the scaling factors are supplied separately in the
  returned object.

  - `asreml="mbf"` indicates the function should put the spline design
    matrices into structures for use with `"mbf"`;

  - `asreml="grp"` indicates the function should add the composite
    spline design matrices (eg. for second-order differencing, matrices
    Xr1:Zc, Xr2:Zc, Zr:Xc1, Zr:Xc2 and Zc:Zr) into the data frame and
    provide a group list structure for each term;

  - `asreml="sepgrp"` indicates the function should generate the
    individual X and Z spline design matrices separately (ie. Xc, Xr, Zc
    and Zr), plus the smooth x smooth interaction term as a whole (ie.
    Zc:Zr), and provide a group list structure for each term.

  - `asreml="own"` indicates the function should generate the composite
    matrix ( Xr:Zc \| Zr:Xc \| Zc:Zr ) as a single set of columns.

- eigenvalues:

  A string. Indicates whether eigenvalues should be included within the
  Z design matrices `eigenvalues="include"`, or whether this scaling
  should be omitted (`eigenvalues="omit"`); default
  `eigenvalues="include"`. If the eigenvalue scaling is omitted from the
  Z design matrices, then it should instead be included in the model as
  a variance structure to obtain the correct TPspline model.

- method:

  A string. Method for forming the penalty; default=`"Lee"` ie the
  penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The
  alternative method is `"Wood"` ie. the method from Wood et al (2012,
  Stat Comp 23, 341-360). This option is a research tool and requires
  further investigation.

- stub:

  A string. Stub to be attached to names in the `mbf` list to avoid
  over-writing structures and general confusion.

## Value

List of length 7, 8 or 9 (according to the `asreml` and `eigenvalues`
parameter settings).

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

2.  `mbflist` = list that can be used in call to asreml (so long as Z
    matrix data frames extracted with right names, eg BcZ\<stub\>.df)

3.  `BcZ.df` = mbf data frame mapping onto smooth part of column spline,
    last column (labelled `TP.col`) gives column index

4.  `BrZ.df` = mbf data frame mapping onto smooth part of row spline,
    last column (labelled `TP.row`) gives row index

5.  `BcrZ.df` = mbf data frame mapping onto smooth x smooth term, last
    column (labelled `TP.CxR`) maps onto col x row combined index

6.  `dim` = list structure, holding dimension values relating to the
    model:

    1.  `"diff.c"` = order of differencing used in column dimension

    2.  `"nbc"` = number of random basis functions in column dimension

    3.  `"nbcn"` = number of nested random basis functions in column
        dimension used in smooth x smooth term

    4.  `"diff.r"` = order of differencing used in column dimension

    5.  `"nbr"` = number of random basis functions in column dimension

    6.  `"nbrn"` = number of nested random basis functions in column
        dimension used in smooth x smooth term

7.  `trace` = list of trace values for ZGZ' for the random TPspline
    terms, where Z is the design matrix and G is the known diagonal
    variance matrix derived from eigenvalues. This can be used to
    rescale the spline design matrix (or equivalently variance
    components).

8.  `grp` = list structure, only added for settings `asreml="grp"`,
    `asreml="sepgrp"` or `asreml="own"`. For `asreml="grp"`, provides
    column indexes for each of the 5 random components of the 2D
    splines. For `asreml="sepgrp"`, provides column indexes for each of
    the X and Z component matrices for the 1D splines, plus the
    composite smooth x smooth interaction term. For `asreml="own"`,
    provides column indexes for the composite random model. Dimensions
    of the components can be derived from the values in the `dim` item.
    The Z terms are scaled by the associated eigenvalues when
    `eigenvalues="include"`, but not when `eigenvalues="omit"`.

9.  `eigen` = list structure, only added for option setting
    `eigenvalues="omit"`. Holds the diagonal elements of the inverse
    variance matrix for the terms Xc:Zr (called `diagr`), Zc:Xr (called
    `diagc`) and Zc:Zr (called `diagcr`).
