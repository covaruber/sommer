\name{tps}
\alias{tps}
\title{Get Tensor Product Spline Mixed Model Incidence Matrices}
\usage{
tps(
  columncoordinates,
  rowcoordinates,
  nsegments=NULL,
  minbound=NULL,
  maxbound=NULL,
  degree = c(3, 3),
  penaltyord = c(2, 2),
  nestorder = c(1, 1),
  asreml = "grp",
  eigenvalues = "include",
  method = "Lee",
  stub = NULL
)
}
\arguments{
\item{columncoordinates}{A string. Gives the name of \code{data} element
holding column locations.}

\item{rowcoordinates}{A string. Gives the name of \code{data} element
holding row locations.}

\item{nsegments}{A list of length 2. Number of segments to split column and
row ranges into, respectively (= number of internal knots + 1). If only
one number is specified, that value is used in both dimensions. If not
specified, (number of unique values - 1) is used in each dimension;
for a grid layout (equal spacing) this gives a knot at each data value.}

\item{minbound}{A list of length 2. The lower bound to be used for column
and row dimensions respectively; default calculated as the minimum value
for each dimension.}

\item{maxbound}{A list of length 2. The upper bound to be used for column
and row dimensions respectively; default calculated as the maximum value
for each dimension.}

\item{degree}{A list of length 2. The degree of polynomial spline to be used
for column and row dimensions respectively; default=3.}

\item{penaltyord}{A list of length 2. The order of differencing for column
and row dimensions, respectively; default=2.}

\item{nestorder}{A list of length 2. The order of nesting for column and row
dimensions, respectively; default=1 (no nesting). A value of 2 generates
a spline with half the number of segments in that dimension, etc. The
number of segments in each direction must be a multiple of the order
of nesting.}

\item{asreml}{A string. Indicates the types of structures to be generated
for use in asreml models; default \code{"mbf"}. The
appropriate eigenvalue scaling is included within the Z matrices unless
setting \code{scaling="none"} is used, and then the scaling factors are
supplied separately in the returned object.
\itemize{
\item {\code{asreml="mbf"}} indicates the function should put the
spline design matrices into structures for use with \code{"mbf"};
\item {\code{asreml="grp"}} indicates the function should add the
composite spline design matrices (eg. for second-order differencing,
matrices Xr1:Zc, Xr2:Zc, Zr:Xc1, Zr:Xc2 and Zc:Zr) into the data frame
and provide a group list structure for each term;
\item {\code{asreml="sepgrp"}} indicates the function should generate the
individual X and Z spline design matrices separately (ie. Xc, Xr, Zc and
Zr), plus the smooth x smooth interaction term as a whole (ie. Zc:Zr),
and provide a group list structure for each term.
\item {\code{asreml="own"}} indicates the function should generate the
composite matrix ( Xr:Zc | Zr:Xc | Zc:Zr ) as a single set of columns.
}}

\item{eigenvalues}{A string. Indicates whether eigenvalues should be
included within the Z design matrices \code{eigenvalues="include"}, or
whether this scaling should be omitted (\code{eigenvalues="omit"});
default \code{eigenvalues="include"}. If the eigenvalue scaling is
omitted from the Z design matrices, then it should instead be included in
the model as a variance structure to obtain the correct TPspline model.}

\item{method}{A string. Method for forming the  penalty; default=\code{"Lee"}
ie the penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The
alternative method is \code{"Wood"} ie. the method from Wood et al (2012,
Stat Comp 23, 341-360).
This option is a research tool and requires further investigation.}

\item{stub}{A string. Stub to be attached to names in the \code{mbf} list to
avoid over-writing structures and general confusion.}
}
\value{
List of length 7, 8 or 9 (according to the \code{asreml} and
    \code{eigenvalues} parameter settings).
    \enumerate{
     \item \code{data} = the input data frame augmented with structures required
    to fit tensor product splines in \code{asreml-R}. This data frame can be used
    to fit the TPS model.

    Added columns:
    \itemize{
    \item \code{TP.col}, \code{TP.row} = column and row coordinates
    \item \code{TP.CxR} = combined index for use with smooth x smooth term
    \item \code{TP.C.n} for n=1:(diff.c) = X parts of column spline for use
    in random model (where diff.c is the order of column differencing)
    \item \code{TP.R.n} for n=1:(diff.r) = X parts of row spline for use in
    random model (where diff.r is the order of row differencing)
    \item \code{TP.CR.n} for n=1:((diff.c*diff.r)) = interaction between the
    two X parts for use in fixed model. The first variate is
    a constant term which should be omitted from the model when the constant
    (1) is present. If all elements are
    included in the model then the constant term should be omitted,
    eg. \code{y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...}
    \item when \code{asreml="grp"} or \code{"sepgrp"}, the spline basis
    functions are also added into the data frame. Column numbers for each
    term are given in the \code{grp} list structure.
    }
\item \code{mbflist} = list that can be used in call to asreml (so long as Z
    matrix data frames extracted with right names, eg BcZ<stub>.df)
\item \code{BcZ.df} = mbf data frame mapping onto smooth part of column
    spline, last column (labelled \code{TP.col}) gives column index
\item \code{BrZ.df} = mbf data frame mapping onto smooth part of row spline,
    last column (labelled \code{TP.row}) gives row index
\item \code{BcrZ.df} = mbf data frame mapping onto smooth x smooth term, last
    column (labelled \code{TP.CxR}) maps onto col x row combined index
 \item \code{dim} = list structure, holding dimension values relating to the
    model:
    \enumerate{
    \item \code{"diff.c"} = order of differencing used in column dimension
    \item \code{"nbc"} = number of random basis functions in column dimension
    \item \code{"nbcn"} = number of nested random basis functions in column dimension
           used in smooth x smooth term
    \item \code{"diff.r"} = order of differencing used in column dimension
    \item \code{"nbr"} = number of random basis functions in column dimension
    \item \code{"nbrn"} = number of nested random basis functions in column dimension
           used in smooth x smooth term
    }
\item \code{trace} = list of trace values for ZGZ' for the random TPspline
   terms, where Z is the design matrix and G
   is the known diagonal variance matrix derived from eigenvalues. This can
   be used to rescale the spline design matrix (or equivalently variance
   components).
\item \code{grp} = list structure, only added for settings
   \code{asreml="grp"},  \code{asreml="sepgrp"} or \code{asreml="own"}.
   For \code{asreml="grp"}, provides column indexes for each of the 5
   random components of the 2D splines.
   For \code{asreml="sepgrp"}, provides column indexes for each of the X and
   Z component matrices for the 1D splines, plus the composite smooth x
   smooth interaction term. For \code{asreml="own"}, provides column indexes
   for the composite random model.
   Dimensions of the components can be derived from the values in the
   \code{dim} item.  The Z terms are scaled by the associated
   eigenvalues when \code{eigenvalues="include"}, but not when
   \code{eigenvalues="omit"}.
\item \code{eigen} = list structure, only added for option setting
   \code{eigenvalues="omit"}. Holds the diagonal elements of the inverse
   variance matrix for the terms Xc:Zr (called \code{diagr}), Zc:Xr
   (called \code{diagc}) and Zc:Zr (called \code{diagcr}).
}
}
\description{
\code{tps} is a wrapper of tpsmmb function from the TPSbits package to avoid version dependencies but if you're using this function for your research please cite the TPSbits package. This function is internally used by the \code{spl2Dmatrices} function to get Tensor-Product P-Spline Mixed Model Bits
    (design matrices) for use with \code{sommer} and its main function \code{mmer}.
}

