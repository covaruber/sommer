\name{spl2Dmats}
\alias{spl2Dmats}
\title{Get Tensor Product Spline Mixed Model Incidence Matrices}
\usage{
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
}
\arguments{
\item{x.coord.name}{A string. Gives the name of \code{data} element
holding column locations.}

\item{y.coord.name}{A string. Gives the name of \code{data} element
holding row locations.}

\item{data}{A dataframe. Holds the dataset to be used for fitting.}

\item{at.name}{name of a variable defining if the 2D spline matrices should be created at different units (e.g., at different environments).}

\item{at.levels}{a vector of names indicating which levels of the at.name variable should be used for fitting the 2D spline function.}

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

\item{method}{A string. Method for forming the  penalty; default=\code{"Lee"}
ie the penalty from Lee, Durban & Eilers (2013, CSDA 61, 22-37). The
alternative method is \code{"Wood"} ie. the method from Wood et al (2012,
Stat Comp 23, 341-360).
This option is a research tool and requires further investigation.}

}
\value{
List of length 7 elements:
\enumerate{
\item \code{data} = the input data frame augmented with structures required
    to fit tensor product splines in \code{asreml-R}. This data frame can be used
    to fit the TPS model.

    Added columns:
    \itemize{
    \item \code{TP.col}, \code{TP.row} = column and row coordinates
    \item \code{TP.CxR} = combined index for use with smooth x smooth term
    \item \code{TP.C.n} for n=1:{diff.c} = X parts of column spline for use
    in random model (where diff.c is the order of column differencing)
    \item \code{TP.R.n} for n=1:{diff.r} = X parts of row spline for use in
    random model (where diff.r is the order of row differencing)
    \item \code{TP.CR.n} for n=1:{(diff.c*diff.r)} = interaction between the
    two X parts for use in fixed model. The first variate is
    a constant term which should be omitted from the model when the constant
    (1) is present. If all elements are
    included in the model then the constant term should be omitted,
    eg. \code{y ~ -1 + TP.CR.1 + TP.CR.2 + TP.CR.3 + TP.CR.4 + other terms...}
    \item when \code{asreml="grp"} or \code{"sepgrp"}, the spline basis
    functions are also added into the data frame. Column numbers for each
    term are given in the \code{grp} list structure.
    }
  
     
\item \code{fR} = Xr1:Zc
\item \code{fC} = Xr2:Zc
\item \code{fR.C} = Zr:Xc1
\item \code{R.fC} = Zr:Xc2 
\item \code{fR.fC} = Zc:Zr
\item \code{all} = Xr1:Zc | Xr2:Zc | Zr:Xc1 | Zr:Xc2 | Zc:Zr

}

}
\description{
\code{spl2Dmats} gets Tensor-Product P-Spline Mixed Model Incidence Matrices
    for use with \code{sommer} and its main function \code{mmer}. We thank Sue Welham for making the TPSbits package available to the community. If you're using this function for your research please cite her TPSbits package :) this is mostly a wrapper of her tpsmmb function to enable the use in sommer.
}
\examples{

data("DT_cpdata")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix

M <- spl2Dmats(x.coord.name = "Col", y.coord.name = "Row", data=DT, nseg =c(14,21))
head(M$data)
# m1g <- mmer(Yield~1+TP.CR.2+TP.CR.3+TP.CR.4,
#             random=~Rowf+Colf+vsr(M$fC)+vsr(M$fR)+
#               vsr(M$fC.R)+vsr(M$C.fR)+vsr(M$fC.fR)+
#               vsr(id,Gu=A),
#             data=M$data, tolpar = 1e-6,
#             iters=30)
# 
# summary(m1g)$varcomp

}
