\name{redmm}
\alias{redmm}
\title{Reduced Model Matrix}
\description{
  `redmm` reduces a model matrix by performing a singular value decomposition or Cholesky on an incidence matrix.
}
\usage{
  redmm(x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE)
}

\arguments{
  \item{x}{as vector with values to form a model matrix or the complete incidence matrix itself for an effect of interest.}
  \item{M}{an optional matrix of features explaining the levels of x. If not provided is assumed that the entire incidence matrix has been provided in x. But if provided, the decomposition occurs in the matrix M. }
  \item{Lam}{a matrix of loadings in case is already available to avoid recomputing it.}
  \item{nPC}{number of principal components to keep from the matrix of loadings to form the model matrix.}
  \item{cholD}{should a Cholesky or a Singular value decomposition should be used. The default is the SVD.}
  \item{returnLam}{should the function return the loading matrix in addition to the incidence matrix. Default is FALSE.}
}
\value{
  \describe{
    
    \item{$S3}{ A list with 3 elements:
    
    1) The model matrix to be used in the mixed modeling.
    
    2) The reduced matrix of loadings (nPC columns).
    
    3) The full matrix of loadings.
    
    }
  }
}
\references{
  
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

data(DT_technow)
DT <- DT_technow
Md <- Md_technow

M <- tcrossprod(Md)
xx = with(DT, redmm(x=dent, M=M, nPC=10))

# ans <- mmec(GY~1,
#             # new model matrix instead of dent
#             random=~vsc(isc(xx$Z)),
#             rcov=~units,
#             data=DT)
# summary(ans)$varcomp
# u = xx$Lam * ans$uList[[1]] # change * for matrix product

}
\seealso{The core functions of the package \code{\link{mmec}} }