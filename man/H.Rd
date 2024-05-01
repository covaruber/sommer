\name{H}
\alias{H}
\title{Two-way id by features table}
\description{
  \code{H} creates a two way id by features table that can be used as the H argument in the \code{\link{rrc}} function that extracts latent covariates.
}
\usage{
  H(timevar=NULL, idvar=NULL, response=NULL, Gu=NULL)
}

\arguments{
  \item{timevar}{vector of the dataset containing the variable to be used to form columns in the wide table.}
  \item{idvar}{vector of the dataset containing the variable to be used to form rows in the wide table.}
  \item{response}{vector of the dataset containing the response variable to be used to fill the cells of the wide table.}
  \item{Gu}{an optional covariance matrix (\strong{not the inverse}) between levels of the idvar in case a sparse (unbalanced) design between timevar and idvar exist.}
  
}
\details{

This is just an aggregate, reshape and imputation of a long format table to a wide format table.

}
\value{
  \describe{
    
    \item{$H}{ two way table of id by features effects.}
    
  }
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

# data(DT_h2)
# DT <- DT_h2
# DT=DT[with(DT, order(Env)), ]
# H0 <- with(DT, H(Env, Name, y) )
# Z <- with(DT, rrc(Env, H0, 2))


}
\seealso{The function \code{\link{vsc}} to know how to use \code{H} in the \code{\link{rrc}} function.}