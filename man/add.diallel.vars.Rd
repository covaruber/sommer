\name{add.diallel.vars}
\alias{add.diallel.vars}
\title{add.diallel.vars}
\description{
  `add.diallel.vars` adds 4 columns to the provided diallel dataset. Specifically, the user provides a dataset with indicator variables for who is the male and female parent and the function returns the same dataset with 4 new dummy variables to allow the model fit of diallel models.
}
\usage{
  add.diallel.vars(df, par1="Par1", par2="Par2",sep.cross="-")
}

\arguments{
  \item{df}{a dataset with the two indicator variables for who is the male and female parent.}
  \item{par1}{the name of the column indicating who is the first parent (e.g. male).}
  \item{par2}{the name of the column indicating who is the second parent (e.g. female).}
  \item{sep.cross}{the character that should be used when creating the column for cross.id. A simple paste of the columns par1 and par2.}
}
\value{
  \describe{
    
    A new data set with the following 4 new dummy variables to allow the fit of complex diallel models:
    
    \item{is.cross}{returns a 0 if is a self and a 1 for a cross.}
    \item{is.self}{returns a 0 if is a cross and a 1 is is a self.}
    \item{cross.type}{returns a -1 for a direct cross, a 0 for a self and a 1 for a reciprocal cross.}
    \item{cross.id}{returns a column psting the par1 and par2 columns.}
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
data(DT_mohring)
DT <- DT_mohring
head(DT)
DT2 <- add.diallel.vars(DT,par1="Par1", par2="Par2")
head(DT2)
## see ?DT_mohring for an example on how to use the data to fit diallel models.

}
\seealso{The \code{\link{overlay}} function and the \code{\link{DT_mohring}} example. }


