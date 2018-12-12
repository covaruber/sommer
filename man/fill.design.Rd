\name{fill.design}
\alias{fill.design}
\title{Filling the design of an experiment}
\description{
  The fill.design function allows the user to fill the missing rows and ranges from an experiment to run specific designs or apply a post blocking. The data frame requires the presence of 2 numeric columns indicating the x and y coordinates, here denominated rows and ranges. This can be effectively used for multi environment trials by using the `by` argument where you specify the column that indicates the environments. 
}
\usage{
  fill.design(x,rows="ROW",ranges="RANGE",by, extra)
}

\arguments{
  \item{x}{a dataframe with 2 obligatory columns; rows and ranges which can have different names and can be matched with the next arguments.}
  \item{rows}{the name of the numeric column that indicates one direction in the field.}
  \item{ranges}{the name of the numeric column that indicates the other direction in the field.}
  \item{by}{optional argument to indicate the name of the column of the dataframe x that indicates the environments so the field is filled by environment.}
  \item{extra}{name of the extra columns to be filled in the dataset based. This are filled based on the rows, ranges information.}
  
}
\value{
  \describe{
    
    \item{$xnew}{ a new dataframe identical to the one provided but with missing rows and ranges filled in.}
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
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
# data(DT_cpdata)
# #### look at the data
# head(DT)
# #### fill the design
# gg <- fill.design(x=DT, rows="Row",ranges="Col")
# head(gg)

}
\seealso{The core functions of the package \code{\link{mmer}}}