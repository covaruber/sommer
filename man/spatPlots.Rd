\name{spatPlots}
\alias{spatPlots}

\title{
Spatial plots
}
\description{
Plot fitted values for all terms in a model to assess spatial fit.
}
\usage{
spatPlots(object, by=NULL, colfunc=NULL,
         row="ROW",range="RANGE", wire=FALSE,
         terms=NULL)

}

\arguments{
  \item{object}{a mixed model fitted using mmer2.}
\item{by}{
A character argument specifying the column name in the data frame to fit a difference variance component for each level of such column. For example, if the model to be fitted should have a variance component for each field where the spatial effect row*range will be estimated, then the column that specifies the different fields should be provided to the `by` argument as a character vector of length one. 
}
\item{colfunc}{
an optional function to creat colors for the spatial plots.
}
\item{row}{
name of the column in the dataset that identifies the spatial position in the x coordinate.
}
\item{range}{
name of the column in the dataset that identifies the spatial position in the y coordinate. 
}
\item{wire}{
a TRUE/FALSE statement indicating if the function should return a wire plot or the default levelplot.
}
\item{terms}{
a character vector specifying the term(s) to be fitted in the spatial plot.
}
}
\details{

The function only takes the asreml model and builds the fitted values for the terms required by building the incidence matrix and getting the BLUPs from the asreml model to build fitted = Zu. Then it makes a levelplot or wireplot to asses the spatial fit.

}
\value{
plots: spatial plots for all fitted values from the mmer2 model terms in the random part.

\describe{
\item{fits}{ A new dataset with the fitted values for all terms in the model.}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient two-dimensional smoothing with P-spline ANOVA mixed models and nested bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Rodriguez-Alvarez, M.X, Boer, M.P., van Eeuwijk, F.A., and Eilers, P.H.C. (2017). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics (to appear). https://doi.org/10.1016/j.spasta.2017.10.003.

}

\examples{
# data("DT_cpdata")
# ### mimic two fields
# aa <- DT; bb <- DT
# aa$FIELD <- "A";bb$FIELD <- "B"
# set.seed(1234)
# aa$Yield <- aa$Yield + rnorm(length(aa$Yield),0,4)
# DT2 <- rbind(aa,bb)
# head(DT2)
# A <- A.mat(GT)
# mix <- mmer(Yield~1,
#             random=~vs(ds(FIELD),id, Gu=A) +
#               vs(ds(FIELD),Rowf) +
#               vs(ds(FIELD),Colf) +
#                 vs(ds(FIELD),spl2D(Row,Col, at=FIELD)),
#             rcov=~vs(ds(FIELD),units),
#             data=DT2)
# summary(mix)$varcomp
# 
# 
# ss <- spatPlots(object=mix, by="FIELD", colfunc=NULL,
#           row="Row",range="Col", wire=FALSE, terms = c("Rowf","Colf","Row"))
}
\seealso{\code{\link{mmer}} -- the core function of the package }
