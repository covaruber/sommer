\name{corImputation}
\alias{corImputation}
\title{Imputing a matrix using correlations}
\description{
  corImputation imputes missing data based on the correlation that exists between row levels.
}
\usage{
  corImputation(wide, Gu=NULL, nearest=10, roundR=FALSE)
}

\arguments{
  \item{wide}{numeric matrix with individuals in rows and time variable in columns (e.g., environments, genetic markers, etc.).}
  \item{Gu}{optional correlation matrix between the individuals or row levels. If NULL it will be computed as the correlation of t(wide).}
  \item{nearest}{integer value describing how many nearest neighbours (the ones showing the highest correlation) should be used to average and return the imputed value.}
  \item{roundR}{a TRUE/FALSE statement describing if the average result should be rounded or not. This may be specifically useful for categorical data in the form of numbers (e.g., -1,0,1).}
}
\value{
  \describe{
    
    \item{$res}{ a list with the imputed matrix and the original matrix.}
    
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

####################################
### imputing genotype data example
####################################
# data(DT_cpdata)
# X <- GT_cpdata
# # add missing data
# v <- sample(1:length(X), 500)
# Xna <- X
# Xna[v]<- NA
# ## impute (can take some time)
# Y <- corImputation(wide=Xna, Gu=NULL, nearest=20, roundR=TRUE) 
# cm <- table(Y$imputed[v],X[v])
# ## calculate accuracy
# sum(diag(cm))/length(v)
####################################
### imputing phenotypic data example
####################################
# data(DT_h2)
# X <- reshape(DT_h2[,c("Name","Env","y")], direction = "wide", idvar = "Name",
#                 timevar = "Env", v.names = "y", sep= "_")
# rownames(X) <- X$Name
# X <- as.matrix(X[,-1])
# head(X)
# # add missing data
# v <- sample(1:length(X), 50)
# Xna <- X
# Xna[v]<- NA
# ## impute
# Y <- corImputation(wide=Xna, Gu=NULL, nearest=20, roundR=TRUE)
# plot(y=Y$imputed[v],x=X[v], xlab="true",ylab="predicted")
# cor(Y$imputed[v],X[v], use = "complete.obs")


}

