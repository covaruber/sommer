\name{H.mat}
\alias{H.mat}

\title{
Combined relationship matrix H
}
\description{
Given a matrix A and a matrix G returns a H matrix with the C++ Armadillo library.
}
\usage{
H.mat(A, G, tau = 1, omega = 1, tolparinv=1e-6)
}

\arguments{
  \item{A}{
Additive relationship matrix based on pedigree.
}
\item{G}{
Additive relationship matrix based on marker data.  
}
\item{tau}{
As described by Martini et al. (2018).
}
\item{omega}{
As described by Martini et al. (2018). 
}
\item{tolparinv}{
Tolerance parameter for matrix inverse used when singularities are encountered in the estimation procedure.
}
}

\details{
See references 
}
\value{

H Matrix with the relationship between the individuals based on pedigree and corrected by molecular information

}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Martini, J. W., Schrauf, M. F., Garcia-Baccino, C. A., Pimentel, E. C., Munilla, S., Rogberg-Munoz, A., ... & Simianer, H. (2018). The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution, 50(1), 16.

}

\examples{
####=========================================####
####random population of 200 lines with 1000 markers
####=========================================####
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}
rownames(M) <- 1:nrow(M)
v <- sample(1:nrow(M),100)
M2 <- M[v,]

A <- A.mat(M) # assume this is a pedigree-based matrix for the sake of example
G <- A.mat(M2)

H <- H.mat(A,G)
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(H[1:15,1:15], col = colfunc(100),Colv = "Rowv")
}
\seealso{The core functions of the package \code{\link{mmer}}}
