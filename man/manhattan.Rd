\name{manhattan}
\alias{manhattan}
\title{Creating a manhattan plot}
\description{
This function was designed to create a manhattan plot using a data frame with columns "Chrom" (Chromosome),  "Position" and "p.val" (significance for the test).  
}
\usage{
manhattan(map, col=NULL, fdr.level=0.05, show.fdr=TRUE, PVCN=NULL, ylim=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{the data frame with 3 columns with names; "Chrom" (Chromosome),  "Position" and "p.val" (significance for the test).}
  \item{col}{colors prefered by the user to be used in the manhattan plot. The default is NULL which will use the red-blue palette.}
  \item{fdr.level}{false discovery rate to be drawn in the plot.}
   \item{show.fdr}{a TRUE/FALSE value indicating if the FDR value should be shown in the manhattan plot or not. By default is TRUE meaning that will be displayed.}
   \item{PVCN}{In case the user wants to provide the name of the column that should be treated as the "p.val" column expected by the program in the 'map' argument.}
    \item{ylim}{the y axis limits for the manhattan plot if the user wants to customize it. By default the plot will reflect the minimum and maximum values found.}
    \item{...}{additional arguments to be passed to the plot function such as pch, cex, etc.}
}


\value{
If all parameters are correctly indicated the program will return:
\describe{

\item{$plot.data}{ a manhattan plot}

}
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}
\examples{
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
set.seed(1234)
pp <- abs(rnorm(500,0,3));pp[23:34] <- abs(rnorm(12,0,20))
geno <- data.frame(Locus=paste("m",1:500, sep="."),Chrom=sort(rep(c(1:5),100)),
                   Position=rep(seq(1,100,1),5),
                   p.val=pp, check.names=FALSE)
geno$Locus <- as.character(geno$Locus)
## look at the data, 5LGs, 100 markers in each
## -log(p.val) value for simulated trait
head(geno)
tail(geno)
manhattan(geno)

}
\seealso{The core functions of the package \code{\link{mmer}} }