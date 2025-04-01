\name{LD.decay}
\alias{LD.decay}
\title{Calculation of linkage disequilibrium decay}
\description{
This function calculates the LD decay based on a marker matrix and a map with distances between markers in cM or base pairs.  
}
\usage{
LD.decay(markers,map,silent=FALSE,unlinked=FALSE,gamma=0.95)
}

\arguments{
  \item{markers}{a numeric matrix of markers (columns) by individuals (rows) in -1, 0, 1 format.}
  \item{map}{a data frame with 3 columns "Locus" (name of markers), "LG" (linkage group or chromosome), and "Position" (in cM or base pairs). }
  \item{silent}{a TRUE/FALSE value statement indicating if the program should or should not display the progress bar. silent=TRUE means that will not be displayed.}
  \item{unlinked}{a TRUE/FALSE value statement indicating if the program should or should not calculate the alpha(see next argument) percentile of interchromosomal LD.}
  \item{gamma}{a percentile value for LD breakage to be used in the calculation of interchromosomal LD extent.}

}

\value{
\describe{
\item{$resp}{ a list with 3 elements; "by.LG", "all.LG", "LDM". The first element (by.LG) is a list with as many elements as chromosomes where each contains a matrix with 3 columns, the distance, the r2 value, and the p-value associated to the chi-square test for disequilibrium. The second element (all.LG) has a big matrix with distance, r2 values and p-values, for each point from all chromosomes in a single data.frame. The third element (LDM) is the matrix of linkage disequilibrium between pairs of markers.

If unlinked is selected the program should return the gamma percentile interchromosomal LD (r2) for each chromosome and average.
}
}

}
\references{

Laido, Giovanni, et al. Linkage disequilibrium and genome-wide association mapping in tetraploid wheat (Triticum turgidum L.). PloS one 9.4 (2014): e95211.

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using 
#### command + shift + C |OR| control + shift + C
####=========================================####
data(DT_cpdata)
#### get the marker matrix
CPgeno <- GT_cpdata; CPgeno[1:5,1:5]
#### get the map
mapCP <- MP_cpdata; head(mapCP)
names(mapCP) <- c("Locus","Position","LG")
#### with example purposes we only do 3 chromosomes
mapCP <- mapCP[which(mapCP$LG <= 3),]
#### run the function
# res <- LD.decay(CPgeno, mapCP)
# names(res)
#### subset only markers with significant LD
# res$all.LG <- res$all.LG[which(res$all.LG$p < .001),]
#### plot the LD decay
# with(res$all.LG, plot(r2~d,col=transp("cadetblue"),
#                     xlim=c(0,55), ylim=c(0,1), 
#                     pch=20,cex=0.5,yaxt="n",
#                     xaxt="n",ylab=expression(r^2),
#                     xlab="Distance in cM")
#                     )
# axis(1, at=seq(0,55,5), labels=seq(0,55,5))
# axis(2,at=seq(0,1,.1), labels=seq(0,1,.1), las=1)

#### if you want to add the loess regression lines
#### this could take a long time!!!
# mod <- loess(r2 ~ d, data=res$all.LG)
# par(new=T)
# lilo <- predict(mod,data.frame(d=1:55))
# plot(lilo, bty="n", xaxt="n", yaxt="n", col="green", 
#      type="l", ylim=c(0,1),ylab="",xlab="",lwd=2)
# res3 <- LD.decay(markers=CPgeno, map=mapCP, 
#                 unlinked = TRUE,gamma = .95)
# abline(h=res3$all.LG, col="red")
}
