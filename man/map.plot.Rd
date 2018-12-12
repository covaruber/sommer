\name{map.plot}
\alias{map.plot}
\title{Creating a genetic map plot}
\description{
This function was designed to create a genetic map plot using a data frame indicating the Linkage Group (LG), Position and marker names (Locus). 
}
\usage{
map.plot(data, trait = NULL, trait.scale = "same", 
        col.chr = NULL, col.trait = NULL, type = "hist", cex = 0.4,
        lwd = 1, cex.axis = 0.4, cex.trait=0.8, jump = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the data frame with 3 columns with names; Locus, LG and Position}
  \item{trait}{if something wants to be plotted next the linkage groups the user must indicate the name of the column containing the values to be ploted, i.e. p-values, LOD scores, X2 segregation distortion values, etc.}
  \item{trait.scale}{is trait is not NULL, this is a character value indicating if the y axis limits for the trait plotted next to the chromosomes should be the same or different for each linkage group. The default value is "same", which means that the same y axis limit is conserved across linkage groups. For giving an individual y axis limit for each linkage group write "diff".}
  \item{col.chr}{a vector with color names for the chromosomes. If NULL they will be drawn in gray-black scale.}
  \item{col.trait}{a vector with color names for the dots, lines or histogram for the trait plotted next to the LG's}
  \item{type}{a character value indicating if the trait should be plotted as scatterplot `dot`, histogram `hist`, line `line` next to the chromosomes.}
  \item{cex}{the cex value determining the size of the cM position labels in the LGs}
  \item{lwd}{the width of the lines in the plot}
  \item{cex.axis}{the cex value for sizing the labels of LGs and traits plotted (top labels)}
  \item{cex.trait}{the cex value for sizing the dots or lines of the trait plotted}
  \item{jump}{a scalar value indicating how often should be drawn a number next to the LG indicating the position. The default is 5 which means every 5 cM a label will be drawn, i.e. 0,5,10,15,... cM.}
}

\value{
If all parameters are correctly indicated the program will return:
\describe{

\item{$plot.data}{ a plot with the LGs and the information used to create a plot}

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
geno <- data.frame(Locus=paste("m",1:500, sep="."),LG=sort(rep(c(1:5),100)),
                   Position=rep(seq(1,100,1),5),
                   X2=rnorm(500,10,4), check.names=FALSE)
geno$Locus <- as.character(geno$Locus)
## look at the data, 5LGs, 100 markers in each
## X2 value for segregation distortion simulated
head(geno)
tail(geno)
table(geno$LG) # 5 LGs, 100 marks
map.plot(geno, trait="X2", type="line")
map.plot(geno, trait="X2", type="hist")
map.plot(geno, trait="X2", type="dot")

}
\seealso{The core functions of the package \code{\link{mmer}}}