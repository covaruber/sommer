\name{DT_polyploid}
\alias{DT_polyploid}
\docType{data}
\title{
Genotypic and Phenotypic data for a potato polyploid population
}
\description{
This dataset contains phenotpic data for 18 traits measured in 187 individuals from a potato diversity panel. In addition contains genotypic data for 221 individuals genotyped with 3522 SNP markers. Please if using this data for your own research make sure you cite Rosyara's (2015) publication (see References).

}
\usage{data("DT_polyploid")}
\format{
  The format is:
 chr "DT_polyploid"
}
\source{
This data was extracted from Rosyara (2016).
}
\references{
If using this data for your own research please cite:

Rosyara Umesh R., Walter S. De Jong, David S. Douches, Jeffrey B. Endelman. Software for genome-wide association studies in autopolyploids and its application to potato. The Plant Genome 2015.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####

data(DT_polyploid)

# ####=========================================####
# ####### convert markers to numeric format
# ####=========================================####
# numo <- atcg1234(data=GT, ploidy=4); 
# numo$M[1:5,1:5]; 
# numo$ref.allele[,1:5]
# 
# ###=========================================####
# ###### plants with both genotypes and phenotypes
# ###=========================================####
# common <- intersect(DT$Name,rownames(numo$M))
# 
# ###=========================================####
# ### get the markers and phenotypes for such inds
# ###=========================================####
# marks <- numo$M[common,]; marks[1:5,1:5]
# DT2 <- DT[match(common,DT$Name),];
# DT2 <- as.data.frame(DT2)
# DT2[1:5,]
# 
# ###=========================================####
# ###### Additive relationship matrix, specify ploidy
# ###=========================================####
# A <- A.mat(marks, ploidy=4)
# D <- D.mat(marks, ploidy=4)
# ###=========================================####
# ### run as mixed model
# ###=========================================####
# ans <- mmer(tuber_shape~1, 
#             random=~vs(Name, Gu=A),
#             data=DT2)
# summary(ans)
# 
# ###=========================================####
# ### run it as GWAS model
# ###=========================================####
# ans2 <- GWAS(tuber_shape~1, 
#              random=~vs(Name,Gu=A),
#              rcov=~units,
#              gTerm = "Name",
#              M=marks, data=DT2)
# summary(ans2)
# plot(ans2$scores[1,])
# plot(ans2$r2m[1,])
}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} }
