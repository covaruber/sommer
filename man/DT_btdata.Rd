\name{DT_btdata}
\alias{DT_btdata}
\docType{data}
\title{
Blue Tit Data for a Quantitative Genetic Experiment
}
\description{

a data frame with 828 rows and 7 columns, with variables tarsus length (tarsus) and colour (back) measured on 828 individuals (animal). The mother of each is also recorded (dam) together with the foster nest (fosternest) in which the chicks were reared. The date on which the first egg in each nest hatched (hatchdate) is recorded together with the sex (sex) of the individuals.

}
\usage{data("DT_btdata")}
\format{
  The format is:
 chr "DT_btdata"
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{

# ####=========================================####
# #### For CRAN time limitations most lines in the 
# #### examples are silenced with one '#' mark, 
# #### remove them and run the examples
# ####=========================================####
# ####=========================================####
# ####=========================================####
# #### EXAMPLE 1
# #### simple example
# ####=========================================####
# ####=========================================####
# data(DT_btdata)
# DT <- DT_btdata
# head(DT)
# mix4 <- mmes(tarsus ~ sex,
#              random = ~ dam + fosternest,
#              rcov=~units,
#              data = DT)
# summary(mix4)$varcomp

# MULTI-TRAIT EXAMPLE
# head(DT)
# DT2 <- stackTrait(DT, traits = c("tarsus","back"))
# head(DT2$long)
# DT2$long=DT2$long[with(DT2$long, order(trait)), ]
# 
# mix3 <- mmes(valueS ~ trait:sex - 1, henderson=TRUE,
#                 random = ~ vsm(usm(trait),ism(dam)) +  
#                            vsm(usm(trait), ism(fosternest)),
#                 rcov= ~ vsm(dsm(trait),ism(units)),
#                 data = DT2$long)
# 
# 
# summary(mix3)$varcomp
# #### calculate the genetic correlation
# cov2cor(mix3$theta[[1]])
# cov2cor(mix3$theta[[2]])

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmes}}  }
