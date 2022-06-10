\name{DT_expdesigns}
\alias{DT_expdesigns}
\docType{data}
\title{
  Data for different experimental designs
}
\description{
  The following data is a list containing data frames for different type of experimental designs relevant in plant breeding:
  
  1) Augmented designs (2 examples)
  
  2) Incomplete block designs (1 example)
  
  3) Split plot design (2 examples)
  
  4) Latin square designs (1 example)
  
  5) North Carolina designs I,II and III
  
  How to fit each is shown at the Examples section. This may help you get introduced to experimental designs relevant to plant breeding. Good luck.
}

\format{
 Different based on the design.
}

\source{
  Datasets and more detail about them can be found in the agricolae package. Here we just show the datasets and how to analyze them using the \code{\link{sommer}} package.
}

\references{
 Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}

\examples{

# #### =================================== ####
# #### ===== Augmented Block Design 1 ==== ####
# #### =================================== ####
# data(DT_expdesigns)
# DT <- DT_expdesigns
# names(DT)
# data1 <- DT$au1
# head(data1)
# ## response variable: "yield"
# ## check indicator: "entryc" ('nc' for all unreplicated, but personal.name for checks)
# ## blocking factor: "block"
# ## treatments, personal names for replicated and non-replicated: "trt"
# ## check no check indicator: "new"
# mix1 <- mmer(yield~entryc,
#              random=~block+trt,
#              rcov=~units, tolpar = 1e-6,
#              data=data1)
# summary(mix1)$varcomp
# 
# mix1b <- mmec(yield~entryc,
#              random=~block+trt,
#              rcov=~units, tolParConv = 1e-6,
#              data=data1)
# summary(mix1b)$varcomp

}

