\name{neMarker}
\alias{neMarker}
\title{Effective population size based on marker matrix}
\description{
  `neMarker` uses a marker matrix to approximate the effective population size (Ne) by discovering how many individuals are needed to sample all possible alleles in a population.
}
\usage{
  neMarker(M, maxNe=100, maxMarker=1000, nSamples=5)
}

\arguments{
  \item{M}{marker matrix coded in a numerical faashion (any allele dosage is fine).}
  \item{maxNe}{maximum number of effective population size to be calculated.}
  \item{maxMarker}{maximum number of markers to use for the analysis.}
  \item{nSamples}{number of individuals to sample for the Ne calculation.}
}
\value{
  \describe{
    
    \item{$S3}{ A vector with allele coverage based on different number of individuals  }
  }
}
\references{
  
  Not based on any theory published yet but in a solid intuition on what is really important for a breeding program when we ask what is the effective population size
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

# data(DT_cpdata) # Madison cranberries
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# M <- GT
# # run the function
# ne <- neMarker(M, maxNe = 30, nSamples = 10)
# ################
# data(DT_technow) # maize
# M <- Md_technow # dent
# M <- (M*2) - 1
# M <- M + 1
# # run the function
# ne <- neMarker(M, maxNe = 60, nSamples = 10)
# ##
# M <- Mf_technow # flint
# M <- (M*2) - 1
# M <- M + 1
# # run the function
# ne <- neMarker(M, maxNe = 60, nSamples = 10)
# ################
# data(DT_wheat) # cimmyt wheat
# M <- GT_wheat + 1
# # run the function
# ne <- neMarker(M, maxNe = 60, nSamples = 10)
# ###############
# data(DT_rice) # Zhao rice
# M <- atcg1234(GT_rice)$M
# # run the function
# ne <- neMarker(M,  maxNe = 60, nSamples = 10)
# ###############
# data(DT_polyploid) # endelman potatoes
# M <- atcg1234(data=GT_polyploid, ploidy=4)$M
# # run the function
# ne <- neMarker(M,  maxNe = 60, nSamples = 10)
# 
# library(ggplot2) #For making plots
# ggplot(ne,aes(x=Ne,y=allelesCovered))+
#   geom_ribbon(aes(x=Ne,ymin=allelesCovered-allelesCoveredSe,
#                   ymax=allelesCovered+allelesCoveredSe),
#                   alpha=0.2,linetype=0)+
#   geom_line(linewidth=1)+
#   guides(alpha=FALSE)+
#   theme_bw()+ 
#   scale_x_continuous("Individual number")+
#   scale_y_continuous("Allele coverage")  + 
#                geom_hline(yintercept = 0.95) + 
#                geom_hline(yintercept = 0.975)

}
\seealso{The core functions of the package \code{\link{mmec}} }