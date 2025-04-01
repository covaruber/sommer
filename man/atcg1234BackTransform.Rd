\name{atcg1234BackTransform}
\alias{atcg1234BackTransform}
\title{
Letter to number converter
}
\description{
This function was designed to help users back transform the numeric marker matrices from the function atcg1234 into letters. 

}
\usage{
atcg1234BackTransform(marks, refs)
}
\arguments{
  \item{marks}{a centered marker matrix coming from atcg1234.}
  \item{refs}{a 2 x m matrix for m markers (columns) and 2 rows where the reference and alternate alleles for each marker are indicated.}
  
}
\value{
\describe{

\item{markers}{ a new marker matrix leter coded according to the reference allele matrix.}

}
}

\author{
Giovanny Covarrubias-Pazaran
}

\examples{
data(DT_polyploid)
genotypes <- GT_polyploid
genotypes[1:5,1:5] # look the original format

# ####=================================================####
# #### convert markers to numeric format polyploid potatoes
# ####=================================================####
# numo <- atcg1234(data=genotypes, ploidy=4)
# numo$M[1:5,1:5]
# numob <- atcg1234BackTransform(marks =  numo$M, refs =  numo$ref.alleles)
# numob[1:4,1:4]
# 
# ####=================================================####
# #### convert markers to numeric format diploid rice lines
# #### single letter code for inbred lines from GBS pipeline
# #### A=AA, T=TT, C=CC, G=GG
# ####=================================================####
# data(DT_rice)
# X <- GT_rice; X[1:5,1:5]; dim(X)
# numo2 <- atcg1234(data=X, ploidy=2)
# numo2$M[1:5,1:5]
# Xb <- atcg1234BackTransform(marks= numo2$M, refs= numo2$ref.alleles)
# Xb[1:4,1:4]


}