\name{DT_rice}
\alias{DT_rice}
\alias{GT_rice}
\alias{GTn_rice}
\title{Rice lines dataset}
\description{

Information from a collection of 413 rice lines.  The DT_rice data set is from 
Rice Diversity Org. Program. The lines are genotyped with 36,901 SNP markers and phenotyped for more than 30 traits. This data set was included in the package to play with it. If using it for your research make sure you cite the original publication from Zhao et al.(2011).

}

\usage{
  data(DT_rice)
}

\format{
 RicePheno contains the phenotypes
 RiceGeno contains genotypes letter code
 RiceGenoN contains the genotypes in numerical code using atcg1234 converter function
}

\source{
  Rice Diversity Organization
  http://www.ricediversity.org/data/index.cfm.
}

\references{

Keyan Zhao, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, M. Rafiqul Islam, Andy Reynolds, Jason Mezey, Anna M. McClung, Carlos D. Bustamante & Susan R. McCouch (2011). Genome-wide association mapping reveals a rich genetic architecture of complex traits in Oryza sativa. Nat Comm 2:467 DOI: 10.1038/ncomms1467, Published Online 13 Sep 2011.

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
data(DT_rice)
# DT <- DT_rice
# GT <- GT_rice
# GTn <- GTn_rice
# head(DT)
# M <- atcg1234(GT)
# A <- A.mat(M$M)
# mix <- mmer(Protein.content~1,
#             random = ~vsr(geno, Gu=A) + geno,
#             rcov=~units,
#             data=DT)
# summary(mix)$varcomp
# 
# Ai <- as(solve(A + diag(1e-6,ncol(A),ncol(A))), Class="dgCMatrix")
# mixb <- mmec(Protein.content~1,
#             random = ~vsc(isc(geno), Gu=Ai) + geno,
#             rcov=~units,
#             data=DT)
# summary(mixb)$varcomp
}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} }
