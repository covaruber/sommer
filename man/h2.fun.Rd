\name{h2.fun}
\alias{h2.fun}
\title{
Obtain heritabilities with three different methods
}
\description{
Obtain heritabilities based on three different methods; Cullis et al. (2006), Oakey at al. (2006), and line-mean h2 (Falconer, 1995).

Cullis at al. (2006): h2= 1 - (PEM.mu/2*Vg)

Oakey at al. (2006): h2= 1 - (tr((0.5*G^-1) * Czz)/m)

Falconer (1995): h2 = Vg / (Vg+(Ve/r*e))

where "PEV.mu" is the average prediction error variance for the genetic term, "Vg" and "Ve" are the genetic and residual variance respectively estimated by REML, "G^-1" is the inverse of A*Vg where "A" is the additive relationship matrix, "Czz" is the prediction error variance for the genetic term, "m" is the number of test lines, "r" is the replicates per environment and "e" the number of environments.
}
\usage{

h2.fun(object, data, gTerm=NULL, eTerm=NULL, md=NULL)
      
}
\arguments{
  \item{object}{a model fitted with the mmer or mmer2 functions.}
  \item{data}{the dataset used to fit the model provided in the object argument.}
  \item{gTerm}{a character vector specifying the genetic terms fitted in the model.}
  \item{eTerm}{a character vector specifying the environment term fitted in the model.}
  \item{md}{a numeric value to multiply the genetic variance in the heritability formulas (see details below). If NULL it will use the mean value of the diagonal from the genomic or relationship matrix.}
}
\details{
Please see the description or go to the canonical papers where methods are explained with more detail.

For each level from the eTerm (environment) the heritability is calculated as:

h2.stdr= Vg/(Vg+Ve/(ne*nr))
h2.cullis= 1 - (PEV/(md*Vg))
h2.oakey= 1 - tr[ ( C22 (Gi*1/md) / m ) ] 

where "Vg" refers to the genotype variance "Ve" the error variance, "ne" number of environments, "nr" number of replicates, "PEV" is the predicted error variance for the genotype (gTerm), "md" is the mean value from the diagonal of the relationship (pedigree or genomic) matrix "G"" and "m" is the number of lines, "Gi" is the inverse of the relationship matrix.
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Oakey, Helena, et al. "Joint modeling of additive and non-additive genetic line effects in single field trials." Theoretical and Applied Genetics 113.5 (2006): 809-819.

Cullis, Brian R., Alison B. Smith, and Neil E. Coombes. "On the design of early generation variety trials with correlated data." Journal of Agricultural, Biological, and Environmental Statistics 11.4 (2006): 381-393.

Falconer, Douglas S., Trudy FC Mackay, and Richard Frankham. "Introduction to quantitative genetics (4th edn)." Trends in Genetics 12.7 (1996): 280.
}
\examples{

data(DT_example)
head(DT)
# ####=========================================####
# #### fit the mixed model (very heavy model)
# ####=========================================####
# ans1 <- mmer(Yield~Env,
#               random=~vs(ds(Env),Name) + vs(ds(Env),Block),
#               rcov=~vs(ds(Env),units),
#               data=DT)
# summary(ans1)

}
\seealso{
\code{\link{sommer}}
}