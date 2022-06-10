\name{DT_sleepstudy}
\alias{DT_sleepstudy}
\docType{data}
\title{
Reaction times in a sleep deprivation study
}
\description{
The average reaction time per day for subjects in a sleep deprivation study. On day 0 the subjects had their normal amount of sleep. Starting that night they were restricted to 3 hours of sleep per night. The observations represent the average reaction time on a series of tests given each day to each subject. Data from sleepstudy to see how lme4 models can be translated in sommer.

}
\usage{data("DT_sleepstudy")}
\format{
  The format is:
 chr "DT_sleepstudy"
}
\source{
These data are from the study described in Belenky et al. (2003), for the sleep deprived group and for the first 10 days of the study, up to the recovery period.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Gregory Belenky et al. (2003) Patterns of performance degradation and restoration during sleep restrictions and subsequent recovery: a sleep dose-response study. Journal of Sleep Research 12, 1-12.

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####
# library(lme4)
data(DT_sleepstudy)
DT <- DT_sleepstudy
head(DT)
##################################
## lme4
# fm1 <- lmer(Reaction ~ Days + (1 | Subject), data=DT)
# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject, 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
## lme4
# fm1 <- lmer(Reaction ~ Days + (Days || Subject), data=DT)
# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject + vsr(Days, Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
## lme4
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=DT)
# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
## sommer
## no equivalence in sommer to find the correlation between the 2 vc
## this is the most similar which is equivalent to (intercept || slope)
fm2 <- mmer(Reaction ~ Days,
            random= ~ Subject + vsr(Days, Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
## lme4
# fm1 <- lmer(Reaction ~ Days + (0 + Days | Subject), data=DT)
# vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
## sommer
fm2 <- mmer(Reaction ~ Days,
            random= ~ vsr(Days, Subject), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} }
