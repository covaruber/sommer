## ----setup, include=FALSE-----------------------------------------------------
library(sommer)

## -----------------------------------------------------------------------------
ar1 <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
round(ar1(5, 0.6), 3)

## -----------------------------------------------------------------------------
set.seed(3)
nRange <- 16; nRow <- 24; nGeno <- 96
g <- setNames(rnorm(nGeno, sd=3), sprintf("G%02d", seq_len(nGeno)))

field <- expand.grid(range=seq_len(nRange), row=seq_len(nRow))
field$geno <- sample(rep(names(g), each=4))
trueField <- as.vector(t(chol(8 * kronecker(ar1(nRow, 0.4), ar1(nRange, 0.7)))) %*%
                         rnorm(nrow(field)))
field$yield <- 50 + g[field$geno] + trueField + rnorm(nrow(field), sd=1)
field <- field[-sample(nrow(field), 15), ]

field$geno <- factor(field$geno)
field$rangef <- factor(field$range, levels=seq_len(nRange))
field$rowf <- factor(field$row, levels=seq_len(nRow))
field$site <- factor("F1")
head(field)

## -----------------------------------------------------------------------------
fits <- list(
  # M0: no spatial model
  iid    = mmes(yield ~ 1,
  random=~geno,
  data=field,
  control = mmes_control(diagnostics = list(verbose = FALSE))),
  # M1: independent random range and row effects (large-scale strips)
  rowcol = mmes(yield ~ 1,
  random=~geno + rangef + rowf,
  data=field,
  control = mmes_control(diagnostics = list(verbose = FALSE))),
  # M2: a smooth 2D surface built from tensor-product P-splines
  spline = mmes(yield ~ 1,
  random=~geno + vsm(ism(spl2Dc(row, range)$Z$`A:all`)),
  data=field,
  control = mmes_control(diagnostics = list(verbose = FALSE))),
  # M3: correlated residuals (R side)
  ar1Res = mmes(yield ~ 1,
  random=~geno,
  rcov=~vsm(ar1m(rangef), ar1m(rowf), ism(units)),
  data=field,
  control = mmes_control(diagnostics = list(verbose = FALSE))),
  # M4: a random AR1 x AR1 spatial field (G side) plus an iid nugget
  ar1Fld = mmes(yield ~ 1,
  random=~geno + vsm(ar1m(rangef), ar1m(rowf), ism(site)),
  rcov=~units,
  data=field,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
)

## -----------------------------------------------------------------------------
compare <- t(sapply(fits, function(m){
  blup <- m$uList[[1]][, 1]
  c(logLik = tail(as.numeric(m$llik), 1),
    AIC = m$AIC,
    genVar = covparams_mmes(m)$estimate[1],
    accuracy = cor(blup[names(g)], g))
}))
knitr::kable(round(compare, 3))

## -----------------------------------------------------------------------------
knitr::kable(covparams_mmes(fits$ar1Fld)[, c("factor", "parameter", "estimate")], digits=3)

## -----------------------------------------------------------------------------
data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
DT$row <- as.numeric(as.character(DT$row))
DT$col <- as.numeric(as.character(DT$col))
DT$R <- as.factor(DT$row)
DT$C <- as.factor(DT$col)

# SPATS MODEL
# m1.SpATS <- SpATS(response = "Y",
#                   spatial = ~ PSANOVA(col, row, nseg = c(14,21), degree = 3, pord = 2),
#                   genotype = "V", fixed = ~ 1,
#                   random = ~ R + C, data = DT,
#                   control = list(tolerance = 1e-04))
# 
# summary(m1.SpATS, which = "variances")
# 
# Spatial analysis of trials with splines 
# 
# Response:                   Y         
# Genotypes (as fixed):       V         
# Spatial:                    ~PSANOVA(col, row, nseg = c(14, 21), degree = 3, pord = 2)
# Fixed:                      ~1        
# Random:                     ~R + C    
# 
# 
# Number of observations:        72
# Number of missing data:        0
# Effective dimension:           17.09
# Deviance:                      483.405
# 
# Variance components:
#                   Variance            SD     log10(lambda)
# R                 1.277e+02     1.130e+01           0.49450
# C                 2.673e-05     5.170e-03           7.17366
# f(col)            4.018e-15     6.339e-08          16.99668
# f(row)            2.291e-10     1.514e-05          12.24059
# f(col):row        1.025e-04     1.012e-02           6.59013
# col:f(row)        8.789e+01     9.375e+00           0.65674
# f(col):f(row)     8.036e-04     2.835e-02           5.69565
# 
# Residual          3.987e+02     1.997e+01 

# SOMMER MODEL
M <- spl2Dmats(x.coord.name = "col", y.coord.name = "row", data=DT, 
               nseg =c(14,21), degree = c(3,3), penaltyord = c(2,2) 
               )
mix <- mmes(Y~V,
  random=~ R + C + vsm(ism(M$fC)) + vsm(ism(M$fR)) +
              vsm(ism(M$fC.R)) + vsm(ism(M$C.fR)) +
              vsm(ism(M$fC.fR)),
  rcov=~units,
  data=M$data,
  control = mmes_control(formulation = "henderson",
  diagnostics = list(verbose = FALSE)))
summary(mix)$varcomp


## -----------------------------------------------------------------------------
simTrial <- function(trial, nRange, nRow, rhoRange, rhoRow, sigma2, g){
  plots <- expand.grid(range=seq_len(nRange), row=seq_len(nRow))
  plots$geno <- sample(rep(names(g), length.out=nrow(plots)))
  K <- sigma2 * kronecker(ar1(nRow, rhoRow), ar1(nRange, rhoRange))
  plots$yield <- 50 + g[plots$geno] + as.vector(t(chol(K)) %*% rnorm(nrow(plots)))
  plots$trial <- trial
  plots[-sample(nrow(plots), round(0.05 * nrow(plots))), ]
}

set.seed(2026)
gMET <- setNames(rnorm(80, sd=3), sprintf("G%02d", 1:80))
MET <- rbind(simTrial("T1", 16, 20, 0.8, 0.2,  9, gMET),
             simTrial("T2", 12, 24, 0.1, 0.7, 16, gMET),
             simTrial("T3", 14, 18, 0.5, 0.5,  4, gMET))
MET$trial <- factor(MET$trial)
MET$geno <- factor(MET$geno)
MET$range <- factor(MET$range, levels=1:16)
MET$row <- factor(MET$row, levels=1:24)
table(MET$trial)

## -----------------------------------------------------------------------------
mShared <- mmes(yield ~ trial,
  random=~geno,
  rcov=~vsm(dsm(trial), ar1m(range), ar1m(row), ism(units)),
  data=MET,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
knitr::kable(covparams_mmes(mShared)[, c("factor", "parameter", "estimate")], digits=3)

## -----------------------------------------------------------------------------
mDsum <- mmes(yield ~ trial,
  random=~geno,
  rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial),
  data=MET,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
knitr::kable(covparams_mmes(mDsum)[, c("factor", "section", "parameter", "estimate")],
             digits=3)

## -----------------------------------------------------------------------------
logLikShared <- tail(as.numeric(mShared$llik), 1)
logLikDsum <- tail(as.numeric(mDsum$llik), 1)
LR <- 2 * (logLikDsum - logLikShared)
c(LR = LR, df = 4, p.value = pchisq(LR, df=4, lower.tail=FALSE))

## -----------------------------------------------------------------------------
MET2 <- MET
MET2$range[MET2$trial == "T3"] <- NA  # T3 has no field coordinates
MET2$row[MET2$trial == "T3"] <- NA

mLevels <- mmes(yield ~ trial,
  random=~geno,
  rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)),
                            by=trial, levels=c("T1", "T2")),
  data=MET2,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
knitr::kable(covparams_mmes(mLevels)[, c("factor", "section", "parameter", "estimate")],
             digits=3)
c(used = nrow(mLevels$y), available = nrow(MET2))

## ----eval=FALSE---------------------------------------------------------------
# random = ~ geno + vsm(dsm(trial), ar1m(range), ar1m(row), ism(site))

