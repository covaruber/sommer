pkgname <- "sommer"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "sommer-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('sommer')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("A.mat")
### * A.mat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: A.mat
### Title: Additive relationship matrix
### Aliases: A.mat

### ** Examples

####=========================================####
#### random population of 200 lines with 1000 markers
####=========================================####
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

A <- A.mat(X)

####=========================================####
#### take a look at the Genomic relationship matrix 
#### (just a small part)
####=========================================####
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(A[1:15,1:15], col = colfunc(100),Colv = "Rowv")
# str(hv)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("A.mat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("APY")
### * APY

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: APY
### Title: Construct an APY Relationship Precision Matrix
### Aliases: APY

### ** Examples

## Marker-derived GBLUP for genotyped animals.
set.seed(2)
ids <- paste0("animal", seq_len(12))
markers <- matrix(sample(c(-1, 0, 1), 12 * 80, replace = TRUE), nrow = 12,
                  dimnames = list(ids, NULL))
G <- A.mat(markers)

## Blend slightly with the identity so this small illustrative G is positive
## definite, then construct APY precision for a six-animal core.
G <- 0.95 * G + 0.05 * diag(nrow(G))
dimnames(G) <- list(ids, ids)
Ginv.apy <- APY(G, core = ids[seq_len(6)])

phenotypes <- data.frame(
  id = factor(rep(ids, each = 3), levels = ids),
  y = rnorm(length(ids) * 3)
)
fit <- mmes(y ~ 1,
  random = ~vsm(ism(id), Gu = Ginv.apy),
  rcov = ~units,
  data = phenotypes,
  control = mmes_control(diagnostics = list(verbose = FALSE)))

## In a single-step evaluation, this APY approximation supplies G^{-1} for
## the genotyped block correction; it is not itself the complete H^{-1}.



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("APY", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("D.mat")
### * D.mat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: D.mat
### Title: Dominance relationship matrix
### Aliases: D.mat

### ** Examples

####=========================================####
#### EXAMPLE 1
####=========================================####
####random population of 200 lines with 1000 markers
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}

D <- D.mat(X)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("D.mat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("DEPRECATED_GWAS")
### * DEPRECATED_GWAS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GWAS
### Title: Genome wide association study analysis
### Aliases: GWAS

### ** Examples


####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
#####========================================####
##### potato example
#####========================================####
# 
# data(DT_polyploid, package="enhancer")
# DT <- DT_polyploid
# GT <- GT_polyploid
# MP <- MP_polyploid
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
# A <- A.mat(marks)
# ###=========================================####
# ### run it as GWAS model
# ###=========================================####
# ans2 <- GWAS(tuber_shape~1,
#              random=~vsr(Name,Gu=A),
#              rcov=~units,
#              gTerm = "u:Name",
#              M=marks, data=DT2)
# plot(ans2$scores[,1])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("DEPRECATED_GWAS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("DEPRECATED_VS")
### * DEPRECATED_VS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vs
### Title: variance structure specification
### Aliases: vs

### ** Examples


data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

## ============================ ##
## example to without structure
## ============================ ##

mix <- mmer(Yield~Env, 
            random= ~ vs(Name),
            rcov=~ vs(units),
            data=DT)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("DEPRECATED_VS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("DEPRECATED_VSR")
### * DEPRECATED_VSR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vsr
### Title: variance structure specification
### Aliases: vsr

### ** Examples


data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

## ============================ ##
## example to without structure
## ============================ ##

mix <- mmer(Yield~Env, 
            random= ~ vsr(Name),
            rcov=~ vsr(units),
            data=DT)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("DEPRECATED_VSR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("E.mat")
### * E.mat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: E.mat
### Title: Epistatic relationship matrix
### Aliases: E.mat

### ** Examples

####=========================================####
####random population of 200 lines with 1000 markers
####=========================================####
X <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  X[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}

E <- E.mat(X, type="A#A") 
# if heterozygote markers are present can be used "A#D" or "D#D"



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("E.mat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("H.mat")
### * H.mat

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: H.mat
### Title: Combined relationship matrix H
### Aliases: H.mat

### ** Examples

####=========================================####
####random population of 200 lines with 1000 markers
####=========================================####
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- sample(c(-1,0,0,1), size=1000, replace=TRUE)
}
rownames(M) <- 1:nrow(M)
v <- sample(1:nrow(M),100)
M2 <- M[v,]

A <- A.mat(M) # assume this is a pedigree-based matrix for the sake of example
G <- A.mat(M2)

H <- H.mat(A,G)
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(H[1:15,1:15], col = colfunc(100),Colv = "Rowv")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("H.mat", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNR")
### * MNR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNR
### Title: Multivariate Newton-Raphson algorithm
### Aliases: MNR

### ** Examples


data(DT_cpdata, package="enhancer")
DT <- DT_cpdata

## No test: 
# response matrix
y <- cbind(imputev(DT$Yield), imputev(DT$Firmness))
# fixed effect incidence matrix
X <- model.matrix(~Rowf, data=DT)
# random effect incidence matrix
Z <- Matrix::sparse.model.matrix(~id - 1, data=DT)
colnames(Z) <- gsub("id","",colnames(Z))
# covariance for random effect
GT <- GT_cpdata
A <- A.mat(GT) # additive relationship matrix
A <- A[colnames(Z),colnames(Z)] # make sure of the order
A <- A + diag(1e-4,nrow(A), nrow(A))
# residual effect incidence matrix (dimensions equal nrow(y))
R1 <- Matrix::Diagonal(n=nrow(y)) 
# weights matrix (dimensions equal nrow(y))
W <- diag(nrow(y)) 
# some other parameters
maxIter=3
stepWeight <- rep(0.9, maxIter)
stepWeight[1:2] <- c(0.5, 0.7)
emWeights <- rep(0,maxIter)
# model fit
res <- MNR( Y=y, # multi-trait response
            X=list(X), Gx=list(diag(2)), # fixed effects
            Z=list(Z), K=list(A), # random effects
            R=list(R1), # residual effects
            Ge=list( (diag(2)*.3)+.15 , diag(2)*0.75 ), # inital vc 
            GeI=list(unsm2(2),diag(2)), # vc constraints
            W=W, isInvW=TRUE, # weights for records
            iters=maxIter,
            tolpar=1e-4, tolparinv=1e-6, # tolerance
            ai=FALSE, pev=FALSE, # algorithm specifics
            verbose=TRUE, 
            retscaled=FALSE, 
            stepweight=stepWeight, # second derivatives weights
            emweight=emWeights # em update weights
)
res$sigma
res$Beta
res$U[[1]]
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("antem")
### * antem

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: antem
### Title: Antedependence covariance structure
### Aliases: antem

### ** Examples

## Not run: 
##D vsm(antem(time, order=2), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("antem", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ar1m")
### * ar1m

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ar1m
### Title: First-order autoregressive covariance structure
### Aliases: ar1m

### ** Examples

## Not run: 
##D vsm(ar1m(time), ism(id))
##D 
##D # Separable row-by-range residual shared across trials:
##D # rcov = ~ vsm(dsm(trial), ar1m(range), ar1m(row), ism(units))
##D # Trial-specific correlations and variances:
##D # rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ar1m", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ar2m")
### * ar2m

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ar2m
### Title: Second-order autoregressive covariance structure
### Aliases: ar2m

### ** Examples

## Not run: 
##D vsm(ar2m(time), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ar2m", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ar3m")
### * ar3m

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ar3m
### Title: Third-order autoregressive covariance structure
### Aliases: ar3m

### ** Examples

## Not run: 
##D vsm(ar3m(time), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ar3m", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("atm")
### * atm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: atm
### Title: Selected-level diagonal covariance structure
### Aliases: atm

### ** Examples

## Not run: 
##D vsm(atm(environment, c("E1","E3")), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("atm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("atr")
### * atr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: atr
### Title: Selected-level diagonal covariance structure for 'mmer' and
###   'vsr'
### Aliases: atr

### ** Examples

# Select levels B and C
A <- atr(factor(c("A","B","A","C")), levs = c("B","C"))
A$Z
A$thetaC



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("atr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("binm")
### * binm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: binm
### Title: Binomial proportion response with its number of trials
### Aliases: binm
### Keywords: models

### ** Examples

set.seed(1)
d <- data.frame(x = rnorm(60), g = factor(rep(1:12, each = 5)))
d$n <- sample(5:10, 60, TRUE)
d$s <- rbinom(60, d$n, plogis(0.3 * d$x))
fit <- mmes(binm(s, trials = n) ~ x,
  random = ~ g,
  data = d,
  family = binomial(),
  control = mmes_control(update = list(blendWeight = 0),
  diagnostics = list(verbose = FALSE)))
fit$b



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("binm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("car")
### * car

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: car
### Title: Proper conditional autoregressive spatial covariance structure
### Aliases: car

### ** Examples

## Not run: 
##D vsm(car(location, W), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("car", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("corImputation")
### * corImputation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: corImputation
### Title: Imputing a matrix using correlations
### Aliases: corImputation

### ** Examples


####################################
### imputing genotype data example
####################################
# data(DT_cpdata, package="enhancer")
# X <- GT_cpdata
# # add missing data
# v <- sample(1:length(X), 500)
# Xna <- X
# Xna[v]<- NA
# ## impute (can take some time)
# Y <- corImputation(wide=Xna, Gu=NULL, nearest=20, roundR=TRUE) 
# cm <- table(Y$imputed[v],X[v])
# ## calculate accuracy
# sum(diag(cm))/length(v)
####################################
### imputing phenotypic data example
####################################
# data(DT_h2, package="enhancer")
# X <- reshape(DT_h2[,c("Name","Env","y")], direction = "wide", idvar = "Name",
#                 timevar = "Env", v.names = "y", sep= "_")
# rownames(X) <- X$Name
# X <- as.matrix(X[,-1])
# head(X)
# # add missing data
# v <- sample(1:length(X), 50)
# Xna <- X
# Xna[v]<- NA
# ## impute
# Y <- corImputation(wide=Xna, Gu=NULL, nearest=20, roundR=TRUE)
# plot(y=Y$imputed[v],x=X[v], xlab="true",ylab="predicted")
# cor(Y$imputed[v],X[v], use = "complete.obs")





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("corImputation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("corgm")
### * corgm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: corgm
### Title: General positive-definite correlation structure
### Aliases: corgm

### ** Examples

data(DT_example, package = "enhancer")
DT <- DT_example
# genetic correlations between environments with a homogeneous variance
m1 <- mmes(Yield ~ Env,
  random = ~ vsm(corgm(Env), ism(Name)),
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
covparams_mmes(m1, 1)

# heterogeneous variances (equivalent to usm, parameterized as corr/var)
m2 <- mmes(Yield ~ Env,
  random = ~ vsm(corgm(Env, variance = "heterogeneous"), ism(Name)),
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
covparams_mmes(m2, 1)
covmatrix_mmes(m2, 1)$correlation



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("corgm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("covm")
### * covm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: covm
### Title: Covariance Between Two Random Effects
### Aliases: covm

### ** Examples


data(DT_ige, package = "enhancer")
DT <- DT_ige

## Correlate two random effects with identity relationship matrices

covRes <- with(
DT,
covm(
vsm(ism(focal)),
vsm(ism(neighbour))
)
)

str(covRes)

## Custom initial variance-covariance matrix

covRes2 <- with(
DT,
covm(
vsm(ism(focal)),
vsm(ism(neighbour)),
theta = matrix(
c(0.5, 0.1,
0.1, 0.3),
2, 2
)
)
)

## Relationship precision matrix

## Ai must have matching row and column names and be marked as inverse.

##

## attr(Ai, "inverse") <- TRUE

##

## covRes3 <- with(

## DT,

## covm(

## vsm(ism(focal), Gu = Ai),

## vsm(ism(neighbour), Gu = Ai)

## )

## )

## See the DT_ige help page for a complete model-fitting example.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("covm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("covmatrix_mmes")
### * covmatrix_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: covmatrix_mmes
### Title: Covariance and correlation matrix of a fitted term with standard
###   errors
### Aliases: covmatrix_mmes
### Keywords: models

### ** Examples

set.seed(17)
n <- 60
DT <- data.frame(id = factor(rep(seq_len(n), each = 3)))
u <- matrix(rnorm(n * 2), n) 
DT$trait_a <- u[DT$id, 1] + rnorm(nrow(DT), sd = 0.7)
DT$trait_b <- u[DT$id, 2] + rnorm(nrow(DT), sd = 1)
DTL <- stackTraits(DT, traits = c("trait_a", "trait_b"), keep = "id")
fit <- mmes(value ~ trait,
  random = ~ vsm(usm(trait), ism(id)),
  rcov = ~ vsm(usm(trait), ism(record)),
  data = DTL,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
G <- covmatrix_mmes(fit, 1)       # genetic
G$covariance
G$correlation
G$correlation.se
R <- covmatrix_mmes(fit, 2)       # residual
R$correlation



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("covmatrix_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("covparams_mmes")
### * covparams_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: covparams_mmes
### Title: Extract native-scale covariance parameters from an mmes fit
### Aliases: covparams_mmes

### ** Examples

## Not run: 
##D fit <- mmes(Yield ~ Env,
##D             random = ~ vsm(dsm(Env), ism(Name)),
##D             rcov = ~ units, data = DT_example)
##D covparams_mmes(fit)
##D fit$covParNative
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("covparams_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("covparams_mmes_se")
### * covparams_mmes_se

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: covparams_mmes_se
### Title: Native-scale covariance parameters and standard errors
### Aliases: covparams_mmes_se

### ** Examples

## Not run: 
##D fit <- mmes(Yield ~ Env,
##D             random = ~ vsm(dsm(Env), ism(Name)),
##D             rcov = ~ units, data = DT_example)
##D covparams_mmes_se(fit)
##D fit$covParNativeSE
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("covparams_mmes_se", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("csm")
### * csm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: csm
### Title: Compound-symmetry covariance structure
### Aliases: csm

### ** Examples

## Not run: 
##D vsm(csm(environment), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("csm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("csr")
### * csr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: csr
### Title: User-specified covariance-component structure for 'mmer' and
###   'vsr'
### Aliases: csr

### ** Examples

# A custom two-level covariance-component pattern
x <- factor(c("A","B","A","B"))
mm <- matrix(c(1, 2,
               2, 1), 2, 2, byrow = TRUE)
C <- csr(x, mm)
C$Z
C$thetaC



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("csr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dfToMatrix")
### * dfToMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dfToMatrix
### Title: data frame to matrix
### Aliases: dfToMatrix

### ** Examples


library(Matrix)
m <- matrix(1:9,3,3)
m <- tcrossprod(m)

mdf <- as.data.frame(as.table(m))
mdf

dfToMatrix(mdf, row = "Var1", column = "Var2", 
            value = "Freq",returnInverse=FALSE )





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dfToMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dsm")
### * dsm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dsm
### Title: Diagonal heterogeneous covariance structure
### Aliases: dsm

### ** Examples

## Not run: 
##D vsm(dsm(environment), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dsm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dsr")
### * dsr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dsr
### Title: Diagonal covariance structure for 'mmer' and 'vsr'
### Aliases: dsr

### ** Examples

# Diagonal covariance structure across levels
D <- dsr(factor(c("A","B","A","C")))
D$Z
D$thetaC



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dsr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dsumm")
### * dsumm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dsumm
### Title: Direct sum of residual covariance structures across sections
### Aliases: dsumm

### ** Examples

## Not run: 
##D 
##D mix <- mmes(yield ~ trial,
##D             rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by = trial),
##D             data = dat)
##D covparams_mmes(mix)
##D 
##D # Spatial structure only for two trials; remaining trials are iid.
##D mix2 <- mmes(yield ~ trial,
##D              rcov = ~ dsumm(vsm(ar1m(range), ar1m(row), ism(units)),
##D                             by = trial, levels = c("T1", "T2")),
##D              data = dat)
##D 
##D # Site-specific unstructured covariance among environments.
##D mix3 <- mmes(y ~ site,
##D              rcov = ~ dsumm(vsm(usm(env), ism(units)), by = site),
##D              data = dat2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dsumm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fam")
### * fam

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fam
### Title: Factor-analytic covariance structure
### Aliases: fam

### ** Examples

## Not run: 
##D vsm(fam(environment, k=2), ism(genotype))
##D vsm(fam(environment, k=2, variance="homogeneous"), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fam", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("familym")
### * familym

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: familym
### Title: Different distribution families for the traits of a multi-trait
###   model
### Aliases: familym
### Keywords: models

### ** Examples

set.seed(12)
ng <- 100; nr <- 4
G <- matrix(c(1, 0.4, 0.4, 0.8), 2)
u <- matrix(rnorm(ng * 2), ng) %*% chol(G)
d <- data.frame(id = factor(rep(paste0("g", 1:ng), each = nr)))
d$disease <- rbinom(nrow(d), 1, plogis(-0.5 + u[as.integer(d$id), 1]))
d$yield <- 10 + u[as.integer(d$id), 2] + rnorm(nrow(d))
L <- stackTraits(d, traits = c("disease", "yield"))

fit <- mmes(value ~ trait,
  random = ~ vsm(usm(trait), ism(id)),
  rcov = ~ vsm(corgm(trait, variance = "heterogeneous"), ism(record)),
  family = familym("trait", disease = binomial(), yield = gaussian()),
  data = L,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
fit$dispersion
covmatrix_mmes(fit, 1)$correlation   # genetic correlation (link scale)
head(fitted(fit))                    # probabilities for disease rows



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("familym", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitted_mmes")
### * fitted_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitted.mmes
### Title: fitted form a LMM fitted with mmes
### Aliases: fitted.mmes fitted.mmer
### Keywords: models

### ** Examples

# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmes(Yield~1,
#               random=~vsm(ism(id),Gu=A)
#                       + Rowf + Colf + spl2Dc(Row,Col),
#               rcov=~units,
#               data=DT)
# 
# ff=fitted(mix1)
# 
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# lattice::wireframe(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted,  
#           aspect=c(61/87,0.4), drape=TRUE,# col.regions = colfunc,
#           light.source=c(10,0,10))
# lattice::levelplot(`u:Row.fitted`~Row*Col, data=ff$dataWithFitted, col.regions = colfunc)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitted_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fixm")
### * fixm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fixm
### Title: fixed indication matrix
### Aliases: fixm

### ** Examples

fixm(4)
fixm(4,2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fixm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_mme_workload")
### * generate_mme_workload

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_mme_workload
### Title: Generate Parameterized MME Workloads for Solver Evaluation
### Aliases: generate_mme_workload

### ** Examples

# Generate a small 1D autoregressive workload
wl <- generate_mme_workload(
  name = "demo_ar1",
  n_ind = 20,
  n_levels = 4,
  structure = "ar1"
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_mme_workload", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ism")
### * ism

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ism
### Title: Identity covariance structure
### Aliases: ism

### ** Examples

## Not run: 
##D vsm(ism(id))
##D vsm(dsm(environment), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ism", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loadings_mmes")
### * loadings_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loadings_mmes
### Title: Extract factor-analytic loadings from a fitted mmes model
### Aliases: loadings_mmes

### ** Examples

## Not run: 
##D mix <- mmes(BLUEs ~ trial,
##D             random = ~ vsm(fam(trial, 2), ism(genotype)),
##D             rcov = ~ units, data = dt)
##D loadings_mmes(mix)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loadings_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mam")
### * mam

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mam
### Title: Moving-average covariance structures of order one or two
### Aliases: mam ma1m ma2m

### ** Examples

## Not run: 
##D vsm(ma1m(time), ism(id))
##D vsm(ma2m(time), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mam", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("maternm")
### * maternm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: maternm
### Title: Matern spatial covariance structure
### Aliases: maternm

### ** Examples

## Not run: 
##D vsm(maternm(cbind(xcoord,ycoord)), ism(genotype))
##D vsm(maternm(time), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("maternm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("meffects_mmes")
### * meffects_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: meffects_mmes
### Title: Marker effects recovered from a fitted relationship-based term
### Aliases: meffects_mmes
### Keywords: models

### ** Examples

set.seed(7)
M <- matrix(sample(c(-1, 0, 1), 40 * 80, TRUE), 40)
rownames(M) <- paste0("g", 1:40); colnames(M) <- paste0("snp", 1:80)
G <- 0.95 * A.mat(M) + 0.05 * diag(40)
Gi <- solve(G); attr(Gi, "inverse") <- TRUE
d <- data.frame(id = factor(rep(rownames(M), each = 2), levels = rownames(M)))
d$y <- rnorm(80)
fit <- mmes(y ~ 1,
  random = ~ vsm(ism(id), Gu = Gi),
  data = d,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
head(meffects_mmes(fit,
  1,
  M,
  blend = 0.05,
  se = TRUE))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("meffects_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("metricm")
### * metricm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: metricm
### Title: Distance-based (metric) spatial correlation structures
### Aliases: metricm
### Keywords: models

### ** Examples

set.seed(1)
co <- cbind(x = runif(30, 0, 10), y = runif(30, 0, 10))
d <- data.frame(x = rep(co[, 1], each = 2), y = rep(co[, 2], each = 2), one = factor(1))
d$yield <- rnorm(60)
fit <- mmes(yield ~ 1,
  random = ~ vsm(metricm(cbind(x, y), model = "power"), ism(one)),
  data = d,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
covparams_mmes(fit)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("metricm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mmer")
### * mmer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mmer
### Title: *m*ixed *m*odel *e*quations for *r* records
### Aliases: mmer

### ** Examples


####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmer(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####

## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmer(Yield~Env,
             random= ~Name + vsr(dsr(Env),Name),
             rcov= ~ vsr(dsr(Env),units),
             data=DT)
summary(ans2)

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmer(Yield~Env,
             random=~ vsr(usr(Env),Name),
             rcov=~vsr(dsr(Env),units), 
             data=DT)
summary(ans3)

## No test: 

####==========================================####
#### Multivariate homogeneous variance models ####
####==========================================####

## Multivariate Compound simmetry (CS) model
DT$EnvName <- paste(DT$Env,DT$Name)
ans4 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(Name, Gtc = unsm(2)) + vsr(EnvName,Gtc = unsm(2)),
              rcov= ~ vsr(units, Gtc = unsm(2)),
              data=DT)
summary(ans4)

####=============================================####
#### Multivariate heterogeneous variance models  ####
####=============================================####

## Multivariate Compound simmetry (CS) + Diagonal (DIAG) model
ans5 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(Name, Gtc = unsm(2)) + vsr(dsr(Env),Name, Gtc = unsm(2)),
              rcov= ~ vsr(dsr(Env),units, Gtc = unsm(2)),
              data=DT)
summary(ans5)

####===========================================####
#### Multivariate unstructured variance models ####
####===========================================####

ans6 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(usr(Env),Name, Gtc = unsm(2)),
              rcov= ~ vsr(dsr(Env),units, Gtc = unsm(2)),
              data=DT)
summary(ans6)

####=========================================####
####=========================================####
#### EXAMPLE SET 2
#### 2 variance components
#### one random effect with variance covariance structure
####=========================================####
####=========================================####

data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
head(DT)
GT[1:4,1:4]
#### create the variance-covariance matrix
A <- A.mat(GT)
#### look at the data and fit the model
mix1 <- mmer(Yield~1,
             random=~vsr(id, Gu=A) + Rowf,
             rcov=~units,
             data=DT)
summary(mix1)$varcomp

#### multi trait example
mix2 <- mmer(cbind(Yield,color)~1,
              random = ~ vsr(id, Gu=A, Gtc = unsm(2)) + # unstructured at trait level
                            vsr(Rowf, Gtc=diag(2)) + # diagonal structure at trait level
                                vsr(Colf, Gtc=diag(2)), # diagonal structure at trait level
              rcov = ~ vsr(units, Gtc = unsm(2)), # unstructured at trait level
              data=DT)
summary(mix2)

## End(No test)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mmer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mmes")
### * mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mmes
### Title: *m*ixed *m*odel *e*quations *s*olver
### Aliases: mmes

### ** Examples


####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)

## Known variances: solve the MME only (e.g. genetic evaluation)
sol1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT, control=mmes_control(task=list(mode="solve")),
             covPar=ans1)
head(sol1$uList[[1]])

## No test: 

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
DT=DT[with(DT, order(Env)), ]
## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ dsumm(vsm(ism(units)), by=Env),
             data=DT)
summary(ans2)

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmes(Yield~Env,
             random=~ vsm(usm(Env),ism(Name)),
             rcov= ~ dsumm(vsm(ism(units)), by=Env),
             data=DT)
summary(ans3)


## End(No test)

## No test: 

## simulated MET
nInds=300
nMarks=1000
nEnvs=3

#random population of nInds lines with nMarks markers
M <- matrix(rep(0,nInds*nMarks),nInds,nMarks)
for (i in 1:nInds) {
  M[i,] <- ifelse(runif(nMarks)<0.5,-1,1)
}
# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))
# get phenotypes
rg <- enhancer::simGECorMat(nEnvs, 1, mu=0.2, v=0.2)
rg <- Matrix::nearPD(rg)$mat
rg # desired genetic correlation
u <- mvrnorm(n = nMarks, rep(0, nEnvs), Sigma=rg) # marker effects
g <- crossprod(t(M),u) # genetic values
h2 <- 0.2  #heritability 
envEffect <- rnorm(nEnvs, mean = abs(rnorm(1))) # environment means
e <- mvrnorm(n = nInds, envEffect, Sigma=(1-h2)/h2*diag(diag(var(g)))) # errors
p = g + e # phenotypes
y <- as.vector(p) # phenotypes as vector
env <- sort(rep(paste0("e",1:nEnvs), nrow(g))) # environment labels
ind <- rep(paste0("i",1:nrow(g)), nEnvs ) # individuals labels
colnames(MMT) <- rownames(MMT) <- paste0("i",1:nrow(g)) # assign dimnames

MMTi <- solve(MMT + diag(1e-4,ncol(MMT),ncol(MMT)))
MMTi[lower.tri(MMTi)] <- t(MMTi)[lower.tri(MMTi)] # fill the lower triangular
MMTi <- as(as(as( MMTi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(MMTi, 'inverse')=TRUE

DT <- data.frame(ind=ind, env=env, y=y)
# fit the models
mixm2 <- mmes(y ~ env,
              random=~ vsm(dsm(env), ism(ind), Gu=MMTi),
              data=DT
)
summary(mixm2)$varcomp

## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mmes_control")
### * mmes_control

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mmes_control
### Title: Modular Numerical Controls for Mixed Models
### Aliases: mmes_control

### ** Examples

ctrl <- mmes_control(
  formulation = "henderson",
  solve = list(method = "pcg", tol = 1e-8),
  score = list(traceProbes = 16L),
  convergence = list(accuracy = "adaptive")
)
print(ctrl)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mmes_control", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ownm")
### * ownm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ownm
### Title: User-defined covariance structure
### Aliases: ownm

### ** Examples

## Not run: 
##D K <- matrix(c(1,.3,.3,1),2,2)
##D vsm(ownm(group, K=K), ism(id))
##D 
##D cfun <- function(p) {
##D   r <- tanh(p[1]); matrix(c(1,r,r,1),2,2)
##D }
##D native <- function(scale, par, factor, absorb_scale=TRUE) {
##D   c(variance=scale, correlation=tanh(par[1]))
##D }
##D vsm(ownm(group, fun=cfun, par=0, native_report=native), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ownm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.monitor")
### * plot.monitor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pmonitor
### Title: plot the change of VC across iterations
### Aliases: pmonitor
### Keywords: models

### ** Examples

data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
head(DT)
m3 <- mmes(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
pmonitor(m3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.monitor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_mmes")
### * plot_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.mmes
### Title: plot form a LMM plot with mmes
### Aliases: plot.mmes plot.mmer
### Keywords: models

### ** Examples

data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
head(DT)
m3 <- mmes(fixed=Y ~ V + N + V:N,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("postPEV")
### * postPEV

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: postPEV
### Title: Post-fit prediction error variances and inverse coefficient
###   matrix
### Aliases: postPEV

### ** Examples

####=========================================####
#### Fit a mixed model without computing PEVs
#### or the complete inverse during estimation
####=========================================####

# mod <- mmes(
#   fixed = Yield ~ 1,
#   random = ~ vsm(ism(Genotype)),
#   rcov = ~ units,
#   data = DT,
#   control = mmes_control(output = list(inverse = 0))
# )

####=========================================####
#### Compute PEVs afterwards using the
#### Takahashi sparse inverse subset
####=========================================####

# mod <- postPEV(mod, mode = 1)
# mod$uPevList

####=========================================####
#### Compute the complete inverse afterwards
####=========================================####

# mod <- postPEV(mod, mode = 2)
# mod$Ci
# mod$uPevList



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("postPEV", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("postVarU")
### * postVarU

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: postVarU
### Title: Post-fit Sampling Covariance of BLUPs
### Aliases: postVarU

### ** Examples

data(DT_example, package="enhancer")
model <- mmes(Yield~Env,
  random=~Name,
  data=DT_example,
  control = mmes_control(output = list(inverse = 0),
  diagnostics = list(verbose = FALSE)))
model <- postVarU(model)
head(model$uVarList[[1]])
full <- postVarU(model, mode=2)
head(diag(full$VarU))

prediction <- predict(model, D="Name", PEV=TRUE, VarU=TRUE)
head(prediction$pvals)
prediction$PEV[1:3, 1:3]
prediction$VarU[1:3, 1:3]
prediction$sampling.vcov[1:3, 1:3]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("postVarU", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict_mmes")
### * predict_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.mmes
### Title: Predict form of a LMM fitted with mmes
### Aliases: predict.mmes predict.mmer
### Keywords: models

### ** Examples


data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
m3 <- mmes(fixed=Y ~ V + N + V:N ,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)

m3 <- postPEV(m3, mode = 2)
 
#############################
## predict means for nitrogen
#############################
Dt <- m3$Dtable; Dt
# first fixed effect just average
Dt[1,"average"] = TRUE
# second fixed effect include
Dt[2,"include"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N")
pp$pvals

#############################
## predict means for variety
#############################

Dt <- m3$Dtable; Dt
# first fixed effect include
Dt[1,"include"] = TRUE
# second fixed effect just average
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="V")
pp$pvals

#############################
## predict means for nitrogen:variety
#############################
# prediction matrix D based on (equivalent to classify in asreml)
Dt <- m3$Dtable; Dt
# first fixed effect include and average
Dt[1,"include"] = TRUE
Dt[1,"average"] = TRUE
# second fixed effect include and average
Dt[2,"include"] = TRUE
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N:V")
pp$pvals

#############################
## select levels through the argument or the Dtable
#############################
Dt <- m3$Dtable
Dt$average[Dt$term %in% c("1", "V", "V:N")] <- TRUE
Dt$include[Dt$term %in% c("N", "V:N")] <- TRUE

p <- predict(m3, D="N", Dtable=Dt,
             levels=list(N=c("0.6", "0.2")))
p$pvals
p$Dtable$levels[[match("N", p$Dtable$term)]]

Dt$levels[[match("N", Dt$term)]] <- c("0.6", "0.2")
p2 <- predict(m3, D="N", Dtable=Dt)
all.equal(p$D, p2$D)

## Covariances of the requested linear combinations
p <- predict(m3, D="N", Dtable=Dt, PEV=TRUE, VarU=TRUE)
p$PEV
p$VarU
p$sampling.vcov
randomColumns <- seq.int(nrow(m3$b)+1L, nrow(m3$bu))
randomD <- p$D[,randomColumns,drop=FALSE]
full <- postVarU(m3, mode=2)
all.equal(unname(p$VarU),
          unname(as.matrix(randomD %*% full$VarU %*% t(randomD))))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("r2")
### * r2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: r2
### Title: Reliability
### Aliases: r2

### ** Examples

####=========================================####
#### Example population
####=========================================####
data(DT_example, package="enhancer")
DT <- DT_example
head(DT)
ans1 <- mmes(Yield~Env, 
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
ans1 <- postPEV(ans1, mode = 1)
rel=r2(ans1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("r2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randef")
### * randef

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randef
### Title: extracting random effects
### Aliases: randef

### ** Examples

# randef(model)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randef", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rrm")
### * rrm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rrm
### Title: Reduced-rank covariance structure
### Aliases: rrm

### ** Examples

## Not run: 
##D vsm(rrm(environment, k=2), ism(genotype))
##D # Former rrm() model: an estimated common specific variance
##D vsm(fam(environment, k=2, variance="homogeneous"), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rrm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sar")
### * sar

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sar
### Title: Simultaneous autoregressive spatial covariance structure
### Aliases: sar

### ** Examples

## Not run: 
##D vsm(sar(location, W), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sar", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("scores_mmes")
### * scores_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: scores_mmes
### Title: Predict factor-analytic scores from a fitted mmes model
### Aliases: scores_mmes

### ** Examples

## Not run: 
##D mix <- mmes(BLUEs ~ trial,
##D             random = ~ vsm(fam(trial, 2), ism(genotype)),
##D             rcov = ~ units, data = dt)
##D scores_mmes(mix)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("scores_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sommer-package")
### * sommer-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sommer-package
### Title: *So*lving *M*ixed *M*odel *E*quations in *R*
### Aliases: sommer
### Keywords: R package

### ** Examples


####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_example, package="enhancer")

## No test: 



DT <- DT_example
DT=DT[with(DT, order(Env)), ]
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
## Compound simmetry (CS) + Diagonal (DIAG) model
ans3 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ vsm(dsm(Env),ism(units)),
             data=DT)
summary(ans3)



## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sommer-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spl2Dc")
### * spl2Dc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spl2Dc
### Title: Two-dimensional penalised tensor-product of marginal B-Spline
###   basis.
### Aliases: spl2Dc spl2Da spl2Db

### ** Examples

## ============================ ##
## example to use spl2Dc() 
## ============================ ## 
data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# A <- A.mat(GT)
## ============================ ##
## mimic 2 fields
## ============================ ## 
# aa <- DT; bb <- DT
# aa$FIELD <- "A";bb$FIELD <- "B"
# set.seed(1234)
# aa$Yield <- aa$Yield + rnorm(length(aa$Yield),0,4)
# DT2 <- rbind(aa,bb)
# head(DT2)
# mix <- mmes(Yield~1, control = mmes_control(formulation = "direct"),
#             random=~
#               vsm(dsm(FIELD),ism(Rowf)) +
#               vsm(dsm(FIELD),ism(Colf)) +
#                 spl2Dc(Row,Col,at.var=FIELD),
#             rcov=~vsm(dsm(FIELD),ism(units)),
#             data=DT2)
# 
# # extract spatial effects
# blup <- mix$uList$`spl2Dc(Row, Col, at.var = FIELD`
# head(blup) # 2 fields
# # recreate the incidence matrices
# xx=with(DT2, spl2Dc(Row,Col,at.var=FIELD))
# # get fitted values Zu for spatial effects and add them to the dataset
# field1 <- xx$Z$`A:all` %*% blup[,1]
# field2 <- xx$Z$`B:all` %*% blup[,2]
# DT2$spat <- field1+field2
# # plots the spatial effects
# lattice::levelplot(spat~Row*Col|FIELD, data=DT2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spl2Dc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spl2Dmats")
### * spl2Dmats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spl2Dmats
### Title: Get Tensor Product Spline Mixed Model Incidence Matrices
### Aliases: spl2Dmats

### ** Examples


data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix

M <- spl2Dmats(x.coord.name = "Col", y.coord.name = "Row", data=DT, nseg =c(14,21))
head(M$data)
# m1g <- mmes(Yield~1+TP.CR.2+TP.CR.3+TP.CR.4,
#             random=~Rowf+Colf+vsm(ism(M$fC))+vsm(ism(M$fR))+
#               vsm(ism(M$fC.R))+vsm(ism(M$C.fR))+vsm(ism(M$fC.fR))+
#               vsm(ism(id),Gu=A),
#             data=M$data, tolpar = 1e-6,
#             iters=30)
# 
# summary(m1g)$varcomp




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spl2Dmats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stackTraits")
### * stackTraits

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stackTraits
### Title: Stack several trait columns into the long format used by
###   multi-trait models
### Aliases: stackTraits
### Keywords: manip

### ** Examples

set.seed(17)
n <- 60
DT <- data.frame(id = factor(rep(seq_len(n), each = 3)))
u <- matrix(rnorm(n * 2), n) 
DT$trait_a <- u[DT$id, 1] + rnorm(nrow(DT), sd = 0.7)
DT$trait_b <- u[DT$id, 2] + rnorm(nrow(DT), sd = 1)
DTL <- stackTraits(DT, traits = c("trait_a", "trait_b"), keep = "id")
head(DTL)
table(DTL$trait)

# bivariate model with the residual paired through the record key
fit <- mmes(value ~ trait,
  random = ~ vsm(usm(trait), ism(id)),
  rcov = ~ vsm(usm(trait), ism(record)),
  data = DTL,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
covmatrix_mmes(fit, 1)$correlation   # genetic correlation



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stackTraits", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("strm")
### * strm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: strm
### Title: Covariance structure spanning several random terms
### Aliases: strm
### Keywords: models

### ** Examples

set.seed(3)
ids <- paste0("i", 1:30)
A <- diag(30) + diag(0.2, 30); dimnames(A) <- list(ids, ids)
Ai <- solve(A); attr(Ai, "inverse") <- TRUE
d <- data.frame(id = factor(sample(ids, 200, TRUE), levels = ids),
                dam = factor(sample(ids[1:15], 200, TRUE), levels = ids))
d$y <- rnorm(200)
fit <- mmes(y ~ 1,
  random = ~ strm(dir = vsm(ism(id)), mat = vsm(ism(dam)), Gu = Ai),
  data = d,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
covparams_mmes(fit, 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("strm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("toeplitzm")
### * toeplitzm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: toeplitzm
### Title: General positive-definite Toeplitz correlation structure
### Aliases: toeplitzm

### ** Examples

## Not run: 
##D vsm(toeplitzm(time), ism(id))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("toeplitzm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("unsm")
### * unsm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: unsm
### Title: unstructured indication matrix
### Aliases: unsm unsm2

### ** Examples

unsm(3)
unsm(3,2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("unsm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("usm")
### * usm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: usm
### Title: Unstructured positive-definite covariance structure
### Aliases: usm

### ** Examples

## Not run: 
##D vsm(usm(environment), ism(genotype))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("usm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("usr")
### * usr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: usr
### Title: Unstructured covariance structure for 'mmer' and 'vsr'
### Aliases: usr

### ** Examples

# Unstructured covariance pattern across three levels
U <- usr(factor(c("A","B","A","C")))
U$Z
U$thetaC



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("usr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vcc")
### * vcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vcc
### Title: Equality and scaling constraints between variance parameters
### Aliases: vcc vcParams
### Keywords: models

### ** Examples

data(DT_yatesoats, package = "enhancer")
DT <- DT_yatesoats
p <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  control = mmes_control(output = list(setup = TRUE)))
p$vcParams
# block and whole-plot variances equal
m1 <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  vcc = data.frame(parameter = c(1, 2), group = 1),
  control = mmes_control(diagnostics = list(verbose = FALSE)))
# block variance twice the whole-plot variance
m2 <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  vcc = data.frame(parameter = c("vsm(ism(B:MP)):sigma2", "vsm(ism(B)):sigma2"),
                            group = 1, scale = c(1, 2)),
  control = mmes_control(diagnostics = list(verbose = FALSE)))
unlist(m2$covPar)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vpredict")
### * vpredict

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vpredict
### Title: vpredict form of a LMM fitted with mmes
### Aliases: vpredict vpredict.mmes
### Keywords: models

### ** Examples


####=========================================####
####=========================================####
#### EXAMPLE 1
#### simple example with univariate models
####=========================================####
####=========================================####
# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix 
# A <- A.mat(GT)
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmes(Yield~1,
#               random=~vsm(ism(id),Gu=A), 
#               data=DT)
# summary(mix1)$varcomp
# #### run the vpredict function
# vpredict(mix1, h2 ~ V1 / ( V1 + V2 ) )





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vpredict", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vsm")
### * vsm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vsm
### Title: Variance-structure model with arbitrary Kronecker covariance
###   factors
### Aliases: vsm

### ** Examples


####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)


str(with(DT, vsm(dsm(Env),ism(Name))))
str(with(DT, vsm(csm(Env),ism(Name))))
str(with(DT, vsm(toeplitzm(Env),ism(Name))))

## Gu in 3-column (row, column, value) format with level names in rowNames
Ainv <- data.frame(Row=c(1,2,2,3), Column=c(1,1,2,3),
                   Ainverse=c(2,-1,2,1))
attr(Ainv, "rowNames") <- c("a","b","c")
id <- factor(c("a","b","c","a"))
str(vsm(ism(id), Gu=Ainv)$Gu)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vsm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("wald_mmes")
### * wald_mmes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: wald_mmes
### Title: Wald tests for fixed effects with optional small-sample
###   denominator df
### Aliases: wald_mmes print.wald.mmes
### Keywords: models

### ** Examples

data(DT_yatesoats, package = "enhancer")
m <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT_yatesoats,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
wald_mmes(m)
wald_mmes(m, ssType = "conditional", denDF = "kr")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("wald_mmes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
