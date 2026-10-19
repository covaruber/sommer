## -----------------------------------------------------------------------------
library(sommer)
data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

ans1 <- mmes(Yield~1,
  random= ~ Name + Env + Env:Name + Env:Block,
  rcov= ~ units,
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
summary(ans1)$varcomp
(n.env <- length(levels(DT$Env)))
vpredict(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

## -----------------------------------------------------------------------------
data(DT_cpdata, package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
DT$idd <-DT$id; DT$ide <-DT$id
### look at the data
A <- A.mat(GT) # additive relationship matrix
D <- D.mat(GT) # dominance relationship matrix
E <- E.mat(GT) # epistatic relationship matrix

Ai <- solve(A + diag(1e-5, nrow(A), nrow(A)))
Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE

Di <- solve(D+ diag(1e-5, nrow(A), nrow(A)))
Di <- as(as(as( Di,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Di, "inverse")=TRUE


# ans.ADE <- mmes(Yield~1, 
#                  random=~vsm(ism(id),Gu=Ai) + vsm(ism(idd),Gu=Di), 
#                  rcov=~units, data=DT,
#                  control=mmes_control(convergence=list(maxIters=10),
#                    diagnostics=list(verbose=FALSE)))
# (summary(ans.ADE)$varcomp)
# vpredict(ans.ADE, h2 ~ (V1) / ( V1+V3) ) # narrow sense
# vpredict(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense

## ----fig.show='hold'----------------------------------------------------------
data(DT_cornhybrids, package="enhancer")
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids
### fit the model
modFD <- mmes(Yield~1,
  random=~ vsm(atm(Location,c("3","4")),ism(GCA2)),
  rcov= ~ vsm(dsm(Location),ism(units)),
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
summary(modFD)

## -----------------------------------------------------------------------------
data(DT_cpdata, package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
### look at the data
A <- A.mat(GT) # additive relationship matrix
Ai <- solve(A+ diag(1e-4, nrow(A),nrow(A))) 
Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE

ans <- mmes(color~1,
  random=~vsm(ism(id),Gu=Ai),
  rcov=~units,
  data=DT,
  control = mmes_control(convergence = list(maxIters = 10),
  diagnostics = list(verbose = FALSE)))
summary(ans)$varcomp
vpredict(ans, h2 ~ (V1) / ( V1+V2) )


## -----------------------------------------------------------------------------
data(DT_cornhybrids, package="enhancer")
DT <- DT_cornhybrids
DTi <- DTi_cornhybrids
GT <- GT_cornhybrids

modFD <- mmes(Yield~Location,
  random=~GCA1+GCA2+SCA,
  rcov=~units,
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
(suma <- summary(modFD)$varcomp)
Vgca <- sum(suma[1:2,"estimate"])
Vsca <- suma[3,"estimate"]
Ve <- suma[4,"estimate"]
Va = 4*Vgca
Vd = 4*Vsca
Vg <- Va + Vd
(H2 <- Vg / (Vg + (Ve)) )
(h2 <- Va / (Vg + (Ve)) )

## -----------------------------------------------------------------------------
data("DT_halfdiallel", package="enhancer")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)
#### model using overlay
modh <- mmes(sugar~1,
  random=~vsm(ism(overlay(femalef,malef)) )
             + genof,
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
summary(modh)$varcomp

## -----------------------------------------------------------------------------

data(DT_wheat, package="enhancer")
DT <- DT_wheat
GT <- apply(GT_wheat,2,as.numeric)
rownames(GT) <- rownames(GT_wheat)

colnames(DT) <- paste0("X",1:ncol(DT))
DT <- as.data.frame(DT);DT$id <- as.factor(rownames(DT))
# select environment 1
K <- A.mat(GT) # additive relationship matrix
colnames(K) <- rownames(K) <- rownames(DT)
Ki <- solve(K+ diag(1e-4, nrow(K),nrow(K))) 
Ki[lower.tri(Ki)] <- t(Ki)[lower.tri(Ki)] # fill the lower triangular
Ki <- as(as(as( Ki,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ki, "inverse")=TRUE


# GBLUP pedigree-based approach
set.seed(12345)
y.trn <- DT
vv <- sample(rownames(DT),round(nrow(DT)/5))
y.trn[vv,"X1"] <- NA
head(y.trn)
## GBLUP with mmes
ans <- mmes(X1~1,
  random=~vsm(ism(id),Gu=Ki),
  rcov=~units,
  data=y.trn,
  control = mmes_control(diagnostics = list(verbose = FALSE))) # kinship based
cor(ans$u[vv,] ,DT[vv,"X1"], use="complete")

## rrBLUP with mmer
ans2 <- mmer(X1~1,
             random=~vsr(list(GT)), 
             rcov=~units, getPEV = TRUE,
             data=y.trn, verbose = FALSE) # kinship based

u <- GT %*% ans2$U$`u:GT`$X1 # BLUPs for individuals
rownames(u) <- rownames(GT)
cor(u[vv,],DT[vv,"X1"]) # same correlation
# the same can be applied in multi-response models in GBLUP or rrBLUP



## ----eval=FALSE---------------------------------------------------------------
# G <- A.mat(GT)
# Ki <- solve(0.99 * G + 0.01 * diag(nrow(G))); attr(Ki, "inverse") <- TRUE
# ans <- mmes(X1 ~ 1,
#   random = ~ vsm(ism(id), Gu = Ki),
#   data = y.trn,
#   control = mmes_control(diagnostics = list(verbose = FALSE)))
# markers <- meffects_mmes(ans,
#   1,
#   GT,
#   se = TRUE,
#   control = mmes_control(update = list(recovery = 0.01)))
# head(markers[order(markers$p.value), ])

## -----------------------------------------------------------------------------
data(DT_ige, package="enhancer")
DT <- DT_ige
Af <- A_ige
An <- A_ige

# Direct genetic effects model
modDGE <- mmes(trait ~ block,
  random = ~ focal,
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
summary(modDGE)$varcomp


## -----------------------------------------------------------------------------
data(DT_ige, package="enhancer")
DT <- DT_ige
A <- A_ige

## Indirect genetic effects model
modIGE <- mmes(trait ~ block,
  random = ~ focal + neighbour,
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(dateWarning = FALSE, verbose = FALSE)))
summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

modIGE <- mmes(trait ~ block,
  random = ~ strm(dge = vsm(ism(focal)), ige = vsm(ism(neighbour))),
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(dateWarning = FALSE, verbose = FALSE)))
summary(modIGE)$varcomp


## -----------------------------------------------------------------------------

Ai <- solve(A_ige + diag(1e-5, nrow(A_ige),nrow(A_ige) ))
attr(Ai, "inverse") <- TRUE
modIGE <- mmes(trait ~ block,
  random = ~ strm(dge = vsm(ism(focal)), ige = vsm(ism(neighbour)), Gu = Ai),
  rcov = ~ units,
  data = DT,
  control = mmes_control(diagnostics = list(dateWarning = FALSE, verbose = FALSE)))
summary(modIGE)$varcomp


## -----------------------------------------------------------------------------
data(DT_technow, package="enhancer")
DT <- DT_technow

Md <- apply(Md_technow,2,as.numeric)
rownames(Md) <- rownames(Md_technow)
Mf <- apply(Mf_technow,2,as.numeric)
rownames(Mf) <- rownames(Mf_technow)

Md <- (Md*2) - 1
Mf <- (Mf*2) - 1
Ad <- A.mat(Md)
Af <- A.mat(Mf)
Adi <- solve(Ad + diag(1e-4,ncol(Ad),ncol(Ad)))
Adi[lower.tri(Adi)] <- t(Adi)[lower.tri(Adi)] # fill the lower triangular
Adi <- as(as(as( Adi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Adi, 'inverse')=TRUE
Afi <- solve(Af + diag(1e-4,ncol(Af),ncol(Af)))
Afi[lower.tri(Afi)] <- t(Afi)[lower.tri(Afi)] # fill the lower triangular
Afi <- as(as(as( Afi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Afi, 'inverse')=TRUE
# RUN THE PREDICTION MODEL
y.trn <- DT
vv1 <- which(!is.na(DT$GY))
vv2 <- sample(vv1, 100)
y.trn[vv2,"GY"] <- NA
anss2 <- mmes(GY~1,
  random=~vsm(ism(dent),Gu=Adi) + vsm(ism(flint),Gu=Afi),
  rcov=~units,
  data=y.trn,
  control = mmes_control(formulation = "henderson",
  convergence = list(maxIters = 15),
  diagnostics = list(verbose = FALSE)))
summary(anss2)$varcomp

# zu1 <- model.matrix(~dent-1,y.trn) %*% anss2$uList$`vsm(ism(dent), Gu = Adi)`
# zu2 <- model.matrix(~flint-1,y.trn) %*% anss2$uList$`vsm(ism(flint), Gu = Afi)`
# u <- zu1+zu2+as.vector(anss2$b)
# cor(u[vv2,], DT$GY[vv2])

## -----------------------------------------------------------------------------
# data(DT_cpdata, package="enhancer")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# traits <- c("color","Yield")
# DT[,traits] <- apply(DT[,traits],2,scale)
# DTL <- reshape(DT[,c("id", traits)],
#                idvar = c("id"),
#                varying = traits,
#                v.names = "value", direction = "long",
#                timevar = "trait", times = traits )
# DTL <- DTL[with(DTL, order(trait)), ]
# head(DTL)
# 
# A <- A.mat(GT) # additive relationship matrix
# # if using mmes=TRUE you need to provide the inverse
# Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
# Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, 'inverse')=TRUE
# #### be patient this model is heavier
# ansm <- mmes( value ~ trait,
#                random=~ vsm(usm(trait), ism(id), Gu=Ai), # Ai if henderson
#                rcov=~ vsm(dsm(trait), ism(units)),
#                data=DTL,
#                control=mmes_control(formulation="henderson"))
# cov2cor(ansm$theta[[1]])


## -----------------------------------------------------------------------------
# DTL <- stackTraits(DT, traits = c("color", "Yield"), keep = "id")
# ansm <- mmes(value ~ trait,
#              random = ~ vsm(usm(trait), ism(id), Gu = Ai),
#              rcov = ~ vsm(usm(trait), ism(record)), data = DTL)
# covmatrix_mmes(ansm, 1)$correlation  # genetic correlation
# covmatrix_mmes(ansm, 2)$correlation  # residual correlation

## -----------------------------------------------------------------------------
data(DT_legendre)
DT <- DT_legendre
head(DT)
DT$SUBJECT <- paste("s",DT$SUBJECT,sep="_")
DT1 <- DT2 <- DT
DT1$TRAIT <- "T1"
DT2$TRAIT <- "T2"
DT2$Y <- sample(DT2$Y)
DTC <- rbind(DT1,DT2)

## -----------------------------------------------------------------------------
# 
# library(orthopolynom)
# 
# Z <- with(DTC, dsm(leg(X,1)) )$Z
# for(i in 1:ncol(Z)){DTC[,colnames(Z)[i]] <- Z[,i]}
# 
# X <- with(DTC, dsm(TRAIT) )$Z
# for(i in 1:ncol(X)){DTC[,colnames(X)[i]] <- X[,i]}
# 
# A <- diag(length(unique(DTC$SUBJECT)))
# rownames(A) <- colnames(A) <- unique(DTC$SUBJECT)
# Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
# Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
# Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(Ai, 'inverse')=TRUE
# ##
# M <- model.matrix(~ T1:leg0 + T1:leg1 + T2:leg0 + T2:leg1 - 1 , data=DTC)
# mRR2b<-mmes(Y ~ Xf,
#             random=~ vsm( usm( M ) ,  ism(SUBJECT) , Gu = Ai),
#             rcov = ~ vsm( dsm(TRAIT), ism(units) ),
#             data=DTC,
#             control=mmes_control(convergence=list(maxIters=10),
#               diagnostics=list(verbose=FALSE)))
# summary(mRR2b)$varcomp

## -----------------------------------------------------------------------------
library(sommer)
data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
M <- GT_cpdata

################
# MARKER MODEL
################
mix.marker <- mmer(Yield~1,
                   random=~Rowf+vsr(list(M)),
                   rcov=~units,data=DT, 
                   verbose = FALSE)


me.marker <- mix.marker$U$`u:M`$Yield

################
# PARTITIONED GBLUP MODEL
################

MMT <-tcrossprod(M) ## MM' = additive relationship matrix 
MMTinv<-solve(MMT + diag(1e-4, nrow(MMT), nrow(MMT))) ## inverse
MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
MMTinv <- as(as(as( MMTinv,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(MMTinv, 'inverse')=TRUE

mix.part <- mmes(Yield~1,
  random=~Rowf+vsm(ism(id), Gu=MMTinv),
  rcov=~units,
  data=DT,
  control = mmes_control(convergence = list(maxIters = 20),
  diagnostics = list(verbose = FALSE)))

#convert BLUPs to marker effects me=M'(M'M)- u
me.part<-MTMMTinv%*%matrix(mix.part$uList$`vsm(ism(id), Gu = MMTinv`,ncol=1)

# compare marker effects between both models
plot(me.marker,me.part)



## -----------------------------------------------------------------------------

data("DT_wheat", package="enhancer")
rownames(GT_wheat) <- rownames(DT_wheat)
GT <- apply(GT_wheat,2,as.numeric)
rownames(GT) <- rownames(GT_wheat)
A <- A.mat(GT)
Ai <- solve(A + diag(1e-5, nrow(A)))
Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
Ai <- as(as(as(Ai, "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse") <- TRUE
isSymmetric(Ai)
# One complete observation per relationship level gives a balanced example.
DTn <- data.frame(
  id=rownames(A),
  y=as.numeric(DT_wheat[,1])
)

model_regular <- mmes(y~1,
  random=~vsm(ism(id), Gu=Ai),
  rcov=~units,
  data=DTn,
  control = mmes_control(diagnostics = list(verbose = FALSE)))

# Henderson MME formulation with eigenbasis random coefficients.
model_rotation_h <- mmes(y~1,
  random=~vsm(ism(id), Gu=Ai, rotation=TRUE),
  rcov=~units,
  data=DTn,
  control = mmes_control(formulation = "henderson",
  diagnostics = list(verbose = FALSE)))

# Lee--van der Werf direct observation-covariance formulation.
model_rotation_d <- mmes(y~1,
  random=~vsm(ism(id), Gu=Ai, rotation=TRUE),
  rcov=~units,
  data=DTn,
  control = mmes_control(formulation = "direct",
  diagnostics = list(verbose = FALSE)))

model_regular$covParNative
model_rotation_h$covParNative
model_rotation_d$covParNative

plot(model_rotation_h$bu[,1], model_regular$bu[,1])
plot(model_rotation_d$bu[,1], model_regular$bu[,1])


## -----------------------------------------------------------------------------

data("DT_wheat", package="enhancer")
rownames(GT_wheat) <- rownames(DT_wheat)
GT <- apply(GT_wheat,2,as.numeric)
rownames(GT) <- rownames(GT_wheat)
A <- A.mat(GT)
Ai <- solve(A + diag(1e-5, nrow(A)))
Ai[lower.tri(Ai)] <- t(Ai)[lower.tri(Ai)] # fill the lower triangular
Ai <- as(as(as(Ai, "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse") <- TRUE
isSymmetric(Ai)
# One complete observation per relationship level gives a balanced example.
data(DT_wheat)
DT <- DT_wheat
GT <- apply(GT_wheat,2,as.numeric)
rownames(GT) <- rownames(GT_wheat)
DT <- data.frame(pheno=as.vector(DT),
                 env=as.factor(paste0("e", sort(rep(1:4,nrow(DT))))),
                 id=rep(rownames(DT),4))

# Henderson MME formulation with eigenbasis random coefficients.
model_rotation_h <- mmes(pheno~1,
  random=~vsm(usm(env),ism(id), Gu=Ai, rotation=TRUE),
  rcov=~vsm(dsm(env),ism(units)),
  data=DT,
  control = mmes_control(formulation = "henderson",
  diagnostics = list(verbose = FALSE)))

# genetic correlation
cov2cor(model_rotation_h$theta[[1]])


## -----------------------------------------------------------------------------

data(DT_expdesigns, package="enhancer")
DT <- DT_expdesigns$car1
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
# lattice::levelplot(yield~male * malef:femalef|set, data=DT, main="NC design I")
##############################
## Expected Mean Square method
##############################
mix1 <- lm(yield~ setf + setf:repf + femalef:malef:setf + malef:setf, data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=2
nfem=2
Vfm <- (ms2-mse)/nrep
Vm <- (ms1-ms2)/(nrep*nfem)

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg
##############################
## REML method
##############################
mix2 <- mmes(yield~ setf + setf:repf,
  random=~femalef:malef:setf + malef:setf,
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
vc <- summary(mix2)$varcomp; vc
Vfm <- vc[1,"estimate"]
Vm <- vc[2,"estimate"]

## Calculate Va and Vd
Va=4*Vm # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm-Vm) # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg


## -----------------------------------------------------------------------------
DT <- DT_expdesigns$car2
DT <- aggregate(yield~set+male+female+rep, data=DT, FUN = mean)
DT$setf <- as.factor(DT$set)
DT$repf <- as.factor(DT$rep)
DT$malef <- as.factor(DT$male)
DT$femalef <- as.factor(DT$female)
#levelplot(yield~male*female|set, data=DT, main="NC desing II")
head(DT)

N=with(DT,table(female, male, set))
nmale=length(which(N[1,,1] > 0))
nfemale=length(which(N[,1,1] > 0))
nrep=table(N[,,1])
nrep=as.numeric(names(nrep[which(names(nrep) !=0)]))

##############################
## Expected Mean Square method
##############################

mix1 <- lm(yield~ setf + setf:repf + 
             femalef:malef:setf + malef:setf + femalef:setf, data=DT)
MS <- anova(mix1); MS
ms1 <- MS["setf:malef","Mean Sq"]
ms2 <- MS["setf:femalef","Mean Sq"]
ms3 <- MS["setf:femalef:malef","Mean Sq"]
mse <- MS["Residuals","Mean Sq"]
nrep=length(unique(DT$rep))
nfem=length(unique(DT$female))
nmal=length(unique(DT$male))
Vfm <- (ms3-mse)/nrep; 
Vf <- (ms2-ms3)/(nrep*nmale); 
Vm <- (ms1-ms3)/(nrep*nfemale); 

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg

##############################
## REML method
##############################

mix2 <- mmes(yield~ setf + setf:repf,
  random=~femalef:malef:setf + malef:setf + femalef:setf,
  data=DT,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
vc <- summary(mix2)$varcomp; vc
Vfm <- vc[1,"estimate"]
Vm <- vc[2,"estimate"]
Vf <- vc[3,"estimate"]

Va=4*Vm; # assuming no inbreeding (4/(1+F))
Va=4*Vf; # assuming no inbreeding (4/(1+F))
Vd=4*(Vfm); # assuming no inbreeding(4/(1+F)^2)
Vg=c(Va,Vd); names(Vg) <- c("Va","Vd"); Vg


## -----------------------------------------------------------------------------
data(DT_cpdata, package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata[,1:200]
MP <- MP_cpdata
#### create the variance-covariance matrix
A <- A.mat(GT) # additive relationship matrix
n <- nrow(DT) # to be used for degrees of freedom
k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)

## -----------------------------------------------------------------------------
###########################
#### Regular GWAS/EMMAX approach
###########################
# mix2 <- GWAS(color~1,
#              random=~vsm(ism(id), Gu=A) + Rowf + Colf,
#              rcov=~units, M=GT, gTerm = "u:id",
#              verbose = FALSE, 
#              data=DT)

## -----------------------------------------------------------------------------
# ###########################
# #### GWAS by RRBLUP approach
# ###########################
# Z <- GT[as.character(DT$id),]
# mixRRBLUP <- mmer(Yield~1,
#               random=~vsr(list(Z)) + Rowf + Colf,
#               rcov=~units, nIters=10,
#               verbose = FALSE,
#               data=DT)
# 
# a <- mixRRBLUP$U$`u:Z`$Yield
# se.a <- sqrt( diag(kronecker(diag(ncol(Z)),mixRRBLUP$sigma$`u:Z`) - mixRRBLUP$PevU$`u:Z`$Yield ) ) # SE of marker effects
# t.stat <- a/se.a # t-statistic
# pvalRRBLUP <- dt(t.stat,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
# ###########################
# #### GWAS by GBLUP approach
# ###########################
# M<- GT
# MMT <-tcrossprod(M) ## MM' = additive relationship matrix
# MMTinv<-solve(MMT + diag(1e-4, ncol(MMT), ncol(MMT))) ## inverse of MM'
# MTMMTinv<-t(M)%*%MMTinv # M' %*% (M'M)-
# MMTinv <- as(as(as( MMTinv,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
# attr(MMTinv, 'inverse')=TRUE
# 
# mixGBLUP <- mmes(Yield~1,
#              random=~vsm(ism(id), Gu=MMTinv) + Rowf + Colf,
#              rcov=~units, data=DT,
#              control=mmes_control(convergence=list(maxIters=25),
#                diagnostics=list(verbose=TRUE), output=list(inverse=2)))
# a.from.g <-MTMMTinv%*%matrix(mixGBLUP$uList$`vsm(ism(id), Gu = MMTinv`,ncol=1)
# start=mixGBLUP$partitions[[1]][1]
# end=mixGBLUP$partitions[[1]][2]
# var.g <- kronecker(MMT,mixGBLUP$theta[[1]]) - mixGBLUP$Ci[start:end,start:end]
# var.a.from.g <- t(M)%*%MMTinv%*% (var.g) %*% t(MMTinv)%*%M
# se.a.from.g <- sqrt(diag(var.a.from.g))
# t.stat.from.g <- a.from.g/se.a.from.g # t-statistic
# pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) # -log10(pval)

## -----------------------------------------------------------------------------
###########################
#### Compare results
###########################
# plot(mix2$scores[,1], main="GWAS")
# plot(-log(pvalRRBLUP), main="GWAS by RRBLUP/SNP-BLUP") 
# plot(-log(pvalGBLUP), main="GWAS by GBLUP")


