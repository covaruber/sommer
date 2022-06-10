# plotMonitor <- function(x,...){
#   mmin <- min(x)
#   mmax <- max(x)
#   plot(x[1,], type="l", ylim=c(mmin,mmax))
#   if(nrow(x) > 1){
#     for(i in 2:nrow(x)){
#       par(new=TRUE)
#       plot(x[i,], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="l",...)
#     }
#   }
#   legend("topright",legend = rownames(x), col=1:(nrow(x)), bty="n", lty=1)
# }
# 
# library(sommer)
# data("DT_cpdata")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# DT$Yield <- imputev(DT$Yield)
# #### create the variance-covariance matrix
# K <- A.mat(GT) # additive relationship matrix
# K <- K + diag(1e-3,nrow(K))
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#              random=~vsr(id,Gu=K)
#              + Rowf,
#              rcov=~units,
#              data=DT)
# summary(mix1)$varcomp
# 
# Ki <- as(solve(K), Class="sparseMatrix")
# mix2 <- mmec(fixed=Yield~1,
#              random=~vsc(ids(id),Gu=Ki) + Rowf,
#              rcov=~units, return.param = F,
#              data=DT)
# 
# mix2$monitor[,15]
# plotMonitor(mix2$monitor)
# pp <- predict.mmec(object=mix2,classify = c("id","(Intercept)"))#, ignore="Rowf")
# head(pp$pvals)
# 
# Z <- list(model.matrix(~id-1, data=DT),model.matrix(~Rowf-1,data=DT)
# )
# Zind <- 1:2
# 
# A <- list(
#   K, 
#   diag(13)
# ) #
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# thetaI <- list(
#   matrix(1000,1,1),
#   matrix(1000,1,1),
#   matrix(1000,1,1)
# );thetaI
# 
# thetaC <- list(
#   matrix(1,1,1),
#   matrix(1,1,1),
#   matrix(1,1,1)
# );thetaC
# 
# X <- model.matrix(~1, data=DT)
# 
# y <- as.matrix(DT$Yield)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# weightInf=rep(1,40) # weights for the information matrix
# weightEmInf=c(seq(.9,.1,-.2),rep(0,36));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# res3$monitor
# summary(mix1)$varcomp
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# # CS + DIAGONAL MODEL
# 
# library(sommer)
# data("DT_example")
# DT <- DT_example
# K <- A_example
# #### look at the data and fit the model
# head(DT)
# DT$Yield <- scale(DT$Yield)
# mix1 <- mmer(Yield~Env,
#              random= ~Name + vsr(dsr(Env),Name),
#              rcov= ~ vsr(dsr(Env),units),
#              data=DT)
# summary(mix1)$varcomp
# 
# mix2 <- mmec(Yield~Env,
#              random= ~ Name + vsc(dsr(Env),ids(Name)),
#              rcov= ~ vsc(dsr(Env),ids(units)),
#              return.param = F,
#              data=DT)
# plotMonitor(mix2$monitor)
# summary(mix1)$varcomp
# mix2$sigma
# 
# zz <- with(DT, vsr(dsr(Env),Name))
# 
# Z <- c(list(model.matrix(~Name-1, data=DT)),zz$Z)
# 
# Zind <- c(1,2,2,2)
# 
# A <- list(diag(41), diag(41))#rep(list(diag(41)),4)
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# thetaI <- list(
#   matrix(10,1,1),
#   diag(10,3,3),
#   diag(10,3,3)
# );thetaI
# 
# thetaC <- list(
#   matrix(1,1,1),
#   diag(1,3,3),
#   diag(1,3,3)
# );thetaC
# 
# X <- model.matrix(~Env, data=DT)
# 
# y <- as.matrix(DT$Yield)
# 
# DTx <- DT; DTx$units <- as.factor(1:nrow(DTx))
# ss <- with(DTx, vsr(dsr(Env),units) )
# 
# S <- ss$Z #list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# 
# ## apply the function
# weightInf=rep(1,40); # weights for the information matrix
# weightEmInf=c(seq(.9,.1,-.1),rep(0,36)); # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# res3$monitor
# summary(mix1)$varcomp
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# # UNSTRUCTURED MODEL
# library(sommer)
# data("DT_example")
# DT <- DT_example
# K <- A_example
# #### look at the data and fit the model
# head(DT)
# DT$Yield <- scale(DT$Yield)
# mix1 <- mmer(Yield~Env,
#              random= ~ vsr(usr(Env),Name),
#              rcov= ~ vsr(dsr(Env),units),
#              data=DT)
# summary(mix1)$varcomp[,1]
# mix1$Beta
# 
# # > summary(mix1)$varcomp[,1]
# # [1] 0.76660279 0.29901471 0.22166682 0.31243923 0.01923234 0.42072739 0.24321083 0.27761802
# # [9] 0.12513727
# 
# mix2 <- mmec(fixed=Yield~Env-1,
#               random= ~ vsc(usr(Env),ids(Name)),
#               rcov= ~ vsc(dsr(Env),ids(units)),
#               return.param = F,
#              # emweight = rep(1,30),
#              nIters = 30,
#               data=DT)
# mix2$theta
# plotMonitor(mix2$monitor, cex=1)
# mix2$sigma
# mix2$monitor[,12]
# plot(summary(mix1)$varcomp[,1],mix2$sigma)
# 
# zz <- with(DT, vsr(dsr(Env),Name))
# 
# Z <- zz$Z
# 
# Zind <- rep(1,length(Z))
# 
# A <- list( diag(41) )#rep(list(diag(41)),4)
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# tt = ((unsm(3)/5) + diag(.8,3,3)) * 10 ;# tt[lower.tri(tt)]=0;tt
# thetaI <- list(
#   tt,
#   diag(10,3,3)
# );thetaI
# 
# ttc= unsm(3);ttc[lower.tri(ttc)]=0;ttc
# thetaC <- list(
#   ttc,
#   diag(1,3,3)
# );thetaC
# 
# X <- model.matrix(~Env, data=DT)
# 
# y <- as.matrix(DT$Yield)
# 
# DTx <- DT; DTx$units <- as.factor(1:nrow(DTx))
# ss <- with(DTx, vsr(dsr(Env),units) )
# 
# S <- ss$Z #list(diag(length(y)))
# 
# Sind <- rep(1,3)
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# 
# ## apply the function
# weightInf=rep(1,50);weightInf # weights for the information matrix
# weightEmInf=rep(1,50)
# weightEmInf=c(seq(1,.1,-.1),rep(0,50));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# summary(mix1)$varcomp
# res3$monitor
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# 
# 
# library(sommer)
# data(DT_technow)
# DT <- DT_technow
# Md <- Md_technow
# Mf <- Mf_technow
# Md <- (Md*2) - 1
# Mf <- (Mf*2) - 1
# Ad <- A.mat(Md)
# Af <- A.mat(Mf)
# DT$GY <- scale(DT$GY)
# ####=========================================####
# ####=========================================####
# mix1 <- mmer(GY~1,
#              random=~vsr(dent,Gu=Ad) + vsr(flint,Gu=Af),
#              rcov=~units,
#              data=DT)
# summary(mix1)$varcomp[,1]
# 
# Adi <- as(solve(Ad + diag(ncol(Ad))*1e-6 ),Class="sparseMatrix")
# Afi <- as(solve(Af + diag(ncol(Af))*1e-6 ),Class="sparseMatrix")
# 
# mix2 <- mmec(GY~1,
#              random=~vsc(ids(dent),Gu=Adi) + vsc(ids(flint),Gu=Afi),
#              rcov=~units, return.param = F,
#              data=DT)
# mix2$monitor
# mix2$sigma
# 
# z1=model.matrix(~dent-1, data=DT); colnames(z1) <- gsub("dent","",colnames(z1))
# z2=model.matrix(~flint-1, data=DT); colnames(z2) <- gsub("flint","",colnames(z2))
# 
# Z <- list(z1,z2)
# 
# Zind <- 1:2
# 
# A <- list(
#   Ad[colnames(z1),colnames(z1)], 
#   Af[colnames(z2),colnames(z2)]
# ) 
# A <- lapply(A, function(x){x + diag(1e-3,ncol(x),ncol(x))})
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# thetaI <- list(
#   matrix(10,1,1),
#   matrix(10,1,1),
#   matrix(10,1,1)
# );thetaI
# 
# thetaC <- list(
#   matrix(1,1,1),
#   matrix(1,1,1),
#   matrix(1,1,1)
# );thetaC
# 
# X <- model.matrix(~1, data=DT)
# 
# y <- as.matrix(DT$GY)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# ## apply the function
# weightInf=rep(1,40);weightInf # weights for the information matrix
# weightEmInf=c(seq(.9,.1,-.2),rep(0,36));weightEmInf # weights for the EM information matrix
# # v2=rep(1,40)
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=FALSE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# summary(mix1)$varcomp
# res3$monitor
# 
# dim(res3$Ci)
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# 
# library(sommer)
# data(DT_ige)
# DT <- DT_ige
# Af <- A_ige
# An <- A_ige
# DT$trait <- scale(imputev(DT$trait))
# ### Direct genetic effects model
# mix1 <- mmer(trait ~ block,
#              random = ~vsr(focal,Gu=Af) + vsr(neighbour,Gu=An),
#              rcov = ~ units,
#              data = DT)
# summary(mix1)$varcomp
# 
# Afi <- as(solve(Af + diag(ncol(Af))*1e-6 ),Class="sparseMatrix")
# Ani <- as(solve(An + diag(ncol(An))*1e-6 ),Class="sparseMatrix")
# 
# mix2 <- mmec(trait ~ block,
#              random = ~vsc(ids(focal),Gu=Afi) + vsc(ids(neighbour),Gu=Ani),
#              rcov = ~ units, return.param = F,
#              data = DT)
# mix2$monitor
# plotMonitor(mix2$monitor)
# 
# g=(unsm(2)*-.15) + (diag(2)*.4)
# g
# mix3 <- mmec(trait ~ block,
#              random = ~vsc(ids(focal),ids(neighbour),Gu=Afi, meN=2, meThetaC = unsm(2)),# + vsc(ids(neighbour),Gu=Ani),
#              rcov = ~ units, return.param = F, emweight = rep(1,20),
#              tolParConv = 1e-6, nIters=40,
#              data = DT)
# mix3$monitor
# mix3$sigma
# mix3$rTermsNames
# plotMonitor(mix3$monitor)
# 
# 
# z1=model.matrix(~focal-1,data=DT); colnames(z1) <- gsub("focal","",colnames(z1))
# z2=model.matrix(~neighbour-1,data=DT);  colnames(z2) <- gsub("neighbour","",colnames(z2))
# 
# Z <- list(z1,z2)
# 
# Zind <- 1:2
# 
# A <- list( 
#   Af[colnames(z1),colnames(z1)],
#   An[colnames(z2),colnames(z2)]
# )
# 
# A <- lapply(A, function(x){x + diag(1e-3,ncol(x),ncol(x))})
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# tt = diag(1)*10000
# thetaI <- list(
#   tt,
#   tt,
#   diag(10000,1,1)
# );thetaI
# 
# ttc=diag(1)
# thetaC <- list(
#   ttc,
#   ttc,
#   diag(1,1,1)
# );thetaC
# 
# X <- model.matrix(~block, data=DT)
# 
# y <- as.matrix(DT$trait)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# weightInf=rep(1,40);weightInf # weights for the information matrix
# weightEmInf=c(seq(.9,.1,-.1),rep(0,36));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# summary(mix1)$varcomp
# res3$monitor
# 
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# # UNSTRUCTURED MODEL
# 
# library(sommer)
# data(DT_ige)
# DT <- DT_ige
# Af <- A_ige
# An <- A_ige
# DT$trait <- imputev(DT$trait)
# ### Direct genetic effects model
# xx <- c(rep(1,4),rep(0,40))
# mix1 <- mmer(trait ~ block,
#              random = ~ gvs(focal, neighbour, Gu=list(Af,Af)), 
#              rcov = ~ units, iters=10,tolpar = 1e-4,
#              data = DT, emupdate = xx)
# summary(mix1)$varcomp
# 
# z1=model.matrix(~focal-1,data=DT); colnames(z1) <- gsub("focal","",colnames(z1))
# z2=model.matrix(~neighbour-1,data=DT);  colnames(z2) <- gsub("neighbour","",colnames(z2))
# 
# Z <- list(z1,z2)
# 
# Zind <- rep(1,2)
# A <- list( 
#   Af
# )
# 
# A <- lapply(A, function(x){x + diag(1e-3,ncol(x),ncol(x))})
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# tt = ((-unsm(2)/5) + diag(.8,2,2)) * 1000 ; #tt[lower.tri(tt)]=0;tt
# # tt = diag(2)*10000
# thetaI <- list(
#   tt,
#   diag(10000,1,1)
# );thetaI
# 
# ttc= unsm(2);ttc[lower.tri(ttc)]=0;ttc
# # ttc=diag(2)
# thetaC <- list(
#   ttc,
#   diag(1,1,1)
# );thetaC
# 
# X <- model.matrix(~block, data=DT)
# 
# y <- as.matrix(DT$trait)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# 
# weightInf=c(rep(1,100));weightInf
# weightEmInf=c(rep(1,4),seq(1,.1,-.2),rep(0.05,3),rep(0.05,30));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# summary(mix1)$varcomp
# res3$monitor
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# # UNSTRUCTURED MODEL
# 
# library(orthopolynom)
# data(DT_legendre)
# DT <- DT_legendre
# DT$Y <- scale(imputev(DT$Y))
# mix1<-mmer(Y~ 1 + Xf
#            , random=~ vsr(usr(leg(X,1)),SUBJECT)
#            , rcov=~vsr(units), tolpar = 1e-6
#            , data=DT)
# summary(mix1)$varcomp
# 
# mix2<-mmec(Y~ 1 + Xf
#            , random=~ vsc(usr(leg(X,1)),ids(SUBJECT))
#            , rcov=~units,return.param = F
#            , data=DT)
# mix2$monitor
# mix2$sigma
# 
# xx=with(DT,vsr(usr(leg(X,1)),SUBJECT))
# 
# Z <- xx$Z[c(1,3)]
# 
# Zind <- rep(1,2)
# 
# A <- list(
#   diag(100)
# ) #
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# tt = ((unsm(2)/5) + diag(.8,2,2))  ; #tt[lower.tri(tt)]=0;tt
# thetaI <- list(
#   tt,
#   diag(1,1,1)
# );thetaI
# 
# ttc= unsm(2);ttc[lower.tri(ttc)]=0;ttc
# # ttc=diag(2)
# thetaC <- list(
#   ttc,
#   diag(1,1,1)
# );thetaC
# 
# X <- model.matrix(~1 + Xf, data=DT)
# 
# y <- as.matrix(DT$Y)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# weightInf=rep(1,40) # weights for the information matrix
# weightEmInf=rep(1,40)
# weightEmInf=c(seq(1,.1,-.2),rep(0,36));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# summary(mix1)$varcomp
# res3$monitor
# 
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# 
# data(DT_cornhybrids)
# DT <- DT_cornhybrids
# DT<-DT[which(!is.na(DT$Yield)),]
# nrow(DT)
# DTi <- DTi_cornhybrids
# GT <- GT_cornhybrids
# hybrid2 <- DT # extract cross data
# A <- GT
# K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)
# K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
# S <- kronecker(K1, K2) ; dim(S)
# rownames(S) <- colnames(S) <- levels(hybrid2$SCA)
# 
# hybrid2$Yield <- scale(hybrid2$Yield)
# 
# mix1 <- mmer(Yield ~ Location,
#              random = ~ vsr(GCA1,Gu=K1) + vsr(GCA2,Gu=K2),# + vsr(SCA,Gu=S),
#              rcov=~units,
#              data=hybrid2)
# summary(mix1)$varcomp
# 
# K1i <- as(solve(K1 + diag(ncol(K1))*1e-6 ),Class="sparseMatrix")
# K2i <- as(solve(K2 + diag(ncol(K2))*1e-6 ),Class="sparseMatrix")
# Si <- as(solve(S + diag(ncol(S))*1e-6 ),Class="sparseMatrix")
# 
# mix2 <- mmec(Yield ~ Location,
#              random = ~ vsc(ids(GCA1),Gu=K1i) + vsc(ids(GCA2),Gu=K2i),# + vsX(ids(SCA),Gu=Si),
#              rcov=~units, return.param = F,
#              data=hybrid2)
# mix2$monitor
# mix2$sigma
# 
# z1 <- model.matrix(~GCA1-1,data=DT);colnames(z1) <- gsub("GCA1","",colnames(z1))
# z2 <- model.matrix(~GCA2-1,data=DT);colnames(z2) <- gsub("GCA2","",colnames(z2))
# z3 <- model.matrix(~SCA-1,data=DT);colnames(z3) <- gsub("SCA","",colnames(z3))
# 
# Z <- list(
#   z1,z2,z3
# )
# 
# Zind <- 1:3
# 
# A <- list(
#   K1[colnames(z1),colnames(z1)], K2[colnames(z2),colnames(z2)], S[colnames(z3),colnames(z3)]
# ) #
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# thetaI <- list(
#   diag(1)*10,
#   diag(1)*10,
#   diag(1)*10,
#   matrix(50,1,1)
# );thetaI
# 
# thetaC <- list(
#   diag(1),
#   diag(1),
#   diag(1),
#   matrix(1,1,1)
# );thetaC
# 
# X <- model.matrix(~Location, data=DT)
# 
# y <- as.matrix(DT$Yield)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# weightInf=rep(1,40) # weights for the information matrix
# weightEmInf=rep(1,40)
# weightEmInf=c(seq(.9,.5,-.1),rep(0,36));weightEmInf # weights for the EM information matrix
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp2(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=TRUE,
#                       nIters=20, tolParConv=1e-5,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# # compare results
# summary(mix1)$varcomp
# res3$monitor
# 
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# ###############################################################################
# 
# library(sommer)
# 
# data(DT_h2)
# DT <- DT_h2
# 
# head(DT)
# length(table(DT$Env))
# nrow(DT)
# DT$y <- scale(DT$y)
# 
# ans1 <- mmer(y~Env,
#              random=~vsr(dsr(Env),Name) + vsr(dsr(Env),Block),
#              rcov=~vsr(dsr(Env),units),
#              data=DT)
# summary(ans1)$varcomp
# 
# ans2 <- mmec(y~Env,
#              random=~vsc(dsr(Env),ids(Name)) + vsc(dsr(Env),ids(Block)),
#              rcov=~vsc(dsr(Env),ids(units)), return.param = F,
#              # emweight = c(1,rep(0,20)),
#              # stepweight = c(1,1,1,rep(0,30)),
#              nIters=30,
#              data=DT)
# 
# plotMonitor(ans2$monitor, cex=.1)
# ans2$sigma
# 
# d=ans2$monitor
# plotMonitor(d)
# plot(d[,(ncol(d)-1)], summary(ans1)$varcomp[,1])
# 
# DT2 <- DT[with(DT, order(Env)), ]
# library(asreml)
# ans3 <- asreml(y~Env,
#               random=~diag(Env):Name + diag(Env):Block,
#               residual=~dsum(~units | Env),
#               
#               data=DT2)
# 
# plot(d[,(ncol(d)-1)], summary(ans3)$varcomp[,1])
# 
# 
# xx=with(DT,vsr(dsr(Env),Name))
# 
# Z <- xx$Z
# 
# Zind <- rep(1,length(Z))
# 
# A <- list(
#   
#   diag(41)
# ) #
# 
# Ai <- lapply(A, function(x){solve(x)})
# 
# thetaI <- list(
#   diag(15)*5,#+diag(rnorm(15)),
#   matrix(10,1,1)
# );thetaI
# 
# thetaC <- list(
#   diag(15),
#   matrix(1,1,1)
# );thetaC
# 
# X <- model.matrix(~Env, data=DT)
# 
# y <- as.matrix(DT$y)
# 
# S <- list(diag(length(y)))
# 
# H <- diag(length(y))
# 
# addScaleParam <- 0
# nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
# nn2 <- sum(nn[1:max(Zind)])
# ff <- diag(nn2)
# thetaF <- cbind(ff,matrix(0,nn2,1))
# thetaF
# 
# ## apply the function
# weightInf=rep(1,40) # weights for the information matrix
# 
# weightEmInf =c(seq(1,.1,-.1),rep(.01,30))
# # weightEmInf =rep(.01,30)
# # weightInfEMv=c(rep(1,8),rep(0,36));weightInfEMv # weights for the EM information matrix
# 
# 
# Z <- lapply(Z, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# S <- lapply(S, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# Ai <- lapply(Ai, function(x){
#   return(as(x, Class = "sparseMatrix"))
# })
# 
# H <- as(H, Class = "sparseMatrix")
# 
# X <- as(X, Class = "sparseMatrix")
# 
# y <- as(y, Class = "sparseMatrix")
# 
# tt=system.time( 
#   expr = res3<-ai_mme_sp(X=X,Z=Z, Zind=Zind, 
#                       Ai=Ai,y=y,
#                       S=S, 
#                       H=H, useH=FALSE,
#                       nIters=15, tolParConv=1e-4,
#                       tolParInv=1e-6,thetaI=thetaI, 
#                       thetaC=thetaC,thetaF=thetaF, 
#                       addScaleParam=addScaleParam, weightEmInf = weightEmInf, 
#                       weightInf = weightInf,
#                       verbose=TRUE
#                       
#   )
# )
# 
# plot(res3$monitor[,ncol(res3$monitor)], summary(ans1)$varcomp[,1] )
# 
# plot(res3$llik[1,])
# 
# #######################################################
# 
# data(DT_yatesoats)
# DT <- DT_yatesoats
# head(DT)
# DT$Y <- scale(DT$Y)
# m3 <- mmer(fixed=Y ~ V + N + V:N,
#            # random = ~ B + B:MP,
#            rcov=~units,
#            data = DT)
# summary(m3)$varcomp
# 
# m4 <- mmec(fixed=Y ~ V + N + V:N,
#            # random = ~ B + B:MP,
#            rcov=~units, return.param = F,
#            data = DT)
# m4$sigma
# plotMonitor(m4$monitor)
# str(m4)
# 
# 
# ####################################################
# 
# library(sommer)
# data("DT_cpdata")
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# DT$Yield <- imputev(DT$Yield)
# #### create the variance-covariance matrix
# K <- A.mat(GT) # additive relationship matrix
# K <- K + diag(1e-3,nrow(K))
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#              random=~vsr(id,Gu=K)
#              + Rowf,
#              rcov=~units,
#              data=DT)
# summary(mix1)$varcomp
# 
# # Ki <- as(solve(K), Class="sparseMatrix")
# m <- GT[,1:300]
# mix2 <- mmec(fixed=Yield~1,
#               random=~vsc(ids(m)) + Rowf,
#               rcov=~units, return.param = F,
#               data=DT)
# mix2$monitor
# summary(mix1)$varcomp
# 
# 
# 
# ?overlay
# 
# data("DT_halfdiallel")
# DT <- DT_halfdiallel
# head(DT)
# DT$femalef <- as.factor(DT$female)
# DT$malef <- as.factor(DT$male)
# DT$genof <- as.factor(DT$geno)
# 
# A <- diag(7); colnames(A) <- rownames(A) <- 1:7;A # if you want to provide a covariance matrix
# #### model using overlay
# modh <- mmer(sugar~1, 
#              random=~vsr(overlay(femalef,malef), Gu=A) 
#              + genof,
#              data=DT)
# 
# Ai <- as(solve(A), Class="sparseMatrix")
# modh <- mmec(sugar~1, 
#              random=~vsc(ids(overlay(femalef,malef)), Gu=Ai) 
#              + genof, #return.param = T,
#              data=DT)
# modh$monitor
# 
# #####################################################
# 
# ?DT_cpdata
# 
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# Ai <- as( solve(A+diag(1e-5,ncol(A),ncol(A))), Class = "sparseMatrix")
# 
# head(DT)
# DT <- DT[,-c(7:8)]
# DT$color <- as.vector(scale(DT$color))
# DT$Yield <- as.vector(scale(DT$Yield))
# head(DT)
# DT2 <- reshape(DT, idvar = "id", varying = list(5:6),
#                v.names = "y", direction = "long", timevar = "trait", times =colnames(DT)[5:6] )
# DT2$trait <- as.factor(DT2$trait)
# head(DT2)
# 
# g=diag(2)*.05 + matrix(.1,2,2);g
# g2 <- diag(2)*.75;g2
# mix1 <- mmec(y~trait-1,
#              random=~vsc(usr(trait, theta = g),ids(id),Gu=Ai),
#              stepweight = ss,
#              emweight =rep(1,30),
#              return.param = F, tolParConv = 1e-6,
#              rcov=~vsc(dsr(trait),ids(units)),
#              data=DT2)
# 
# mix1$theta
# mix1$thetaC
# mix1$sigma
# mix1$b
# plot(mix1$u)
# plot(mix1$llik[1,])
# 
# # > summary(mix2)$varcomp
# # VarComp  VarCompSE    Zratio Constraint
# # u:id.color-color    0.7557346 0.15068202  5.015426   Positive
# # u:id.color-Yield    0.0848963 0.07171699  1.183768   Unconstr
# # u:id.Yield-Yield    0.1394332 0.07012037  1.988484   Positive
# # u:units.color-color 0.3779553 0.04195655  9.008254   Positive
# # u:units.Yield-Yield 0.8808832 0.07522796 11.709519   Positive
# # # 
# 
