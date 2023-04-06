
gg=with(DT, rrc(Env, Name, y, nPC = 2))
timevar=DT$Env; idvar=DT$Name; response=DT$y; 
Gu=NULL; nPC=1; returnLam=FALSE; cholD=TRUE

library(sommer)
data(DT_h2)
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]

ans1b <- mmec(y~Env,
              random=~vsc( usc( rrc(Env, Name, y, nPC = 3) ) , isc(Name)) + vsc(dsc( rrc(Env, Block, y, nPC = 3)  ),isc(Block)),
              rcov=~units,#vsc(dsc(Env),isc(units)), 
              # we recommend giving more iterations to these models
              nIters = 50, 
              # we recommend giving more EM iterations at the beggining for usc models
              emWeight = c(rep(1,10),logspace(10,1,.05), rep(.05,80)),
              # verbose=FALSE,
              data=DT)
dim(ans1b$Ci)
summary(ans1b)$varcomp
vc <- diag(ans1b$theta[[1]])
vc/sum(vc)
vc2 <- diag(ans1b$theta[[2]])
vc2/sum(vc2)
myr2 <- r2(ans1b)
apply(myr2[[1]],2,mean)

Lam=with(DT, rrc(Env, Name, y, returnLam = TRUE, nPC = 3))$Lam # extract loadings
score.mat <- ans1b$uList[[1]]; # extract factor scores
BLUP = score.mat %*% t(Lam) # BLUPs for all environments

library(asreml)

ans3c <- asreml(y~Env,
                random=~rr(Env,2):Name + diag(Env):Name,# + vsc(dsc(Env),isc(Block)),
                residual=~units, #dsum(~units | Env), 
                maxit = 50,
                # emWeight = rep(1,50),
                data=DT)
summary(ans3c)$varcomp
source("~/Downloads/BCS.FieldTrialData 2/R/gammasRr4.R")
oo=gammasRr4(model=ans3c,gxeTerm = "Env:Name", pedigree = NULL)

pp <- predict(ans3c, classify = "Env:Name")
head(pp$pvals)
wide <- reshape(pp$pvals[,c("Env","Name","predicted.value")], direction = "wide", idvar = "Name",
                timevar = "Env", v.names = "predicted.value", sep= "_")
head(wide)
colnames(wide) <- gsub("predicted.value_","", colnames(wide))

lattice::levelplot(cor(wide[,(colnames(BLUP))[-c(5)]]))
lattice::levelplot(cor(BLUP[,-5]))

lattice::levelplot(cor(wide[,(colnames(BLUP))]))
lattice::levelplot(oo$Cmat)
lattice::levelplot(cor(BLUP))
vc <- ans1b$theta[[1]];vc
G = Lam  %*% vc %*% t(Lam) # change space for matrix product
lattice::levelplot(cov2cor(G[-5,-5]))
lattice::levelplot(cov2cor(G))


for(i in 1:ncol(BLUP)){
  
  plot(BLUP[,i], wide[,(colnames(BLUP))[i]], main=i)
  
}

main <- apply(BLUP,1,mean,na.rm=TRUE)
main2 <- apply(wide[,-1],1,mean,na.rm=TRUE)
plot(main,main2)

qq=predict(ans1b,D="Name")
plot(qq$pvals$predicted.value,main2)
