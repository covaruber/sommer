



data(DT_rice)
DT <- DT_rice
GT <- GT_rice
GTn <- GTn_rice
head(DT)
M <- atcg1234(GT)



### univariate model
A <- A.mat(M$M)
A <- A + diag(1e-4, ncol(A), ncol(A))

### multi-trait model
DT2 <- stackTrait(data=DT, traits = colnames(DT)[3:20])$long
DT2$trait <- as.factor(DT2$trait)
head(DT2)
table(DT2$trait)
system.time(
  mix <- lmebreed(valueS ~ (0+trait|geno),
                  relmat = list(geno=A),
                  control = lmerControl(
                    check.nobs.vs.nlev = "ignore",
                    check.nobs.vs.rankZ = "ignore",
                    check.nobs.vs.nRE="ignore"
                  ), rotation = TRUE,
                  data=DT2)
)
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
lattice::levelplot(cov2cor(vc$geno))

colfunc <- colorRampPalette(c("yellow","orange","red"))
hv <- heatmap(cov2cor(vc$geno), col = colfunc(100),Colv = "Rowv")
str(hv)

