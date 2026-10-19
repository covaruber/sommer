## ----setup, message=FALSE-----------------------------------------------------
library(sommer)
library(Matrix)
set.seed(4)

## ----pedigree-----------------------------------------------------------------
simulate_pedigree <- function(nGen, perGen, nSires){
  n <- nGen * perGen
  gen <- rep(seq_len(nGen), each=perGen)
  male <- rep(c(TRUE, FALSE), length.out=n)
  sire <- dam <- rep(NA_integer_, n)
  for(g in 2:nGen){
    prev <- which(gen == g - 1L)
    cur <- which(gen == g)
    sires <- sample(prev[male[prev]], nSires)
    sire[cur] <- sample(sires, length(cur), replace=TRUE)
    dam[cur] <- sample(prev[!male[prev]], length(cur), replace=TRUE)
  }
  data.frame(index=seq_len(n), id=paste0("A", seq_len(n)), sire=sire, dam=dam,
             gen=gen)
}

ped <- simulate_pedigree(nGen=6, perGen=3000, nSires=60)
head(ped[ped$gen == 2, ])

## ----ainv---------------------------------------------------------------------
pedigree_ainv <- function(ped){
  n <- nrow(ped)
  hasSire <- !is.na(ped$sire)
  hasDam <- !is.na(ped$dam)
  P <- sparseMatrix(i=c(which(hasSire), which(hasDam)),
                    j=c(ped$sire[hasSire], ped$dam[hasDam]),
                    x=0.5, dims=c(n, n))
  Tm <- Diagonal(n) - P
  dinv <- 1 / (1 - 0.25 * (hasSire + hasDam))
  Ainv <- forceSymmetric(crossprod(Tm, Diagonal(x=dinv) %*% Tm))
  dimnames(Ainv) <- list(ped$id, ped$id)
  attr(Ainv, "inverse") <- TRUE
  Ainv
}

## ----truth--------------------------------------------------------------------
G0 <- matrix(c(0.30, 0.23,
               0.23, 0.50), 2, dimnames=list(c("y1","y2"), c("y1","y2")))
R0 <- matrix(c(0.70, 0.24,
               0.24, 0.90), 2, dimnames=list(c("y1","y2"), c("y1","y2")))
cov2cor(G0)[1, 2]   # genetic correlation
diag(G0) / (diag(G0) + diag(R0))   # heritabilities

## ----phenotypes---------------------------------------------------------------
n <- nrow(ped)
bv <- matrix(0, n, 2, dimnames=list(ped$id, c("y1", "y2")))
dvar <- 1 - 0.25 * ((!is.na(ped$sire)) + (!is.na(ped$dam)))
mendelian <- matrix(rnorm(n * 2), n) %*% chol(G0)
for(g in sort(unique(ped$gen))){
  cur <- which(ped$gen == g)
  pa <- if(g == 1) 0 else 0.5 * (bv[ped$sire[cur], ] + bv[ped$dam[cur], ])
  bv[cur, ] <- pa + sqrt(dvar[cur]) * mendelian[cur, ]
}

recorded <- ped[ped$gen < 6, ]
nHerd <- 300
recorded$herd <- factor(sample(seq_len(nHerd), nrow(recorded), replace=TRUE))
herdEffect <- matrix(rnorm(nHerd * 2, sd=c(1, 1.5)), nHerd, byrow=TRUE)
errors <- matrix(rnorm(nrow(recorded) * 2), ncol=2) %*% chol(R0)
recorded$y1 <- 10 + herdEffect[recorded$herd, 1] + bv[recorded$index, 1] + errors[, 1]
recorded$y2 <- 20 + herdEffect[recorded$herd, 2] + bv[recorded$index, 2] + errors[, 2]
recorded$y2[sample(nrow(recorded), round(0.3 * nrow(recorded)))] <- NA
dim(recorded)

## ----long---------------------------------------------------------------------
pheno <- recorded[, c("id", "herd", "gen", "y1", "y2")]
long <- stackTraits(pheno, traits=c("y1", "y2"))
head(long)

## ----subset-------------------------------------------------------------------
trace_pedigree <- function(ped, index){
  keep <- rep(FALSE, nrow(ped))
  todo <- index
  while(length(todo)){
    keep[todo] <- TRUE
    parents <- c(ped$sire[todo], ped$dam[todo])
    todo <- unique(parents[!is.na(parents) & !keep[parents]])
  }
  sub <- ped[keep, ]
  sub$sire <- match(sub$sire, sub$index)
  sub$dam <- match(sub$dam, sub$index)
  sub
}

pilotHerds <- levels(long$herd)[1:60]
pilot <- droplevels(long[long$herd %in% pilotHerds, ])
pilotPed <- trace_pedigree(ped, recorded$index[recorded$herd %in% pilotHerds])
c(records=nrow(pilot), animals=length(unique(pilot$id)), pedigree=nrow(pilotPed))

## ----reml---------------------------------------------------------------------
Ainv <- pedigree_ainv(pilotPed)
timeReml <- system.time(
  reml <- mmes(value ~ trait + trait:herd,
  random = ~ vsm(usm(trait), ism(id), Gu=Ainv),
  rcov = ~ vsm(usm(trait), ism(record)),
  data=pilot,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
)
timeReml[["elapsed"]]

## ----remlEstimates------------------------------------------------------------
G0hat <- covmatrix_mmes(reml, 1)
R0hat <- covmatrix_mmes(reml, 2)
round(G0hat$covariance, 3)
round(G0hat$covariance.se, 3)
round(R0hat$covariance, 3)

## ----evaluation---------------------------------------------------------------
Ainv <- pedigree_ainv(ped)
timeEval <- system.time(
  ebv <- mmes(value ~ trait + trait:herd,
  random = ~ vsm(usm(trait), ism(id), Gu=Ainv),
  rcov = ~ vsm(usm(trait), ism(record)),
  data=long,
  covPar=reml,
  control = mmes_control(task = list(mode = "solve"),
  diagnostics = list(verbose = FALSE)))
)
timeEval[["elapsed"]]
ebv

## ----convergence, fig.width=6, fig.height=3.5---------------------------------
ebv$pcg[c("iterations", "relres", "converged", "setupSeconds", "solveSeconds")]
plot(log10(ebv$pcg$history), type="l", xlab="PCG iteration",
     ylab="log10 relative residual")
abline(h=log10(1e-8), lty=2)

## ----accuracy-----------------------------------------------------------------
u <- ebv$uList[[1]][ped$id, ]
accuracy <- function(rows) diag(cor(u[rows, ], bv[rows, ]))
rbind(recorded = accuracy(ped$gen < 6),
      candidates = accuracy(ped$gen == 6))

## ----ranking------------------------------------------------------------------
candidates <- ped$id[ped$gen == 6]
index <- u[candidates, "y1"] + u[candidates, "y2"]
top <- head(sort(index, decreasing=TRUE), 5)
data.frame(id=names(top), index=round(top, 3),
           trueIndex=round(rowSums(bv[names(top), ]), 3))

## ----truePars-----------------------------------------------------------------
ebvTrue <- mmes(value ~ trait + trait:herd,
  random = ~ vsm(usm(trait), ism(id), Gu=Ainv),
  rcov = ~ vsm(usm(trait), ism(record)),
  data=long,
  covPar=list(G0, R0),
  control = mmes_control(task = list(mode = "solve"),
  diagnostics = list(verbose = FALSE)))
uTrue <- ebvTrue$uList[[1]][ped$id, ]
diag(cor(u, uTrue))
diag(cor(uTrue[ped$gen == 6, ], bv[ped$gen == 6, ]))
ebvTrue$vcParams[, c("term", "parameter", "value")]

## ----singleTrait--------------------------------------------------------------
single <- mmes(y1 ~ herd,
  random = ~ vsm(ism(id), Gu=Ainv, sigma2=0.30, fixedSigma2=TRUE),
  rcov = ~ vsm(ism(units), sigma2=0.70, fixedSigma2=TRUE),
  data=pheno,
  control = mmes_control(task = list(mode = "solve"),
  diagnostics = list(verbose = FALSE)))
cor(single$uList[[1]][candidates, 1], bv[candidates, "y1"])

