varu_fixture <- function(){
  set.seed(410)
  ids <- paste0("g", seq_len(7))
  relationship <- outer(seq_along(ids), seq_along(ids), function(first, second) 0.55^abs(first-second))
  dimnames(relationship) <- list(ids, ids)
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  data <- expand.grid(id=factor(ids, levels=ids), env=factor(c("e1", "e2")),
                      replicate=seq_len(2), KEEP.OUT.ATTRS=FALSE)
  data$y <- 3+(data$env == "e2")+rep(rnorm(length(ids)),4)+rnorm(nrow(data), sd=0.8)
  list(data=data, precision=precision, relationship=relationship)
}

varu_dense_reference <- function(model){
  X <- as.matrix(model$W[,seq_len(nrow(model$b)),drop=FALSE])
  Z <- as.matrix(model$W[,-seq_len(nrow(model$b)),drop=FALSE])
  blocks <- lapply(seq_along(model$uList), function(index){
    Sigma <- covmatrix_mmes(model, index, se=FALSE)$covariance
    kronecker(Sigma, solve(as.matrix(model$randomPrecision[[index]])))
  })
  G <- as.matrix(do.call(Matrix::bdiag, blocks))
  residual <- tail(model$theta,1)[[1]][1,1]
  V <- Z %*% G %*% t(Z)+diag(residual,nrow(X))
  Vi <- solve(V)
  B <- solve(crossprod(X, Vi %*% X))
  P <- Vi-Vi %*% X %*% B %*% t(X) %*% Vi
  S <- G %*% t(Z) %*% P %*% Z %*% G
  list(G=G, S=S, Q=G-S, B=B)
}

test_that("postVarU matches dense reference and preserves PEV outputs", {
  fixture <- varu_fixture()
  for(REML in c(TRUE,FALSE)){
    model <- mmes(y~env, random=~vsm(usm(env), ism(id), Gu=fixture$precision),
                  data=fixture$data, REML=REML, verbose=FALSE, computeCi=0)
    reference <- varu_dense_reference(model)
    original <- postPEV(model,2)
    diagonal <- postVarU(original)
    full <- postVarU(original,2)
    expect_equal(unname(full$VarU), reference$S, tolerance=1e-7)
    expect_equal(as.numeric(diagonal$uVarList[[1]]), diag(reference$S), tolerance=1e-7)
    expect_equal(diagonal$uVarList, full$uVarList, tolerance=1e-7)
    expect_equal(dimnames(full$uVarList[[1]]), dimnames(model$uList[[1]]))
    expect_equal(full$Ci, original$Ci)
    expect_equal(full$uPevList, original$uPevList)
    expect_true(is.null(diagonal$VarU))
    cleared <- postVarU(full,0)
    expect_null(cleared$VarU)
    expect_null(cleared$uVarList)
    expect_equal(cleared$uPevList, original$uPevList)
    legacy <- model
    legacy$randomPrecision <- NULL
    expect_equal(postVarU(legacy)$uVarList, postVarU(model)$uVarList, tolerance=1e-7)
  }
})

test_that("prediction distinguishes full PEV, random VarU and sampling covariance", {
  fixture <- varu_fixture()
  model <- mmes(y~env, random=~vsm(usm(env), ism(id), Gu=fixture$precision),
                data=fixture$data, verbose=FALSE, computeCi=2)
  reference <- varu_dense_reference(model)
  nFixed <- nrow(model$b)
  nRandom <- nrow(model$u)
  D <- matrix(0,3,nFixed+nRandom)
  D[,1] <- 1
  D[1,nFixed+1] <- 1
  D[2,nFixed+c(2,8)] <- c(0.5,0.5)
  D[3,seq_len(nFixed)] <- c(1,1)
  prediction <- predict(model,D=D,sed=TRUE,pairwise=TRUE)
  Db <- D[,seq_len(nFixed),drop=FALSE]
  Du <- D[,-seq_len(nFixed),drop=FALSE]
  expectedU <- Du %*% reference$S %*% t(Du)
  expect_equal(unname(prediction$VarU), expectedU, tolerance=1e-7)
  expect_equal(unname(prediction$sampling.vcov),
               Db %*% reference$B %*% t(Db)+expectedU, tolerance=1e-7)
  expect_equal(unname(prediction$PEV), as.matrix(D %*% model$Ci %*% t(D)), tolerance=1e-7)
  expect_equal(prediction$PEV, prediction$vcov)
  expect_equal(prediction$pvals$pev, unname(diag(prediction$PEV)))
  expect_equal(prediction$pvals$var.u, unname(diag(prediction$VarU)))
  expect_equal(prediction$pvals$std.error^2, unname(diag(prediction$PEV)))
  expect_equal(prediction$VarU[3,3], 0)
  for(PEV in c(TRUE,FALSE)) for(VarU in c(TRUE,FALSE)){
    toggled <- predict(model,D=D,PEV=PEV,VarU=VarU,sed=TRUE,pairwise=TRUE)
    expect_equal(toggled$vcov,prediction$vcov)
    expect_equal(toggled$sed,prediction$sed)
    expect_equal(toggled$pairwise,prediction$pairwise)
    expect_identical(!is.null(toggled$PEV),PEV)
    expect_identical(!is.null(toggled$VarU),VarU)
    expect_identical(!is.null(toggled$sampling.vcov),VarU)
  }
})

test_that("joint VarU retains sampling covariance between independent prior terms", {
  fixture <- varu_fixture()
  model <- mmes(y~env,random=~id+vsm(usm(env),ism(id),Gu=fixture$precision),
                data=fixture$data,verbose=FALSE,nIters=5,computeCi=2)
  reference <- varu_dense_reference(model)
  full <- postVarU(model,2)
  expect_equal(unname(full$VarU),reference$S,tolerance=1e-7)
  expect_gt(max(abs(full$VarU[seq_len(7),8:21])),1e-8)
  expect_equal(sum(lengths(full$uVarList)),nrow(full$VarU))
})

test_that("VarU maps rotation and factor-score effects into public coordinates", {
  fixture <- varu_fixture()
  data <- fixture$data[fixture$data$replicate == 1,]
  ordinary <- mmes(y~env,random=~vsm(ism(id),Gu=fixture$precision),data=data,verbose=FALSE)
  rotated <- mmes(y~env,random=~vsm(ism(id),Gu=fixture$precision,rotation=TRUE),data=data,verbose=FALSE)
  expect_equal(postVarU(rotated)$uVarList[[1]],postVarU(ordinary)$uVarList[[1]],tolerance=1e-6)
  expect_equal(unname(postVarU(rotated,2)$VarU),unname(postVarU(ordinary,2)$VarU),tolerance=1e-6)
  expect_equal(predict(rotated,D="id")$VarU,predict(ordinary,D="id")$VarU,tolerance=1e-6)
  data <- fixture$data
  shape <- fam(data$env,1L,fixed=rep(TRUE,3L))
  fit <- function(mode){
        mmes(y~env,random=~vsm(shape,ism(id),Gu=fixture$precision,sigma2=1,fixedSigma2=TRUE),
          rcov=~vsm(ism(units),sigma2=1,fixedSigma2=TRUE),data=data,
          verbose=FALSE,factorScoreAugmentation=mode)
  }
  marginal <- fit("none")
  augmented <- fit("fixed-shape")
  for(model in list(marginal,augmented)){
    reference <- varu_dense_reference(model)
    expect_equal(as.numeric(postVarU(model)$uVarList[[1]]),diag(reference$S),tolerance=1e-7)
    expect_equal(unname(postVarU(model,2)$VarU),reference$S,tolerance=1e-7)
  }
})

test_that("VarU validates flags and unsupported model types", {
  fixture <- varu_fixture()
  model <- mmes(y~env,random=~id,data=fixture$data,verbose=FALSE)
  expect_error(postVarU(model,3),"mode")
  expect_error(predict(model,D="id",VarU=NA),"TRUE or FALSE")
  expect_error(predict(model,D="id",PEV=1),"TRUE or FALSE")
  direct <- mmes(y~env,random=~id,data=fixture$data,verbose=FALSE,henderson=FALSE)
  expect_error(postVarU(direct),"Henderson")
  class(model) <- c("mmes.glmm",class(model))
  expect_error(postVarU(model),"Gaussian")
  expect_error(predict(model,D="id"),"Gaussian")
  expect_silent(predict(model,D="id",VarU=FALSE))
})