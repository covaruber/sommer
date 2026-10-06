test_that("diagonal VarU batches agree with full covariance beyond 128 effects", {
  set.seed(512)
  data <- data.frame(id=factor(rep(seq_len(135),each=2)))
  data$y <- rep(rnorm(135),each=2)+rnorm(nrow(data),sd=0.5)
  model <- mmes(y~1,random=~id,data=data,verbose=FALSE,computeCi=0)
  diagonal <- postVarU(model)
  full <- postVarU(model,2)
  expect_equal(as.numeric(diagonal$uVarList[[1]]),unname(diag(full$VarU)),tolerance=1e-9)
  expect_false(isTRUE(diagonal$CiComputed))
  expect_equal(diagonal$Ci,model$Ci)
})

test_that("fixed-only models have zero random VarU", {
  data <- data.frame(y=c(2,4,3,5,7,6),group=factor(rep(c("a","b"),each=3)))
  model <- mmes(y~group,data=data,verbose=FALSE)
  full <- postVarU(model,2)
  expect_length(full$uVarList,0L)
  expect_equal(dim(full$VarU),c(0L,0L))
  prediction <- predict(model,D="group")
  expect_equal(unname(prediction$VarU),matrix(0,2,2))
  expect_equal(prediction$sampling.vcov,prediction$PEV)
})

test_that("saved precision makes VarU independent of later external Gu changes", {
  set.seed(513)
  ids <- paste0("g",seq_len(8))
  relationship <- outer(seq_along(ids),seq_along(ids),function(first,second) 0.5^abs(first-second))
  dimnames(relationship) <- list(ids,ids)
  precision <- solve(relationship)
  attr(precision,"inverse") <- TRUE
  data <- data.frame(id=factor(rep(ids[1:7],each=3),levels=ids),y=rnorm(21))
  model <- mmes(y~1,random=~vsm(ism(id),Gu=precision),data=data,verbose=FALSE)
  original <- postVarU(model,2)
  precision[,] <- 0
  expect_equal(postVarU(model,2)$VarU,original$VarU)
  prediction <- predict(model,D="id",levels=setNames(list(c("g1","g8")),names(model$uList)))
  expect_equal(unname(prediction$VarU),unname(original$VarU[c(1,8),c(1,8)]))
})

test_that("Dtable averages project PEV and VarU through the final linear combination", {
  set.seed(514)
  data <- expand.grid(id=factor(seq_len(6)),env=factor(c("a","b","c")),replicate=1:2)
  data$y <- rep(rnorm(6),6)+rnorm(nrow(data))
  model <- mmes(y~env,random=~vsm(usm(env),ism(id)),data=data,verbose=FALSE,computeCi=2)
  Dt <- model$Dtable
  Dt$average[Dt$type == "fixed"] <- TRUE
  Dt$include[Dt$type == "random"] <- TRUE
  Dt$average[Dt$type == "random"] <- TRUE
  prediction <- predict(model,D="id",Dtable=Dt)
  full <- postVarU(model,2)
  randomD <- as.matrix(prediction$D[,-seq_len(nrow(model$b)),drop=FALSE])
  expect_equal(unname(prediction$VarU),unname(randomD %*% full$VarU %*% t(randomD)),tolerance=1e-8)
  expect_equal(unname(prediction$PEV),unname(as.matrix(prediction$D %*% model$Ci %*% t(prediction$D))),tolerance=1e-8)
  expect_equal(unique(as.numeric(randomD)),c(1/3,0))
  explicit <- predict(model,D=prediction$D)
  expect_equal(explicit$PEV,prediction$PEV)
  expect_equal(explicit$VarU,prediction$VarU)
})