.direct_inference_control <- function(formulation="direct", ...){
  mmes_control(formulation=formulation, ...,
    diagnostics=list(verbose=FALSE, dateWarning=FALSE),
    convergence=list(maxIters=100L, likelihoodTol=1e-8, parameterTol=1e-8))
}

.direct_inference_data <- function(){
  set.seed(140)
  data <- expand.grid(group=factor(letters[1:6]), visit=factor(1:4), rep=1:2)
  data$x <- rnorm(nrow(data))
  data$y <- 3 + 0.2 * data$x + rnorm(6L)[data$group] + rnorm(nrow(data))
  data
}

.direct_inference_reference <- function(X, Z, G, R){
  V <- Z %*% G %*% t(Z) + R
  Vi <- solve(V)
  Phi <- solve(crossprod(X, Vi %*% X))
  cross <- -Phi %*% crossprod(X, Vi %*% Z) %*% G
  P <- Vi - Vi %*% X %*% Phi %*% t(X) %*% Vi
  rbind(cbind(Phi, cross), cbind(t(cross), G - G %*% t(Z) %*% P %*% Z %*% G))
}

test_that("sparse identity detection checks stored entries and dimensions", {
  expect_true(sommer:::isIdentity_spmat(Matrix::Diagonal(10000L)))
  expect_true(sommer:::isIdentity_spmat(Matrix::Matrix(diag(1, 4L), sparse=TRUE)))
  expect_false(sommer:::isIdentity_spmat(Matrix::Diagonal(4L, x=c(1, 1, 1, 0))))
  expect_false(sommer:::isIdentity_spmat(
    Matrix::sparseMatrix(i=1L, j=2L, x=1, dims=c(4L, 4L))))
  expect_false(sommer:::isIdentity_spmat(Matrix::Matrix(0, 3L, 4L, sparse=TRUE)))
})

test_that("direct post-fit covariance matches independent joint error algebra", {
  data <- .direct_inference_data()
  K <- 0.25^abs(outer(1:6, 1:6, "-"))
  A <- solve(K)
  dimnames(A) <- list(levels(data$group), levels(data$group))
  attr(A, "inverse") <- TRUE
  fit <- mmes(y~x, random=~vsm(ism(group), Gu=A, sigma2=0.7, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=0.4, fixedSigma2=TRUE), data=data,
    control=.direct_inference_control())
  X <- as.matrix(fit$W[, 1:2])
  Z <- as.matrix(fit$W[, -(1:2)])
  reference <- .direct_inference_reference(X, Z, 0.7*K, diag(0.4, nrow(data)))
  expect_null(fit$C)
  expect_equal(dim(fit$Ci), c(0L, 0L))
  expect_true(fit$convergence)
  full <- postPEV(fit, 2L)
  diagonal <- postPEV(fit)
  expect_equal(as.matrix(full$Ci), unname(reference), tolerance=1e-9)
  expect_equal(as.numeric(diag(diagonal$Ci)), unname(diag(reference)), tolerance=1e-9)
  expect_false(diagonal$CiComputed)
  expect_true(full$CiComputed)
  expect_identical(full$CiMode, 2L)
  expect_equal(full$uPevList, diagonal$uPevList)
  expect_equal(as.numeric(full$uPevList[[1L]]), unname(diag(reference)[-(1:2)]), tolerance=1e-9)
  D <- rbind(c(1, 0, rep(1/6, 6)), c(0, 0.5, 1, -1, 0, 0, 0, 0))
  result <- predict(fit, D=D, sed=TRUE, pairwise=TRUE, df="residual")
  expected <- D %*% reference %*% t(D)
  expect_equal(result$vcov, expected, tolerance=1e-9)
  expect_equal(result$PEV, expected, tolerance=1e-9)
  expect_equal(result$pvals$predicted.value, as.numeric(D %*% fit$bu), tolerance=1e-10)
  expect_equal(result$pvals$std.error, sqrt(diag(expected)), tolerance=1e-9)
  expect_equal(result$pairwise$SED, sqrt(expected[1,1]+expected[2,2]-2*expected[1,2]))
  expect_equal(result$pairwise$df, nrow(data)-ncol(X))
  randomD <- D[, -(1:2), drop=FALSE]
  prior <- randomD %*% (0.7*K) %*% t(randomD)
  expect_equal(unname(result$VarU), prior - randomD %*% reference[-(1:2), -(1:2)] %*% t(randomD),
    tolerance=1e-9)
  expect_null(predict(fit, D=D, PEV=FALSE, VarU=FALSE)$PEV)
  expect_equal(predict(full, D=D)$vcov, predict(fit, D=D)$vcov)
  cleared <- postPEV(full, 0L)
  expect_equal(dim(cleared$Ci), c(0L, 0L))
  expect_length(cleared$uPevList[[1L]], 0L)
  expect_false(cleared$CiComputed)
  expect_equal(cleared$bu, fit$bu)
  expect_equal(cleared$covPar, fit$covPar)
})

test_that("direct predictions match Henderson with weighted structured residuals and filtering", {
  data <- .direct_inference_data()
  data$y[c(4, 17)] <- NA_real_
  data$x[29] <- NA_real_
  W <- Matrix::bandSparse(nrow(data), k=c(-1, 0, 1),
    diagonals=list(rep(-0.15, nrow(data)-1L), rep(1.4, nrow(data)),
      rep(-0.15, nrow(data)-1L)))
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(y~x, random=~vsm(ism(group), sigma2=0.7, fixedSigma2=TRUE),
      rcov=~vsm(ar1m(visit, rho=0.2, fixed=TRUE), ism(group), ism(rep),
        sigma2=0.4, fixedSigma2=TRUE), W=W, data=data,
      control=.direct_inference_control(formulation)))
  expect_equal(fits[[1]]$obsInfo, fits[[2]]$obsInfo)
  for(classify in c("group", "x")){
    direct <- predict(fits[[1]], D=classify, sed=TRUE, pairwise=TRUE)
    henderson <- predict(fits[[2]], D=classify, sed=TRUE, pairwise=TRUE)
    for(name in c("D", "vcov", "PEV", "VarU", "sampling.vcov", "pvals", "sed", "pairwise")){
      expect_equal(direct[[name]], henderson[[name]], tolerance=2e-7, info=name)
    }
  }
  for(mode in 1:2){
    direct <- postPEV(fits[[1]], mode)
    henderson <- postPEV(fits[[2]], mode)
    expect_equal(direct$uPevList, henderson$uPevList, tolerance=2e-7)
    if(mode == 2) expect_equal(as.matrix(direct$Ci), as.matrix(henderson$Ci), tolerance=2e-7)
    expect_equal(postVarU(fits[[1]], mode)$uVarList,
      postVarU(fits[[2]], mode)$uVarList, tolerance=2e-7)
  }
})

test_that("direct inference covers structured random effects and cross-term error covariance", {
  data <- .direct_inference_data()
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(y~visit, random=~vsm(usm(visit, theta=0.3^abs(outer(1:4, 1:4, "-")),
      fixed=rep.int(TRUE, 9L)),
      ism(group), sigma2=0.5, fixedSigma2=TRUE) +
      vsm(ism(rep), sigma2=0.2, fixedSigma2=TRUE), data=data,
      rcov=~vsm(ism(units), sigma2=0.4, fixedSigma2=TRUE),
      control=.direct_inference_control(formulation)))
  full <- lapply(fits, postPEV, mode=2L)
  expect_equal(as.matrix(full[[1]]$Ci), as.matrix(full[[2]]$Ci), tolerance=1e-8)
  expect_equal(full[[1]]$uPevList, full[[2]]$uPevList, tolerance=1e-8)
  for(classify in c("group", "visit:group")){
    direct <- predict(fits[[1]], D=classify, sed=TRUE)
    henderson <- predict(fits[[2]], D=classify, sed=TRUE)
    expect_equal(direct$D, henderson$D)
    expect_equal(direct$vcov, henderson$vcov, tolerance=1e-8)
    expect_equal(direct$VarU, henderson$VarU, tolerance=1e-8)
  }
})

test_that("direct fixed-only and reduced-rank fits preserve public estimability", {
  data <- .direct_inference_data()
  fit <- mmes(y~x, data=data, control=.direct_inference_control())
  reference <- lm(y~x, data=data)
  expect_equal(postPEV(fit, 2L)$Ci, unname(vcov(reference)), tolerance=1e-8)
  expect_length(postPEV(fit)$uPevList, 0L)
  expect_equal(predict(fit, D=as.matrix(fit$W))$pvals$predicted.value,
    as.numeric(fitted(reference)), tolerance=1e-8)

  data$duplicate <- data$x
  fits <- lapply(c("direct", "henderson"), function(formulation){
    expect_warning(result <- mmes(y~x+duplicate,
      random=~vsm(ism(group), sigma2=0.5, fixedSigma2=TRUE), data=data,
      control=.direct_inference_control(formulation,
        representation=list(fixedRank="drop"))), "rank reduction")
    result
  })
  expect_equal(predict(fits[[1]], D=fits[[1]]$W)$vcov,
    predict(fits[[2]], D=fits[[2]]$W)$vcov, tolerance=1e-7)
  expect_equal(as.matrix(postPEV(fits[[1]], 2L)$Ci),
    as.matrix(postPEV(fits[[2]], 2L)$Ci), tolerance=1e-7)
  nonestimable <- matrix(0, 1L, nrow(fits[[1]]$bu))
  nonestimable[1L, 2L] <- 1
  expect_error(predict(fits[[1]], D=nonestimable), "not estimable")
})

test_that("direct PQL inference uses the final working covariance without refitting", {
  data <- .direct_inference_data()
  data$exposure <- seq(0.8, 1.2, length.out=nrow(data))
  data$count <- rpois(nrow(data), data$exposure * exp(0.4+0.2*data$x))
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(count~x+offset(log(exposure)),
      random=~vsm(ism(group), sigma2=0.2, fixedSigma2=TRUE),
      data=data, family=poisson(),
      control=.direct_inference_control(formulation,
        task=list(pql=list(maxit=40, tol=1e-9)))))
  X <- as.matrix(fits[[1]]$W[, 1:2])
  Z <- as.matrix(fits[[1]]$W[, -(1:2)])
  expected <- .direct_inference_reference(X, Z, diag(0.2, ncol(Z)),
    diag(1/fits[[1]]$workingWeights))
  local_mocked_bindings(.mmes_fit=function(...) stop("Unexpected refit"),
    .package="sommer")
  expect_equal(postPEV(fits[[1]], 2L)$Ci, unname(expected), tolerance=1e-8)
  expect_equal(postPEV(fits[[1]])$uPevList, postPEV(fits[[2]])$uPevList, tolerance=1e-7)
  direct <- predict(fits[[1]], D="group", VarU=FALSE, pairwise=TRUE, sed=TRUE)
  henderson <- predict(fits[[2]], D="group", VarU=FALSE, pairwise=TRUE, sed=TRUE)
  expect_equal(direct$vcov, henderson$vcov, tolerance=1e-7)
  expect_equal(direct$pvals, henderson$pvals, tolerance=1e-7)
  expect_equal(direct$pairwise, henderson$pairwise, tolerance=1e-7)
  expect_error(predict(fits[[1]], D="group"), "only for Gaussian")
  expect_error(predict(fits[[1]], D="group", VarU=FALSE, pairwise=TRUE,
    df="satterthwaite"), "not available for PQL")
  expect_true(all(is.na(fits[[1]]$llik)))
})

test_that("direct post-fit operations honor memory controls and serialize without native state", {
  data <- .direct_inference_data()
  fit <- mmes(y~x, random=~group, data=data, control=.direct_inference_control())
  D <- cbind(Matrix::Diagonal(2L), Matrix::Matrix(0, 2L, nrow(fit$u), sparse=TRUE))
  result <- predict(fit, D=D, VarU=FALSE)
  restored <- unserialize(serialize(fit, NULL))
  restored$dataOriginal <- restored$inputArgs <- restored$args <- NULL
  expect_equal(predict(restored, D=D, VarU=FALSE), result)
  expect_equal(postPEV(restored)$uPevList, postPEV(fit)$uPevList)
  expect_null(restored$C)
  expect_error(predict(fit, D=D, maxMemoryMB=1e-5), "maxMemoryMB")
  expect_error(predict(fit, D=D, maxMemoryMB=0), "maxMemoryMB must be")
  fit$control$execution$denseMemoryMB <- 1e-5
  expect_error(postPEV(fit), "workspace exceeds")
  expect_error(postPEV(fit, 2L), "exceeds execution")
  expect_equal(dim(postPEV(fit, 0L)$Ci), c(0L, 0L))
  expect_error(postPEV(fit, NA_integer_), "mode")
  omitted <- mmes(y~x, data=data, control=.direct_inference_control(
    output=list(coefficientMatrix=FALSE)))
  expect_null(omitted$directInference)
  expect_error(predict(omitted, D=matrix(c(1,0), 1L)), "not retained")
  expect_error(postPEV(omitted), "not retained")
  legacy <- restored
  legacy$directInference <- NULL
  expect_error(postPEV(legacy), "refit")
  expect_error(predict(legacy, D=D), "refit")
  invalid <- D
  invalid[1, 1] <- NA_real_
  expect_error(predict(restored, D=invalid), "must be finite")
  invalid[1, ncol(invalid)] <- Inf
  expect_error(predict(restored, D=invalid), "must be finite")
  expect_error(predict(restored, D=matrix(0, 1L, 3L)), "one column per fitted")
})

test_that("direct diagonal PEV is panelled when coefficients outnumber observations", {
  data <- .direct_inference_data()
  labels <- c(letters[1:6], paste0("extra", 1:114))
  A <- Matrix::Diagonal(120L)
  dimnames(A) <- list(labels, labels)
  attr(A, "inverse") <- TRUE
  fit <- mmes(y~x, random=~vsm(ism(group), Gu=A, sigma2=0.7, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=0.4, fixedSigma2=TRUE), data=data,
    control=.direct_inference_control(output=list(panelSize=4L),
      execution=list(denseMemoryMB=0.7)))
  expect_gt(nrow(fit$bu), nrow(fit$W))
  diagonal <- postPEV(fit)
  unseen <- match("extra1", rownames(diagonal$uPevList[[1]]))
  expect_equal(diagonal$uPevList[[1]][unseen, 1], 0.7, tolerance=1e-10)
  expect_error(postPEV(fit, 2L), "Full direct coefficient covariance exceeds")
  selection <- setNames(list(c("extra1", "a")), names(fit$uList))
  prediction <- predict(fit, D="group", levels=selection, sed=TRUE, pairwise="extra1")
  expect_equal(as.character(prediction$pvals[[1]]), c("extra1", "a"))
  expect_equal(dim(prediction$vcov), c(2L, 2L))
  expect_equal(prediction$VarU[1,1], 0, tolerance=1e-10)
  expect_true(all(is.finite(prediction$pvals$std.error)))
  rules <- prediction$Dtable
  expect_equal(predict(fit, D="group", levels=selection, Dtable=rules)$D, prediction$D)
  expect_equal(predict(fit, D=prediction$D)$vcov, prediction$vcov)
})

test_that("direct post-fit inference supports singular reduced-rank random covariance", {
  data <- .direct_inference_data()
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(y~visit, random=~vsm(rrm(visit, 1L,
      loadings=matrix(c(1, 0.7, -0.5, 1.2), 4L), fixed=rep.int(TRUE, 4L)),
      ism(group), sigma2=0.6, fixedSigma2=TRUE),
      rcov=~vsm(ism(units), sigma2=0.4, fixedSigma2=TRUE), data=data,
      control=.direct_inference_control(formulation)))
  expect_equal(as.matrix(postPEV(fits[[1]], 2L)$Ci),
    as.matrix(postPEV(fits[[2]], 2L)$Ci), tolerance=1e-7)
  expect_equal(postPEV(fits[[1]])$uPevList, postPEV(fits[[2]])$uPevList,
    tolerance=1e-7)
  expect_equal(predict(fits[[1]], D="group")$vcov,
    predict(fits[[2]], D="group")$vcov, tolerance=1e-7)
  expect_equal(postVarU(fits[[1]], 2L)$VarU, postVarU(fits[[2]], 2L)$VarU,
    tolerance=1e-7)
})

test_that("direct observed classify predictions preserve retained-cell weighting", {
  data <- .direct_inference_data()
  data <- data[!(data$group == "a" & data$visit == "4"), ]
  fits <- lapply(c("direct", "henderson"), function(formulation){
    expect_warning(fit <- mmes(y~visit+visit:group, data=data,
      control=.direct_inference_control(formulation,
        representation=list(fixedRank="drop"))), "rank reduction")
    fit
  })
  direct <- predict(fits[[1]], D="visit:group", observed=TRUE, sed=TRUE)
  henderson <- predict(fits[[2]], D="visit:group", observed=TRUE, sed=TRUE)
  expect_equal(direct$D, henderson$D)
  expect_equal(direct$pvals, henderson$pvals, tolerance=1e-7)
  expect_equal(direct$vcov, henderson$vcov, tolerance=1e-7)
  expect_equal(direct$observed, henderson$observed)
})

test_that("direct mixed scalar PQL PEV respects final trait-specific dispersion", {
  data <- .direct_inference_data()
  data$count <- rpois(nrow(data), exp(0.5+0.1*data$x))
  data$positive <- rgamma(nrow(data), shape=5, scale=exp(0.5+0.1*data$x)/5)
  data <- rbind(transform(data, trait="count", response=count),
    transform(data, trait="positive", response=positive))
  data$trait <- factor(data$trait, levels=c("count", "positive"))
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(response~trait+trait:x-1,
      random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
      rcov=~vsm(dsm(trait), ism(units)), data=data,
      family=familym("trait", count=poisson(), positive=Gamma("log")),
      control=.direct_inference_control(formulation,
        task=list(pql=list(maxit=40, tol=1e-9)))))
  expect_true(fits[[1]]$pqlConverged)
  expect_equal(as.matrix(postPEV(fits[[1]], 2L)$Ci),
    as.matrix(postPEV(fits[[2]], 2L)$Ci), tolerance=2e-6)
  expect_equal(predict(fits[[1]], D="group", VarU=FALSE)$vcov,
    predict(fits[[2]], D="group", VarU=FALSE)$vcov, tolerance=2e-6)
})

test_that("direct fixed-effect comparisons support small-sample degrees of freedom", {
  data <- .direct_inference_data()
  fits <- lapply(c("direct", "henderson"), function(formulation)
    mmes(y~visit, random=~group, data=data,
      control=.direct_inference_control(formulation)))
  for(df in c("satterthwaite", "kr")){
    direct <- predict(fits[[1]], D="visit", pairwise=TRUE, df=df)
    henderson <- predict(fits[[2]], D="visit", pairwise=TRUE, df=df)
    expect_true(all(is.finite(direct$pairwise$df)))
    expect_equal(direct$pairwise, henderson$pairwise, tolerance=2e-4)
  }
})

test_that("deferred direct PEV is not limited to 5000 random coefficients", {
  data <- .direct_inference_data()[1:12, ]
  labels <- c(letters[1:6], paste0("extra", seq_len(5095L)))
  A <- Matrix::Diagonal(length(labels))
  dimnames(A) <- list(labels, labels)
  attr(A, "inverse") <- TRUE
  fit <- mmes(y~x, random=~vsm(ism(group), Gu=A, sigma2=0.7, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=0.4, fixedSigma2=TRUE), data=data,
    control=.direct_inference_control(output=list(panelSize=16L),
      execution=list(denseMemoryMB=50)))
  result <- postPEV(fit)
  expect_length(result$uPevList[[1L]], 5101L)
  expect_s4_class(result$Ci, "diagonalMatrix")
  expect_equal(as.numeric(result$uPevList[[1L]][-(1:6), 1L]), rep(0.7, 5095L),
    tolerance=1e-10)
  expect_error(postPEV(fit, 2L), "Full direct coefficient covariance exceeds")
})
