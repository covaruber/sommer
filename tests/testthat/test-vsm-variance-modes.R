test_that("correlation variance modes create consistent covariance descriptors", {
  levels3 <- factor(rep(letters[1:3], each=2))
  levels4 <- factor(rep(1:4, each=2))
  levels5 <- factor(rep(1:5, each=2))

  structures <- list(
    csm(levels3, variance="heterogeneous"),
    ar1m(levels4, variance="heterogeneous"),
    ar2m(levels4, variance="heterogeneous"),
    ar3m(levels5, variance="heterogeneous")
  )

  for(structure in structures){
    descriptor <- structure$covFactor
    expect_true(length(descriptor$par) == length(descriptor$free))
    expect_true(length(descriptor$par) == length(descriptor$par_names))
  }

  for(structure in structures[2:4]){
    descriptor <- structure$covFactor
    covariance <- descriptor$evaluator$fun(descriptor$par)
    expect_true(isSymmetric(covariance))
    expect_true(all(is.finite(covariance)))
    expect_gt(min(eigen(covariance, symmetric=TRUE, only.values=TRUE)$values), 0)
  }
})

test_that("native covariance reporting returns model-scale parameters", {
  DT_example <- get("DT_example", envir=asNamespace("sommer"))
  envFixed <- rep(FALSE, nlevels(DT_example$Env) - 1L)
  envFixed[1L] <- TRUE

  fitDiag <- mmes(
    Yield ~ Env,
    random=~vsm(dsm(Env, fixed=envFixed), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  diagNative <- covparams_mmes(fitDiag, 1L)
  expect_equal(nrow(diagNative), nlevels(DT_example$Env))
  expect_true(all(grepl("^variance\\[", diagNative$parameter)))
  expect_false(any(grepl("ratio", diagNative$parameter)))
  expect_equal(diagNative$estimate, diag(fitDiag$theta[[1]]), tolerance=1e-8)
  expect_equal(fitDiag$covParNative, covparams_mmes(fitDiag))

  diagSE <- covparams_mmes_se(fitDiag, 1L)
  expect_equal(diagSE$estimate, diagNative$estimate)
  expect_true(all(is.finite(diagSE$StdError)))
  expect_equal(fitDiag$covParNativeSE, covparams_mmes_se(fitDiag))

  reported <- as.numeric(fitDiag$covPar[[1]])
  scale <- reported[1L]
  ratios <- c(1, reported[-1L])
  localCovariance <- fitDiag$theta_se[seq_along(reported), seq_along(reported),
                                     drop=FALSE]
  free <- fitDiag$covStruct[[1]]$free
  localCovariance[!free,] <- 0
  localCovariance[,!free] <- 0
  analyticJacobian <- matrix(0, nlevels(DT_example$Env), length(reported))
  analyticJacobian[,1L] <- ratios
  if(length(reported) > 1L){
    analyticJacobian[cbind(2:length(ratios), 2:length(reported))] <- scale
  }
  analyticSE <- sqrt(pmax(diag(analyticJacobian %*% localCovariance %*%
                                 t(analyticJacobian)), 0))
  expect_equal(diagSE$StdError, analyticSE, tolerance=1e-5)

  fitUs <- mmes(
    Yield ~ Env,
    random=~vsm(usm(Env), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  usNative <- covparams_mmes(fitUs, 1L)
  expected <- fitUs$theta[[1]][lower.tri(fitUs$theta[[1]], diag=TRUE)]
  expect_equal(usNative$estimate, expected, tolerance=1e-8)
  expect_false(any(grepl("chol", usNative$parameter)))
  usSE <- covparams_mmes_se(fitUs, 1L)
  expect_true(all(is.finite(usSE$StdError)))
  expect_equal(usSE$parameter, usNative$parameter)

  fitCorrelation <- mmes(
    Yield ~ Env,
    random=~vsm(csm(Env, rho=0.2, fixed=TRUE), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  correlationSE <- covparams_mmes_se(fitCorrelation, 1L)
  expect_equal(correlationSE$StdError[correlationSE$parameter == "rho"], 0)
  expect_true(is.na(correlationSE$Zratio[correlationSE$parameter == "rho"]))
})

test_that("ownm native reporting callback is dispatched generically", {
  group <- factor(rep(c("a", "b"), each=4))
  id <- factor(rep(seq_len(4), 2))
  y <- seq_along(group)
  covariance <- function(par){
    rho <- tanh(par[1])
    matrix(c(1, rho, rho, 1), 2, 2)
  }
  reporter <- function(scale, par, factor, absorb_scale=TRUE){
    c(variance=scale, correlation=tanh(par[1]))
  }

  fit <- mmes(
    y ~ 1,
    random=~vsm(ownm(group, fun=covariance, par=0.1,
                     par_names="rho_work", native_report=reporter), ism(id)),
    rcov=~units,
    data=data.frame(y, group, id),
    verbose=FALSE,
    nIters=2
  )
  native <- covparams_mmes(fit, 1L)
  expect_equal(native$parameter, c("variance", "correlation"))
  expect_true(all(is.finite(native$estimate)))
  nativeSE <- covparams_mmes_se(fit, 1L)
  expect_true(all(is.finite(nativeSE$StdError)))
})
