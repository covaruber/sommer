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
  expect_equal(fitDiag$covParNative, covparams_mmes_se(fitDiag))

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

test_that("all covariance-factor families provide finite native reports", {
  x3 <- factor(rep(letters[1:3], each=2))
  x4 <- factor(rep(letters[1:4], each=2))
  coords <- matrix(c(0,0, 1,0, 0,1), ncol=2, byrow=TRUE)
  W <- matrix(c(0,1,1, 1,0,1, 1,1,0), 3, 3)
  dimnames(W) <- list(letters[1:3], letters[1:3])

  structures <- list(
    ism(x3), dsm(x3), usm(x3),
    ar1m(x3), ar1m(x3, variance="heterogeneous"),
    ar2m(x4), ar2m(x4, variance="heterogeneous"),
    ar3m(x4), ar3m(x4, variance="heterogeneous"),
    csm(x3), csm(x3, variance="heterogeneous"),
    mam(x3), corgm(x3), fam(x3, 1), antem(x3), rrcm(x3, 1),
    maternm(coords), toeplitzm(x3), sar(x3, W), car(x3, W)
  )

  for(structure in structures){
    factor <- structure$covFactor
    reporter <- factor$native_report$fun
    natural <- factor$par
    if(length(natural)){
      transform <- factor$report$transform
      for(k in seq_along(natural)){
        natural[k] <- switch(transform[k],
          identity=natural[k], exp=exp(natural[k]), tanh=tanh(natural[k]),
          bounded_logit={
            p <- stats::plogis(natural[k])
            factor$report$lower[k] +
              (factor$report$upper[k] - factor$report$lower[k]) * p
          })
      }
    }
    values <- reporter(scale=2, par=natural, factor=factor, absorb_scale=TRUE)
    expect_true(length(values) > 0L, info=factor$model)
    expect_true(all(is.finite(values)), info=factor$model)
    expect_true(all(nzchar(names(values))), info=factor$model)

    factor$par_start <- 2L
    factor$par_end <- length(natural) + 1L
    fitted <- structure(list(
      covStruct=stats::setNames(list(list(
        factors=list(factor), free=c(TRUE, factor$free)
      )), "term"),
      covPar=stats::setNames(list(c(2, natural)), "term"),
      theta_se=diag(length(natural) + 1L)
    ), class="mmes")
    publicValues <- covparams_mmes(fitted)
    publicSE <- covparams_mmes_se(fitted)
    expect_equal(publicValues$estimate, as.numeric(values), info=factor$model)
    expect_true(all(is.finite(publicSE$StdError)), info=factor$model)
  }
})

test_that("heterogeneous AR1 reports rho and environment variances", {
  DT_example <- get("DT_example", envir=asNamespace("sommer"))
  fit <- mmes(
    Yield ~ Env,
    random=~vsm(ar1m(Env, variance="heterogeneous"), ism(Name)),
    rcov=~vsm(dsm(Env), ism(units)),
    data=DT_example,
    verbose=FALSE,
    nIters=5
  )
  native <- covparams_mmes(fit, 1L)
  varianceRows <- grepl("^variance\\[", native$parameter)
  expect_equal(native$parameter[1L], "rho")
  expect_equal(native$estimate[varianceRows], diag(fit$theta[[1]]), tolerance=1e-8)
  expect_false(any(grepl("ratio|sigma2", native$parameter)))
  nativeSE <- covparams_mmes_se(fit, 1L)
  expect_true(all(is.finite(nativeSE$StdError)))
})

