.direct_diallel <- function(blendWeight=0, maxIters=30L, setup=FALSE,
                            stepWeight=NULL, likelihoodTol=1e-4,
                            parameterTol=1e-4){
  data(DT_mohring, envir=environment())
  dataset <- add.diallel.vars(DT_mohring, par1="Par1", par2="Par2")
  mmes(Ftime~1, random=~Block+vsm(ism(overlay(Par1, Par2)))+cross.id+
    vsm(dsm(cross.type), ism(cross.id)), data=dataset,
    control=mmes_control(formulation="direct",
      update=list(blendWeight=blendWeight),
      proposal=list(stepWeight=stepWeight),
      convergence=list(maxIters=maxIters, likelihoodTol=likelihoodTol,
        parameterTol=parameterTol),
      output=list(setup=setup),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
}

.direct_scalar_likelihood <- function(setup, variances, reml=TRUE){
  y <- as.numeric(setup$yvar)
  X <- as.matrix(setup$X)
  kernels <- c(lapply(setup$Z, function(z) as.matrix(tcrossprod(z))),
    list(diag(length(y))))
  V <- Reduce(`+`, Map(`*`, kernels, variances))
  Vi <- chol2inv(chol(V))
  VX <- Vi %*% X
  Q <- crossprod(X, VX)
  P <- Vi - VX %*% solve(Q, t(VX))
  -0.5 * (as.numeric(determinant(V, logarithm=TRUE)$modulus) +
    if(reml) as.numeric(determinant(Q, logarithm=TRUE)$modulus) else 0) -
    0.5 * (sum(y * (P %*% y)) + (length(y)-if(reml) ncol(X) else 0) * log(2*pi))
}

test_that("direct diallel fitting reaches the variance boundary without freezing other terms", {
  fit <- .direct_diallel()
  expect_true(fit$convergence)
  expect_lt(max(abs(fit$covParNative$estimate -
    c(1e-6, 1276.73089, 1110.99090, 66.02295, 418.47949))), 0.15)
  expect_equal(fit$covParNative$estimate[1L], 1e-6, tolerance=1e-10)
  expect_true(all(diff(as.numeric(fit$llik)) > 0))
  expect_identical(fit$engineDiagnostics$direct$stopReason, "converged")
  expect_lte(fit$engineDiagnostics$direct$finalProjectedScoreNorm,
    fit$engineDiagnostics$direct$scoreTolerance)
  expect_true(all(is.finite(fit$theta_se)))
  expect_gt(fit$covParNative$StdError[1L], 1)
  expect_lt(fit$covParNative$StdError[1L], 100)

  setup <- .direct_diallel(setup=TRUE)
  variances <- fit$covParNative$estimate
  expect_equal(tail(as.numeric(fit$llik), 1L),
    .direct_scalar_likelihood(setup, variances), tolerance=1e-8)
  expect_equal(as.numeric(tail(t(fit$monitor), 1L)), variances, tolerance=1e-10)
  expect_equal(fit$AIC, -2*tail(as.numeric(fit$llik), 1L)+2*ncol(setup$X))
})

test_that("direct stopping never treats tiny updates or iteration limits as convergence", {
  limited <- .direct_diallel(maxIters=1L)
  setup <- .direct_diallel(setup=TRUE)
  initial <- vapply(setup$covStruct, function(cs) exp(cs$par[1L]), numeric(1))
  expect_false(limited$convergence)
  expect_equal(limited$covParNative$estimate, unname(initial), tolerance=1e-10)
  expect_equal(as.numeric(limited$llik), .direct_scalar_likelihood(setup, initial),
    tolerance=1e-8)
  expect_identical(limited$engineDiagnostics$direct$stopReason, "iteration-limit")

  tiny <- .direct_diallel(stepWeight=1e-7, maxIters=4L)
  expect_false(tiny$convergence)
  expect_lt(max(diff(as.numeric(tiny$llik))), 1e-4)
  expect_gt(tiny$engineDiagnostics$direct$finalProjectedScoreNorm,
    tiny$engineDiagnostics$direct$scoreTolerance)

  expect_warning(stalled <- .direct_diallel(stepWeight=1e-20),
    "Direct-engine optimization stalled.*convergence was not certified")
  expect_false(stalled$convergence)
  expect_identical(stalled$engineDiagnostics$direct$stopReason, "line-search-stalled")
  expect_equal(stalled$covPar, limited$covPar)
  expect_equal(stalled$llik, limited$llik)
})

test_that("direct likelihood line search rejects decreasing covariance updates", {
  fit <- .direct_diallel(stepWeight=100, maxIters=100L,
    likelihoodTol=1e-8, parameterTol=1e-8)
  expect_true(fit$convergence)
  expect_true(all(diff(as.numeric(fit$llik)) > 0))
  expect_lt(max(abs(fit$covParNative$estimate[-1L] -
    c(1276.7585, 1111.0853, 66.02665, 418.46485))), 0.02)

  dataset <- data.frame(y=sin(seq_len(24L)))
  backtracked <- mmes(y~1, rcov=~vsm(ism(units), sigma2=0.4), data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", proposal=list(stepWeight=100), convergence=list(maxIters=100L, likelihoodTol=1e-8, parameterTol=1e-8), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(backtracked$convergence)
  expect_gt(backtracked$engineDiagnostics$direct$lineSearchHalvings, 0L)
  expect_true(all(diff(as.numeric(backtracked$llik)) > 0))
  expect_equal(as.numeric(backtracked$covPar[[1L]]), var(dataset$y), tolerance=5e-5)
})

test_that("direct active-set updates retain fixed and tied covariance coordinates", {
  data(DT_example)
  fixed <- mmes(Yield~Env, random=~vsm(ism(Name), sigma2=2, fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=1, fixedSigma2=TRUE), data=DT_example, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(fixed$convergence)
  expect_equal(unname(unlist(fixed$covPar)), c(2, 1))
  expect_equal(fixed$theta_se, matrix(0, 2L, 2L))

  dataset <- DT_example
  dataset$Name2 <- dataset$Name
  tied <- mmes(Yield~Env, random=~Name+Name2, data=dataset, vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)), control=mmes_control(formulation=if(FALSE) "henderson" else "direct", convergence=list(maxIters=100L, likelihoodTol=1e-8, parameterTol=1e-8), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(tied$convergence)
  expect_equal(tied$covPar[[1L]], 2*tied$covPar[[2L]], tolerance=1e-8, ignore_attr=TRUE)
  expect_equal(tied$monitor[1L, ], 2*tied$monitor[2L, ], tolerance=1e-8)
})
