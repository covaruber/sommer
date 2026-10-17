.recovery_fit <- function(...){
  data(DT_example, envir=environment())
  mmes(Yield~Env, random=~Name, data=DT_example,
    verbose=FALSE, dateWarning=FALSE, ...)
}

test_that("healthy default AI skips all EM and PX work and reuses accepted fits", {
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
    for(method in c("em", "pxem")){
      fit <- .recovery_fit(solver=solver, blend=method, nIters=30L, REML=reml)
      diagnostics <- fit$engineDiagnostics
      expect_equal(diagnostics$optimizer, "AI-first-recovery")
      expect_equal(diagnostics$recovery$likelihoodMode, if(reml) "REML" else "ML")
      expect_true(fit$convergence)
      expect_equal(diagnostics$recovery$stopReason, "converged")
      expect_equal(diagnostics$recovery$steps, 0)
      expect_lte(diagnostics$recovery$finalProjectedScoreNorm,
        diagnostics$recovery$scoreTolerance)
      expect_equal(diagnostics$blend$emEvaluations, 0)
      expect_equal(diagnostics$blend$qEvaluations, 0)
      expect_equal(diagnostics$blend$momentPanelSolves, 0)
      expect_true(all(diagnostics$blend$weights == 0))
      expect_gt(diagnostics$blend$acceptedAIFitReuses, 0)
      expect_equal(diagnostics$coefficientEvaluations, ncol(fit$llik))
      expect_null(fit$inputArgs$blendWeight)
    }
  }
  adaptive <- .recovery_fit(nIters=3L)
  explicit <- .recovery_fit(nIters=3L, blendWeight=0, stepWeight=1)
  for(field in c("llik", "theta", "bu", "InfMat")){
    expect_equal(adaptive[[field]], explicit[[field]], tolerance=1e-10)
  }
  expect_null(explicit$engineDiagnostics$recovery)
})

test_that("negligible gains with a non-small score trigger one recovery then retry AI", {
  for(reml in c(TRUE, FALSE)) for(method in c("em", "pxem")){
    fit <- .recovery_fit(nIters=6L, stepWeight=1e-7, blend=method, REML=reml)
    diagnostics <- fit$engineDiagnostics
    expect_false(fit$convergence)
    expect_equal(diagnostics$recovery$iterations, 4L)
    expect_equal(diagnostics$recovery$reasons, "negligible-gains")
    expect_equal(diagnostics$recovery$steps, 1L)
    expect_equal(diagnostics$recovery$acceptedSteps, 1L)
    expect_equal(as.numeric(diagnostics$blend$weights), c(0, 0, 0, 1, 0, 0))
    expect_equal(diagnostics$blend$emEvaluations, 1L)
    expect_gt(diagnostics$recovery$finalProjectedScoreNorm,
      diagnostics$recovery$scoreTolerance)
    gains <- diff(as.numeric(fit$llik))
    expect_true(all(gains >= -1e-10))
    expect_gt(gains[4L], 1e-10)
    expect_equal(diagnostics$recovery$stopReason, "iteration-limit")
  }
})

test_that("singular AI information invokes exact EM recovery", {
  data(DT_example)
  dataset <- DT_example
  dataset$Name2 <- dataset$Name
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
    fit <- mmes(Yield~Env, random=~Name+Name2, data=dataset, solver=solver,
      nIters=4L, REML=reml, verbose=FALSE, dateWarning=FALSE)
    diagnostics <- fit$engineDiagnostics
    expect_equal(diagnostics$recovery$reasons[1L], "information-solve")
    expect_gt(diagnostics$blend$emEvaluations, 0)
    expect_true(all(diff(as.numeric(fit$llik)) > 1e-10))
    expect_true(all(is.finite(fit$InfMat)))
    expect_false(fit$convergence)
  }
})

test_that("failed recovery returns the accepted point and warns without false convergence", {
  for(reml in c(TRUE, FALSE)){
  expect_warning(fit <- .recovery_fit(nIters=8L, stepWeight=1e-12, REML=reml),
    "AI-first optimization stalled")
  diagnostics <- fit$engineDiagnostics
  expect_false(fit$convergence)
  expect_equal(diagnostics$recovery$stopReason, "recovery-stalled")
  expect_equal(diagnostics$recovery$iterations, 4L)
  expect_equal(diagnostics$recovery$acceptedSteps, 0L)
  expect_true(diagnostics$blend$lineSearchStalled)
  expect_true(all(diff(as.numeric(fit$llik)) >= -1e-10))
  expect_equal(ncol(fit$llik), 4L)
  expect_true(all(is.finite(fit$bu)))
  expect_true(all(is.finite(fit$dLu)))
  }
})

test_that("a boundary solution does not itself trigger recovery", {
  dataset <- expand.grid(id=factor(seq_len(10L)), replicate=seq_len(4L))
  dataset$y <- rep(c(-1, 1, -1, 1), each=10L)
  for(reml in c(TRUE, FALSE)){
  fit <- mmes(y~1, random=~vsm(ism(id), sigma2=1e-6), data=dataset, nIters=80L,
    tolParInv=1e-6, REML=reml, verbose=FALSE, dateWarning=FALSE)
  expect_true(fit$convergence)
  expect_equal(fit$engineDiagnostics$recovery$steps, 0L)
  expect_lte(fit$theta[[1L]][1L], 1e-4)
  expect_lte(fit$engineDiagnostics$recovery$finalProjectedScoreNorm,
    fit$engineDiagnostics$recovery$scoreTolerance)
  }
})

test_that("ML recovery with no random coefficients uses deterministic residual moments", {
  data(DT_example)
  response <- DT_example$Yield
  expected <- sum(residuals(lm(Yield~Env, data=DT_example))^2) / length(response)
  for(solver in c("ldlt", "cholmod")) for(method in c("em", "pxem")){
    fit <- mmes(Yield~Env, data=DT_example, REML=FALSE, solver=solver,
      nIters=7L, blend=method, stepWeight=1e-7, verbose=FALSE, dateWarning=FALSE)
    diagnostics <- fit$engineDiagnostics
    expect_equal(diagnostics$recovery$acceptedSteps, 1L)
    expect_equal(diagnostics$recovery$reasons, "negligible-gains")
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-10))
    expect_equal(diagnostics$DsymbolicAnalyses, 0)
    em <- mmes(Yield~Env, data=DT_example, REML=FALSE, solver=solver,
      nIters=2L, blend=method, blendWeight=1, verbose=FALSE, dateWarning=FALSE)
    expect_equal(as.numeric(em$theta[[1L]]), expected, tolerance=1e-10)
    expect_equal(diagnostics$blend$momentPanelSolves, 0L)
  }
})
