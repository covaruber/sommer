.recovery_fit <- function(solver="auto", blend="em", nIters=30L, REML=TRUE,
                          stepWeight=NULL, blendWeight=NULL, tolParInv=1e-6){
  data(DT_example, envir=environment())
  mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(
    solve=list(method=solver),
    task=list(objective=if(REML) "reml" else "ml"),
    convergence=list(maxIters=nIters),
    execution=list(inverseTol=tolParInv),
    update=list(recovery=blend, blendWeight=blendWeight),
    proposal=list(stepWeight=stepWeight),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
}

test_that("regularized AI solves in trust units and preserves fixed and tied tangent directions", {
  helper <- .native_ai_recovery()
  information <- diag(c(1e-30, 2, 3, 4))
  score <- c(-1e-8, 3, 2, 1)
  tangent <- cbind(c(1, 0, 0, 0), c(0, 2, 1, 0))
  caps <- c(1, 0.5, 2, 1)
  candidate <- helper$regularized_ai_probe(information, score, tangent, caps)
  expect_true(candidate$reliable)
  expect_false(candidate$ordinaryReliable)
  expect_true(all(is.finite(candidate$direction)))
  expect_gt(sum(candidate$direction * score), 0)
  expect_lte(max(abs(candidate$direction) / caps), 1)
  expect_equal(candidate$direction[2L], 2 * candidate$direction[3L])
  expect_equal(candidate$direction[4L], 0)
  expect_gt(candidate$ridge, 0)
  rescaled <- helper$regularized_ai_probe(information * 100, score * 100, tangent, caps)
  expect_equal(rescaled$direction, candidate$direction, tolerance=1e-10)
  expect_equal(rescaled$ridge, candidate$ridge * 100, tolerance=1e-10)
  indefinite <- helper$regularized_ai_probe(diag(c(-1, 2)), c(1, 1), diag(2), c(1, 1))
  expect_true(indefinite$reliable)
  expect_gt(sum(indefinite$direction), 0)
  empty <- helper$regularized_ai_probe(information, score, matrix(0, 4, 0), caps)
  flat <- helper$regularized_ai_probe(matrix(0, 4, 4), score, diag(4), caps)
  expect_false(empty$reliable)
  expect_false(flat$reliable)
  healthy <- helper$regularized_ai_probe(diag(4), score, diag(4), rep(1, 4))
  expect_true(healthy$ordinaryReliable)
  threshold <- sqrt(.Machine$double.eps)
  for(ratio in c(0.5, 2)){
    probe <- helper$regularized_ai_probe(diag(c(1, ratio * threshold)),
      c(1, 1), diag(2), c(1, 1))
    expect_identical(probe$ordinaryReliable, ratio > 1)
  }
  units <- helper$regularized_ai_probe(diag(c(1e16, 1)), c(1, 1),
    diag(2), c(1e-8, 1))
  expect_true(units$ordinaryReliable)
})

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
      expect_equal(diagnostics$recovery$regularizedAttempts, 0L)
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

test_that("singular AI information tries regularized AI before computing EM moments", {
  data(DT_example)
  dataset <- DT_example
  dataset$Name2 <- dataset$Name
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
    fit <- mmes(Yield~Env, random=~Name+Name2, data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=4L), task=list(objective=if(reml) "reml" else "ml"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    diagnostics <- fit$engineDiagnostics
    expect_gt(diagnostics$recovery$regularizedAttempts, 0L)
    expect_gt(diagnostics$recovery$regularizedAcceptedSteps, 0L)
    expect_equal(diagnostics$recovery$steps, 0L)
    expect_equal(diagnostics$blend$emEvaluations, 0L)
    expect_equal(diagnostics$blend$qEvaluations, 0L)
    expect_equal(diagnostics$blend$momentPanelSolves, 0L)
    expect_true(all(diagnostics$recovery$regularizedRidges > 0))
    expect_true(all(diff(as.numeric(fit$llik)) > 1e-10))
    expect_true(all(is.finite(fit$InfMat)))
    if(fit$convergence){
      expect_lte(diagnostics$recovery$finalProjectedScoreNorm,
        diagnostics$recovery$scoreTolerance)
    }
  }
})

test_that("failed regularized line search still falls back to EM without accepting a negligible step", {
  data(DT_example)
  dataset <- DT_example
  dataset$Name2 <- dataset$Name
  for(reml in c(TRUE, FALSE)) for(method in c("em", "pxem")){
    expect_warning(fit <- mmes(Yield~Env, random=~Name+Name2, data=dataset, control=mmes_control(convergence=list(maxIters=4L), proposal=list(stepWeight=1e-12), task=list(objective=if(reml) "reml" else "ml"), update=list(recovery=method), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "AI-first optimization stalled")
    recovery <- fit$engineDiagnostics$recovery
    expect_gt(recovery$regularizedAttempts, 0L)
    expect_equal(recovery$regularizedAcceptedSteps, 0L)
    expect_equal(recovery$reasons, "information-solve")
    expect_equal(recovery$steps, 1L)
    expect_equal(recovery$acceptedSteps, 0L)
    expect_gt(fit$engineDiagnostics$blend$emEvaluations, 0L)
    expect_equal(recovery$stopReason, "recovery-stalled")
    expect_false(fit$convergence)
    initial <- mmes(Yield~Env, random=~Name+Name2, data=dataset, control=mmes_control(convergence=list(maxIters=1L), task=list(objective=if(reml) "reml" else "ml"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    for(field in c("theta", "bu", "llik")){
      expect_equal(fit[[field]], initial[[field]], tolerance=1e-9)
    }
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
  fit <- mmes(y~1, random=~vsm(ism(id), sigma2=1e-6), data=dataset, control=mmes_control(convergence=list(maxIters=80L), execution=list(inverseTol=1e-6), task=list(objective=if(reml) "reml" else "ml"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
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
    fit <- mmes(Yield~Env, data=DT_example, control=mmes_control(task=list(objective=if(FALSE) "reml" else "ml"), solve=list(method=solver), convergence=list(maxIters=7L), update=list(recovery=method), proposal=list(stepWeight=1e-7), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    diagnostics <- fit$engineDiagnostics
    expect_equal(diagnostics$recovery$acceptedSteps, 1L)
    expect_equal(diagnostics$recovery$reasons, "negligible-gains")
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-10))
    expect_equal(diagnostics$DsymbolicAnalyses, 0)
    em <- mmes(Yield~Env, data=DT_example, control=mmes_control(task=list(objective=if(FALSE) "reml" else "ml"), solve=list(method=solver), convergence=list(maxIters=2L), update=list(recovery=method, blendWeight=1), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(as.numeric(em$theta[[1L]]), expected, tolerance=1e-10)
    expect_equal(diagnostics$blend$momentPanelSolves, 0L)
  }
})
