test_that("profiling is off by default and does not change the fit", {
  set.seed(7L)
  data <- expand.grid(id = factor(seq_len(6L)), level = factor(seq_len(4L)), rep = seq_len(2L))
  data$response <- rnorm(nrow(data)) + as.integer(data$level) / 4
  fitWith <- function(mode){
    mmes(response ~ level, random=~vsm(ar1m(level), ism(id)), rcov=~units, data=data,
      control=mmes_control(solve=list(method="cholmod"),
        convergence=list(maxIters=2L), diagnostics=list(profile=mode,
          verbose=FALSE, dateWarning=FALSE)))
  }
  off <- fitWith("off")
  on <- fitWith("force")
  expect_null(off$engineDiagnostics$profile)
  expect_equal(on$control$diagnostics$profile, "force")
  for(component in c("llik", "theta", "bu", "InfMat")){
    expect_identical(on[[component]], off[[component]])
  }
})

test_that("phase profile reports consistent seconds, iterations and counters", {
  set.seed(3L)
  data <- expand.grid(id = factor(seq_len(6L)), level = factor(seq_len(4L)), rep = seq_len(2L))
  data$response <- rnorm(nrow(data)) + as.integer(data$level) / 4
  for(solver in c("ldlt", "cholmod")){
    fit <- mmes(response ~ level, random=~vsm(ar2m(level), ism(id)), rcov=~units, data=data,
      control=mmes_control(solve=list(method=solver), convergence=list(maxIters=2L),
        diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)))
    p <- fit$engineDiagnostics$profile
    expect_identical(names(p$seconds), sommer:::.mme_phase_names)
    expect_true(all(p$seconds >= 0))
    expect_equal(sum(p$seconds) / sum(p$seconds), 1)
    byIteration <- p$secondsByIteration
    expect_equal(nrow(byIteration), p$counters$iterations)
    expect_equal(unname(colSums(byIteration)[c("setup", "final")]), c(0, 0))
    expect_equal(unname(colSums(byIteration)),
                 unname(p$seconds * as.numeric(!names(p$seconds) %in% c("setup", "final"))),
                 tolerance = 1e-9)
    expect_gt(p$counters$nnzC, 0)
    expect_gt(p$counters$freeCovarianceParameters, 0)
    profiledSeconds <- c(
      p$counters$traceDerivativePreparationSeconds,
      p$counters$traceRandomPrecisionDerivativeSeconds,
      p$counters$traceFactorCovarianceDerivativeSeconds,
      p$counters$traceFactorPrecisionDerivativeSeconds,
      p$counters$traceRandomPrecisionOtherSeconds,
      p$counters$traceResidualPrecisionDerivativeSeconds,
      p$counters$traceGridDerivativeBasisSeconds,
      p$counters$traceResidualDerivativeBasisSeconds,
      p$counters$traceDerivativeOtherSeconds,
      p$counters$traceSelectedInverseSeconds,
      p$counters$sensitivitySolveSeconds,
      p$counters$sensitivitySolveReachEntries,
      p$counters$sensitivitySolveReachPossibleEntries,
      p$counters$scoreTraceFallbackSeconds,
      p$counters$aiAssemblyWallSeconds,
      p$counters$aiRandomFirstTermSeconds,
      p$counters$aiRandomCurvatureSeconds,
      p$counters$aiResidualBlockSeconds
    )
    expect_true(all(is.finite(profiledSeconds)))
    expect_true(all(profiledSeconds >= 0))
    expect_gt(p$counters$sensitivityRHSNonzeros, 0)
    expect_gt(p$counters$sensitivitySolveSeconds, 0)
    expect_lte(p$counters$sensitivitySolveSeconds, p$seconds[["sensitivities"]] + 1e-9)
    expect_lte(p$counters$sensitivitySolveReachEntries,
               p$counters$sensitivitySolveReachPossibleEntries)
    expect_equal(
      p$counters$selectedInverseLookupRequests,
      p$counters$selectedInverseLookupHits + p$counters$selectedInverseLookupMisses
    )
    expect_lte(p$counters$selectedInverseLookupHits, p$counters$selectedInverseLookupRequests)
    expect_lte(p$counters$selectedInverseLookupMisses, p$counters$selectedInverseLookupRequests)
    expect_lte(p$counters$selectedInverseDemandUniqueEntries,
               p$counters$selectedInverseLookupHits)
    expect_true(p$counters$selectedInverseDemandTrackingComplete)
    expect_gte(p$counters$selectedInverseDemandClosureEntries,
               p$counters$selectedInverseDemandUniqueEntries)
    expect_gt(p$counters$traceDerivativePreparationSeconds, 0)
    expect_gt(p$counters$traceSelectedInverseSeconds, 0)
    expect_gt(p$counters$aiAssemblyWallSeconds, 0)
    expect_equal(
      sum(unlist(p$counters[c(
        "traceFactorCovarianceDerivativeSeconds",
        "traceFactorPrecisionDerivativeSeconds",
        "traceRandomPrecisionOtherSeconds",
        "traceResidualPrecisionDerivativeSeconds",
        "traceGridDerivativeBasisSeconds",
        "traceResidualDerivativeBasisSeconds",
        "traceDerivativeOtherSeconds"
      )])),
      p$counters$traceDerivativePreparationSeconds,
      tolerance = 1e-2
    )
    expect_equal(
      p$counters$traceRandomPrecisionDerivativeSeconds,
      p$counters$traceFactorCovarianceDerivativeSeconds +
        p$counters$traceFactorPrecisionDerivativeSeconds +
        p$counters$traceRandomPrecisionOtherSeconds,
      tolerance = 1e-2
    )
    expect_lte(p$counters$traceDerivativePreparationSeconds, p$seconds[["tracePrep"]] + 1e-9)
    expect_lte(p$counters$traceSelectedInverseSeconds, p$seconds[["tracePrep"]] + 1e-9)
    expect_lte(p$counters$scoreTraceFallbackSeconds, p$seconds[["scoreTraces"]] + 1e-9)
    expect_lte(p$counters$aiAssemblyWallSeconds, p$seconds[["aiAssembly"]] + 1e-9)
    expect_lte(p$counters$aiRandomFirstTermSeconds, p$seconds[["aiAssembly"]] + 1e-9)
    expect_lte(p$counters$aiRandomCurvatureSeconds, p$seconds[["aiAssembly"]] + 1e-9)
    expect_lte(p$counters$aiResidualBlockSeconds, p$seconds[["aiAssembly"]] + 1e-9)
  }
})

test_that("exact trace panels preserve AR2 fits, scores and uncertainty", {
  skip_if(is.null(sommer:::.mmes_blas_controller()),
    "No supported BLAS controller available")
  fitWith <- function(workload, mode, threads, inverseMode, reml = TRUE, iterations = 2L, budget = 1000,
                      policy = "speed"){
    mmes(workload$fixed, random=workload$random, rcov=workload$rcov, data=workload$data,
      control=mmes_control(execution=list(blasThreads=threads,
          denseMemoryMB=budget, memoryPolicy=policy),
        solve=list(blockSchur="off"), score=list(traceBatch=mode),
        convergence=list(maxIters=iterations),
        output=list(inverse=inverseMode),
        task=list(objective=if(reml) "reml" else "ml"),
        update=list(blendWeight=0),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se", "uPevList")
  for(relationship in c("none", "sparse", "dense")){
    workload <- generate_mme_workload(structure = "ar2", n_levels = 9L, n_ind = 4L,
                                     relationship = relationship, missing_rate = 0.1, seed = 303L)
    for(inverseMode in 0:2){
      reml <- inverseMode != 1L
      threads <- if(inverseMode == 0L) 1L else 4L
      reference <- fitWith(workload, "off", 1L, inverseMode, reml)
      candidate <- fitWith(workload, "force", threads, inverseMode, reml)
      for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-8) }
      if(inverseMode == 2L){
        expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance = 1e-8)
        product <- as.matrix(candidate$C %*% candidate$Ci)
        expect_lt(max(abs(product / mean(diag(product)) - diag(nrow(product)))), 1e-8)
      }
      counters <- candidate$engineDiagnostics$profile$counters
      expect_equal(counters$traceFallbackSolves, reference$engineDiagnostics$profile$counters$traceFallbackSolves)
      expect_true(all(counters$traceBatchSolves <= counters$traceFallbackSolves))
      expect_true(all(counters$traceBatchWorkspaceBytes > 0))
      expect_equal(candidate$control$score$traceBatch, "force")
    }
  }
  workload <- generate_mme_workload(structure = "ar2", n_levels = 12L, n_ind = 4L, seed = 7L)
  for(pacf in list(c(0, 0), c(0.94, 0.85))){
    shape <- ar2m(workload$data$level, pacf = pacf, variance = "heterogeneous",
                  fixed = c(TRUE, TRUE, rep(FALSE, 11L)),
                  values = seq(0.5, 1.5, length.out = 12L))
    workload$random <- ~vsm(shape, ism(id))
    reference <- fitWith(workload, "off", 1L, 0L)
    candidate <- fitWith(workload, "force", 4L, 0L)
    for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7) }
  }
  workload <- generate_mme_workload(structure = "ar2", n_levels = 9L, n_ind = 4L, seed = 303L)
  reference <- fitWith(workload, "off", 1L, 0L, iterations = 8L)
  candidate <- fitWith(workload, "force", 4L, 0L, iterations = 8L)
  for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7) }
  workload$data$batch <- factor(workload$data$rep)
  workload$random <- ~vsm(ar2m(level), ism(id)) + batch
  reference <- fitWith(workload, "off", 1L, 1L)
  candidate <- fitWith(workload, "force", 4L, 1L)
  for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-8) }
  expect_error(fitWith(workload, "force", 1L, 0L, budget = 0.001), "denseMemoryMB is 0")
  limited <- fitWith(workload, "force", 1L, 0L, budget = 0.001, policy = "limit")
  reference <- fitWith(workload, "off", 1L, 0L, budget = 0.001, policy = "limit")
  for(field in fields){ expect_equal(limited[[field]], reference[[field]], tolerance = 1e-10) }
  expect_equal(sum(limited$engineDiagnostics$profile$counters$traceBatchSolves), 0)
})
