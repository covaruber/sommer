test_that("profiling is off by default and does not change the fit", {
  set.seed(7L)
  data <- expand.grid(id = factor(seq_len(6L)), level = factor(seq_len(4L)), rep = seq_len(2L))
  data$response <- rnorm(nrow(data)) + as.integer(data$level) / 4
  fitWith <- function(mode){
    mmes(response ~ level, random = ~vsm(ar1m(level), ism(id)), rcov = ~units,
         data = data, nIters = 2L, solver = "cholmod", engine = mmes_control(solver = "cholmod",
           experimental = list(profile = mode)),
         verbose = FALSE, dateWarning = FALSE)
  }
  off <- fitWith("off")
  on <- fitWith("force")
  expect_null(off$engineDiagnostics$profile)
  expect_equal(on$engineDiagnostics$control$profile, "force")
  for(component in c("llik", "theta", "bu", "InfMat")){
    expect_identical(on[[component]], off[[component]])
  }
})

test_that("phase profile reports consistent seconds, iterations and counters", {
  set.seed(3L)
  data <- expand.grid(id = factor(seq_len(6L)), level = factor(seq_len(4L)), rep = seq_len(2L))
  data$response <- rnorm(nrow(data)) + as.integer(data$level) / 4
  for(solver in c("ldlt", "cholmod")){
    fit <- mmes(response ~ level, random = ~vsm(ar2m(level), ism(id)), rcov = ~units,
      data = data, nIters = 2L, solver = solver, verbose = FALSE, dateWarning = FALSE,
      engine = mmes_control(solver = solver, experimental = list(profile = "force")))
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
  }
})

test_that("exact trace panels preserve AR2 fits, scores and uncertainty", {
  fitWith <- function(workload, mode, threads, inverseMode, reml = TRUE, iterations = 2L, budget = 1000,
                      policy = "speed"){
    mmes(workload$fixed, random = workload$random, rcov = workload$rcov, data = workload$data,
         nIters = iterations, nThreads = threads, computeCi = inverseMode, getPEV = inverseMode > 0L,
         REML = reml, blendWeight = 0, verbose = FALSE, dateWarning = FALSE,
         engine = mmes_control(solver = "ldlt", denseMemoryMB = budget, memoryPolicy = policy,
           experimental = list(profile = "force", blockSchur = "off", traceBatch = mode)))
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
      expect_gt(counters$traceBatchSolves, 0)
      expect_lte(counters$traceBatchSolves, counters$traceFallbackSolves)
      expect_gt(counters$traceBatchWorkspaceBytes, 0)
      expect_equal(candidate$engineDiagnostics$control$traceBatch, "force")
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
  expect_error(fitWith(workload, "force", 1L, 0L, budget = 0.001), "Batched exact trace solves needs about")
  limited <- fitWith(workload, "force", 1L, 0L, budget = 0.001, policy = "limit")
  reference <- fitWith(workload, "off", 1L, 0L, budget = 0.001, policy = "limit")
  for(field in fields){ expect_equal(limited[[field]], reference[[field]], tolerance = 1e-10) }
  expect_equal(limited$engineDiagnostics$profile$counters$traceBatchSolves, 0)
})
