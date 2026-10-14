test_that("profiling is off by default and does not change the fit", {
  wl <- generate_mme_workload(n_ind = 12L, n_levels = 4L, structure = "ar1", seed = 7L)
  fitWith <- function(mode){
    mmes(wl$fixed, random = wl$random, rcov = wl$rcov, data = wl$data, nIters = 3L,
         solver = "cholmod", engine = mmes_control(solver = "cholmod",
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
  wl <- generate_mme_workload(n_ind = 12L, n_levels = 4L, structure = "ar2", seed = 3L)
  for(solver in c("ldlt", "cholmod")){
    p <- profile_mme_workload(wl, mmes_control(solver = solver), nIters = 3L)
    expect_s3_class(p, "sommer_mme_profile")
    expect_identical(p$table$phase, sommer:::.mme_phase_names)
    expect_true(all(p$table$seconds >= 0))
    expect_equal(sum(p$table$share), 1)
    byIteration <- p$secondsByIteration
    expect_equal(nrow(byIteration), p$counters$iterations)
    expect_equal(unname(colSums(byIteration)[c("setup", "final")]), c(0, 0))
    expect_equal(unname(colSums(byIteration)),
                 p$table$seconds * (!p$table$phase %in% c("setup", "final")),
                 tolerance = 1e-9)
    expect_gt(p$counters$nnzC, 0)
    expect_gt(p$counters$freeCovarianceParameters, 0)
    expect_output(print(p), "Phase Profile")
  }
  expect_error(profile_mme_workload(list()), "mme_workload")
})

test_that("benchmark results carry phase columns only when profiling is on", {
  wl <- generate_mme_workload(n_ind = 10L, n_levels = 3L, structure = "ar1")
  bench <- run_mme_benchmark(wl, list(
    plain = mmes_control(solver = "ldlt"),
    profiled = mmes_control(solver = "ldlt", experimental = list(profile = "force"))),
    replications = 1L, nIters = 2L, verbose = FALSE)
  expect_true(is.na(bench$results$phase_assembly[1L]))
  expect_false(is.na(bench$results$phase_assembly[2L]))
  expect_true(is.na(bench$results$dominantPhase[1L]))
  expect_true(nzchar(bench$results$dominantPhase[2L]))
})

test_that("exact trace panels preserve AR2 fits, scores and uncertainty", {
  fitWith <- function(workload, mode, threads, inverseMode, reml = TRUE, iterations = 2L, budget = 1000,
                      policy = "speed"){
    mmes(workload$fixed, random = workload$random, rcov = workload$rcov, data = workload$data,
         nIters = iterations, nThreads = threads, computeCi = inverseMode, getPEV = inverseMode > 0L,
         REML = reml, verbose = FALSE, dateWarning = FALSE,
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
