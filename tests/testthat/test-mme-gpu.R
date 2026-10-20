gpuEngineArguments <- function(workload, solver = "cholmod", inverseMode = 0L,
                               reml = TRUE){
  setup <- mmes(workload$fixed, random=workload$random, rcov=workload$rcov,
    data=workload$data,
    control=mmes_control(solve=list(method=solver,
        blockSchur=if(solver == "cholmod") "force" else "off"),
      convergence=list(maxIters=2L),
      output=list(inverse=inverseMode, setup=TRUE),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE),
      task=list(objective=if(reml) "reml" else "ml")))
  list(X = setup$X, ZI = setup$Z, Zind = setup$Zind, AiI = setup$Ai,
    y0 = setup$yvar, H = setup$W, useH = setup$useH,
    residualBlockI = setup$residualBlock, residualIndexI = setup$residualIndex,
    nIters = setup$nIters, tolParConvLL = setup$tolParConvLL,
    tolParConvNorm = setup$tolParConvNorm, tolParInv = setup$tolParInv,
    covStructI = setup$covStruct, weightEmInf = setup$blendWeight,
    weightInf = setup$stepWeight, verbose = FALSE, computeCi = inverseMode,
    solver = solver, reml = reml, responsePrepared = setup$responsePrepared,
    preparedMean = setup$preparedMean, preparedSd = setup$preparedSd,
    preparedIntercept = setup$preparedIntercept)
}

test_that("unified solver preserves CPU results with GPU disabled", {
  for(structure in c("ar1", "ar2")){
    workload <- generate_mme_workload(n_ind = 6L, n_levels = 4L,
      structure = structure, relationship = "dense", seed = 37L)
    for(solver in c("ldlt", "cholmod", "pcg")){
      arguments <- gpuEngineArguments(workload, solver)
      reference <- do.call(sommer:::ai_mme_sp2, arguments)
      candidate <- do.call(sommer:::ai_mme_sp2, c(arguments, list(gpu = FALSE)))
      for(field in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se", "C")){
        expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-10)
      }
      expect_false(candidate$gpuDiagnostics$enabled)
      expect_equal(candidate$gpuDiagnostics$choleskyCalls, 0)
      expect_equal(candidate$gpuDiagnostics$triangularSolveCalls, 0)
      expect_equal(candidate$gpuDiagnostics$productCalls, 0)
    }
  }
})

test_that("GPU requests either use CUDA or report an explicit backend error", {
  workload <- generate_mme_workload(n_ind = 6L, n_levels = 4L,
    structure = "ar1", relationship = "dense", seed = 40L)
  result <- tryCatch(do.call(sommer:::ai_mme_sp2,
    c(gpuEngineArguments(workload), list(gpu = TRUE, gpuMinSize = 1L))),
    error = identity)
  if(inherits(result, "error")){
    expect_match(conditionMessage(result), "CUDA|cuBLAS|cuSOLVER|GPU")
  }else{
    expect_true(result$gpuDiagnostics$enabled)
    expect_gt(result$gpuDiagnostics$choleskyCalls, 0)
  }
})

test_that("CUDA kernels match CPU fits when GPU testing is enabled", {
  skip_if(Sys.getenv("SOMMER_TEST_CUDA") != "true",
    "Set SOMMER_TEST_CUDA=true on a CUDA build with an NVIDIA GPU")
  workload <- generate_mme_workload(n_ind = 8L, n_levels = 4L,
    structure = "ar1", relationship = "dense", seed = 39L)
  for(inverseMode in c(0L, 2L)){
    arguments <- gpuEngineArguments(workload, inverseMode = inverseMode)
    reference <- do.call(sommer:::ai_mme_sp2, arguments)
    candidate <- do.call(sommer:::ai_mme_sp2,
      c(arguments, list(gpu = TRUE, gpuMinSize = 1L)))
    for(field in c("llik", "theta", "bu", "InfMat", "Ci", "uPevList")){
      expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7)
    }
    expect_gt(candidate$gpuDiagnostics$choleskyCalls, 0)
    expect_gt(candidate$gpuDiagnostics$triangularSolveCalls, 0)
    expect_gt(candidate$gpuDiagnostics$productCalls, 0)
  }
})

test_that("GPU controls honor global defaults and explicit per-fit overrides", {
  old <- options(sommer.mme.gpu = TRUE, sommer.mme.gpuMinSize = 512L)
  on.exit(options(old), add = TRUE)
  expect_true(mmes_control()$execution$gpu)
  expect_identical(mmes_control()$execution$gpuMinSize, 512L)
  expect_false(mmes_control(execution=list(gpu=FALSE))$execution$gpu)
  expect_error(mmes_control(execution=list(gpu=NA)), "TRUE/FALSE")
  expect_error(mmes_control(execution=list(gpuMinSize=1.5)), "positive integer")
  expect_error(mmes_control(gpu=TRUE), "unused argument")
  options(sommer.mme.gpu = FALSE)
  workload <- generate_mme_workload(n_ind = 6L, n_levels = 4L, structure = "ar1")
  fit <- mmes(workload$fixed, random=workload$random, rcov=workload$rcov, data=workload$data,
    control=mmes_control(execution=list(gpu=FALSE),
      convergence=list(maxIters=2L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(fit$gpuDiagnostics$enabled)
  expect_identical(fit$gpuDiagnostics$minimumSize, 512L)
  expect_error(mmes(workload$fixed, data=workload$data, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", execution=list(gpu=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "execution\\$gpu")
  if(!sommer:::.mmes_gpu_status()$available){
    expect_error(mmes(workload$fixed, data=workload$data, control=mmes_control(execution=list(gpu=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "gpu=TRUE cannot run")
  }
})
