test_that("native PCG primitives pass strict standalone checks", {
  names <- c("test-mme-pcg", "test-mme-pcg-estimators", "test-mme-sufficient")
  for(name in names){
    .native_run(name, c("mme_pcg.h", "mme_pcg_estimators.h", "mme_sufficient.h"),
      eigen=TRUE, strict=TRUE)
  }
})

.pcg_module_values <- function(settings){
  mapping <- c(tol="solve.tol", maxIters="solve.maxIters",
    preconditioner="solve.preconditioner", nystromRank="solve.nystromRank",
    warmStart="solve.warmStart", recycleRank="solve.recycleRank",
    recycleMethod="solve.recycleMethod", recycleState="solve.recycleState",
    exportRecycle="solve.exportRecycle", rhsMethod="solve.rhsMethod",
    blockSize="solve.blockSize", residualRefresh="solve.residualRefresh",
    logDetScaling="logdet.scaling", lanczosSteps="logdet.lanczosSteps",
    traceProbes="score.traceProbes", traceEstimator="score.estimator",
    traceRank="score.rank", accuracy="convergence.accuracy",
    maxTraceProbes="convergence.maxTraceProbes",
    maxLanczosSteps="convergence.maxLanczosSteps",
    logDetTol="convergence.logDetTol", scoreTol="convergence.scoreTol",
    seed="execution.seed", compression="representation.compression")
  unknown <- setdiff(names(settings), names(mapping))
  if(length(unknown)) stop("Unknown PCG settings: ", paste(unknown, collapse=", "))
  modules <- list()
  for(name in names(settings)){
    path <- strsplit(unname(mapping[name]), ".", fixed=TRUE)[[1L]]
    modules[[path[1L]]][[path[2L]]] <- settings[[name]]
  }
  modules
}

.pcg_numerical_fit <- function(workload, pcg=list(), matrixFree="force", iterations=1L,
                               solver="pcg", ...){
  extra <- list(...)
  unknown <- setdiff(names(extra), "denseMemoryMB")
  if(length(unknown)) stop("Unknown MME controls: ", paste(unknown, collapse=", "))
  modules <- .pcg_module_values(pcg)
  merge <- function(module, defaults=list()){
    values <- modules[[module]]
    if(is.null(values)) values <- list()
    utils::modifyList(defaults, values, keep.null=TRUE)
  }
  execution <- modules$execution
  if(is.null(execution)) execution <- list()
  execution <- utils::modifyList(execution, extra, keep.null=TRUE)
  control <- mmes_control(
    solve=merge("solve", c(list(method=solver), if(solver == "pcg") list(tol=1e-10))),
    logdet=merge("logdet", if(solver == "pcg") list(lanczosSteps=40L) else list()),
    score=merge("score", if(solver == "pcg") list(traceProbes=8L) else list()),
    convergence=merge("convergence", list(maxIters=iterations)),
    representation=merge("representation", if(solver == "pcg") list(matrixFree=matrixFree) else list()),
    execution=execution,
    update=list(blendWeight=0),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  suppressMessages(mmes(workload$fixed, random=workload$random,
    rcov=workload$rcov, data=workload$data, control=control))
}

test_that("PCG operator preconditioners and block RHS agree with assembled solves", {
  cases <- expand.grid(structure=c("ar1", "us", "fa"),
    relationship=c("none", "sparse", "dense"), stringsAsFactors=FALSE)
  for(case in seq_len(nrow(cases))){
    relationship <- cases$relationship[case]
    workload <- generate_mme_workload(n_ind=8L, n_levels=3L,
      structure=cases$structure[case], relationship=relationship, seed=331L)
    reference <- .pcg_numerical_fit(workload, solver="ldlt")
    variants <- list(
      list(), list(preconditioner="block-jacobi"),
      list(preconditioner="nystrom", nystromRank=6L),
      list(rhsMethod="block", blockSize=4L, residualRefresh=3L),
      list(logDetScaling="diagonal"))
    for(pcg in variants){
      operator <- .pcg_numerical_fit(workload, pcg)
      assembled <- .pcg_numerical_fit(workload, pcg, matrixFree="off")
      expect_true(operator$pcgMatrixFree)
      expect_false(assembled$pcgMatrixFree)
      expect_equal(operator$bu, reference$bu, tolerance=1e-7)
      expect_equal(assembled$bu, reference$bu, tolerance=1e-7)
      expect_equal(operator$dLu, assembled$dLu, tolerance=1e-7)
      expect_equal(operator$llik, assembled$llik, tolerance=1e-7)
      expect_lte(operator$pcgSolveMaxRelativeResidual, 1e-10)
      diagnostics <- operator$engineDiagnostics$pcg
      expect_gt(diagnostics$operatorCalls, 0)
      expect_gte(diagnostics$operatorVectorEquivalent, diagnostics$operatorCalls)
      expect_length(diagnostics$logDetProbeEstimates, 8L)
      expect_equal(ncol(diagnostics$scoreProbeEstimates), 8L)
      expect_true(is.finite(diagnostics$scoreMaxSE))
    }
  }
})

test_that("full-rank trace control variates recover direct score and information", {
  for(relationship in c("sparse", "dense")){
    workload <- generate_mme_workload(n_ind=8L, n_levels=3L,
      structure="us", relationship=relationship, seed=381L)
    reference <- .pcg_numerical_fit(workload, solver="ldlt")
    candidate <- .pcg_numerical_fit(workload,
      list(traceEstimator="control-variate", traceRank=1000L))
    expect_equal(candidate$bu, reference$bu, tolerance=1e-7)
    expect_equal(candidate$dLu, reference$dLu, tolerance=1e-7)
    expect_equal(candidate$InfMat, reference$InfMat, tolerance=1e-7)
    expect_lte(candidate$engineDiagnostics$pcg$scoreMaxSE, 1e-7)
  }
})

test_that("PCG seeds reproduce probes without modifying the R random stream", {
  workload <- generate_mme_workload(n_ind=6L, n_levels=3L,
    structure="ar1", relationship="none", seed=491L)
  set.seed(193L)
  before <- .Random.seed
  first <- .pcg_numerical_fit(workload, list(seed=51L))
  expect_identical(.Random.seed, before)
  second <- .pcg_numerical_fit(workload, list(seed=51L))
  third <- .pcg_numerical_fit(workload, list(seed=52L))
  expect_identical(first$llik, second$llik)
  expect_identical(first$engineDiagnostics$pcg$logDetProbeEstimates,
    second$engineDiagnostics$pcg$logDetProbeEstimates)
  expect_false(identical(first$engineDiagnostics$pcg$logDetProbeEstimates,
    third$engineDiagnostics$pcg$logDetProbeEstimates))
})

test_that("PCG controls preserve weighted incomplete residual score contractions", {
  dataset <- expand.grid(id=factor(seq_len(4L)), environment=factor(seq_len(3L)),
    replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 5
  dataset <- dataset[-c(1L, 7L), ]
  weights <- diag(nrow(dataset)) + matrix(0.02, nrow(dataset), nrow(dataset))
  fit <- function(solver, pcg=list()){
    modules <- .pcg_module_values(pcg)
    merge <- function(module, defaults=list()){
      values <- modules[[module]]
      if(is.null(values)) values <- list()
      utils::modifyList(defaults, values, keep.null=TRUE)
    }
    control <- mmes_control(
      solve=merge("solve", c(list(method=solver), if(solver == "pcg") list(tol=1e-10))),
      logdet=merge("logdet", if(solver == "pcg") list(lanczosSteps=40L) else list()),
      score=merge("score", if(solver == "pcg") list(traceProbes=8L) else list()),
      convergence=list(maxIters=1L),       update=list(blendWeight=0),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))
    suppressMessages(mmes(response~environment, random=~id,
      rcov=~vsm(usm(environment), ism(record)), W=weights, data=dataset,
      control=control))
  }
  reference <- fit("ldlt")
  for(preconditioner in c("diagonal", "block-jacobi", "nystrom")){
    candidate <- fit("pcg", list(preconditioner=preconditioner,
      rhsMethod="block", traceEstimator="control-variate", traceRank=1000L))
    for(component in c("bu", "dLu", "InfMat")){
      expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
    }
    expect_lte(candidate$engineDiagnostics$pcg$scoreMaxSE, 1e-7)
  }
})

test_that("adaptive PCG calibrates before fitting and refuses false certification", {
  workload <- generate_mme_workload(n_ind=6L, n_levels=3L,
    structure="ar1", relationship="none", seed=591L)
  settings <- list(accuracy="adaptive", traceProbes=2L, lanczosSteps=2L,
    maxTraceProbes=4L, maxLanczosSteps=8L, logDetTol=1e-15, scoreTol=1e-15)
  expect_warning(candidate <- .pcg_numerical_fit(workload, settings),
    "final empirical accuracy checks failed")
  expect_false(candidate$convergence)
  diagnostics <- candidate$engineDiagnostics$pcg
  expect_false(diagnostics$empiricalAccuracyPassed)
  expect_gt(diagnostics$accuracyRefinements, 0)
  expect_equal(diagnostics$traceProbes, 4L)
  expect_equal(diagnostics$lanczosSteps, 8L)
  loose <- .pcg_numerical_fit(workload, list(accuracy="adaptive",
    logDetTol=1e6, scoreTol=1e6))
  expect_true(loose$engineDiagnostics$pcg$empiricalAccuracyPassed)
  expect_equal(loose$engineDiagnostics$pcg$accuracyRefinements, 0L)
})

test_that("PCG fails explicitly on exhausted solve and memory budgets", {
  workload <- generate_mme_workload(n_ind=6L, n_levels=3L,
    structure="ar1", relationship="dense", seed=621L)
  expect_error(.pcg_numerical_fit(workload, list(maxIters=1L)),
    "failed to converge")
  expect_error(.pcg_numerical_fit(workload, denseMemoryMB=1e-6),
    "PCG.*needs about")
})
