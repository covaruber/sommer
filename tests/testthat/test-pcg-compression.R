.compression_fixture <- function() {
  dataset <- expand.grid(id=factor(seq_len(8L)), environment=factor(seq_len(3L)),
    replicate=seq_len(4L))
  dataset$response <- sin(seq_len(nrow(dataset)) * 0.7) +
    as.integer(dataset$id) / 4 + as.integer(dataset$environment) / 3
  dataset$x <- as.integer(dataset$id) / 8
  dataset
}

.compression_fit <- function(dataset, method="off", pcg=list(), weights=NULL,
                              iterations=1L, solver="pcg", rcov=~vsm(dsm(environment), ism(units)),
                              inverse=0L) {
  solveNames <- c("preconditioner", "nystromRank", "warmStart", "recycleRank",
    "recycleMethod", "rhsMethod", "blockSize", "residualRefresh")
  logdetNames <- "logDetScaling"
  scoreNames <- c("traceEstimator", "traceRank")
  unknown <- setdiff(names(pcg), c(solveNames, logdetNames, scoreNames))
  if(length(unknown)) stop("Unknown PCG settings: ", paste(unknown, collapse=", "))
  solve <- c(list(method=solver, tol=1e-11),
    pcg[names(pcg) %in% solveNames])
  logdet <- c(list(lanczosSteps=40L),
    if("logDetScaling" %in% names(pcg))
      list(scaling=pcg$logDetScaling) else list())
  score <- c(list(traceProbes=12L),
    if("traceEstimator" %in% names(pcg))
      list(estimator=pcg$traceEstimator) else list(),
    if("traceRank" %in% names(pcg))
      list(rank=pcg$traceRank) else list())
  control <- mmes_control(
    representation=list(compression=method),
    solve=solve, logdet=logdet, score=score,
    convergence=list(maxIters=iterations),
    output=list(inverse=inverse),
    update=list(blendWeight=0),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  arguments <- list(fixed=response~environment+x, random=~id, rcov=rcov,
    data=dataset, control=control)
  if(!is.null(weights)) arguments$W <- weights
  suppressMessages(do.call(mmes, arguments))
}

test_that("compression preserves weighted likelihood, score, AI and original records", {
  dataset <- .compression_fixture()
  for(weights in list(NULL, Matrix::Diagonal(x=seq(0.5, 2, length.out=nrow(dataset))))){
    for(estimator in c("hutchinson", "control-variate")){
      controls <- list(traceEstimator=estimator, traceRank=1000L)
      reference <- .compression_fit(dataset, pcg=controls, weights=weights)
      for(method in c("groups", "classes")){
        candidate <- .compression_fit(dataset, method, pcg=controls, weights=weights)
        for(field in c("bu", "C", "llik", "dLu", "InfMat", "theta", "theta_se",
                      "monitor", "yWorking")){
          expect_equal(candidate[[field]], reference[[field]], tolerance=1e-7, info=paste(method, estimator, field))
        }
        expect_equal(candidate$W, reference$W)
        expect_equal(candidate$y, reference$y)
        expect_equal(residuals(candidate), residuals(reference), tolerance=1e-7)
        expect_equal(fitted(candidate), fitted(reference), tolerance=1e-7)
        diagnostics <- candidate$engineDiagnostics$pcg$compression
        expect_equal(diagnostics$observations, nrow(dataset))
        expect_equal(diagnostics$groups, 24)
        expect_gt(diagnostics$retainedBytes, 0)
        expect_gte(diagnostics$setupSeconds, 0)
        expect_lte(candidate$pcgSolveMaxRelativeResidual, 1e-11)
      }
    }
  }
})

test_that("exact grouping keeps distinct covariates and residual variance classes", {
  dataset <- .compression_fixture()
  dataset$x <- seq_len(nrow(dataset)) / nrow(dataset)
  reference <- .compression_fit(dataset)
  candidate <- .compression_fit(dataset, "groups")
  expect_equal(candidate$engineDiagnostics$pcg$compression$groups, nrow(dataset))
  expect_equal(candidate$bu, reference$bu, tolerance=1e-7)
  expect_equal(candidate$dLu, reference$dLu, tolerance=1e-7)
  fit <- function(method) suppressMessages(mmes(response~environment+x,
    random=~id, rcov=~vsm(dsm(environment), ism(units)),
    data=.compression_fixture(), control=mmes_control(
      solve=list(method="pcg"),
      representation=list(compression=method),
      score=list(traceProbes=12L), logdet=list(lanczosSteps=40L),
      convergence=list(maxIters=1L),       update=list(blendWeight=0),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
  reference <- fit("off")
  candidate <- fit("groups")
  expect_equal(candidate$engineDiagnostics$pcg$compression$groups, 24)
  expect_equal(candidate$dLu, reference$dLu, tolerance=1e-7)
  expect_equal(candidate$InfMat, reference$InfMat, tolerance=1e-7)
})

test_that("compression explicitly rejects changing weights and correlated residuals", {
  dataset <- .compression_fixture()
  expect_error(.compression_fit(dataset, "groups", rcov=~vsm(usm(environment), ism(units))),
    "independent residuals")
  denseWeights <- diag(nrow(dataset)) + matrix(0.01, nrow(dataset), nrow(dataset))
  expect_error(.compression_fit(dataset, "groups", weights=denseWeights),
    "independent residuals")
  expect_error(suppressMessages(mmes(response~environment+x, random=~id, data=dataset,
    family=poisson(), control=mmes_control(
      solve=list(method="pcg"),
      representation=list(compression="groups"),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))), "Gaussian")
})

test_that("Krylov recycling works without observation compression across residual models", {
  for(structure in c("ar1", "us")){
    workload <- generate_mme_workload(n_ind=8L, n_levels=3L,
      structure=structure, relationship="sparse", seed=771L)
    fit <- function(pcg) {
      unknown <- setdiff(names(pcg), c("recycleRank", "recycleMethod", "warmStart"))
      if(length(unknown)) stop("Unknown PCG settings: ", paste(unknown, collapse=", "))
      control <- mmes_control(solve=c(list(method="pcg"), pcg),
        score=list(traceProbes=8L), logdet=list(lanczosSteps=40L),
        convergence=list(maxIters=4L),         update=list(blendWeight=0),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE))
      suppressMessages(mmes(workload$fixed, random=workload$random,
        rcov=workload$rcov, data=workload$data, control=control))
    }
    reference <- fit(list())
    candidate <- fit(list(recycleRank=4L, recycleMethod="krylov", warmStart=TRUE))
    for(field in c("bu", "llik", "dLu", "theta", "InfMat")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-5)
    }
    expect_gt(candidate$engineDiagnostics$pcg$balancedPreconditionerApplications, 0)
    expect_true(candidate$pcgMatrixFree)
  }
})
