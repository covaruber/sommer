test_that("PCG controls validate all scalar inputs before coercion", {
  control <- mmes_control(solver="pcg")
  expect_false(control$pcg$warmStart)
  expect_identical(control$pcg$recycleRank, 0L)
  expect_identical(control$pcg$rhsMethod, "independent")
  expect_identical(control$pcg$logDetScaling, "none")
  expect_identical(control$pcg$traceEstimator, "hutchinson")
  expect_identical(control$pcg$accuracy, "fixed")
  expect_identical(control$pcg$maxTraceProbes, control$pcg$traceProbes)
  expect_identical(control$pcg$maxLanczosSteps, control$pcg$lanczosSteps)
  for(name in c("tol", "logDetTol", "scoreTol")){
    for(value in list(NA_real_, Inf, NaN, 0, -1, numeric(), c(1, 2), "1", 1i)){
      expect_error(mmes_control(pcg=setNames(list(value), name)), paste0("pcg\\$", name))
    }
  }
  for(name in c("maxIters", "traceProbes", "lanczosSteps", "nystromRank",
               "recycleRank", "blockSize", "residualRefresh", "traceRank",
               "maxTraceProbes", "maxLanczosSteps", "seed")){
    for(value in list(NA_real_, Inf, -1, 1.5, numeric(), c(1, 2), "2",
                     .Machine$integer.max + 1)){
      expect_error(mmes_control(pcg=setNames(list(value), name)), paste0("pcg\\$", name))
    }
  }
  for(name in c("preconditioner", "rhsMethod", "logDetScaling", "traceEstimator", "accuracy")){
    for(value in list(NA_character_, "", c("diagonal", "none"), FALSE)){
      expect_error(mmes_control(pcg=setNames(list(value), name)), paste0("pcg\\$", name))
    }
  }
  for(value in list(NA, 1, logical(), c(TRUE, FALSE))){
    expect_error(mmes_control(pcg=list(warmStart=value)), "warmStart")
  }
  expect_error(mmes_control(pcg=list(traceProbes=16L, maxTraceProbes=8L)), "accuracy limits")
  expect_error(mmes_control(pcg=list(lanczosSteps=40L, maxLanczosSteps=20L)), "accuracy limits")
  expect_error(mmes_control(pcg=list(unknownSetting=1)), "Unknown pcg")
  expect_error(mmes_control(pcg=list(1)), "unique")
  expect_error(mmes_control(pcg=list(tol=1e-8, tol=1e-6)), "unique")
})

test_that("PCG controls survive direct, engine and legacy merging", {
  pcg <- list(preconditioner="block-jacobi", warmStart=TRUE, recycleRank=4L,
    rhsMethod="block", blockSize=4L, residualRefresh=20L,
    logDetScaling="diagonal", traceEstimator="control-variate", traceRank=3L,
    seed=43L)
  direct <- sommer:::.merge_mmes_control(solver="pcg", pcg=pcg)
  engine <- sommer:::.merge_mmes_control(
    engine=mmes_control(solver="pcg"), pcg=pcg)
  expect_identical(direct$pcg, engine$pcg)
  expect_identical(direct$pcg$preconditioner, "block-jacobi")
  overridden <- sommer:::.merge_mmes_control(
    engine=mmes_control(solver="pcg"), pcg=list(traceProbes=32L))
  expect_identical(overridden$pcg$maxTraceProbes, 32L)
  expect_error(sommer:::.merge_mmes_control(
    engine=mmes_control(solver="pcg"), pcg=list(seed=NA)), "seed")
  expect_error(sommer:::.merge_mmes_control(pcgMaxIters=2.5), "maxIters")
  expect_error(sommer:::.merge_mmes_control(
    engine=mmes_control(), pcgTraceProbes=2.5), "traceProbes")
  expect_error(sommer:::.merge_mmes_control(pcgTol="1e-8"), "tol")
  expect_identical(sommer:::.merge_mmes_control(
    pcgTraceProbes=2L, pcg=list(maxTraceProbes=4L))$pcg$maxTraceProbes, 4L)
})

test_that("non-PCG setup does not forward experimental PCG state", {
  set.seed(19)
  dataset <- data.frame(id=factor(rep(1:8, each=3)), y=rnorm(24))
  pcg <- list(preconditioner="block-jacobi", warmStart=TRUE, recycleRank=3L,
    rhsMethod="block", logDetScaling="diagonal")
  for(solver in c("ldlt", "cholmod", "auto")){
    plain <- suppressMessages(mmes(y~1, random=~id, data=dataset,
      solver=solver, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE))
    configured <- suppressMessages(mmes(y~1, random=~id, data=dataset,
      engine=mmes_control(solver=solver, pcg=pcg), returnParam=TRUE,
      verbose=FALSE, dateWarning=FALSE))
    expect_null(attr(configured$covStruct, "pcgControl"))
    for(name in c("X", "Z", "Ai", "yvar", "W", "residualBlock", "residualIndex")){
      expect_identical(configured[[name]], plain[[name]], info=paste(solver, name))
    }
    fit <- function(settings){
      suppressMessages(mmes(y~1, random=~id, data=dataset, nIters=2L,
        engine=mmes_control(solver=solver, pcg=settings), blendWeight=0,
        verbose=FALSE, dateWarning=FALSE))
    }
    reference <- fit(list())
    candidate <- fit(pcg)
    expect_null(reference$engineDiagnostics$pcg)
    expect_null(candidate$engineDiagnostics$pcg)
    for(name in c("llik", "theta", "covPar", "bu", "C", "InfMat", "theta_se")){
      expect_identical(candidate[[name]], reference[[name]], info=paste(solver, name))
    }
  }
})

test_that("known-variance PCG verifies residuals and supports compatible warm starts", {
  set.seed(29)
  dataset <- data.frame(id=factor(rep(1:12, each=3)), x=rnorm(36), y=rnorm(36))
  base <- list(fixed=y~x,
    random=~vsm(ism(id), sigma2=0.7, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=1.2, fixedSigma2=TRUE),
    data=dataset, solveOnly=TRUE, verbose=FALSE, dateWarning=FALSE)
  first <- suppressMessages(do.call(mmes, c(base, list(pcgTol=1e-11))))
  expect_lte(first$pcg$trueRelativeResidual, 1e-11)
  expect_gt(first$pcg$operatorApplications, 0)
  warm <- suppressMessages(do.call(mmes, c(base, list(covPar=first,
    engine=mmes_control(pcg=list(tol=1e-11, warmStart=TRUE, residualRefresh=2L))))))
  expect_true(warm$pcg$warmStarted)
  expect_true(warm$convergence)
  expect_lte(warm$pcg$trueRelativeResidual, 1e-11)
  expect_equal(warm$bu, first$bu, tolerance=1e-9)
  expect_lte(warm$pcg$iterations, first$pcg$iterations)
  expect_error(suppressMessages(do.call(mmes, c(base, list(
    engine=mmes_control(pcg=list(preconditioner="nystrom")))))), "solveOnly")
  expect_error(suppressMessages(do.call(mmes, c(base, list(
    engine=mmes_control(pcg=list(warmStart=TRUE)))))), "warmStart")
})
