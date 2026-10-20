test_that("PCG controls validate all scalar inputs before coercion", {
  control <- mmes_control(solve=list(method="pcg"))
  expect_false(control$solve$warmStart)
  expect_identical(control$solve$recycleRank, 0L)
  expect_identical(control$solve$recycleMethod, "coarse")
  expect_identical(control$representation$compression, "off")
  expect_identical(control$solve$rhsMethod, "independent")
  expect_identical(control$logdet$scaling, "none")
  expect_identical(control$score$estimator, "hutchinson")
  expect_identical(control$convergence$accuracy, "fixed")
  expect_identical(control$convergence$maxTraceProbes, control$score$traceProbes)
  expect_identical(control$convergence$maxLanczosSteps, control$logdet$lanczosSteps)
  for(name in c("tol")){
    for(value in list(NA_real_, Inf, NaN, 0, -1, numeric(), c(1, 2), "1", 1i)){
      expect_error(mmes_control(solve=setNames(list(value), name)),
        paste0(name, " must be one positive finite number"))
    }
  }
  for(entry in list(
      list(module="solve", name="maxIters"),
      list(module="solve", name="nystromRank"),
      list(module="solve", name="recycleRank"),
      list(module="solve", name="blockSize"),
      list(module="solve", name="residualRefresh"),
      list(module="score", name="traceProbes"),
      list(module="score", name="rank"),
      list(module="logdet", name="lanczosSteps"),
      list(module="convergence", name="maxTraceProbes"),
      list(module="convergence", name="maxLanczosSteps"),
      list(module="execution", name="seed"))){
    name <- entry$name
    for(value in list(NA_real_, Inf, -1, 1.5, numeric(), c(1, 2), "2",
                     .Machine$integer.max + 1)){
      expect_error(do.call(mmes_control,
        setNames(list(setNames(list(value), name)), entry$module)),
        paste0(entry$module, "\\$", name))
    }
  }
  for(entry in list(
      list(module="solve", name="preconditioner"),
      list(module="solve", name="rhsMethod"),
      list(module="solve", name="recycleMethod"),
      list(module="logdet", name="scaling"),
      list(module="score", name="estimator"),
      list(module="convergence", name="accuracy"),
      list(module="representation", name="compression"))){
    name <- entry$name
    for(value in list(NA_character_, "", c("diagonal", "none"), FALSE)){
      expect_error(do.call(mmes_control,
        setNames(list(setNames(list(value), name)), entry$module)),
        paste0(entry$module, "\\$", name))
    }
  }
  for(value in list(NA, 1, logical(), c(TRUE, FALSE))){
    expect_error(mmes_control(solve=list(warmStart=value)), "warmStart")
  }
  expect_error(mmes_control(score=list(traceProbes=16L),
    convergence=list(maxTraceProbes=8L)), "must not be smaller than")
  expect_error(mmes_control(logdet=list(lanczosSteps=40L),
    convergence=list(maxLanczosSteps=20L)), "must not be smaller than")
  expect_error(mmes_control(solve=list(unknownSetting=1)), "Unknown solve")
  expect_error(mmes_control(solve=list(1)), "unique")
  expect_error(mmes_control(solve=list(tol=1e-8, tol=1e-6)), "unique")
})

test_that("non-PCG solvers reject PCG-only settings instead of ignoring them", {
  set.seed(19)
  dataset <- data.frame(id=factor(rep(1:8, each=3)), y=rnorm(24))
  pcg <- list(
    solve=list(preconditioner="block-jacobi", warmStart=TRUE, recycleRank=3L,
      rhsMethod="block", recycleMethod="krylov"),
    logdet=list(scaling="diagonal"))
  for(solver in c("ldlt", "cholmod", "auto")){
    plain <- suppressMessages(mmes(y~1, random=~id, data=dataset, control=mmes_control(solve=list(method=solver), output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
    expect_null(attr(plain$covStruct, "pcgControl"))
    configuredControl <- utils::modifyList(list(solve=list(method=solver)),
      pcg, keep.null=TRUE)
    expect_error(do.call(mmes_control, configuredControl),
      "solve\\$preconditioner, solve\\$warmStart, solve\\$recycleRank, solve\\$recycleMethod, solve\\$rhsMethod, logdet\\$scaling")
  }
})

test_that("known-variance PCG verifies residuals and supports compatible warm starts", {
  set.seed(29)
  dataset <- data.frame(id=factor(rep(1:12, each=3)), x=rnorm(36), y=rnorm(36))
  base <- list(fixed=y~x,
    random=~vsm(ism(id), sigma2=0.7, fixedSigma2=TRUE),
    rcov=~vsm(ism(units), sigma2=1.2, fixedSigma2=TRUE),
    data=dataset)
  solveControl <- mmes_control(task=list(mode="solve"),
    solve=list(tol=1e-11), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  first <- suppressMessages(do.call(mmes, c(base, list(control=solveControl))))
  expect_lte(first$pcg$trueRelativeResidual, 1e-11)
  expect_gt(first$pcg$operatorApplications, 0)
  warm <- suppressMessages(do.call(mmes, c(base, list(covPar=first,
    control=mmes_control(task=list(mode="solve"),
      solve=list(tol=1e-11, warmStart=TRUE, residualRefresh=2L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))))))
  expect_true(warm$pcg$warmStarted)
  expect_true(warm$convergence)
  expect_lte(warm$pcg$trueRelativeResidual, 1e-11)
  expect_equal(warm$bu, first$bu, tolerance=1e-9)
  expect_lte(warm$pcg$iterations, first$pcg$iterations)
  expect_error(suppressMessages(do.call(mmes, c(base, list(
    control=mmes_control(task=list(mode="solve"),
      solve=list(preconditioner="nystrom")))))), "task\\$mode='solve'")
  expect_error(suppressMessages(do.call(mmes, c(base, list(
    control=mmes_control(task=list(mode="solve"),
      solve=list(warmStart=TRUE)))))), "warmStart")
})
