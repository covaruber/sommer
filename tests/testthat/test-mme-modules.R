test_that("modular controls preserve the established default numerical bundle", {
  control <- mmes_control()
  expect_s3_class(control, "sommer_mmes_control")
  expect_output(print(control), "modular numerical controls")
  engine <- sommer:::.mmes_compile_control(control)
  expect_equal(control$formulation, "henderson")
  expect_equal(control$task$objective, "reml")
  expect_equal(control$update$blendWeight, NULL)
  expect_equal(control$curvature$method, "exact-ai")
  expect_equal(control$proposal$extrapolation, "none")
  expect_equal(control$convergence$maxIters, 30L)
  expect_equal(engine$solver, "auto")
  expect_equal(engine$pcg, sommer:::.mme_pcg_control())
  expect_equal(engine$acceleration, "none")
  expect_equal(engine$compression, "off")
  expect_equal(engine$covarianceInformation, "full")
  expect_true(engine$retainCoefficientMatrix)
})

test_that("old controls and misspelled module entries fail explicitly", {
  expect_error(mmes_control(solver="ldlt"), "unused argument")
  expect_error(mmes_control(experimental=list()), "unused argument")
  expect_error(mmes_control(acceleration="aitken"), "unused argument")
  expect_error(mmes_control(solve=list(sovler="ldlt")), "Unknown solve settings")
  expect_error(mmes_control(score=list(traceProbes=NA)), "traceProbes")
  expect_error(mmes_control(output=list(inverse=1.2)), "output\\$inverse")
  expect_error(mmes_control(proposal=list(extrapolation="aitken"),
    curvature=list(method="hybrid-ai")), "not a supported combination")
  expect_error(mmes_control(solve=list(tol=0)), "tol")
  expect_error(mmes_control(task=list(pql=list(typo=TRUE))), "Unknown task\\$pql")
  expect_error(mmes_control(solve=setNames(list(1,2), c("tol","tol"))), "unique")
  expect_error(mmes(y~1, solver="ldlt"), "unused argument")
  expect_error(mmes(y~1, control=list()), "created by mmes_control")
})

test_that("module-owned estimator budgets and outputs reach private transport", {
  control <- mmes_control(
    solve=list(method="pcg", tol=1e-10, rhsMethod="block"),
    logdet=list(scaling="diagonal", lanczosSteps=30),
    score=list(traceProbes=20, estimator="control-variate", rank=8),
    convergence=list(accuracy="adaptive", maxTraceProbes=40, maxLanczosSteps=60),
    output=list(covariance="none", inverse=2),
    execution=list(seed=9, denseMemoryMB=500))
  engine <- sommer:::.mmes_compile_control(control)
  expect_equal(engine$pcg$traceProbes, 20L)
  expect_equal(engine$pcg$maxTraceProbes, 40L)
  expect_equal(engine$pcg$lanczosSteps, 30L)
  expect_equal(engine$pcg$seed, 9L)
  expect_equal(engine$pcg$traceEstimator, "control-variate")
  expect_equal(engine$covarianceInformation, "none")
  expect_equal(control$output$inverse, 2L)
})

test_that("modular fits retain defaults, public calls and rebuildable inference inputs", {
  data(DT_example)
  control <- mmes_control(convergence=list(maxIters=5),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  fit <- mmes(Yield~Env, random=~Name, data=DT_example, control=control)
  expect_identical(fit$control, control)
  expect_null(fit$engineControl)
  expect_equal(fit$executionPlan$solver$selected, "ldlt")
  expect_equal(as.character(fit$call[[1L]]), "mmes")
  inputs <- sommer:::.mmes_rebuild_inputs(fit)
  expect_equal(nrow(inputs$X), nrow(fit$W))
  setupControl <- control
  setupControl$output$setup <- TRUE
  setup <- mmes(Yield~Env, random=~Name, data=DT_example, control=setupControl)
  expect_equal(setup$executionPlan$status, "assembled")
  expect_equal(setup$blendWeight, rep(0, 5))
  expect_true(isTRUE(attr(setup$covStruct, "adaptiveRecovery")))
  expect_equal(setup$stepWeight, rep(1,5))
})

test_that("control mutation is revalidated before fitting", {
  control <- mmes_control()
  control$solve$typo <- 1
  expect_error(mmes(y~1, control=control), "Unknown solve settings")
  control <- mmes_control()
  control$convergence$maxIters <- 0
  expect_error(mmes(y~1, control=control), "convergence\\$maxIters")
})

test_that("incompatible or unused module settings stop before fitting", {
  conflicts <- list(
    list(list(formulation="direct", solve=list(method="pcg")), "solve\\$method"),
    list(list(formulation="direct", output=list(inverse=1L)), "output\\$inverse"),
    list(list(formulation="direct", execution=list(gpu=TRUE)), "execution\\$gpu"),
    list(list(formulation="direct", curvature=list(method="hybrid-ai")), "curvature\\$method"),
    list(list(formulation="direct", representation=list(compression="groups")), "representation\\$compression"),
    list(list(task=list(mode="solve"), output=list(inverse=2L)), "output\\$inverse"),
    list(list(task=list(mode="solve"), solve=list(method="ldlt")), "solve\\$method"),
    list(list(task=list(mode="solve"), formulation="direct"), "formulation"),
    list(list(solve=list(method="ldlt"), score=list(traceProbes=50L)), "score\\$traceProbes"),
    list(list(logdet=list(lanczosSteps=40L)), "logdet\\$lanczosSteps"),
    list(list(solve=list(method="pcg"), proposal=list(extrapolation="aitken")), "proposal\\$extrapolation"),
    list(list(solve=list(method="pcg"), update=list(blendWeight=0.5)), "update\\$blendWeight"),
    list(list(solve=list(method="pcg"), task=list(objective="ml")), "task\\$objective"),
    list(list(solve=list(method="pcg"), update=list(recovery="pxem")), "update\\$recovery"),
    list(list(representation=list(equations="saddle-point"), proposal=list(extrapolation="aitken")), "representation\\$equations"),
    list(list(representation=list(factorScores="profile"), output=list(inverse=2L)), "output\\$inverse"),
    list(list(update=list(blendWeight=0.3, commonScale="force")), "update\\$commonScale"),
    list(list(output=list(coefficientMatrix=FALSE)), "output\\$coefficientMatrix"))
  for(case in conflicts){
    error <- tryCatch(do.call(mmes_control, case[[1L]]), error=conditionMessage)
    expect_type(error, "character")
    expect_match(error, "Incompatible mmes_control\\(\\) settings", info=case[[2L]])
    expect_match(error, paste0("Offending: [^.]*", case[[2L]]), info=case[[2L]])
  }
  expect_match(tryCatch(mmes_control(formulation="direct", solve=list(method="pcg"),
    curvature=list(method="hybrid-ai")), error=conditionMessage),
    "Offending: solve\\$method, curvature\\$method")
})

test_that("compatible settings, defaults and round trips remain accepted", {
  valid <- list(list(), list(formulation="direct", output=list(inverse=2L),
      update=list(recovery="pxem", blendWeight=0.3), proposal=list(stepWeight=0.8)),
    list(solve=list(method="pcg", tol=1e-9), score=list(traceProbes=4L),
      logdet=list(lanczosSteps=30L), execution=list(seed=3L)),
    list(task=list(mode="solve"), solve=list(method="pcg", warmStart=TRUE)),
    list(curvature=list(method="hybrid-ai"), update=list(blendWeight=0)),
    list(representation=list(equations="saddle-point"),
      output=list(coefficientMatrix=FALSE)),
    list(solve=list(method="ldlt"), curvature=list(reduced="off"),
      update=list(commonScale="off")))
  for(case in valid){
    control <- do.call(mmes_control, case)
    expect_s3_class(control, "sommer_mmes_control")
    expect_s3_class(do.call(mmes_control, unclass(control)), "sommer_mmes_control")
  }
  dataset <- data.frame(y=c(1.2, 2.1, 2.9, 4.2, 5.1, 5.8),
    id=factor(rep(seq_len(3L), 2L)))
  expect_error(mmes(y~1, random=~id, data=dataset, control=mmes_control(
    task=list(pql=list(maxit=3L)), diagnostics=list(verbose=FALSE))),
    "Offending: task\\$pql")
})
