test_that("mmes_control exposes strict module-specific controls", {
  expect_true("mmes_control" %in% getNamespaceExports("sommer"))
  control <- mmes_control()
  expect_s3_class(control, "sommer_mmes_control")
  expect_output(print(control), "sommer modular numerical controls")
  expect_identical(control$formulation, "henderson")
  expect_identical(control$solve$method, "auto")
  expect_identical(control$task$objective, "reml")
  expect_identical(control$convergence$maxIters, 30L)
  expect_identical(control$output$inverse, 0L)
  expect_false("pev" %in% names(control$output))
  expect_null(control$update$blendWeight)
  expect_identical(mmes_control(update=list(blendWeight=0))$update$blendWeight, 0)

  # Test custom options and three-state normalizations
  custom <- mmes_control(
    solve=list(method="cholmod", blockBanded="force"),
    execution=list(denseMemoryMB=500),
    score=list(covarianceChecks="off"),
    convergence=list(maxIters=4L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)
  )
  expect_identical(custom$solve$method, "cholmod")
  expect_identical(custom$solve$blockBanded, "force")
  expect_identical(custom$execution$denseMemoryMB, 500)
  expect_identical(custom$score$covarianceChecks, "off")
  expect_identical(custom$convergence$maxIters, 4L)

  # Validation errors
  expect_error(mmes_control(solve=list(unknownSetting=1)), "Unknown solve settings")
  expect_error(mmes_control(output=list(pev=FALSE)), "Unknown output settings")
  expect_error(mmes_control(convergence=list(maxIters=0)), "positive integer")
  expect_error(mmes_control(update=list(blendWeight=NA_real_)), "blendWeight")
})

test_that("fit revalidates canonical controls and exposes public control schema", {
  data(DT_example)
  control <- mmes_control(
    solve=list(method="ldlt"),
    convergence=list(maxIters=2L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)
  )
  fit <- suppressMessages(mmes(Yield~Env, random=~Name, data=DT_example,
    control=control))
  expect_s3_class(fit$control, "sommer_mmes_control")
  expect_identical(fit$control$solve$method, "ldlt")
  expect_identical(fit$control$convergence$maxIters, 2L)
  expect_identical(fit$CiMode, 0L)
  expect_true(all(vapply(fit$uPevList, ncol, integer(1)) == 0L))

  invalid <- control
  invalid$solve$method <- "unknown"
  expect_error(mmes(Yield~Env, random=~Name, data=DT_example,
    control=invalid), "solve\\$method")

  expect_true("..." %in% names(formals(mmes)))
  expect_false("engine" %in% names(formals(mmes)))
  expect_true(".internal" %in% names(formals(mmes)))

  legacyControl <- mmes_control(output=list(inverse=2L),
    convergence=list(maxIters=2L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  expect_warning(
    legacyFit <- mmes(Yield~Env, random=~Name, data=DT_example, control=legacyControl,
      getPEV=FALSE),
    "deprecated"
  )
  expect_identical(legacyFit$control$output$inverse, 0L)
})

test_that("benchmarking framework generates workloads and checks solver parity", {
  workload <- generate_mme_workload(
    name="test_ar1_workload",
    n_ind=20L,
    n_levels=4L,
    structure="ar1",
    seed=123L
  )
  expect_s3_class(workload, "mme_workload")
  expect_equal(workload$dimensions$n_ind, 20L)

  # Run benchmark with two engine configurations
  engines <- list(
    ldlt_base=mmes_control(solve=list(method="ldlt")),
    cholmod=mmes_control(solve=list(method="cholmod"))
  )
  benchmark <- run_mme_benchmark(
    workloads=workload,
    engines=engines,
    replications=1L,
    verifyParity=TRUE,
    tolerance=1e-4,
    nIters=3L,
    verbose=FALSE
  )
  expect_s3_class(benchmark, "sommer_mme_benchmark")
  expect_equal(nrow(benchmark$results), 2L)
  expect_true(all(benchmark$results$parity_ok))
})
