test_that("mmes_control constructs valid solver and experimental configuration", {
  ctrl <- mmes_control()
  expect_s3_class(ctrl, "sommer_mme_control")
  expect_equal(ctrl$solver, "auto")
  expect_equal(ctrl$denseMemoryMB, 1000)
  expect_equal(ctrl$experimental$blockBanded, "auto")
  expect_equal(ctrl$experimental$factorSchurMultiFactor, "auto")
  expect_equal(ctrl$experimental$dFirstML, "auto")
  expect_equal(ctrl$experimental$matrixFreePcg, "auto")
  expect_equal(ctrl$experimental$cholmodFinalSelected, "auto")
  expect_equal(ctrl$experimental$pdCertificates, "auto")
  expect_equal(ctrl$experimental$panelCi2, 64L)
  expect_equal(ctrl$pcg$tol, 1e-8)
  expect_equal(ctrl$pcg$preconditioner, "diagonal")

  # Test custom options and three-state normalizations
  custom <- mmes_control(
    solver = "cholmod",
    denseMemoryMB = 500,
    experimental = list(
      blockBanded = "force",
      factorSchurMultiFactor = FALSE,
      dFirstML = "off",
      panelCi2 = 16
    ),
    pcg = list(tol = 1e-6, preconditioner = "nystrom")
  )
  expect_equal(custom$solver, "cholmod")
  expect_equal(custom$denseMemoryMB, 500)
  expect_equal(custom$experimental$blockBanded, "force")
  expect_equal(custom$experimental$factorSchurMultiFactor, "off")
  expect_equal(custom$experimental$dFirstML, "off")
  expect_equal(custom$experimental$panelCi2, 16L)
  expect_equal(custom$pcg$tol, 1e-6)
  expect_equal(custom$pcg$preconditioner, "nystrom")

  # Direct dot arguments for experimental options
  dot_ctrl <- mmes_control(blockBanded = TRUE, factorSchurMultiFactor = "force")
  expect_equal(dot_ctrl$experimental$blockBanded, "force")
  expect_equal(dot_ctrl$experimental$factorSchurMultiFactor, "force")

  # Validation errors
  expect_error(mmes_control(denseMemoryMB = -10), "denseMemoryMB must be a positive finite number")
  expect_error(mmes_control(experimental = list(blockBanded = "invalid")), "Experimental toggle must be")
})

test_that("mmes seamlessly accepts engine argument and reports control diagnostics", {
  data(DT_example)
  ctrl <- mmes_control(
    solver = "ldlt",
    denseMemoryMB = 800,
    experimental = list(pdCertificates = "off")
  )
  fit <- mmes(Yield ~ Env, random = ~ Name, data = DT_example,
              nIters = 2L, engine = ctrl, verbose = FALSE, dateWarning = FALSE)

  expect_s3_class(fit$engineControl, "sommer_mme_control")
  expect_equal(fit$engineControl$solver, "ldlt")
  expect_equal(fit$engineControl$denseMemoryMB, 800)
  expect_equal(fit$engineDiagnostics$control$solver, "ldlt")
  expect_equal(fit$engineDiagnostics$control$denseMemoryMB, 800)
  expect_equal(fit$engineDiagnostics$control$pdCertificates, "off")

  # returnParam includes engine
  param <- mmes(Yield ~ Env, random = ~ Name, data = DT_example,
                returnParam = TRUE, engine = ctrl, verbose = FALSE, dateWarning = FALSE)
  expect_s3_class(param$engine, "sommer_mme_control")
  expect_equal(param$engine$denseMemoryMB, 800)
})

test_that("benchmarking framework generates workloads and checks solver parity", {
  wl <- generate_mme_workload(
    name = "test_ar1_workload",
    n_ind = 20L,
    n_levels = 4L,
    structure = "ar1",
    seed = 123L
  )
  expect_s3_class(wl, "mme_workload")
  expect_equal(wl$dimensions$n_ind, 20L)
  expect_equal(wl$dimensions$n_levels, 4L)
  expect_equal(wl$dimensions$structure, "ar1")

  # Run benchmark with two engine configurations
  engines <- list(
    ldlt_base = mmes_control(solver = "ldlt"),
    cholmod_auto = mmes_control(solver = "cholmod")
  )
  bench <- run_mme_benchmark(
    workloads = wl,
    engines = engines,
    replications = 1L,
    verifyParity = TRUE,
    tolerance = 1e-4,
    nIters = 3L,
    verbose = FALSE
  )
  expect_s3_class(bench, "sommer_mme_benchmark")
  expect_equal(nrow(bench$results), 2L)
  expect_true(all(bench$results$parity_ok))
  expect_true(is.numeric(bench$results$elapsed_sec))

  # Crossover analysis
  co <- analyze_mme_crossovers(bench, baseline_engine = "ldlt_base", candidate_engine = "cholmod_auto")
  expect_s3_class(co, "sommer_mme_crossover")
  expect_true("ar1" %in% names(co$summaries))
  expect_equal(co$summaries$ar1$total_cases, 1L)
})

test_that("mmes formals exclude internal recursion hooks while dots support backward compatibility", {
  m_formals <- names(formals(mmes))
  excluded_hooks <- c(".pqlInner", ".pqlFixedDispersion", ".pqlWorkingPrecision",
                      ".pqlBaseW", ".pqlBaseFactor", ".pqlStart", ".factorScoreParameters",
                      "acceleration", "factorScoreAugmentation",
                      "pcgTol", "pcgMaxIters", "pcgTraceProbes", "pcgLanczosSteps",
                      "pcgPreconditioner", "pcgNystromRank")
  for(hook in excluded_hooks){
    expect_false(hook %in% m_formals, info = paste("Hook should not be in formals:", hook))
  }
  expect_true("engine" %in% m_formals)
  expect_true("..." %in% m_formals)
  expect_true(".internal" %in% m_formals)

  # Passing legacy arguments via dots is cleanly forwarded to engineControl
  data(DT_example)
  fit_accel <- mmes(Yield ~ Env, random = ~ Name, data = DT_example,
                    nIters = 2L, acceleration = "aitken", verbose = FALSE, dateWarning = FALSE)
  expect_equal(fit_accel$engineControl$acceleration, "aitken")

  # Passing internal hooks via dots or .internal works cleanly
  fit_internal <- mmes(Yield ~ Env, random = ~ Name, data = DT_example,
                       nIters = 2L, .internal = list(pqlInner = FALSE),
                       verbose = FALSE, dateWarning = FALSE)
  expect_s3_class(fit_internal, "mmes")
})

test_that("mme_workload_grid expands multiple problem dimensions correctly", {
  grid <- mme_workload_grid(
    structures = c("ar1", "us"),
    n_ind = c(20L),
    n_levels = c(3L, 5L),
    relationships = c("none")
  )
  expect_equal(length(grid), 4L)
  expect_true(all(vapply(grid, function(w) inherits(w, "mme_workload"), logical(1))))
})

test_that("calibrated crossover rules activate multi-factor latent Schur under auto mode", {
  set.seed(42L)
  dataset8 <- expand.grid(
    id = factor(seq_len(15L)),
    environment = factor(seq_len(8L)),
    replicate = factor(seq_len(2L))
  )
  dataset8$response <- stats::rnorm(nrow(dataset8))
  M8 <- matrix(stats::rnorm(15L * 10L), 15L, 10L)
  G8 <- stats::cor(t(M8)) + diag(0.1, 15L)
  dimnames(G8) <- list(levels(dataset8$id), levels(dataset8$id))
  precision8 <- solve(G8)
  attr(precision8, "inverse") <- TRUE

  fit_auto8 <- mmes(
    response ~ environment,
    random = ~ vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu = precision8),
    data = dataset8,
    solver = "cholmod",
    nIters = 2L,
    computeCi = 0L,
    getPEV = FALSE,
    verbose = FALSE,
    dateWarning = FALSE,
    engine = mmes_control(solver = "cholmod", experimental = list(factorSchurMultiFactor = "auto"))
  )
  expect_true(fit_auto8$engineDiagnostics$factorSchurActive)

  dataset4 <- expand.grid(
    id = factor(seq_len(15L)),
    environment = factor(seq_len(4L)),
    replicate = factor(seq_len(2L))
  )
  dataset4$response <- stats::rnorm(nrow(dataset4))
  fit_auto4 <- mmes(
    response ~ environment,
    random = ~ vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu = precision8),
    data = dataset4,
    solver = "cholmod",
    nIters = 2L,
    computeCi = 0L,
    getPEV = FALSE,
    verbose = FALSE,
    dateWarning = FALSE,
    engine = mmes_control(solver = "cholmod", experimental = list(factorSchurMultiFactor = "auto"))
  )
  expect_false(fit_auto4$engineDiagnostics$factorSchurActive)
})
