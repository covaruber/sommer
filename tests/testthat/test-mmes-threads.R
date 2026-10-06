test_that("BLAS controls are opt-in and restored", {
  configured <- 4L
  updates <- integer()
  controller <- list(provider="mock", get=function() configured,
    set=function(value){ configured <<- value; updates <<- c(updates, value) })
  scope <- sommer:::.mmes_blas_scope(controller=controller)
  expect_identical(scope$info$blasThreadsConfigured, 4L)
  scope$restore()
  expect_length(updates, 0L)
  scope <- sommer:::.mmes_blas_scope(2L, controller)
  expect_identical(configured, 2L)
  scope$restore()
  expect_identical(configured, 4L)
  expect_identical(updates, c(2L, 4L))
})

test_that("BLAS controls reject invalid and unsupported requests", {
  for(value in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), TRUE, "all", .Machine$integer.max + 1)){
    expect_error(sommer:::.mmes_blas_scope(value, NULL), "blasThreads must")
  }
  expect_error(sommer:::.mmes_blas_scope(2L, NULL), "Runtime BLAS control is unavailable")
  expect_true(is.na(sommer:::.mmes_blas_scope(controller=NULL)$info$blasThreadsConfigured))
})

test_that("BLAS restoration also handles setter errors and ignored settings", {
  configured <- 4L
  controller <- list(provider="mock", get=function() configured,
    set=function(value){
      configured <<- value
      if(value == 2L) stop("setter failed")
    })
  expect_error(sommer:::.mmes_blas_scope(2L, controller), "setter failed")
  expect_identical(configured, 4L)
  controller$set <- function(value) invisible(NULL)
  expect_error(sommer:::.mmes_blas_scope(2L, controller), "did not accept")
  expect_identical(configured, 4L)
})

test_that("CPU measurements are labelled separately from configured threads", {
  scope <- sommer:::.mmes_blas_scope(controller=NULL)
  expect_output(sommer:::.mmes_thread_report(scope), "configured: unknown")
  result <- sommer:::.mmes_thread_finish(list(), scope)
  expect_gte(result$threadingDiagnostics$cpuSeconds, 0)
  expect_gte(result$threadingDiagnostics$elapsedSeconds, 0)
  expect_output(sommer:::.mmes_thread_finish(list(), scope, TRUE), "not an active-thread count")
})

test_that("mmes restores real BLAS settings after success, early return, and error", {
  controller <- sommer:::.mmes_blas_controller()
  skip_if(is.null(controller), "No supported BLAS controller available")
  previous <- controller$get()
  on.exit(controller$set(previous), add=TRUE)
  data(DT_example)
  fit <- function(...) mmes(Yield~Env, random=~Name, rcov=~units,
    data=DT_example, nIters=2, verbose=FALSE, dateWarning=FALSE, blasThreads=1L, ...)
  result <- fit(solver="ldlt")
  expect_equal(controller$get(), previous)
  expect_identical(result$threadingDiagnostics$blasThreadsConfigured, 1L)
  expect_gt(result$factorizationDiagnostics$timings$ldlt$calls, 0L)
  expect_gte(result$factorizationDiagnostics$timings$ldlt$elapsedSeconds, 0)
  setup <- fit(returnParam=TRUE)
  expect_identical(setup$blasThreads, 1L)
  expect_equal(controller$get(), previous)
  expect_error(fit(solver="invalid"), "solver must")
  expect_equal(controller$get(), previous)
})

test_that("thread control preserves estimates and reports the actual numerical path", {
  controller <- sommer:::.mmes_blas_controller()
  skip_if(is.null(controller), "No supported BLAS controller available")
  previous <- controller$get()
  on.exit(controller$set(previous), add=TRUE)
  probe <- try(sommer:::.mmes_blas_scope(2L, controller), silent=TRUE)
  skip_if(inherits(probe, "try-error"), "The backend does not support two BLAS threads")
  probe$restore()
  data(DT_example)
  fit <- function(threads) mmes(Yield~Env, random=~Name, rcov=~units,
    data=DT_example, nIters=4, verbose=FALSE, dateWarning=FALSE,
    solver="cholmod", blasThreads=threads)
  serial <- fit(1L)
  threaded <- fit(2L)
  expect_equal(threaded$theta, serial$theta, tolerance=1e-8)
  expect_equal(threaded$bu, serial$bu, tolerance=1e-8)
  expect_equal(threaded$llik, serial$llik, tolerance=1e-8)
  expect_identical(threaded$threadingDiagnostics$numericalPath, "block Schur")
  expect_gt(threaded$threadingDiagnostics$factorizations$blockSchur$calls, 0L)
  old <- options(sommer.mme.denseMemoryMB=1e-6)
  on.exit(options(old), add=TRUE)
  supernodal <- fit(1L)
  expect_false(supernodal$engineDiagnostics$blockSchurActive)
  expect_gt(supernodal$threadingDiagnostics$factorizations$cholmod$calls, 0L)
  expect_equal(supernodal$theta, serial$theta, tolerance=1e-8)
  expect_equal(supernodal$bu, serial$bu, tolerance=1e-8)
  expect_equal(controller$get(), previous)
})

test_that("auto BLAS sizing uses allocation-aware CPU availability", {
  skip_if_not_installed("parallelly")
  configured <- 1L
  controller <- list(provider="mock", get=function() configured,
    set=function(value) configured <<- value)
  scope <- sommer:::.mmes_blas_scope("auto", controller)
  expect_equal(configured, as.integer(parallelly::availableCores()))
  scope$restore()
  expect_identical(configured, 1L)
})

test_that("the scaling runner records both paths and numerical agreement", {
  controller <- sommer:::.mmes_blas_controller()
  skip_if(is.null(controller), "No supported BLAS controller available")
  previous <- controller$get()
  on.exit(controller$set(previous), add=TRUE)
  benchmarkEnvironment <- new.env(parent=asNamespace("sommer"))
  sys.source(test_path("..", "..", "benchmarks", "thread-scaling.R"),
    envir=benchmarkEnvironment)
  models <- benchmarkEnvironment$benchmark_models(ids=8L, envs=2L, iterations=2L)
  output <- tempfile(fileext=".csv")
  stem <- sub("\\.csv$", "", output)
  on.exit(unlink(c(output, paste0(stem, ".summary.csv"), paste0(stem, ".runtime.rds"))), add=TRUE)
  expect_warning(
    result <- withCallingHandlers(benchmarkEnvironment$run_thread_scaling(
      1L, models, repeats=1L, output=output), warning=function(condition){
        if(grepl("^Set OMP_NUM_THREADS", conditionMessage(condition))){
          invokeRestart("muffleWarning")
        }
      }), NA)
  expect_equal(nrow(result$results), 2L)
  expect_true(all(result$results$numericalAgreement))
  expect_setequal(result$results$numericalPath, c("cholmod", "block Schur"))
  expect_true(all(result$results$factorCalls > 0L))
  expect_true(all(result$summary$fitSpeedup == 1))
  expect_true(file.exists(output))
  expect_true(file.exists(paste0(stem, ".summary.csv")))
  expect_true(file.exists(paste0(stem, ".runtime.rds")))
  expect_equal(nrow(utils::read.csv(output)), 2L)
  expect_equal(controller$get(), previous)
})