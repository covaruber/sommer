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

test_that("mmes restores real BLAS settings after success, early return, and error", {
  controller <- sommer:::.mmes_blas_controller()
  skip_if(is.null(controller), "No supported BLAS controller available")
  previous <- controller$get()
  on.exit(controller$set(previous), add=TRUE)
  data(DT_example)
  fit <- function(method="ldlt", setup=FALSE){
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
      control=mmes_control(solve=list(method=method),
        convergence=list(maxIters=2), output=list(setup=setup),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE),
        execution=list(blasThreads=1L)))
  }
  result <- fit()
  expect_equal(controller$get(), previous)
  expect_identical(result$threadingDiagnostics$blasThreadsConfigured, 1L)
  expect_gt(result$factorizationDiagnostics$timings$ldlt$calls, 0L)
  expect_gte(result$factorizationDiagnostics$timings$ldlt$elapsedSeconds, 0)
  setup <- fit(setup=TRUE)
  expect_identical(setup$blasThreads, 1L)
  expect_equal(controller$get(), previous)
  expect_error(fit(method="invalid"), "solve\\$method")
  expect_equal(controller$get(), previous)
})

test_that("parallel AI contractions preserve serial fits and cache accounting", {
  environment <- c(R_LIBS=paste(.libPaths(), collapse=.Platform$path.sep),
    OMP_NUM_THREADS="1", OMP_THREAD_LIMIT="2",
    OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
  previous <- Sys.getenv(names(environment), unset=NA_character_)
  on.exit({
    Sys.unsetenv(names(previous)[is.na(previous)])
    do.call(Sys.setenv, as.list(previous[!is.na(previous)]))
  }, add=TRUE)
  do.call(Sys.setenv, as.list(environment))
  script <- tempfile(fileext=".R")
  outputs <- c(tempfile(fileext=".rds"), tempfile(fileext=".rds"))
  on.exit(unlink(c(script, outputs)), add=TRUE)
  writeLines(deparse(substitute({
    arguments <- commandArgs(TRUE)
    library(sommer)
    dataset <- expand.grid(trial=factor(seq_len(6L)),
      genotype=factor(seq_len(12L)), block=factor(seq_len(3L)))
    dataset$y <- 100 + c(45, -30, -12, 52, 18, -20)[dataset$trial] +
      sin(as.integer(dataset$genotype) / 3) * 8 +
      cos(seq_len(nrow(dataset)) / 5) * 3
    models <- list(
      ~vsm(fam(trial, 2), ism(genotype)),
      ~vsm(usm(trial), ism(genotype)),
      ~vsm(fam(trial, 2), ism(genotype)) + vsm(usm(block), ism(genotype))
    )
    results <- list()
    for(repetition in seq_len(2L)){
      for(model in models){
        fit <- mmes(y~trial+block, random=model, data=dataset, control=mmes_control(convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
        results[[length(results) + 1L]] <- c(
          fit[c("theta", "bu", "llik", "InfMat", "dLu")],
          list(contractionDigitEntries=fit$engineDiagnostics$contractionDigitEntries)
        )
      }
    }
    saveRDS(results, arguments[1L])
  })), script)
  for(index in seq_along(outputs)){
    threads <- c(1L, 2L)[index]
    Sys.setenv(OMP_NUM_THREADS=threads)
    output <- system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script), shQuote(outputs[index])),
      stdout=TRUE, stderr=TRUE)
    expect_null(attr(output, "status"), info=paste(output, collapse="\n"))
    expect_true(file.exists(outputs[index]), info=paste(output, collapse="\n"))
    if(!file.exists(outputs[index])) return(invisible(NULL))
  }
  serial <- readRDS(outputs[1L])
  threaded <- readRDS(outputs[2L])
  expect_true(any(vapply(serial, function(fit) fit$contractionDigitEntries > 0, logical(1))))
  expect_equal(threaded, serial, tolerance=1e-8)
})
