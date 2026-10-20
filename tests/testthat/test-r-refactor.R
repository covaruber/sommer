refactor_fit <- function(){
  set.seed(814)
  data <- expand.grid(id=factor(paste0("g", seq_len(8))),
                      env=factor(c("a", "b")), rep=seq_len(2))
  data$y <- 2 + as.integer(data$env) + rnorm(nrow(data))
  mmes(y~env, random=~id, data=data,
       control=mmes_control(convergence=list(maxIters=4),
                            diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
}

test_that("operation-local covariance solves agree with independent solves", {
  fit <- refactor_fit()
  set.seed(91)
  D <- Matrix::Matrix(matrix(rnorm(3*nrow(fit$bu)), 3), sparse=TRUE)
  createFactor <- sommer:::mmes_vcov_factor_cpp
  factorizations <- 0L
  local_mocked_bindings(mmes_vcov_factor_cpp=function(model){
    factorizations <<- factorizations + 1L
    createFactor(model)
  }, .package="sommer")
  solve <- sommer:::.mmes_vcov_solver(fit)
  expected <- sommer:::predict_mmes_vcov_cpp(fit, D)
  expect_equal(solve(D), expected, tolerance=1e-12)
  expect_equal(solve(D[1,,drop=FALSE]), expected[1,1,drop=FALSE],
               tolerance=1e-12)
  expect_identical(factorizations, 1L)
  expect_equal(sommer:::.mmes_henderson_vcov(fit, list(D, 2*D)),
               list(expected, 4*expected), tolerance=1e-12)
  expect_error(solve(D[,-1,drop=FALSE]), "one column")
  factor <- sommer:::mmes_vcov_factor_cpp(fit)
  expired <- unserialize(serialize(factor, NULL))
  expect_error(sommer:::mmes_vcov_solve_cpp(expired, D), "invalid")
  expect_false(any(vapply(fit, inherits, logical(1),
                          what="sommer_prediction_factor")))
})

test_that("benchmark rows retain failures and actual convergence status", {
  calls <- 0L
  local_mocked_bindings(mmes=function(..., control){
    if(control$diagnostics$dateWarning) stop("unexpected date warning")
    if(isTRUE(control$diagnostics$fail)) stop("deliberate benchmark failure")
    calls <<- calls + 1L
    list(llik=c(-10, -2), monitor=matrix(1, 1, 2),
         theta=list(matrix(1)), bu=matrix(1), convergence=TRUE,
         engineDiagnostics=list())
  }, .package="sommer")
  workload <- generate_mme_workload(n_ind=4L, n_levels=2L)
  baseline <- mmes_control()
  failed <- baseline
  failed$diagnostics$fail <- TRUE
  result <- run_mme_benchmark(workload,
    engines=list(baseline=baseline, failed=failed), replications=1L,
    verbose=FALSE)$results
  expect_equal(nrow(result), 2L)
  expect_identical(result$converged, c(TRUE, FALSE))
  expect_true(result$parity_ok[1])
  expect_true(is.na(result$parity_ok[2]))
  expect_match(result$error[2], "deliberate benchmark failure")
  expect_identical(calls, 2L)
})

test_that("benchmark parity compares final finite estimates, not histories", {
  reference <- list(llik=c(-10,-2), theta=list(1), bu=matrix(1))
  fit <- reference
  fit$llik <- c(-100, -20, -2)
  expect_equal(sommer:::.mme_benchmark_difference(fit, reference), 0)
  fit$llik <- NA_real_
  expect_identical(sommer:::.mme_benchmark_difference(fit, reference), Inf)
  workloads <- mme_workload_grid(structures="ar1", n_ind=4L, n_levels=2L,
                                 relationships="none", computeCi=c(0L, 1L, 2L))
  expect_length(workloads, 3L)
  expect_identical(anyDuplicated(names(workloads)), 0L)
})

test_that("timed benchmark failures retain a complete result row", {
  calls <- 0L
  local_mocked_bindings(mmes=function(..., control){
    calls <<- calls + 1L
    if(calls == 2L) stop("timed run failed")
    list(llik=-2, theta=list(1), bu=matrix(1), monitor=matrix(1),
         convergence=TRUE, engineDiagnostics=list(
           profile=list(seconds=c(assembly=0.1), counters=list(nEffects=2L))))
  }, .package="sommer")
  workload <- generate_mme_workload(n_ind=4L, n_levels=2L)
  result <- run_mme_benchmark(workload,
    engines=list(baseline=mmes_control(), candidate=mmes_control()),
    replications=1L, verbose=FALSE)$results
  expect_equal(nrow(result), 2L)
  expect_match(result$error[1L], "timed run failed")
  expect_true(is.na(result$parity_ok[2L]))
  expect_equal(result$nEffects[2L], 2L)
  expect_true(all(is.na(result$freeParams)))
})

test_that("legacy exports and covariance interfaces survive consolidation", {
  expect_true(all(vapply(c("mmer", "vs", "vsr", "atr", "csr", "dsr", "usr",
                           "unsm", "unsm2", "fixm", "GWAS"),
                         function(name) is.function(getExportedValue("sommer", name)),
                         logical(1))))
  expect_identical(formals(vs), formals(vsr))
  expect_identical(body(vs), body(vsr))
  expect_equal(unsm(2), matrix(c(1,2,2,1), 2))
  expect_equal(unsm2(2), matrix(c(1,0,2,1), 2))
  expect_equal(fixm(2), matrix(3, 2, 2))
  expect_true("fixedSigma2" %in% names(formals(covm)))
})
