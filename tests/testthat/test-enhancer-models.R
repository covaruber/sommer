test_that("documented enhancer model families retain numerical solver parity", {
  skip_on_cran()
  for(arguments in .enhancer_regression_models()){
    fit <- function(method) do.call(sommer:::.mmes_dispatch, c(arguments, list(control=sommer:::.mmes_control_internal(
      solve=list(method=method), convergence=list(maxIters=5L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))))
    reference <- fit("ldlt")
    candidate <- fit("cholmod")
    expect_equal(candidate$llik, reference$llik, tolerance=1e-6)
    expect_equal(candidate$bu, reference$bu, tolerance=1e-6)
    expect_equal(candidate$covPar, reference$covPar, tolerance=1e-6)
    expect_true(all(is.finite(fitted(candidate))))
  }
})
