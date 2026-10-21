test_that("documented enhancer model families retain numerical solver parity", {
  for(arguments in .enhancer_regression_models()){
    fit <- function(method) do.call(mmes, c(arguments, list(control=mmes_control(
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
