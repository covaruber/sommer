test_that("mmes accepts namespace-qualified covariance calls", {
  data("DT_example", package="enhancer")

  fit <- sommer::mmes(Yield ~ Env, random=~sommer::vsm(sommer::dsm(Env), sommer::ism(Name)), rcov=~sommer::vsm(sommer::dsm(Env), sommer::ism(units)), data=DT_example, control=mmes_control(solve=list(method="ldlt"), convergence=list(maxIters=2), diagnostics=list(verbose=FALSE)))

  expect_s3_class(fit, "mmes")
  expect_identical(fit$control$solve$method, "ldlt")
  expect_length(fit$covStruct, 2L)
})
