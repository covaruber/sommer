test_that("mmes accepts namespace-qualified covariance calls", {
  data("DT_example", package="enhancer")

  fit <- sommer::mmes(
    Yield ~ Env,
    random=~sommer::vsm(sommer::dsm(Env), sommer::ism(Name)),
    rcov=~sommer::vsm(sommer::dsm(Env), sommer::ism(units)),
    nIters=2,
    verbose=FALSE,
    data=DT_example
  )

  expect_s3_class(fit, "mmes")
  expect_identical(fit$solver, "ldlt")
  expect_length(fit$covStruct, 2L)
})
