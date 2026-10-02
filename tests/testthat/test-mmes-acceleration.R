test_that("the default EM schedule uses the validated twenty-step taper", {
  data(DT_example)
  setup <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
    nIters=30, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
  expected <- c(exp(seq(log(1), log(0.03), length.out=20L)), rep(0.03, 10L))
  expect_equal(setup$emWeight, expected)
})

test_that("guarded Aitken acceleration preserves the converged deterministic fit", {
  data(DT_example)
  fit <- function(acceleration){
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
      nIters=80, tolParConvLL=1e-10, tolParConvNorm=1e-10,
      acceleration=acceleration, verbose=FALSE, dateWarning=FALSE)
  }
  plain <- fit("none")
  accelerated <- fit("aitken")
  expect_equal(plain$engineDiagnostics$acceleratedProposals, 0L)
  expect_gt(accelerated$engineDiagnostics$acceleratedProposals, 0L)
  expect_true(all(diff(as.numeric(accelerated$llik)) >= -1e-8))
  expect_equal(unlist(accelerated$theta), unlist(plain$theta), tolerance=1e-5)
  expect_equal(as.numeric(accelerated$bu), as.numeric(plain$bu), tolerance=1e-5)
  expect_gte(accelerated$engineDiagnostics$coefficientEvaluations,
    ncol(accelerated$llik) + accelerated$engineDiagnostics$lineSearchHalvings)
})

test_that("acceleration rejects unsupported solver modes", {
  data(DT_example)
  expect_error(mmes(Yield~Env, random=~Name, data=DT_example,
    acceleration="aitken", solver="pcg"), "deterministic")
  expect_error(mmes(Yield~Env, random=~Name, data=DT_example,
    acceleration="aitken", henderson=FALSE), "deterministic")
})