test_that("the default EM schedule uses a thirteen-step taper", {
  data(DT_example)
  for(nIters in c(1L, 5L, 13L, 30L)){
    setup <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
      nIters=nIters, returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)
    taperIters <- min(nIters, 13L)
    expected <- if(nIters == 1L) 1 else
      c(exp(seq(log(1), log(0.03), length.out=taperIters)),
        rep(0.03, nIters - taperIters))
    expect_equal(setup$emWeight, expected)
  }
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