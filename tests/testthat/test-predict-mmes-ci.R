test_that("predict.mmes gives exact SEs without requiring the full Ci inverse", {
  data(DT_example)
  DT <- DT_example

  m0 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, data=DT)
  expect_identical(m0$CiMode, 0L)
  expect_false(isTRUE(m0$CiComputed))

  m1 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, computeCi=1, data=DT)
  expect_identical(m1$CiMode, 1L)
  expect_false(isTRUE(m1$CiComputed))

  m2 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, computeCi=2, data=DT)
  expect_identical(m2$CiMode, 2L)
  expect_true(isTRUE(m2$CiComputed))

  # predict_mmes_vcov_cpp() solves C %*% X = t(D) for the requested rows
  # instead of using Ci, so all three computeCi modes must agree exactly
  # with the full-inverse reference, for both a fixed-effect (Env) and a
  # random-effect (Name) classify.
  p0 <- predict(m0, D="Env")
  p1 <- predict(m1, D="Env")
  p2 <- predict(m2, D="Env")
  expect_true(all(p2$pvals$std.error > 0))
  expect_equal(p0$pvals$predicted.value, p2$pvals$predicted.value)
  expect_equal(p0$pvals$std.error, p2$pvals$std.error)
  expect_equal(p1$pvals$predicted.value, p2$pvals$predicted.value)
  expect_equal(p1$pvals$std.error, p2$pvals$std.error)

  p0n <- predict(m0, D="Name")
  p2n <- predict(m2, D="Name")
  expect_equal(p0n$pvals$std.error, p2n$pvals$std.error)
})
