test_that("predict.mmes requires the complete coefficient-matrix inverse", {
  data(DT_example)
  DT <- DT_example

  m0 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, data=DT)
  expect_identical(m0$CiMode, 0L)
  expect_false(isTRUE(m0$CiComputed))
  expect_error(predict(m0, D="Env"), "complete coefficient-matrix inverse")

  m1 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, computeCi=1, data=DT)
  expect_identical(m1$CiMode, 1L)
  expect_false(isTRUE(m1$CiComputed))
  expect_error(predict(m1, D="Env"), "complete coefficient-matrix inverse")

  m2 <- mmes(Yield~Env, random=~Name, rcov=~units, nIters=3,
             verbose=FALSE, computeCi=2, data=DT)
  expect_identical(m2$CiMode, 2L)
  expect_true(isTRUE(m2$CiComputed))
  p2 <- predict(m2, D="Env")
  expect_true(all(p2$pvals$std.error > 0))

  m1p <- postPEV(m1, mode=2)
  expect_identical(m1p$CiMode, 2L)
  expect_true(isTRUE(m1p$CiComputed))
  p1p <- predict(m1p, D="Env")
  expect_equal(p1p$pvals$predicted.value, p2$pvals$predicted.value)
  expect_equal(p1p$pvals$std.error, p2$pvals$std.error)
})
