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

  envColumns <- as.vector(m0$partitionsX[["Env"]])
  interceptColumn <- as.vector(m0$partitionsX[["1"]])
  expect_equal(unname(as.matrix(p0n$D[, envColumns, drop=FALSE])),
               matrix(1 / (length(envColumns) + 1L),
                      nrow=nrow(p0n$D), ncol=length(envColumns)))
  expect_equal(as.numeric(p0n$D[, interceptColumn]), rep(1, nrow(p0n$D)))
})

test_that("predict.mmes averages fixed interactions over their full level space", {
  data(DT_yatesoats)

  model <- mmes(Y ~ B + B:MP, random=~V, rcov=~units, nIters=1,
                verbose=FALSE, data=DT_yatesoats)
  prediction <- predict(model, D="V")

  bColumns <- as.vector(model$partitionsX[["B"]])
  interactionColumns <- as.vector(model$partitionsX[["B:MP"]])
  randomRange <- model$partitions[[1L]][1L, ]
  randomColumns <- seq.int(randomRange[1L], randomRange[2L])

  expect_equal(unique(as.numeric(prediction$D[, bColumns])), 1 / 6)
  expect_equal(unique(as.numeric(prediction$D[, interactionColumns])), 1 / 18)
  expect_equal(unname(as.matrix(prediction$D[, randomColumns])), diag(3))
})
