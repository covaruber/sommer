test_that("starting-value least squares matches dense QR for full-rank designs", {
  set.seed(5)
  X <- Matrix::sparse.model.matrix(~ g + x, data.frame(g = factor(sample(letters[1:6], 60, TRUE)), x = rnorm(60)))
  Y <- matrix(rnorm(60), ncol = 1)
  fast <- sommer:::.mmes_start_lsq(X, Y)
  qrX <- qr(as.matrix(X))
  expect_equal(fast$rank, qrX$rank)
  expect_equal(as.numeric(fast$beta), as.numeric(qr.coef(qrX, Y)), tolerance = 1e-8)
  expect_equal(as.numeric(fast$resid), as.numeric(Y - as.matrix(X) %*% qr.coef(qrX, Y)), tolerance = 1e-8)
})

test_that("starting-value least squares flags aliased columns like QR", {
  set.seed(6)
  g <- factor(sample(letters[1:4], 40, TRUE))
  X <- cbind(Matrix::sparse.model.matrix(~ g - 1), all = 1)
  fast <- sommer:::.mmes_start_lsq(X, matrix(rnorm(40), ncol = 1))
  expect_lt(fast$rank, ncol(X))
  expect_true(anyNA(fast$beta))
})
