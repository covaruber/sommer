set.seed(11)
n <- 60; m <- 300
X <- sapply(runif(m, 0.1, 0.9), function(p) rbinom(n, 2, p)) - 1
X <- X[, apply(X, 2, var) > 0]
p <- colMeans(X + 1) / 2; q <- 1 - p; pq2 <- 2 * p * q

test_that("D.mat nishio uses the Vitezica orthogonal coding", {
  W <- matrix(NA_real_, nrow(X), ncol(X))
  for (k in seq_len(ncol(X))) {
    W[, k] <- c(-2 * p[k]^2, pq2[k], -2 * q[k]^2)[X[, k] + 2]
  }
  Dref <- tcrossprod(W) / sum(pq2^2)
  expect_equal(unname(D.mat(X, nishio = TRUE)), Dref, tolerance = 1e-10)
})

test_that("D.mat Su et al. (2012) is unchanged", {
  M <- sweep(1 - abs(X), 2, pq2)
  Dref <- tcrossprod(M) / sum(pq2 * (1 - pq2))
  expect_equal(unname(D.mat(X, nishio = FALSE)), Dref, tolerance = 1e-10)
})

test_that("D.mat mean-imputes missing markers", {
  Xna <- X
  Xna[cbind(c(1, 5, 9), c(2, 7, 11))] <- NA
  res <- suppressMessages(capture.output(out <- D.mat(Xna, return.imputed = TRUE)))
  expect_false(anyNA(out$X))
  expect_equal(out$X[1, 2], mean(Xna[, 2], na.rm = TRUE))
  expect_equal(out$D, D.mat(out$X))
})
