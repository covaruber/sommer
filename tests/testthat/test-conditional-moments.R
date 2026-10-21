test_that("conditional moments reuse symbolic analysis without reusing numeric values", {
  H <- Matrix::bandSparse(9L, k=c(-1L, 0L, 1L),
    diagonals=list(rep(-0.2, 8L), rep(2, 9L), rep(-0.2, 8L)))
  pattern <- H
  pattern@x[] <- 1
  rhs <- seq_len(9L) / 7
  first <- sommer:::mme_conditional_moments_cpp(H, pattern, rhs, NULL,
    TRUE, 3L, 10)
  expected <- solve(as.matrix(H))
  entries <- summary(pattern)
  expect_equal(as.numeric(first$mean), as.numeric(expected %*% rhs), tolerance=1e-12)
  expect_equal(first$covariance[cbind(entries$i, entries$j)],
    expected[cbind(entries$i, entries$j)], tolerance=1e-12)
  expect_equal(first$logDet, as.numeric(determinant(H, logarithm=TRUE)$modulus))
  changed <- H + Matrix::Diagonal(9L)
  second <- sommer:::mme_conditional_moments_cpp(changed, pattern, rhs,
    first$workspace, TRUE, 2L, 10)
  expect_equal(second$symbolicAnalyses, 1L)
  expect_equal(second$numericFactorizations, 2L)
  expect_equal(as.numeric(second$mean), as.numeric(solve(changed, rhs)), tolerance=1e-12)
  expect_equal(second$panelColumns, 2L)
  expect_error(sommer:::mme_conditional_moments_cpp(H, pattern, rhs, NULL,
    TRUE, 9L, 1e-6), "exceed")
})

test_that("conditional workspace invalidates changed patterns and admits mean-only solves", {
  H <- as(Matrix::Diagonal(4L, 2), "generalMatrix")
  pattern <- H
  rhs <- seq_len(4L)
  first <- sommer:::mme_conditional_moments_cpp(H, pattern, rhs, NULL,
    FALSE, 4L, 1e-6)
  expect_null(first$covariance)
  expect_equal(first$panelBytes, 0)
  changed <- H
  changed[1L,2L] <- changed[2L,1L] <- .1
  expanded <- changed
  expanded@x[] <- 1
  second <- sommer:::mme_conditional_moments_cpp(changed, expanded, rhs,
    first$workspace, TRUE, 2L, 10)
  expect_equal(second$symbolicAnalyses, 2L)
  expect_equal(as.numeric(second$mean), as.numeric(solve(changed,rhs)), tolerance=1e-12)
  expect_error(sommer:::mme_conditional_moments_cpp(changed, pattern, rhs, NULL,
    TRUE, 2L, 10), "omits precision")
  expect_error(sommer:::mme_conditional_moments_cpp(-H, pattern, rhs, NULL,
    TRUE, 2L, 10), "not positive definite")
})
