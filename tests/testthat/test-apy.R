apy_relationship <- function(){
  coreCovariance <- matrix(c(2, 0.4, 0.4, 1.5), 2L)
  regression <- matrix(c(0.2, -0.1, 0.4, 0.3, -0.2, 0.5), 3L, 2L)
  conditional <- matrix(c(0.6, 0.08, 0.03, 0.08, 0.5, -0.04,
                          0.03, -0.04, 0.7), 3L)
  G <- rbind(
    cbind(coreCovariance, coreCovariance %*% t(regression)),
    cbind(regression %*% coreCovariance,
          regression %*% coreCovariance %*% t(regression) + conditional)
  )
  dimnames(G) <- list(paste0("id", seq_len(5L)), paste0("id", seq_len(5L)))
  list(G=G, regression=regression, conditional=conditional)
}

test_that("APY precision matches the core and cross-block construction", {
  example <- apy_relationship()
  result <- APY(example$G, core=c("id2", "id1"), return.details=TRUE)
  Q <- as.matrix(result$Gu)
  core <- result$core
  nonCore <- result$nonCore
  B <- result$B
  conditionalVariance <- result$conditionalVariance

  expected <- matrix(0, 5L, 5L)
  expected[core, core] <- solve(example$G[core, core, drop=FALSE]) +
    crossprod(B, B / conditionalVariance)
  expected[core, nonCore] <- -t(B / conditionalVariance)
  expected[nonCore, core] <- t(expected[core, nonCore, drop=FALSE])
  diag(expected)[nonCore] <- 1 / conditionalVariance
  dimnames(expected) <- dimnames(example$G)
  expect_equal(Q, expected, tolerance=1e-12)
  expect_true(inherits(result$Gu, "CsparseMatrix"))
  expect_true(isTRUE(attr(result$Gu, "inverse")))
  expect_identical(rownames(result$Gu), rownames(example$G))
  expect_identical(attr(result$Gu, "APY")$core, c("id2", "id1"))
  reconstructed <- solve(Q)
  expect_equal(reconstructed[core, core], example$G[core, core], tolerance=1e-12)
  expect_equal(reconstructed[core, nonCore], example$G[core, nonCore], tolerance=1e-12)
  expect_equal(diag(reconstructed), diag(example$G), tolerance=1e-12)
  expect_equal(unname(reconstructed[nonCore, nonCore] -
    reconstructed[nonCore, core, drop=FALSE] %*%
      solve(reconstructed[core, core, drop=FALSE]) %*%
      reconstructed[core, nonCore, drop=FALSE]),
    unname(diag(result$conditionalVariance)), tolerance=1e-12)
})

test_that("APY is exact when non-core conditional covariance is diagonal", {
  Gcc <- matrix(c(2, 0.4, 0.4, 1.5), 2L)
  B <- matrix(c(0.2, -0.1, 0.4, 0.3, -0.2, 0.5), 3L, 2L)
  residual <- diag(c(0.6, 0.5, 0.7))
  G <- rbind(cbind(Gcc, Gcc %*% t(B)),
             cbind(B %*% Gcc, B %*% Gcc %*% t(B) + residual))
  dimnames(G) <- list(paste0("g", 1:5), paste0("g", 1:5))
  precision <- APY(G, core=1:2)
  expect_equal(solve(as.matrix(precision)), G, tolerance=1e-12)
})

test_that("APY precision can be used directly as Gu", {
  example <- apy_relationship()
  precision <- APY(example$G, core=1:2)
  set.seed(19)
  data <- data.frame(id=factor(rep(rownames(example$G), each=3L),
                               levels=rownames(example$G)), y=rnorm(15L))
  fit <- mmes(y~1, random=~vsm(ism(id), Gu=precision), rcov=~units,
              data=data, nIters=3, verbose=FALSE, dateWarning=FALSE)
  expect_s3_class(fit, "mmes")
  expect_true(all(is.finite(fit$theta[[1L]])))
  expect_true(all(is.finite(fit$uList[[1L]])))
})

test_that("APY validates the relationship, core and conditional variances", {
  example <- apy_relationship()
  expect_error(APY(example$G, core=character()), "core")
  expect_error(APY(example$G, core=c(1, 1)), "unique")
  expect_error(APY(example$G, core=6), "indices")
  nonsymmetric <- example$G
  nonsymmetric[1L, 2L] <- nonsymmetric[1L, 2L] + 0.01
  expect_error(APY(nonsymmetric, core=1:2, tol=1e-8), "symmetric")
  expect_error(APY(example$G, core=1:2, tol=0), "tol")

  singular <- matrix(c(1, 1, 1, 1), 2L,
                     dimnames=list(c("a", "b"), c("a", "b")))
  expect_error(APY(singular, core="a"), "conditional variances")
})