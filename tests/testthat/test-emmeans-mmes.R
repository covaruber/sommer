test_that("emmeans obtains mmes fixed-effect means and covariance", {
  skip_if_not_installed("emmeans")
  data(DT_example)

  fit <- mmes(Yield ~ Env, random=~Name, rcov=~units, data=DT_example,
              nIters=3, verbose=FALSE)
  recovered <- emmeans::recover_data(fit)
  grid <- emmeans::ref_grid(fit)
  means <- emmeans::emmeans(fit, ~ Env)
  summaryMeans <- summary(means)
  fixedContrast <- as.matrix(grid@linfct)
  D <- cbind(Matrix::Diagonal(ncol(fixedContrast)),
             Matrix::Matrix(0, nrow=ncol(fixedContrast),
                            ncol=nrow(fit$bu) - ncol(fixedContrast), sparse=TRUE))
  fixedVcov <- sommer:::predict_mmes_vcov_cpp(fit, D)
  expectedVcov <- fixedContrast %*% fixedVcov %*% t(fixedContrast)

  expect_true(is.data.frame(recovered))
  expect_true("Env" %in% names(recovered))
  expect_equal(summaryMeans$emmean, as.vector(fixedContrast %*% fit$b))
  expect_equal(summaryMeans$SE, unname(sqrt(diag(expectedVcov))))
  expect_true(all(is.infinite(summaryMeans$df)))
})

test_that("emmeans support rejects multivariate mmes fits", {
  skip_if_not_installed("emmeans")
  data(DT_example)
  fit <- mmes(Yield ~ Env, random=~Name, rcov=~units, data=DT_example,
              nIters=1, verbose=FALSE)
  fit$y <- cbind(fit$y, fit$y)

  expect_error(emmeans::emmeans(fit, ~ Env), "univariate response")
})