test_that("loadings_mmes and scores_mmes reconstruct the fitted FA/RR covariance", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(6L))),
    genotype=factor(paste0("G", seq_len(24L)))
  )

  trial_effect <- c(45, -30, -12, 52, 18, -20)
  genotype_effect <- sin(seq_len(24L) / 3) * 8
  interaction <- cos(seq_len(nrow(data)) / 5) * 3
  data$BLUEs <- 100 + trial_effect[data$trial] +
    genotype_effect[data$genotype] + interaction

  fit_fa <- mmes(
    BLUEs ~ trial,
    random=~vsm(fam(trial, 2), ism(genotype)),
    rcov=~units,
    data=data,
    nIters=20,
    verbose=FALSE
  )

  term_fa <- names(fit_fa$uList)[1]
  fa <- loadings_mmes(fit_fa, term_fa, varianceScale = FALSE, rotation = FALSE)

  expect_equal(dim(fa$loadings), c(6L, 2L))
  expect_length(fa$specific, 6L)

  reconstructed <- fa$sigma2 * (
    fa$loadings %*% t(fa$loadings) + diag(fa$specific)
  )
  expect_equal(unname(reconstructed), unname(fit_fa$theta[[term_fa]]), tolerance=1e-6)

  scores_fa <- scores_mmes(fit_fa, term_fa, varianceScale = FALSE, rotation = FALSE)
  expect_equal(dim(scores_fa), c(24L, 2L))
  expect_true(all(is.finite(scores_fa)))

  scores_fa_bartlett <- scores_mmes(fit_fa, term_fa, method="bartlett", varianceScale = FALSE, rotation = FALSE)
  expect_equal(dim(scores_fa_bartlett), c(24L, 2L))
  expect_true(all(is.finite(scores_fa_bartlett)))

  fit_rr <- mmes(
    BLUEs ~ trial,
    random=~vsm(rrcm(trial, 1), ism(genotype)),
    rcov=~units,
    data=data,
    nIters=20,
    verbose=FALSE
  )

  term_rr <- names(fit_rr$uList)[1]
  rr <- loadings_mmes(fit_rr, term_rr, varianceScale = FALSE, rotation = FALSE)

  expect_equal(dim(rr$loadings), c(6L, 1L))
  expect_true(all(rr$specific > 0))

  reconstructed_rr <- rr$sigma2 * (
    rr$loadings %*% t(rr$loadings) + diag(rr$specific)
  )
  expect_equal(unname(reconstructed_rr), unname(fit_rr$theta[[term_rr]]), tolerance=1e-6)

  scores_rr <- scores_mmes(fit_rr, term_rr, varianceScale = FALSE, rotation = FALSE)
  expect_equal(dim(scores_rr), c(24L, 1L))
  expect_true(all(is.finite(scores_rr)))
})

test_that("loadings_mmes gives an informative error for unsupported terms", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(4L))),
    genotype=factor(paste0("G", seq_len(10L)))
  )
  data$BLUEs <- rnorm(nrow(data))

  fit <- mmes(
    BLUEs ~ trial,
    random=~vsm(ism(genotype)),
    rcov=~units,
    data=data,
    nIters=2,
    verbose=FALSE
  )

  expect_error(loadings_mmes(fit, varianceScale = FALSE, rotation = FALSE), "No fam\\(\\)/rrcm\\(\\) covariance term")
})
