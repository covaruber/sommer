test_that("rank-two factor-analytic updates remain finite and positive definite", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(6L))),
    genotype=factor(paste0("G", seq_len(24L)))
  )

  trial_effect <- c(45, -30, -12, 52, 18, -20)
  genotype_effect <- sin(seq_len(24L) / 3) * 8
  interaction <- cos(seq_len(nrow(data)) / 5) * 3
  data$BLUEs <- 100 + trial_effect[data$trial] +
    genotype_effect[data$genotype] + interaction

  fit <- mmes(
    BLUEs ~ trial,
    random=~vsm(fam(trial, 2), ism(genotype)),
    rcov=~units,
    data=data,
    nIters=20,
    verbose=FALSE
  )

  expect_s3_class(fit, "mmes")
  expect_true(all(is.finite(fit$monitor)))
  expect_true(all(vapply(fit$covPar, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(vapply(
    fit$theta,
    function(x) min(eigen(x, symmetric=TRUE, only.values=TRUE)$values) > 0,
    logical(1)
  )))
})
