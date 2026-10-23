test_that("legacy multivariate AI retains CRAN 4.4.87 estimates", {
  set.seed(815)
  legacyLevels <- sprintf("line%02d", seq_len(24L))
  legacyKinship <- toeplitz(0.6^(seq_along(legacyLevels) - 1L))
  dimnames(legacyKinship) <- list(legacyLevels, legacyLevels)
  legacyData <- data.frame(genotype=factor(rep(legacyLevels, each=2L)))
  legacyGenetic <- t(chol(legacyKinship)) %*% matrix(rnorm(48L), 24L, 2L) %*%
    chol(matrix(c(1.2, 0.4, 0.4, 0.8), 2L))
  legacyResidual <- matrix(rnorm(96L), 48L, 2L) %*%
    chol(matrix(c(0.5, 0.15, 0.15, 0.7), 2L))
  legacyData$yield <- 3 + legacyGenetic[legacyData$genotype, 1L] +
    legacyResidual[, 1L]
  legacyData$height <- 8 + legacyGenetic[legacyData$genotype, 2L] +
    legacyResidual[, 2L]

  legacyFit <- suppressMessages(mmer(
    cbind(yield, height) ~ 1,
    random=~vsr(genotype, Gu=legacyKinship, Gtc=unsm(2L)),
    rcov=~vsr(units, Gtc=unsm(2L)),
    data=legacyData, method="AI", nIters=8L,
    tolParConvLL=1e-10, dateWarning=FALSE, verbose=FALSE
  ))

  expect_equal(unname(unlist(legacyFit$sigma)),
    c(0.870056043197485, 0.0887502084261768, 0.0887502084261768,
      0.392046234562069, 0.680945760698956, 0.188121084047432,
      0.188121084047432, 0.497428809879275), tolerance=1e-7)
  expect_equal(unname(legacyFit$Beta$Estimate),
    c(2.62174876871535, 7.90364716499132), tolerance=1e-7)
})
