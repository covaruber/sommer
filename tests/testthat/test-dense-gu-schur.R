test_that("multiple dense-Gu effects use exact Schur blocks", {
  set.seed(17)
  nInd <- 18
  ids <- paste0("id", seq_len(nInd))
  markers <- matrix(sample(c(-1, 0, 1), nInd * 40, replace=TRUE), nrow=nInd)
  G1 <- tcrossprod(scale(markers)) / ncol(markers) + diag(0.2, nInd)
  G2 <- 0.7 * G1 + 0.3 * diag(nInd)
  Ginv1 <- solve(G1)
  Ginv2 <- solve(G2)
  dimnames(Ginv1) <- dimnames(Ginv2) <- list(ids, ids)
  attr(Ginv1, "inverse") <- TRUE
  attr(Ginv2, "inverse") <- TRUE

  dat <- expand.grid(id=ids, env=paste0("E", 1:3), KEEP.OUT.ATTRS=FALSE)
  dat$id <- factor(dat$id, levels=ids)
  dat$env <- factor(dat$env)
  dat$y <- rnorm(nrow(dat)) + as.numeric(G1 %*% rnorm(nInd))[as.integer(dat$id)]
  random <- ~vsm(dsm(env), ism(id), Gu=Ginv1) +
    vsm(dsm(env), ism(id), Gu=Ginv2)

  cholmodOutput <- capture.output(cholmodFit <- mmes(
    y~env, random=random, rcov=~units, data=dat, solver="cholmod",
    nIters=8, verbose=TRUE, dateWarning=FALSE
  ))
  ldltFit <- suppressMessages(mmes(
    y~env, random=random, rcov=~units, data=dat, solver="ldlt",
    nIters=8, verbose=FALSE, dateWarning=FALSE
  ))

  expect_true(any(grepl("Dense block-Schur engine active \\(3 random-effect blocks, 3 border effects\\)", cholmodOutput)))
  expect_false(any(grepl("Materialized cross-group inverse block", cholmodOutput)))
  expect_equal(unlist(cholmodFit$theta), unlist(ldltFit$theta), tolerance=1e-8)
  expect_equal(cholmodFit$u, ldltFit$u, tolerance=1e-8)
  expect_true(cholmodFit$engineDiagnostics$blockSchurActive)
  expect_equal(cholmodFit$engineDiagnostics$blockSchurGroups, 3L)

  mlSchur <- suppressMessages(mmes(y~env, random=random, rcov=~units,
    data=dat, solver="cholmod", REML=FALSE, nIters=8,
    verbose=FALSE, dateWarning=FALSE))
  mlLdlt <- suppressMessages(mmes(y~env, random=random, rcov=~units,
    data=dat, solver="ldlt", REML=FALSE, nIters=8,
    verbose=FALSE, dateWarning=FALSE))
  expect_true(mlSchur$engineDiagnostics$blockSchurActive)
  expect_equal(unlist(mlSchur$theta), unlist(mlLdlt$theta), tolerance=1e-8)
  expect_equal(as.numeric(mlSchur$llik), as.numeric(mlLdlt$llik), tolerance=1e-8)
  expect_equal(mlSchur$bu, mlLdlt$bu, tolerance=1e-8)
})