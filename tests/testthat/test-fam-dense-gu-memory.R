test_that("factor analytic random effects with dense Gu preserve exact fits", {
  set.seed(15)
  n <- 18L
  data <- expand.grid(envf=factor(seq_len(4L)), id=factor(seq_len(n)),
                      rep=factor(seq_len(2L)))
  data$repf <- interaction(data$envf, data$rep, drop=TRUE)
  data$y <- rnorm(nrow(data))
  markers <- matrix(rnorm(n * 55L), n)
  relationship <- solve(tcrossprod(markers) / 55 + diag(0.5, n))
  dimnames(relationship) <- list(levels(data$id), levels(data$id))
  attr(relationship, "inverse") <- TRUE

  fit <- suppressMessages(mmes(
    y~envf,
    random=~vsm(fam(envf, 1), ism(id), Gu=relationship) +
      vsm(dsm(envf), ism(repf)),
    rcov=~units, data=data, solver="cholmod", nIters=2,
    verbose=FALSE, dateWarning=FALSE
  ))

  expect_equal(as.numeric(fit$llik), c(-82.92354132851109, -81.10458105294998),
               tolerance=1e-7)
  expect_equal(head(as.numeric(fit$bu), 4L),
               c(-0.01976577434524176, 0.31310560257702447,
                 -0.14463758346058719, -0.10188982610961937), tolerance=1e-7)
})
test_that("unused factor levels do not enter covariance or fixed designs", {
  set.seed(3)
  n <- 12L
  data <- expand.grid(envf=factor(1:3, levels=1:8), id=factor(seq_len(n)))
  data$y <- rnorm(nrow(data))
  relationship <- diag(n)
  dimnames(relationship) <- list(levels(data$id), levels(data$id))
  attr(relationship, "inverse") <- TRUE

  expect_equal(ncol(sommer:::rrm(data$envf, 2)$Z), 3L)
  p <- suppressMessages(mmes(y~envf, random=~vsm(rrm(envf, 2), ism(id), Gu=relationship),
                             data=data, returnParam=TRUE, dateWarning=FALSE))
  expect_equal(ncol(p$X), 3L)
  expect_equal(p$covStruct[[1]]$dim, 3L)
})
