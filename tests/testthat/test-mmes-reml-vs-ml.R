test_that("REML=FALSE gives ML estimates matching a dense brute-force optimizer", {
  data(DT_example)
  DT <- DT_example
  DT <- DT[!is.na(DT$Yield) & !is.na(DT$Env) & !is.na(DT$Name), ]

  X <- model.matrix(~Env, data=DT)
  Z <- model.matrix(~Name-1, data=DT)
  y <- DT$Yield
  n <- length(y)

  negLogLikML <- function(par){
    su2 <- exp(par[1]); se2 <- exp(par[2])
    V <- su2 * tcrossprod(Z) + se2 * diag(n)
    cholV <- chol(V)
    logdetV <- 2 * sum(log(diag(cholV)))
    Vi <- chol2inv(cholV)
    XtVi <- t(X) %*% Vi
    bhat <- solve(XtVi %*% X, XtVi %*% y)
    resid <- y - X %*% bhat
    quad <- as.numeric(t(resid) %*% Vi %*% resid)
    0.5 * (n * log(2 * pi) + logdetV + quad)
  }
  opt <- optim(c(log(5), log(8)), negLogLikML, method="BFGS",
               control=list(reltol=1e-12, maxit=500))
  refVarcomp <- exp(opt$par)

  fitReml <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                  nIters=40, tolParConvLL=1e-10, verbose=FALSE)
  fitMl <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                nIters=40, tolParConvLL=1e-10, verbose=FALSE, REML=FALSE)

  expect_true(fitReml$REML)
  expect_false(fitMl$REML)

  mlVarcomp <- c(fitMl$theta[[1]][1,1], fitMl$theta[[2]][1,1])
  expect_equal(mlVarcomp, refVarcomp, tolerance=1e-3)

  # ML estimates are biased downward relative to REML for this classical
  # one-way random-effects model.
  remlVarcomp <- c(fitReml$theta[[1]][1,1], fitReml$theta[[2]][1,1])
  expect_true(all(mlVarcomp < remlVarcomp))
})

test_that("REML=FALSE with no random effects recovers the classical SSE/n MLE", {
  data(DT_example)
  DT <- DT_example
  DT <- DT[!is.na(DT$Yield) & !is.na(DT$Env), ]
  n <- nrow(DT)
  p <- length(unique(DT$Env))
  sse <- sum(residuals(lm(Yield~Env, data=DT))^2)

  fitReml <- mmes(Yield~Env, rcov=~units, data=DT, nIters=20, verbose=FALSE)
  fitMl <- mmes(Yield~Env, rcov=~units, data=DT, nIters=20, verbose=FALSE, REML=FALSE)

  expect_equal(as.numeric(fitReml$theta[[1]]), sse / (n - p), tolerance=1e-4)
  expect_equal(as.numeric(fitMl$theta[[1]]), sse / n, tolerance=1e-4)
})

test_that("REML=FALSE validates its inputs and the solver restriction", {
  data(DT_example)
  DT <- DT_example

  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, nIters=2,
         verbose=FALSE, REML=NA),
    "REML must be a single TRUE/FALSE value"
  )
  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, nIters=2,
         verbose=FALSE, REML=FALSE, solver="pcg"),
    "REML=FALSE"
  )
  # cholmod is a supported REML=FALSE backend (unlike pcg above)
  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, nIters=2,
         verbose=FALSE, REML=FALSE, solver="cholmod"),
    NA
  )
})

test_that("REML=FALSE with solver='cholmod' matches solver='ldlt' exactly", {
  data(DT_example)
  DT <- DT_example
  DT <- DT[!is.na(DT$Yield) & !is.na(DT$Env) & !is.na(DT$Name), ]

  fLdlt <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                nIters=30, tolParConvLL=1e-10, verbose=FALSE,
                REML=FALSE, solver="ldlt")
  fChol <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                nIters=30, tolParConvLL=1e-10, verbose=FALSE,
                REML=FALSE, solver="cholmod")

  expect_equal(as.numeric(fLdlt$theta[[1]]), as.numeric(fChol$theta[[1]]),
               tolerance=1e-6)
  expect_equal(as.numeric(fLdlt$theta[[2]]), as.numeric(fChol$theta[[2]]),
               tolerance=1e-6)
  expect_equal(tail(fLdlt$llik, 1), tail(fChol$llik, 1), tolerance=1e-6)
})

test_that("REML=FALSE with a dense Gu (solver='cholmod') matches solver='ldlt'", {
  set.seed(1)
  nInd <- 40
  M <- matrix(sample(c(-1, 0, 1), nInd * 150, replace=TRUE), nrow=nInd)
  Gu <- tcrossprod(scale(M)) / ncol(M)
  Gu <- Gu + diag(1e-3, nInd)
  rownames(Gu) <- colnames(Gu) <- paste0("id", 1:nInd)
  GuInv <- solve(Gu)
  attr(GuInv, "inverse") <- TRUE

  DTg <- data.frame(id=factor(rep(paste0("id", 1:nInd), 3)))
  set.seed(2)
  u <- MASS::mvrnorm(1, mu=rep(0, nInd), Sigma=5 * Gu)
  DTg$y <- 10 + u[as.integer(DTg$id)] + rnorm(nrow(DTg), sd=2)

  fLdlt <- mmes(y~1, random=~vsm(ism(id), Gu=GuInv), rcov=~units, data=DTg,
                nIters=25, verbose=FALSE, REML=FALSE, solver="ldlt")
  fAuto <- mmes(y~1, random=~vsm(ism(id), Gu=GuInv), rcov=~units, data=DTg,
                nIters=25, verbose=FALSE, REML=FALSE)

  expect_equal(as.numeric(fLdlt$theta[[1]]), as.numeric(fAuto$theta[[1]]),
               tolerance=1e-6)
  expect_equal(as.numeric(fLdlt$theta[[2]]), as.numeric(fAuto$theta[[2]]),
               tolerance=1e-6)
})
