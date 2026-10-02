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
  expect_equal(fitMl$engineDiagnostics$DsymbolicAnalyses, 1L)
  expect_equal(fitMl$engineDiagnostics$DselectedTopologyBuilds, 1L)
  expect_gt(fitMl$engineDiagnostics$DselectedTopologyReuses, 0L)

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

test_that("matrix-free PCG handles multi-factor Kronecker random precision", {
  set.seed(17)
  dat <- expand.grid(
    environment=factor(seq_len(3)),
    trait=factor(seq_len(3)),
    id=factor(seq_len(8))
  )
  dat$y <- rnorm(nrow(dat))

  fit <- function(solver, nIters=5){
    mmes(
      y ~ 1,
      random=~vsm(csm(environment, rho=0.15), dsm(trait), ism(id)),
      rcov=~units,
      data=dat,
      nIters=nIters,
      verbose=FALSE,
      getPEV=FALSE,
      computeCi=0,
      solver=solver,
      pcgTol=1e-10,
      pcgTraceProbes=24,
      pcgLanczosSteps=30
    )
  }

  pcg1 <- fit("pcg")
  pcg2 <- fit("pcg")
  pcgOne <- fit("pcg", nIters=1)
  ldltOne <- fit("ldlt", nIters=1)

  expect_equal(pcg1$covPar, pcg2$covPar, tolerance=1e-12)
  expect_equal(pcg1$bu, pcg2$bu, tolerance=1e-12)
  expect_true(pcg1$pcgMatrixFree)
  expect_false(ldltOne$pcgMatrixFree)
  expect_equal(pcgOne$bu, ldltOne$bu, tolerance=1e-7)
  expect_true(inherits(pcg1$C, "sparseMatrix"))
  expect_equal(dim(pcg1$C), c(length(pcg1$bu), length(pcg1$bu)))
  expect_gt(Matrix::nnzero(pcg1$C), Matrix::nnzero(Matrix::crossprod(pcg1$W)))

  pcgLogLik <- mmes(
    y ~ 1,
    random=~vsm(csm(environment, rho=0.15), dsm(trait), ism(id)),
    rcov=~units,
    data=dat,
    nIters=1,
    verbose=FALSE,
    getPEV=FALSE,
    computeCi=0,
    solver="pcg",
    pcgTol=1e-10,
    pcgTraceProbes=128,
    pcgLanczosSteps=72
  )
  exactLogLik <- fit("ldlt", nIters=1)
  expect_equal(as.numeric(pcgLogLik$llik), as.numeric(exactLogLik$llik),
               tolerance=0.2)
})

test_that("Nyström PCG preconditioning preserves dense-Gu fits and reduces iterations", {
  set.seed(11)
  nId <- 40L
  ids <- paste0("i", seq_len(nId))
  markers <- matrix(sample(c(-1, 0, 1), nId * 120L, replace=TRUE), nrow=nId)
  relationship <- tcrossprod(scale(markers)) / ncol(markers) + diag(0.2, nId)
  precision <- solve(relationship)
  dimnames(precision) <- list(ids, ids)
  attr(precision, "inverse") <- TRUE
  data <- data.frame(id=factor(rep(ids, each=3L), levels=ids), y=rnorm(nId * 3L))
  fit <- function(preconditioner){
    mmes(y~1, random=~vsm(ism(id), Gu=precision), rcov=~units,
      data=data, nIters=3, solver="pcg", computeCi=2,
      pcgPreconditioner=preconditioner, pcgNystromRank=16L,
      pcgTraceProbes=32L, pcgLanczosSteps=30L,
      verbose=FALSE, dateWarning=FALSE)
  }
  diagonal <- fit("diagonal")
  nystrom <- fit("nystrom")
  expect_false(diagonal$pcgMatrixFree)
  expect_false(nystrom$pcgMatrixFree)
  expect_identical(nystrom$engineDiagnostics$pcgPreconditioner, "nystrom")
  expect_equal(nystrom$engineDiagnostics$pcgNystromRank, 16L)
  expect_lt(nystrom$engineDiagnostics$pcgBatchIterations,
            diagonal$engineDiagnostics$pcgBatchIterations)
  expect_equal(as.numeric(nystrom$llik), as.numeric(diagonal$llik), tolerance=1e-7)
  expect_equal(as.numeric(nystrom$bu), as.numeric(diagonal$bu), tolerance=1e-7)
})
