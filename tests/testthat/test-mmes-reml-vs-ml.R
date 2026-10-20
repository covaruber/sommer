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

  fitReml <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=40, likelihoodTol=1e-10), diagnostics=list(verbose=FALSE)))
  fitMl <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=40, likelihoodTol=1e-10), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml")))

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

  fitReml <- mmes(Yield~Env, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=20), diagnostics=list(verbose=FALSE)))
  fitMl <- mmes(Yield~Env, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=20, likelihoodTol=1e-10, parameterTol=1e-10), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml")))

  expect_equal(as.numeric(fitReml$theta[[1]]), sse / (n - p), tolerance=1e-4)
  expect_equal(as.numeric(fitMl$theta[[1]]), sse / n, tolerance=1e-4)
})

test_that("REML=FALSE validates its inputs and the solver restriction", {
  data(DT_example)
  DT <- DT_example

  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=2), diagnostics=list(verbose=FALSE), task=list(objective=NA_character_))),
    "task\\$objective"
  )
  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=2), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml"), solve=list(method="pcg"))),
    "task\\$objective"
  )
  # cholmod is a supported REML=FALSE backend (unlike pcg above)
  expect_error(
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT, control=mmes_control(convergence=list(maxIters=2), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml"), solve=list(method="cholmod"))),
    NA
  )
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

  fLdlt <- mmes(y~1, random=~vsm(ism(id), Gu=GuInv), rcov=~units, data=DTg, control=mmes_control(convergence=list(maxIters=25), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml"), solve=list(method="ldlt")))
  fAuto <- mmes(y~1, random=~vsm(ism(id), Gu=GuInv), rcov=~units, data=DTg, control=mmes_control(convergence=list(maxIters=25), diagnostics=list(verbose=FALSE), task=list(objective=if(FALSE) "reml" else "ml")))

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
    mmes(y ~ 1, random=~vsm(csm(environment, rho=0.15), dsm(trait), ism(id)), rcov=~units, data=dat, control=mmes_control(convergence=list(maxIters=nIters), diagnostics=list(verbose=FALSE), output=list(inverse=0), solve=c(list(method=solver), if(solver == "pcg") list(tol=1e-10)), score=if(solver == "pcg") list(traceProbes=24) else list(), logdet=if(solver == "pcg") list(lanczosSteps=30) else list()))
  }

  pcg1 <- fit("pcg")
  pcg2 <- fit("pcg")
  pcgOne <- fit("pcg", nIters=1)
  ldltOne <- fit("ldlt", nIters=1)

  expect_equal(pcg1$covPar, pcg2$covPar, tolerance=1e-12)
  expect_equal(pcg1$bu, pcg2$bu, tolerance=1e-12)
  expect_true(is.finite(pcg1$pcgLogDetProbeStandardError))
  expect_lte(pcg1$pcgSolveMaxRelativeResidual, 1.01e-10)
  expect_true(pcg1$pcgMatrixFree)
  expect_false(ldltOne$pcgMatrixFree)
  expect_equal(pcgOne$bu, ldltOne$bu, tolerance=1e-7)
  expect_true(inherits(pcg1$C, "sparseMatrix"))
  expect_equal(dim(pcg1$C), c(length(pcg1$bu), length(pcg1$bu)))
  expect_gt(Matrix::nnzero(pcg1$C), Matrix::nnzero(Matrix::crossprod(pcg1$W)))

  pcgLogLik <- mmes(y ~ 1, random=~vsm(csm(environment, rho=0.15), dsm(trait), ism(id)), rcov=~units, data=dat, control=mmes_control(convergence=list(maxIters=1), diagnostics=list(verbose=FALSE), output=list(inverse=0), solve=list(method="pcg", tol=1e-10), score=list(traceProbes=128), logdet=list(lanczosSteps=72)))
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
    mmes(y~1, random=~vsm(ism(id), Gu=precision), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=3), solve=list(method="pcg", preconditioner=preconditioner, nystromRank=16L), output=list(inverse=2), score=list(traceProbes=32L), logdet=list(lanczosSteps=30L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
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

test_that("PCG contracts correlated complete and ragged residual derivatives through probes", {
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(seq_len(4L)),
    replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 4
  for(missing in c(FALSE, TRUE)){
    active <- if(missing) dataset[-c(2L, 7L, 15L), ] else dataset
    fitModel <- function(solver){
      mmes(response~environment, random=~id, rcov=~vsm(ar1m(environment), ism(record)), data=active, control=mmes_control(solve=c(list(method=solver), if(solver == "pcg") list(tol=1e-10)), convergence=list(maxIters=1L), score=if(solver == "pcg") list(traceProbes=24L) else list(), logdet=if(solver == "pcg") list(lanczosSteps=30L) else list(), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    reference <- fitModel("ldlt")
    candidate <- fitModel("pcg")
    repeated <- fitModel("pcg")
    expect_gt(candidate$engineDiagnostics$pcgResidualTraceContractions, 0L)
    expect_equal(candidate$bu, reference$bu, tolerance=1e-7)
    expect_equal(candidate$InfMat, reference$InfMat, tolerance=1e-7)
    expect_equal(candidate$dLu, repeated$dLu, tolerance=1e-12)
  }
})

test_that("D-first ML matches joint factorization and full uncertainty for both direct backends", {
  data(DT_example)
  for(solver in c("ldlt", "cholmod")){
    fitModel <- function(mode){
      mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(solve=list(method=solver), task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=4L), output=list(inverse=mode), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    joint <- fitModel(1L)
    reduced <- fitModel(0L)
    full <- fitModel(2L)
    expect_false(joint$engineDiagnostics$dFirstMLActive)
    expect_true(reduced$engineDiagnostics$dFirstMLActive)
    expect_equal(reduced$engineDiagnostics$CsymbolicAnalyses, 0L)
    expect_true(full$engineDiagnostics$dFirstMLActive)
    if(solver == "cholmod"){
      expect_gt(reduced$engineDiagnostics$cholmodDSelectedBuilds, 0L)
    }
    for(component in c("llik", "theta", "bu", "InfMat", "dLu", "theta_se")){
      expect_equal(reduced[[component]], joint[[component]], tolerance=1e-7)
      expect_equal(full[[component]], joint[[component]], tolerance=1e-7)
    }
    expect_equal(full$uPevList, joint$uPevList, tolerance=1e-7)
    expect_equal(as.matrix(full$C) %*% as.matrix(full$Ci),
      diag(nrow(full$C)) * full$Cscale, tolerance=1e-7)
  }
})
