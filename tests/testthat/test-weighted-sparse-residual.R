test_that("sparse non-diagonal weights preserve diagonal-residual fits", {
  set.seed(7)
  n <- 60L
  data <- data.frame(id=factor(rep(seq_len(n / 2L), each=2L)), y=rnorm(n))
  weights <- Matrix::bdiag(replicate(n / 2L,
                                    matrix(c(1.4, -0.2, -0.2, 1.4), 2L),
                                    simplify=FALSE))
  sparse <- suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data,
                                  W=weights, nIters=3, solver="ldlt",
                                  verbose=FALSE, dateWarning=FALSE))
  expect_equal(unname(unlist(sparse$theta)),
               c(0.1532013323059318, 1.0711724582035496), tolerance=1e-8)
  expect_equal(head(as.numeric(sparse$bu), 5L),
               c(0.21835699494394373, 0.09423111073305052,
                 -0.19829464766940430, -0.30130806010703842,
                 0.02750949952965763), tolerance=1e-8)
})

test_that("weights formula declares independent blocks without changing the fit", {
  set.seed(18)
  data <- data.frame(
    trial=factor(rep(seq_len(30L), each=2L)),
    id=factor(rep(seq_len(15L), each=4L)),
    y=rnorm(60L)
  )
  W <- Matrix::bdiag(replicate(30L,
    matrix(c(1.4, -0.2, -0.2, 1.4), 2L), simplify=FALSE))
  plain <- suppressMessages(mmes(y~1, random=~id, rcov=~units,
    data=data, W=W, nIters=4, solver="ldlt", verbose=FALSE,
    dateWarning=FALSE))
  grouped <- suppressMessages(mmes(y~1, random=~id, rcov=~units,
    data=data, W=W, weights=~trial, nIters=4, solver="ldlt",
    verbose=FALSE, dateWarning=FALSE))
  expect_equal(unlist(grouped$theta), unlist(plain$theta), tolerance=1e-10)
  expect_equal(as.numeric(grouped$bu), as.numeric(plain$bu), tolerance=1e-10)
  expect_equal(grouped$engineDiagnostics$weightBlockCount, 30L)
})

test_that("weights formula rejects W cross-block entries and unsupported use", {
  data <- data.frame(trial=factor(rep(1:2, each=2L)), y=rnorm(4L))
  W <- diag(4L)
  W[1L, 3L] <- W[3L, 1L] <- 0.1
  expect_error(mmes(y~1, rcov=~units, data=data, W=W, weights=~trial),
               "across groups")
  expect_error(mmes(y~1, rcov=~units, data=data, weights=~trial),
               "provide W")
  expect_error(mmes(y~1, rcov=~units, data=data, W=diag(4L),
                   weights=stats::as.formula("trial ~ y")),
               "one-sided formula")
})

test_that("weight groups do not restrict cross-group residual covariance", {
  set.seed(24)
  data <- expand.grid(id=factor(seq_len(12L)), trial=factor(seq_len(3L)))
  data$y <- rnorm(nrow(data))
  diagonalWeights <- 1 + rep(c(0.1, 0.4, 0.8), each=12L)
  W <- Matrix::Diagonal(nrow(data), x=diagonalWeights)
  fit <- function(declareBlocks){
    args <- list(fixed=y~1, rcov=~vsm(usm(trial), ism(id)), data=data,
      W=W, nIters=5, solver="ldlt", verbose=FALSE, dateWarning=FALSE)
    if(declareBlocks) args$weights <- ~trial
    suppressMessages(do.call(mmes, args))
  }
  ungrouped <- fit(FALSE)
  grouped <- fit(TRUE)
  expect_equal(unlist(grouped$theta), unlist(ungrouped$theta), tolerance=1e-9)
  expect_equal(as.numeric(grouped$bu), as.numeric(ungrouped$bu), tolerance=1e-9)
})

test_that("weight block formulas work when group rows are interleaved", {
  set.seed(27)
  data <- expand.grid(id=factor(seq_len(12L)), trial=factor(c("A", "B")))
  data$y <- rnorm(nrow(data))
  W <- Matrix::Matrix(0, nrow(data), nrow(data), sparse=TRUE)
  for(trial in levels(data$trial)){
    rows <- which(data$trial == trial)
    block <- Matrix::bdiag(replicate(6L,
      matrix(c(1.3, -0.15, -0.15, 1.3), 2L), simplify=FALSE))
    W[rows, rows] <- block
  }
  fit <- function(declareBlocks){
    args <- list(fixed=y~1, random=~id, rcov=~units, data=data,
      W=W, nIters=3, solver="ldlt", verbose=FALSE, dateWarning=FALSE)
    if(declareBlocks) args$weights <- ~trial
    suppressMessages(do.call(mmes, args))
  }
  expect_equal(unlist(fit(TRUE)$theta), unlist(fit(FALSE)$theta), tolerance=1e-9)
})

test_that("PCG retains off-diagonal weighted residual derivative traces", {
  set.seed(7)
  data <- data.frame(id=factor(rep(seq_len(30L), each=2L)), y=rnorm(60L))
  weights <- Matrix::bdiag(replicate(30L,
    matrix(c(1.4, -0.2, -0.2, 1.4), 2L), simplify=FALSE))
  fit <- function(solver){
    suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data,
      W=weights, nIters=3, solver=solver, pcgTraceProbes=512,
      pcgLanczosSteps=40, computeCi=if(solver == "pcg") 2 else 0,
      verbose=FALSE, dateWarning=FALSE))
  }
  exact <- fit("ldlt")
  iterative <- fit("pcg")
  expect_gt(iterative$engineDiagnostics$pcgBatchSolves, 0L)
  expect_gt(iterative$engineDiagnostics$pcgBatchIterations, 0L)
  expect_equal(unlist(iterative$theta), unlist(exact$theta), tolerance=0.01)
  expect_equal(as.numeric(iterative$bu), as.numeric(exact$bu), tolerance=0.01)
})

test_that("batched irregular weighted residual assembly preserves REML and ML", {
  set.seed(19)
  data <- expand.grid(env=factor(1:2), id=factor(1:40))
  data$y <- rnorm(nrow(data))
  weights <- Matrix::bdiag(replicate(40L,
    matrix(c(1.4, -0.2, -0.2, 1.4), 2L), simplify=FALSE))
  reference <- list(
    c(0.117592930734105, 1.03075080183896, -0.467457634944157, 0.929889253815438),
    c(0.112735152698391, 1.01003588929121, -0.451567555673451, 0.927177476757712)
  )
  likelihood <- c(-54.9015401677679, -51.2584732333244)
  for(index in seq_along(reference)){
    fit <- suppressMessages(mmes(y~env, random=~id,
      rcov=~vsm(usm(env), ism(units)), data=data, W=weights,
      nIters=3, solver="ldlt", REML=index == 1L,
      verbose=FALSE, dateWarning=FALSE))
    expect_equal(unname(unlist(fit$covPar)), reference[[index]], tolerance=1e-8)
    expect_equal(tail(as.numeric(fit$llik), 1L), likelihood[index], tolerance=1e-8)
  }
})