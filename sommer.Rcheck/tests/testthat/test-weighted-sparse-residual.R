test_that("CHOLMOD selected residual traces match LDLT for weighted incomplete blocks", {
  dataset <- expand.grid(id=factor(seq_len(4L)), environment=factor(seq_len(3L)),
    replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 5
  dataset <- dataset[-c(1L, 7L), ]
  weights <- diag(nrow(dataset)) + matrix(0.02, nrow(dataset), nrow(dataset))
  fitModel <- function(solver){
    mmes(response~environment, random=~id, rcov=~vsm(usm(environment), ism(record)), W=weights, data=dataset, control=mmes_control(solve=list(method=solver), convergence=list(maxIters=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  reference <- fitModel("ldlt")
  candidate <- fitModel("cholmod")
  expect_gt(candidate$engineDiagnostics$cholmodRSelectedBuilds, 0L)
  for(component in c("llik", "theta", "bu", "InfMat", "dLu")){
    expect_equal(candidate[[component]], reference[[component]], tolerance=1e-7)
  }
})

test_that("sparse non-diagonal weights preserve diagonal-residual fits", {
  set.seed(7)
  n <- 60L
  data <- data.frame(id=factor(rep(seq_len(n / 2L), each=2L)), y=rnorm(n))
  weights <- Matrix::bdiag(replicate(n / 2L,
                                    matrix(c(1.4, -0.2, -0.2, 1.4), 2L),
                                    simplify=FALSE))
  sparse <- suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data, W=weights, control=mmes_control(convergence=list(maxIters=3), solve=list(method="ldlt"), update=list(blendWeight=20), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
  expect_equal(unname(unlist(sparse$theta)),
               c(0.141923343971598, 1.24385475608586), tolerance=1e-8)
  expect_equal(head(as.numeric(sparse$bu), 5L),
               c(0.218356994943944, 0.0792772728361476,
                 -0.16682663255206, -0.253492515401901,
                 0.0231439286116111), tolerance=1e-8)
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
  plain <- suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data, W=W, control=mmes_control(convergence=list(maxIters=4), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
  grouped <- suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data, W=W, weights=~trial, control=mmes_control(convergence=list(maxIters=4), solve=list(method="ldlt"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
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
      W=W, control=mmes_control(convergence=list(maxIters=5L),
        solve=list(method="ldlt"), diagnostics=list(verbose=FALSE,
          dateWarning=FALSE)))
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
      W=W, control=mmes_control(convergence=list(maxIters=3L),
        solve=list(method="ldlt"), diagnostics=list(verbose=FALSE,
          dateWarning=FALSE)))
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
    suppressMessages(mmes(y~1, random=~id, rcov=~units, data=data, W=weights, control=mmes_control(convergence=list(maxIters=3), solve=list(method=solver), score=if(solver == "pcg") list(traceProbes=512) else list(), update=list(blendWeight=0), logdet=if(solver == "pcg") list(lanczosSteps=40) else list(), output=list(inverse=if(solver == "pcg") 2 else 0), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
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
    c(0.223291894770803, 1.30307957166026, -0.829268180618696, 0.483128144582169),
    c(0.220896144065077, 0.65129496259526, -0.683012787315524, 0.731017895495925)
  )
  likelihood <- c(-52.34982835241, -55.4451301835064)
  for(index in seq_along(reference)){
    fit <- suppressMessages(mmes(y~env, random=~id, rcov=~vsm(usm(env), ism(units)), data=data, W=weights, control=mmes_control(convergence=list(maxIters=3), solve=list(method="ldlt"), task=list(objective=if(index == 1L) "reml" else "ml"), update=list(blendWeight=if(index == 1L) 20 else 0), diagnostics=list(verbose=FALSE, dateWarning=FALSE))))
    expect_equal(unname(unlist(fit$covPar)), reference[[index]], tolerance=1e-8)
    expect_equal(tail(as.numeric(fit$llik), 1L), likelihood[index], tolerance=1e-8)
  }
})
