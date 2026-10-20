test_that("rrm has exactly zero specifics and an identified loading scale", {
  environment <- factor(letters[1:5])
  loading <- matrix(c(2, 0.4, -0.6, 1.2, 0.8, 0, 1, 0.3, -0.2, 0.7), 5, 2)
  factor <- rrm(environment, 2L, loadings=loading)$covFactor
  covariance <- factor$evaluator$fun(factor$par)
  expect_equal(covariance, tcrossprod(loading) / loading[1L, 1L]^2)
  expect_equal(sum(eigen(covariance, symmetric=TRUE)$values > 1e-10), 2L)
  expect_identical(factor$par[1L], 0)
  expect_false(factor$free[1L])
  expect_equal(sum(factor$free), 8L)
  expect_true(factor$representation$loadings$zero_specific)
  rescaled <- rrm(environment, 2L, loadings=loading * 3)$covFactor
  expect_equal(rescaled$par, factor$par)
  expect_error(rrm(environment, NA), "k in rrm")
  expect_error(rrm(environment, 1.5), "k in rrm")
  expect_error(rrm(environment, 1, fixed=rep(NA, 5)), "logical vector")
})

test_that("homogeneous fam preserves the former RR covariance family", {
  environment <- factor(letters[1:5])
  loading <- matrix(c(0.6, 0.2, -0.3, 0.1, 0.4, 0, 0.7, 0.2, -0.1, 0.3), 5, 2)
  homogeneous <- fam(environment, 2L, loadings=loading, specific=1,
                     variance="homogeneous")$covFactor
  heterogeneous <- fam(environment, 2L, loadings=loading, specific=rep(1, 5))$covFactor
  expect_length(homogeneous$par, 9L)
  expect_length(heterogeneous$par, 13L)
  expect_false(homogeneous$representation$loadings$heterogeneous_specific)
  expect_false(homogeneous$representation$loadings$zero_specific)
  natural <- homogeneous$par
  diagonal <- homogeneous$report$transform == "exp"
  natural[diagonal] <- exp(natural[diagonal])
  native <- homogeneous$native_report$fun(2, natural, homogeneous)
  reported <- matrix(0, 5, 2)
  reported[cbind(homogeneous$fa_row, homogeneous$fa_col)] <- native[seq_len(9L)]
  expected <- 2 * (tcrossprod(loading) + diag(5)) / (1 + loading[1, 1]^2)
  expect_equal(unname(tcrossprod(reported) + diag(native[10L], 5)), expected)
  expect_error(fam(environment, variance="homogeneous", specific=1:5),
               "single value or q equal")
  expect_error(fam(environment, variance="homogeneous", specific=0), "positive")
})

.rr_test_data <- function(){
  set.seed(4821)
  n <- 12L
  dataset <- expand.grid(id=factor(seq_len(n)), environment=factor(letters[1:4]),
                         replicate=seq_len(2L))
  relationship <- 0.35^abs(outer(seq_len(n), seq_len(n), "-"))
  dimnames(relationship) <- rep(list(levels(dataset$id)), 2L)
  loading <- c(1.2, 0.7, -0.5, 1.6)
  scores <- drop(t(chol(relationship)) %*% rnorm(n))
  dataset$response <- c(2, -1, 0, 1)[dataset$environment] +
    loading[dataset$environment] * scores[dataset$id] + rnorm(nrow(dataset), sd=0.6)
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  list(data=dataset, relationship=relationship, precision=precision)
}

test_that("zero-specific RR matches independent observation-space REML", {
  example <- .rr_test_data()
  dataset <- example$data
  X <- model.matrix(~environment, dataset)
  relationship <- example$relationship[as.integer(dataset$id), as.integer(dataset$id)]
  objective <- function(parameters){
    loading <- c(exp(parameters[1L]), parameters[2:4])
    V <- tcrossprod(loading[dataset$environment]) * relationship
    diag(V) <- diag(V) + exp(parameters[5L])
    cholesky <- chol(V)
    Vi <- chol2inv(cholesky)
    information <- crossprod(X, Vi %*% X)
    beta <- solve(information, crossprod(X, Vi %*% dataset$response))
    residual <- dataset$response - drop(X %*% beta)
    as.numeric(0.5 * (2 * sum(log(diag(cholesky))) +
      2 * sum(log(diag(chol(information)))) + crossprod(residual, Vi %*% residual)))
  }
  optimum <- optim(c(log(1.2), 0.7, -0.5, 1.6, log(0.36)), objective,
                   method="BFGS", control=list(reltol=1e-12, maxit=500))
  expect_identical(optimum$convergence, 0L)
  expectedLoading <- c(exp(optimum$par[1L]), optimum$par[2:4])
  precision <- example$precision
  for(henderson in c(TRUE, FALSE)){
    fit <- mmes(response~environment, random=~vsm(rrm(environment, 1L), ism(id), Gu=precision), data=dataset, control=mmes_control(formulation=if(henderson) "henderson" else "direct", convergence=list(maxIters=150, likelihoodTol=1e-10, parameterTol=1e-8), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(unname(fit$theta[[1L]]), tcrossprod(expectedLoading), tolerance=2e-3)
    expect_equal(unname(fit$theta[[2L]][1L, 1L]), exp(optimum$par[5L]), tolerance=2e-3)
    normalized <- loadings_mmes(fit, varianceScale=FALSE, rotation=FALSE)
    scaled <- loadings_mmes(fit, varianceScale=TRUE, rotation=FALSE)
    expect_equal(unname(scaled$specific), rep(0, 4L))
    expect_equal(unname(tcrossprod(scaled$loadings)), unname(fit$theta[[1L]]), tolerance=1e-8)
    expect_equal(unname(normalized$sigma2 * tcrossprod(normalized$loadings)),
                 unname(fit$theta[[1L]]), tolerance=1e-8)
    native <- covparams_mmes_se(fit)
    specifics <- grepl("^specific_variance\\[", native$parameter)
    expect_equal(native$estimate[specifics], rep(0, 4L))
    expect_true(all(is.finite(native$StdError[!specifics])))
    expect_true(all(is.finite(unlist(fit$uPevList))))
    expect_true(all(unlist(fit$uPevList) >= -1e-8))
    if(henderson){
      expect_true(fit$engineDiagnostics$reducedRankActive)
      expect_lt(fit$engineDiagnostics$latentCoefficients,
                fit$engineDiagnostics$originalCoefficients)
      expect_equal(nrow(fit$C), fit$engineDiagnostics$latentCoefficients)
    }
    for(rotation in c(TRUE, FALSE)) for(varianceScale in c(TRUE, FALSE)){
      info <- loadings_mmes(fit, rotation=rotation, varianceScale=varianceScale)
      scores <- scores_mmes(fit, rotation=rotation, varianceScale=varianceScale)
      expect_equal(unname(scores %*% t(info$loadings)), unname(fit$uList[[1L]]),
                   tolerance=1e-8)
    }
    expect_error(scores_mmes(fit, method="bartlett"), "undefined for rrm")
  }
})

test_that("rank-two RR agrees across engines with extra terms and residual heterogeneity", {
  set.seed(1904)
  dataset <- expand.grid(id=factor(seq_len(9L)), environment=factor(letters[1:5]),
                         replicate=factor(seq_len(2L)))
  loading <- matrix(c(1, 0.5, -0.7, 0.4, 1.2, 0, 0.9, 0.2, -0.6, 0.5), 5, 2)
  scores <- matrix(rnorm(18), 9, 2)
  genetic <- scores %*% t(loading)
  dataset$response <- genetic[cbind(as.integer(dataset$id), as.integer(dataset$environment))] +
    rnorm(nrow(dataset), sd=0.6)
  dataset <- dataset[-c(3, 15, 72), ]
  shape <- rrm(dataset$environment, 2L, loadings=loading, fixed=rep(TRUE, 9L))
  fits <- lapply(c(TRUE, FALSE), function(henderson){
    mmes(response~environment, random=~vsm(shape, ism(id)) + vsm(ism(replicate)), rcov=~vsm(dsm(environment), ism(units)), data=dataset, control=mmes_control(formulation=if(henderson) "henderson" else "direct", convergence=list(maxIters=100, likelihoodTol=1e-10, parameterTol=1e-8), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  })
  expect_equal(unname(fits[[1L]]$theta), unname(fits[[2L]]$theta), tolerance=2e-3)
  expect_equal(fits[[1L]]$bu, fits[[2L]]$bu, tolerance=2e-3)
  expect_equal(fits[[1L]]$uPevList, fits[[2L]]$uPevList, tolerance=2e-3)
  expect_identical(fits[[1L]]$control$formulation, "henderson")
  expect_identical(fits[[2L]]$control$formulation, "direct")
})

test_that("singular RR misuse fails explicitly", {
  example <- .rr_test_data()
  dataset <- example$data
  expect_error(mmes(response~1, rcov=~vsm(rrm(environment), ism(units)), data=dataset, control=mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "only supported for random effects")
  expect_error(mmes(response~1, random=~vsm(rrm(environment), ism(id)), data=dataset, control=mmes_control(task=list(mode=if(TRUE) "solve" else "fit"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "task\\$mode='solve' does not support")
  expect_error(mmes(response~1, random=~vsm(rrm(environment), ism(id)), data=dataset, control=mmes_control(representation=list(factorScores="fixed-shape"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "not representation\\$factorScores")
})

test_that("latent RR prediction and post-fit PEV use the original effect space", {
  example <- .rr_test_data()
  dataset <- example$data
  precision <- example$precision
  shape <- rrm(dataset$environment, 1L, fixed=rep(TRUE, 4L))
  fit <- mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision, sigma2=0.8, fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=0.5, fixedSigma2=TRUE), data=dataset, control=mmes_control(output=list(inverse=0L), convergence=list(maxIters=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  full <- postPEV(fit, mode=2L)
  selected <- postPEV(fit, mode=1L)
  expect_equal(diag(selected$Ci), diag(full$Ci), tolerance=1e-8)
  expect_equal(selected$uPevList, full$uPevList, tolerance=1e-8)
  D <- diag(nrow(fit$bu))[seq_len(8L), , drop=FALSE]
  prediction <- predict(fit, D=D)
  expect_equal(prediction$pvals$predicted.value, drop(D %*% fit$bu), tolerance=1e-8)
  expect_equal(prediction$pvals$std.error^2, diag(as.matrix(full$Ci))[seq_len(8L)],
               tolerance=1e-8)
  expect_error(mmes(response~environment, random=~vsm(shape, ism(id)), data=dataset, control=mmes_control(solve=list(method="pcg"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "requires solver='ldlt' or solver='cholmod'")
})
