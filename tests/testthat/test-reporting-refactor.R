reportingRefactorFit <- function(){
  structure(list(
    covStruct=list(a=list(factors=list(), free=TRUE),
                   b=list(factors=list(), free=TRUE),
                   residual=list(factors=list(), free=TRUE)),
    covPar=list(a=1.7, b=2, residual=1),
    theta_se=matrix(c(1, 0.1, 0.2, 0.1, 1, 0.3, 0.2, 0.3, 1), 3L)
  ), class="mmes")
}

test_that("native covariance SE propagates selected terms jointly", {
  fit <- reportingRefactorFit()
  fit$covStruct[[1L]]$free[] <- FALSE
  selected <- c(3L, 1L, 2L)
  reported <- covparams_mmes_se(fit, selected)
  expected <- fit$theta_se[selected, selected]
  expected[2L,] <- expected[,2L] <- 0
  expect_identical(reported$term, names(fit$covStruct)[selected])
  expect_equal(attr(reported, "vcov"), expected, tolerance=1e-7)
  expect_equal(reported$StdError[2L], 0)
  expect_equal(attr(covparams_mmes_se(fit, 3L), "vcov"),
               expected[1L, 1L, drop=FALSE])
  empty <- covparams_mmes_se(fit, character())
  expect_equal(nrow(empty), 0L)
  expect_equal(dim(attr(empty, "vcov")), c(0L, 0L))
})

test_that("vpredict uses the full native covariance quadratic form", {
  fit <- reportingRefactorFit()
  fit$covParNative <- covparams_mmes(fit)
  covariance <- fit$theta_se
  attr(fit$covParNative, "vcov") <- covariance

  result <- vpredict(fit, ~V1/(V1+V2))
  estimate <- fit$covParNative$estimate
  gradient <- c(estimate[2L], -estimate[1L], 0) /
    (estimate[1L] + estimate[2L])^2
  expectedSE <- sqrt(abs(drop(crossprod(gradient, covariance %*% gradient))))
  expect_s3_class(result, "vpredict.mmes")
  expect_named(result, c("Estimate", "SE"))
  expect_equal(as.numeric(result$Estimate), estimate[1L] / sum(estimate[1:2]))
  expect_equal(result$SE, expectedSE)
})

test_that("unavailable uncertainty does not contaminate independent reports", {
  fit <- structure(list(
    covStruct=list(a=list(factors=list(), free=TRUE),
                   b=list(factors=list(), free=TRUE)),
    covPar=list(a=1, b=2),
    theta_se=matrix(c(NA, NA, NA, 0.5), 2, 2)
  ), class="mmes")
  reported <- covparams_mmes_se(fit)
  expect_true(is.na(reported$StdError[1]))
  expect_equal(reported$StdError[2], sqrt(0.5))
  expect_true(all(is.na(attr(reported, "vcov")[1,])))
  fit$covParNative <- reported
  expect_equal(vpredict(fit, ~V2)$SE, sqrt(0.5))
  expect_true(is.na(vpredict(fit, ~V1)$SE))
})

test_that("covmatrix SE agrees with direct Jacobian covariance propagation", {
  set.seed(213)
  data <- expand.grid(id=factor(seq_len(8L)), env=factor(letters[1:3]),
                      replicate=seq_len(2L))
  data$y <- rnorm(nrow(data)) + as.integer(data$id) / 5
  fit <- suppressMessages(mmes(y~env, random=~vsm(usm(env), ism(id)),
    data=data, control=mmes_control(convergence=list(maxIters=3L),
      diagnostics=list(verbose=FALSE))))
  term <- 1L
  covarianceResult <- covmatrix_mmes(fit, term, se=TRUE)
  descriptor <- fit$covStruct[[term]]
  natural <- as.numeric(fit$covPar[[term]])
  counts <- vapply(fit$covPar, length, integer(1))
  start <- sum(counts[seq_len(term - 1L)]) + 1L
  global <- start + seq_along(natural) - 1L
  parameterCovariance <- fit$theta_se[global, global, drop=FALSE]
  free <- as.logical(descriptor$free)
  parameterCovariance[!free,] <- 0
  parameterCovariance[,!free] <- 0
  evaluate <- function(values){
    sommer:::evaluate_covstruct_cpp(
      descriptor, sommer:::.natural_to_working(descriptor, values))
  }
  dimension <- nrow(covarianceResult$covariance)
  expect_equal(dimension, 3L)
  jacobian <- matrix(0, dimension^2, length(natural))
  for(k in which(free)){
    step <- 1e-5 * max(abs(natural[k]), 1e-3)
    plus <- natural
    minus <- natural
    plus[k] <- plus[k] + step
    minus[k] <- minus[k] - step
    jacobian[, k] <- as.vector(evaluate(plus) - evaluate(minus)) / (2 * step)
  }
  expected <- matrix(
    sqrt(pmax(diag(jacobian %*% parameterCovariance %*% t(jacobian)), 0)),
    dimension, dimension
  )
  expect_equal(unname(covarianceResult$covariance.se), unname(expected),
               tolerance=1e-8)
})

test_that("native reporting callbacks must preserve their output shape", {
  factor <- list(par_start=2L, par_end=1L, native_report=list(
    backend="R", fun=function(scale, par, factor, absorb_scale){
      stats::setNames(scale, if(scale > 1) "changed" else "variance")
    }))
  fit <- structure(list(covStruct=list(a=list(factors=list(factor), free=TRUE)),
    covPar=list(a=1), theta_se=matrix(1)), class="mmes")
  expect_error(covparams_mmes_se(fit), "changed output shape")
})

test_that("random-effect outputs share the original response scale", {
  set.seed(214)
  dataset <- data.frame(id=factor(rep(seq_len(8L), each=3L)))
  dataset$y <- 20 + 7 * (rep(rnorm(8L), each=3L) +
    rnorm(nrow(dataset), sd=0.5))
  X <- matrix(1, nrow(dataset), 1L)
  Z <- model.matrix(~id-1, dataset)
  W <- cbind(X, Z)
  precision <- crossprod(W) / 12.25 + diag(c(0, rep(1 / 49, 8L)))
  expected <- solve(precision, crossprod(W, dataset$y) / 12.25)
  for(method in c("ldlt", "cholmod", "pcg", "direct")){
    control <- mmes_control(
      formulation=if(method == "direct") "direct" else "henderson",
      solve=if(method == "direct") list() else list(method=method),
      convergence=list(maxIters=1L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))
    fit <- mmes(y~1, random=~vsm(ism(id), sigma2=49, fixedSigma2=TRUE),
      rcov=~vsm(ism(units), sigma2=12.25, fixedSigma2=TRUE),
      data=dataset, control=control)
    expect_equal(as.numeric(fit$bu), as.numeric(expected), tolerance=1e-6)
    expect_equal(as.numeric(fit$u), as.numeric(expected[-1L]), tolerance=1e-6)
    expect_equal(as.numeric(fit$uList[[1L]]), as.numeric(fit$u), tolerance=1e-8)
    expect_equal(as.numeric(fit$bu), c(as.numeric(fit$b), as.numeric(fit$u)),
      tolerance=1e-8)
  }
})
