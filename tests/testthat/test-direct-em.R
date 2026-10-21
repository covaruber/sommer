.direct_em_data <- function(){
  set.seed(982)
  dataset <- expand.grid(id=factor(seq_len(6L)), trait=factor(seq_len(3L)),
                         replicate=factor(seq_len(2L)))
  effects <- matrix(rnorm(18L), 6L)
  dataset$y <- as.numeric(dataset$trait) / 2 +
    effects[cbind(as.integer(dataset$id), as.integer(dataset$trait))] +
    rnorm(nrow(dataset), sd=0.7)
  dataset
}

.direct_em_reference <- function(setup, response){
  X <- as.matrix(setup$X)
  Z <- lapply(setup$Z, as.matrix)
  D <- do.call(cbind, c(list(X), Z))
  theta <- lapply(setup$covStruct, function(cs){
    sommer:::evaluate_covstruct_cpp(cs, cs$par)
  })
  R0 <- matrix(0, length(response), length(response))
  for(rows in split(seq_along(response), setup$residualBlock)){
    index <- setup$residualIndex[rows]
    R0[rows, rows] <- theta[[length(theta)]][index, index, drop=FALSE]
  }
  U <- if(setup$useH) chol(as.matrix(setup$W)) else diag(length(response))
  Uinv <- solve(U)
  R <- t(Uinv) %*% R0 %*% Uinv
  Ri <- solve(R)
  prior <- matrix(0, ncol(D), ncol(D))
  indices <- list()
  first <- ncol(X) + 1L
  for(term in seq_along(setup$Ai)){
    selected <- which(setup$Zind == term)
    width <- ncol(Z[[selected[1L]]])
    index <- first:(first + length(selected) * width - 1L)
    prior[index, index] <- kronecker(solve(theta[[term]]), as.matrix(setup$Ai[[term]]))
    indices[[term]] <- matrix(index, width)
    first <- max(index) + 1L
  }
  information <- crossprod(D, Ri %*% D) + prior
  posterior <- solve(information)
  mean <- drop(posterior %*% crossprod(D, Ri %*% response))
  if(!setup$REML){
    random <- seq.int(ncol(X) + 1L, ncol(D))
    posterior[,] <- 0
    posterior[random, random] <- solve(information[random, random, drop=FALSE])
    randomMean <- drop(D[, random, drop=FALSE] %*% mean[random])
    mean[seq_len(ncol(X))] <- solve(crossprod(X, Ri %*% X),
      crossprod(X, Ri %*% (response - randomMean)))
  }
  updated <- lapply(seq_along(indices), function(term){
    index <- indices[[term]]
    means <- matrix(mean[index], nrow(index))
    precision <- as.matrix(setup$Ai[[term]])
    moment <- crossprod(means, precision %*% means)
    for(a in seq_len(ncol(index))) for(b in seq_len(ncol(index))){
      moment[a, b] <- moment[a, b] +
        sum(precision * posterior[index[, a], index[, b], drop=FALSE])
    }
    moment / nrow(index)
  })
  error <- drop(t(U) %*% (response - D %*% mean))
  transformed <- t(U) %*% D
  residualMoment <- tcrossprod(error) + transformed %*% posterior %*% t(transformed)
  updated[[length(updated) + 1L]] <- matrix(sum(diag(solve(R0, residualMoment))) /
    length(response) * theta[[length(theta)]][1L, 1L], 1L)
  updated
}

test_that("direct one-step EM matches independent ML and REML posterior moments", {
  dataset <- .direct_em_data()
  relationship <- 0.3^abs(outer(seq_len(6L), seq_len(6L), "-"))
  dimnames(relationship) <- rep(list(levels(dataset$id)), 2L)
  relationship <- solve(relationship)
  attr(relationship, "inverse") <- TRUE
  for(reml in c(TRUE, FALSE)) for(shape in c("scalar", "diagonal", "us", "weighted")){
    random <- switch(shape,
      scalar=~vsm(ism(id), Gu=relationship, sigma2=0.8) +
        vsm(ism(replicate), sigma2=0.3),
      diagonal=~vsm(dsm(trait), ism(id), Gu=relationship, sigma2=0.8),
      us=~vsm(usm(trait), ism(id), Gu=relationship, sigma2=0.8),
      weighted=~vsm(ism(id), Gu=relationship, sigma2=0.8))
    arguments <- list(fixed=y~trait, random=random,
      rcov=~vsm(ism(units), sigma2=0.6), data=dataset,
      control=mmes_control(formulation="direct",
        task=list(objective=if(reml) "reml" else "ml"),
        update=list(blendWeight=1), proposal=list(stepWeight=1),
        convergence=list(maxIters=2L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    if(shape == "weighted"){
      arguments$W <- Matrix::Matrix(
        0.2^abs(outer(seq_len(nrow(dataset)), seq_len(nrow(dataset)), "-")), sparse=TRUE)
    }
    setupControl <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(output=list(setup=TRUE)),
      keep.null=TRUE))
    setupArguments <- arguments
    setupArguments$control <- setupControl
    setup <- do.call(mmes, setupArguments)
    reference <- .direct_em_reference(setup, dataset$y)
    if(shape == "diagonal") reference[[1L]] <- diag(diag(reference[[1L]]))
    fit <- do.call(mmes, arguments)
    for(term in seq_along(reference)){
      expect_equal(unname(fit$theta[[term]]), reference[[term]], tolerance=1e-8,
        info=paste(reml, shape, term))
    }
    expect_equal(fit$engineDiagnostics$direct$aiUpdateSolves, 0L)
    expect_equal(fit$engineDiagnostics$blend$closedFormSteps, 1L)
    expect_gte(fit$engineDiagnostics$blend$qGains[1L], 0)
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-9))
  }
})

test_that("direct EM preserves fixed parameters and covariance ties", {
  dataset <- .direct_em_data()
  dataset$id2 <- dataset$id
  for(reml in c(TRUE, FALSE)){
    fit <- mmes(y~trait, random=~vsm(ism(id), sigma2=0.8, fixedSigma2=TRUE), data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=1), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(unname(fit$covPar[[1L]]), 0.8, ignore_attr=TRUE)
    tied <- mmes(y~trait, random=~id+id2, data=dataset, vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)), control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=1), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(tied$monitor[1L, ], 2*tied$monitor[2L, ], tolerance=1e-8)
    expect_true(all(diff(as.numeric(tied$llik)) >= -1e-8))
  }
})

test_that("direct EM retains fixed tiny variances and admits starts below the free scale floor", {
  dataset <- .direct_em_data()
  for(reml in c(TRUE, FALSE)) for(fixed in c(TRUE, FALSE)){
    fit <- mmes(y~trait, random=~vsm(ism(id), sigma2=1e-10, fixedSigma2=fixed), data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=1), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_true(all(diff(as.numeric(fit$llik)) >= -1e-8))
    if(fixed){
      expect_equal(as.numeric(fit$covPar[[1L]]) / 1e-10, 1, tolerance=1e-14)
    } else {
      expect_true(all(fit$monitor[1L, ] >= 1e-10 * (1-1e-10)))
    }
  }
})

test_that("direct blending and PX-GEM use true EM proposals", {
  dataset <- .direct_em_data()
  for(reml in c(TRUE, FALSE)) for(method in c("em", "pxem")){
    for(weight in list(0, 1, 0.4, c(1, 0.4, 0), 3)){
      fit <- mmes(y~trait, random=~id, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(recovery=method, blendWeight=weight), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      expect_true(all(diff(as.numeric(fit$llik)) >= -1e-8))
      if(identical(weight, 0)){
        expect_equal(fit$engineDiagnostics$blend$eStepEvaluations, 0L)
        expect_equal(fit$engineDiagnostics$blend$qEvaluations, 0L)
      } else {
        expect_gt(fit$engineDiagnostics$blend$eStepEvaluations, 0L)
      }
      if(identical(weight, 1)) expect_equal(fit$engineDiagnostics$direct$aiUpdateSolves, 0L)
      if(method == "pxem" && !identical(weight, 0)){
        expect_true(all(na.omit(fit$engineDiagnostics$blend$pxQGains) >= -1e-7))
      }
    }
  }
})

test_that("direct AI-first recovery is lazy, resumes AI, and reports its trigger", {
  dataset <- .direct_em_data()
  healthy <- mmes(y~trait, random=~id, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(healthy$engineDiagnostics$blend$eStepEvaluations, 0L)
  expect_equal(healthy$engineDiagnostics$blend$qEvaluations, 0L)
  for(reml in c(TRUE, FALSE)) for(method in c("em", "pxem")){
    fit <- mmes(y~trait, random=~id, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(recovery=method), proposal=list(stepWeight=1e-7), convergence=list(maxIters=6L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    recovery <- fit$engineDiagnostics$recovery
    expect_identical(recovery$automatic, TRUE)
    expect_equal(recovery$iterations, 4L)
    expect_equal(recovery$reasons, "negligible-gains")
    expect_equal(recovery$accepted, 1L)
    expect_equal(fit$engineDiagnostics$blend$eStepEvaluations, 1L)
    expect_equal(fit$engineDiagnostics$blend$weights, c(0, 0, 0, 1, 0, 0))
    expect_false(fit$convergence)
    setup <- mmes(y~trait, random=~id, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_true(isTRUE(attr(setup$covStruct, "adaptiveRecovery")))
  }
})

test_that("direct exhausted AI searches and unreliable information recover explicitly", {
  dataset <- .direct_em_data()
  for(method in c("em", "pxem")){
    recovered <- mmes(y~trait, random=~id, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", update=list(recovery=method), proposal=list(stepWeight=1e-20), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(recovered$engineDiagnostics$recovery$iterations, c(1L, 2L))
    expect_true(all(recovered$engineDiagnostics$recovery$reasons == "line-search"))
    expect_equal(recovered$engineDiagnostics$recovery$accepted, 2L)
    expect_true(all(diff(as.numeric(recovered$llik)) > 0))
  }
  dataset$id2 <- dataset$id
  regularized <- mmes(y~trait, random=~id+id2, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_gt(regularized$engineDiagnostics$recovery$regularizedAttempts, 0L)
  expect_gt(regularized$engineDiagnostics$recovery$regularizedAcceptedSteps, 0L)
  expect_true(all(diff(as.numeric(regularized$llik)) > 0))
})

test_that("direct residual-only EM returns analytical ML and REML estimates", {
  dataset <- data.frame(x=seq(-1, 1, length.out=24L))
  dataset$y <- 3 + 2*dataset$x + sin(seq_len(nrow(dataset)))
  errors <- residuals(lm(y~x, data=dataset))
  for(reml in c(TRUE, FALSE)){
    fit <- mmes(y~x, data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", task=list(objective=if(reml) "reml" else "ml"), update=list(blendWeight=1), convergence=list(maxIters=100L, likelihoodTol=1e-10, parameterTol=1e-10), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expected <- sum(errors^2) / (nrow(dataset)-if(reml) 2L else 0L)
    expect_true(fit$convergence)
    expect_equal(as.numeric(fit$covPar[[1L]]), expected, tolerance=1e-7)
    expect_equal(fit$engineDiagnostics$direct$aiUpdateSolves, 0L)
  }
})

test_that("direct EM keeps evaluated states, PEV, and original-scale inference consistent", {
  dataset <- .direct_em_data()
  fitModel <- function(iterations=1L, inverse=0L, setup=FALSE, denseMemoryMB=NULL){
    execution <- if(is.null(denseMemoryMB)) list() else
      list(denseMemoryMB=denseMemoryMB)
    mmes(y~trait, random=~vsm(ism(id), sigma2=0.8),
      rcov=~vsm(ism(units), sigma2=0.6), data=dataset,
      control=mmes_control(formulation="direct", execution=execution,
        update=list(blendWeight=1),
        convergence=list(maxIters=iterations), output=list(inverse=inverse,
          setup=setup), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  initial <- fitModel()
  expect_false(initial$convergence)
  setup <- fitModel(setup=TRUE)
  expectedInitial <- vapply(setup$covStruct, function(cs) exp(cs$par[1L]), numeric(1))
  expect_equal(unname(unlist(initial$covPar)), unname(expectedInitial))
  expect_equal(length(initial$llik), 1L)
  expect_equal(initial$engineDiagnostics$blend$eStepEvaluations, 0L)
  fitted <- fitModel(iterations=4L, inverse=2L)
  known <- mmes(y~trait, random=~vsm(ism(id), sigma2=as.numeric(fitted$covPar[[1L]]), fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=as.numeric(fitted$covPar[[2L]]), fixedSigma2=TRUE), data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", update=list(blendWeight=0), convergence=list(maxIters=1L), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(tail(as.numeric(fitted$llik), 1L), as.numeric(known$llik), tolerance=1e-9)
  expect_equal(fitted$bu, known$bu, tolerance=1e-9)
  expect_equal(fitted$VarBeta, known$VarBeta, tolerance=1e-9)
  expect_equal(unname(fitted$uPevList), unname(known$uPevList), tolerance=1e-9)
  expect_equal(as.numeric(fitted$monitor[, ncol(fitted$monitor)]),
    as.numeric(unlist(fitted$covPar)), tolerance=1e-9)
  expect_true(all(is.finite(fitted$theta_se)))
  expect_equal(fitted$AIC, -2*tail(as.numeric(fitted$llik), 1L) +
    2*ncol(model.matrix(~trait, dataset)))
  expect_error(fitModel(iterations=2L, denseMemoryMB=0.01),
    "Exact direct EM moments exceed denseMemoryMB")
})

test_that("direct EM evaluates residual callback failures explicitly", {
  dataset <- .direct_em_data()
  broken <- function(par) stop("intentional covariance callback failure")
  expect_error(mmes(y~trait, random=~vsm(ownm(trait, fun=broken, par=0.2), ism(id)), data=dataset, control=mmes_control(formulation=if(FALSE) "henderson" else "direct", update=list(blendWeight=1), convergence=list(maxIters=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "intentional covariance callback failure")
})
