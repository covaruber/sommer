.rr_ai_fixture <- function(){
  set.seed(376)
  dataset <- expand.grid(id=factor(seq_len(7L)), environment=factor(letters[1:4]),
                         replicate=factor(seq_len(2L)))
  loading <- matrix(c(1, 0.6, -0.4, 1.1, 0, 0.8, 0.3, -0.5), 4L, 2L)
  scores <- matrix(rnorm(14L), 7L, 2L)
  dataset$response <- c(1, -0.5, 0.3, 0.8)[dataset$environment] +
    rowSums(scores[dataset$id, ] * loading[dataset$environment, ]) +
    rnorm(nrow(dataset), sd=0.7)
  dataset <- dataset[-c(3L, 17L), ]
  relationship <- 0.25^abs(outer(seq_len(7L), seq_len(7L), "-"))
  dimnames(relationship) <- rep(list(levels(dataset$id)), 2L)
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  list(data=dataset, precision=precision, loading=loading)
}

.rr_ai_dense_reference <- function(setup, working, reml, residualCovariance=NULL){
  counts <- vapply(setup$covStruct, function(structure) length(structure$par), integer(1))
  offsets <- cumsum(c(0L, counts))
  evaluate <- function(parameters){
    covariance <- matrix(0, nrow(setup$X), nrow(setup$X))
    for(term in seq_along(counts)){
      local <- parameters[offsets[term] + seq_len(counts[term])]
      if(term < length(counts)){
        shape <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[term]], local)
        designs <- setup$Z[which(setup$Zind == term)]
        relationship <- solve(as.matrix(setup$Ai[[term]]))
        for(first in seq_along(designs)) for(second in seq_along(designs)){
          covariance <- covariance + shape[first, second] *
            as.matrix(designs[[first]] %*% relationship %*% t(designs[[second]]))
        }
      }else{
        if(is.null(residualCovariance)){
          shape <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[term]], local)
          sameBlock <- outer(setup$residualBlock, setup$residualBlock, "==")
          residual <- shape[setup$residualIndex, setup$residualIndex, drop=FALSE] * sameBlock
        }else{
          residual <- residualCovariance(local)
        }
        if(setup$useH){
          inverseTransform <- solve(chol(as.matrix(setup$W)))
          residual <- t(inverseTransform) %*% residual %*% inverseTransform
        }
        covariance <- covariance + residual
      }
    }
    covariance
  }
  covariance <- evaluate(working)
  inverse <- solve(covariance)
  X <- as.matrix(setup$X)
  y <- as.numeric(setup$yvar)
  informationFixed <- crossprod(X, inverse %*% X)
  projection <- inverse - inverse %*% X %*% solve(informationFixed, t(X) %*% inverse)
  projected <- projection %*% y
  derivatives <- lapply(seq_along(working), function(parameter){
    h <- 1e-5 * max(1, abs(working[parameter]))
    plus <- minus <- working
    plus[parameter] <- plus[parameter] + h
    minus[parameter] <- minus[parameter] - h
    (evaluate(plus) - evaluate(minus)) / (2 * h)
  })
  traceProjection <- if(reml) projection else inverse
  score <- vapply(derivatives, function(derivative)
    0.5 * (drop(crossprod(projected, derivative %*% projected)) -
      sum(traceProjection * derivative)), numeric(1))
  vectors <- do.call(cbind, lapply(derivatives, function(derivative) derivative %*% projected))
  averageInformation <- 0.5 * crossprod(vectors, projection %*% vectors)
  likelihood <- function(parameters){
    V <- evaluate(parameters)
    Vi <- solve(V)
    information <- crossprod(X, Vi %*% X)
    residual <- y - X %*% solve(information, crossprod(X, Vi %*% y))
    -0.5 * (as.numeric(determinant(V, logarithm=TRUE)$modulus) +
      if(reml) as.numeric(determinant(information, logarithm=TRUE)$modulus) else 0) -
      0.5 * drop(crossprod(residual, Vi %*% residual))
  }
  numericalScore <- vapply(seq_along(working), function(parameter){
    h <- 1e-5 * max(1, abs(working[parameter]))
    plus <- minus <- working
    plus[parameter] <- plus[parameter] + h
    minus[parameter] <- minus[parameter] - h
    (likelihood(plus) - likelihood(minus)) / (2 * h)
  }, numeric(1))
  list(score=score, information=averageInformation, numericalScore=numericalScore,
       likelihood=likelihood, covariance=evaluate)
}

test_that("latent RR analytic scores and AI match independent dense formulas", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  precision <- fixture$precision
  shape <- rrm(dataset$environment, 2L, loadings=fixture$loading)
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
    fitModel <- function(returnParam, blendWeight=0){
      mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision, sigma2=0.9) +
          vsm(ism(replicate), sigma2=0.2), rcov=~vsm(dsm(environment, values=c(0.6, 0.8, 0.5, 0.9)), ism(units),
                  sigma2=0.6), data=dataset, control=mmes_control(task=list(objective=if(reml) "reml" else "ml"), solve=list(method=solver), convergence=list(maxIters=1L), update=list(blendWeight=blendWeight), proposal=list(stepWeight=1), output=list(inverse=0L, setup=returnParam), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    setup <- fitModel(TRUE)
    fit <- fitModel(FALSE)
    working <- unlist(fit$covParWorking, use.names=FALSE)
    reference <- .rr_ai_dense_reference(setup, working, reml)
    expect_equal(reference$score, reference$numericalScore, tolerance=1e-6)
    diagnostics <- fit$engineDiagnostics
    expect_match(diagnostics$optimizer, "AI")
    expect_identical(diagnostics$parameterInformation, "average")
    expect_identical(diagnostics$traceMethod, "shared-covariance-contractions")
    expect_equal(diagnostics$randomTraceDirections, 3)
    expect_equal(as.numeric(diagnostics$scoreWorking), reference$score, tolerance=2e-6)
    expect_equal(as.matrix(diagnostics$averageInformationWorking),
                 reference$information, tolerance=2e-6)
    automatic <- fitModel(FALSE, blendWeight=NULL)
    free <- unlist(lapply(setup$covStruct, function(term) term$free), use.names=FALSE)
    expectedNorm <- max(abs(reference$score[free]) /
      sqrt(pmax(setup$tolParInv, abs(diag(reference$information)[free]))))
    expect_equal(automatic$engineDiagnostics$recovery$finalProjectedScoreNorm,
      expectedNorm, tolerance=2e-6)
    expect_equal(automatic$engineDiagnostics$recovery$steps, 0L)
    expect_equal(automatic$engineDiagnostics$blend$emEvaluations, 0L)
    expect_gt(automatic$engineDiagnostics$blend$latentAIEvaluations, 0L)
  }
})

test_that("RR trace caches select pattern entries and restore general lookups on demand", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  precision <- fixture$precision
  for(zeroLoading in c(FALSE, TRUE)){
    loading <- matrix(c(1, 0.6, if(zeroLoading) 0 else -0.4, 1.1), 4L)
    shape <- rrm(dataset$environment, 1L, loadings=loading)
    fitModel <- function(solver="cholmod", setup=FALSE, inverse=0L){
      mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision, sigma2=0.9) +
          vsm(dsm(environment), ism(id), Gu=precision, sigma2=0.3),
        data=dataset, control=mmes_control(
          solve=list(method=solver, blockSchur="off"),
          output=list(setup=setup, inverse=inverse),
          convergence=list(maxIters=1L), update=list(blendWeight=0),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    fit <- fitModel()
    diagnostics <- fit$engineDiagnostics
    expect_gt(diagnostics$patternInverseCacheBuilds, 0L)
    expect_gt(diagnostics$sparseRandomTracePanels, 0)
    expect_equal(diagnostics$denseRandomTracePanels, 0)
    if(zeroLoading){
      expect_gt(diagnostics$patternInverseCacheFallbacks, 0L)
    }else{
      expect_equal(diagnostics$patternInverseCacheFallbacks, 0L)
    }
    setup <- fitModel(setup=TRUE)
    reference <- .rr_ai_dense_reference(setup, unlist(fit$covParWorking), TRUE)
    expect_equal(as.numeric(diagnostics$scoreWorking), reference$score, tolerance=2e-6)
    expect_equal(diagnostics$averageInformationWorking, reference$information, tolerance=2e-6)
    ldlt <- fitModel(solver="ldlt")
    for(component in c("llik", "bu", "InfMat", "dLu")){
      expect_equal(fit[[component]], ldlt[[component]], tolerance=1e-7)
    }
    full <- fitModel(inverse=2L)
    ldltFull <- fitModel(solver="ldlt", inverse=2L)
    expect_equal(as.matrix(full$Ci), as.matrix(ldltFull$Ci), tolerance=1e-7)
  }
})

test_that("latent RR updates honor step weights and retain an iteration history", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  precision <- fixture$precision
  shape <- rrm(dataset$environment, 1L, loadings=matrix(c(1, 0.6, -0.4, 1.1), 4L))
  fitModel <- function(step, iterations, em=0){
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision, sigma2=0.9), rcov=~vsm(ism(units), sigma2=0.6), data=dataset, control=mmes_control(convergence=list(maxIters=iterations, likelihoodTol=1e-12, parameterTol=1e-10), update=list(blendWeight=rep(em, iterations)), proposal=list(stepWeight=step), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  full <- fitModel(1, 2L)
  half <- fitModel(0.5, 2L)
  initial <- fitModel(1, 1L)
  freeInitial <- unlist(lapply(initial$covStruct, `[[`, "free"), use.names=FALSE)
  direction <- solve(initial$engineDiagnostics$averageInformationWorking[freeInitial, freeInitial],
                     initial$engineDiagnostics$scoreWorking[freeInitial])
  trustScale <- min(1, 2.5 / max(abs(direction)))
  for(fit in list(full, half)){
    progress <- fit$engineDiagnostics$iterationProgress
    expect_equal(unname(progress[2L, "trustAlpha"]), trustScale, tolerance=1e-8)
    expect_equal(as.numeric(fit$monitor[freeInitial, 2L] - fit$monitor[freeInitial, 1L]),
                 as.numeric(direction * trustScale * progress[2L, "lineAlpha"]),
                 tolerance=1e-8)
  }
  expect_gt(ncol(full$monitor), 1L)
  expect_equal(ncol(full$monitor), length(full$llik))
  expect_equal(ncol(half$monitor), length(half$llik))
  expect_gt(max(abs(unlist(full$covPar) - unlist(half$covPar))), 1e-7)
  stabilized <- fitModel(1, 2L, em=1)
  expect_gt(max(abs(unlist(full$covPar) - unlist(stabilized$covPar))), 1e-7)
  converged <- fitModel(1, 100L)
  expect_gt(length(converged$llik), 2L)
  expect_true(all(diff(as.numeric(converged$llik)) >= -1e-7))
  expect_true(all(is.finite(converged$theta_se)))
  expect_equal(ncol(converged$normMonitor), ncol(converged$monitor))
  expect_gt(diff(range(converged$normMonitor[2L, ])), 1e-7)
  free <- unlist(lapply(converged$covStruct, `[[`, "free"), use.names=FALSE)
  information <- converged$engineDiagnostics$averageInformationWorking
  workingCovariance <- matrix(0, length(free), length(free))
  workingCovariance[free, free] <- solve(information[free, free, drop=FALSE])
  natural <- unlist(converged$covPar, use.names=FALSE)
  jacobian <- rep(1, length(natural))
  kinds <- converged$vcParams$kind
  jacobian[kinds %in% c("scale", "exp")] <- natural[kinds %in% c("scale", "exp")]
  expected <- diag(jacobian) %*% workingCovariance %*% diag(jacobian)
  expect_equal(as.matrix(converged$theta_se), expected, tolerance=2e-5)
})

test_that("latent AI differentiates covariance-product axes and correlated residuals", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  dataset$record <- interaction(dataset$id, dataset$replicate, drop=TRUE)
  precision <- fixture$precision
  shape <- rrm(dataset$environment, 1L, loadings=matrix(c(1, 0.6, -0.4, 1.1), 4L))
  for(position in 1:2){
    random <- if(position == 1L){
      ~vsm(shape, ar1m(replicate, rho=0.2), ism(id), Gu=precision, sigma2=0.9)
    }else{
      ~vsm(ar1m(replicate, rho=0.2), shape, ism(id), Gu=precision, sigma2=0.9)
    }
    fitModel <- function(returnParam){
      mmes(response~environment, random=random, rcov=~vsm(ar1m(environment, rho=0.25), ism(record), sigma2=0.6), data=dataset, control=mmes_control(output=list(setup=returnParam), convergence=list(maxIters=1L), proposal=list(stepWeight=1), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    setup <- fitModel(TRUE)
    fit <- fitModel(FALSE)
    reference <- .rr_ai_dense_reference(setup,
      unlist(fit$covParWorking, use.names=FALSE), TRUE)
    expect_equal(reference$score, reference$numericalScore, tolerance=1e-6)
    expect_equal(as.numeric(fit$engineDiagnostics$scoreWorking), reference$score,
                 tolerance=2e-6)
    expect_equal(as.matrix(fit$engineDiagnostics$averageInformationWorking),
                 reference$information, tolerance=2e-6)
  }
})

test_that("latent AI scores and information honor residual weights", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  precision <- fixture$precision
  shape <- rrm(dataset$environment, 1L, loadings=matrix(c(1, 0.6, -0.4, 1.1), 4L))
  weights <- diag(seq(0.8, 1.2, length.out=nrow(dataset)))
  fitModel <- function(returnParam, weights, heterogeneous){
    residual <- if(heterogeneous)
      ~vsm(dsm(environment, values=c(0.6, 0.8, 0.5, 0.9)), ism(units), sigma2=0.6)
    else ~vsm(ism(units), sigma2=0.6)
    mmes(response~environment, random=~vsm(shape, ism(id), Gu=precision, sigma2=0.9), rcov=residual, W=weights, data=dataset, control=mmes_control(output=list(setup=returnParam), convergence=list(maxIters=1L), update=list(blendWeight=0), proposal=list(stepWeight=1), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  correlated <- weights + 0.01 * 0.5^abs(outer(seq_len(nrow(dataset)),
                                             seq_len(nrow(dataset)), "-"))
  for(weightMatrix in list(weights, correlated)) for(heterogeneous in c(FALSE, TRUE)){
    setup <- fitModel(TRUE, weightMatrix, heterogeneous)
    fit <- fitModel(FALSE, weightMatrix, heterogeneous)
    reference <- .rr_ai_dense_reference(setup,
      unlist(fit$covParWorking, use.names=FALSE), TRUE)
    expect_equal(reference$score, reference$numericalScore, tolerance=1e-6)
    expect_equal(as.numeric(fit$engineDiagnostics$scoreWorking), reference$score,
                 tolerance=2e-6)
    expect_equal(as.matrix(fit$engineDiagnostics$averageInformationWorking),
                 reference$information, tolerance=2e-6)
    if(identical(weightMatrix, weights)){
      expect_gt(fit$engineDiagnostics$sparseRandomTracePanels, 0)
      expect_equal(fit$engineDiagnostics$denseRandomTracePanels, 0)
    }else{
      expect_equal(fit$engineDiagnostics$sparseRandomTracePanels, 0)
      expect_gt(fit$engineDiagnostics$denseRandomTracePanels, 0)
    }
  }
})

test_that("latent AI matches independently optimized ML", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  precision <- fixture$precision
  fitModel <- function(returnParam){
    mmes(response~environment, random=~vsm(rrm(environment, 1L), ism(id), Gu=precision), rcov=~units, data=dataset, control=mmes_control(task=list(objective=if(FALSE) "reml" else "ml"), convergence=list(maxIters=150L, likelihoodTol=1e-10, parameterTol=1e-8), output=list(setup=returnParam, inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  setup <- fitModel(TRUE)
  fit <- fitModel(FALSE)
  working <- unlist(fit$covParWorking, use.names=FALSE)
  reference <- .rr_ai_dense_reference(setup, working, FALSE)
  free <- unlist(lapply(setup$covStruct, `[[`, "free"), use.names=FALSE)
  initial <- unlist(lapply(setup$covStruct, `[[`, "par"), use.names=FALSE)
  objective <- function(parameters){
    current <- initial
    current[free] <- parameters
    -reference$likelihood(current)
  }
  optimum <- optim(initial[free], objective, method="BFGS",
                   control=list(reltol=1e-12, maxit=500L))
  expected <- initial
  expected[free] <- optimum$par
  expect_identical(optimum$convergence, 0L)
  expect_equal(working[free], expected[free], tolerance=2e-4)
  expect_lt(max(abs(as.numeric(fit$engineDiagnostics$scoreWorking) -
                      reference$score)), 1e-8)
  expect_equal(as.matrix(fit$engineDiagnostics$averageInformationWorking),
               reference$information, tolerance=2e-6)
})

test_that("diagonal residual contractions are exact across layouts, solvers and weights", {
  fixture <- .rr_ai_fixture()
  dataset <- fixture$data
  dataset <- dataset[rev(seq_len(nrow(dataset))), ]
  precision <- fixture$precision
  shape <- rrm(dataset$environment, 1L, loadings=fixture$loading[, 1L, drop=FALSE])
  diagonal <- diag(seq(0.7, 1.4, length.out=nrow(dataset)))
  correlated <- diagonal + 0.02 * 0.4^abs(outer(seq_len(nrow(dataset)),
                                               seq_len(nrow(dataset)), "-"))
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(residual in list(~dsumm(vsm(ism(units), sigma2=0.6), by=environment),
                        ~vsm(dsm(environment), ism(units), sigma2=0.6))){
      for(weights in list(NULL, diagonal, correlated)){
        fitModel <- function(setup, blend=0){
          arguments <- list(fixed=response~environment,
            random=~vsm(shape, ism(id), Gu=precision, sigma2=0.9) +
              vsm(csm(replicate, rho=0.2), sigma2=0.2),
            rcov=residual, data=dataset,
            control=mmes_control(solve=list(method=solver),
              task=list(objective=if(reml) "reml" else "ml"),
              convergence=list(maxIters=1L), output=list(setup=setup),
              update=list(blendWeight=blend),
              diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
          if(!is.null(weights)) arguments$W <- weights
          do.call(mmes, arguments)
        }
        setup <- fitModel(TRUE)
        fit <- fitModel(FALSE)
        diagnostics <- fit$engineDiagnostics
        reference <- .rr_ai_dense_reference(setup,
          unlist(fit$covParWorking, use.names=FALSE), reml)
        expect_equal(reference$score, reference$numericalScore, tolerance=1e-6)
        expect_equal(as.numeric(diagnostics$scoreWorking), reference$score, tolerance=2e-6)
        expect_equal(as.matrix(diagnostics$averageInformationWorking),
          reference$information, tolerance=2e-6)
        eligible <- is.null(weights) || identical(weights, diagonal)
        expect_equal(diagnostics$diagonalResidualTraceCalls, as.numeric(eligible))
        if(eligible){
          expect_lte(diagnostics$diagonalResidualTraceSolveColumns,
            diagnostics$latentCoefficients)
          if(reml && solver == "cholmod" && diagnostics$blockSchurActive){
            expect_equal(diagnostics$diagonalResidualTraceSolveColumns, 0)
          }
        }
        automatic <- fitModel(FALSE, NULL)
        expect_equal(as.numeric(automatic$engineDiagnostics$scoreWorking),
          reference$score, tolerance=2e-6)
        expect_gte(automatic$engineDiagnostics$diagonalResidualTraceCalls,
          as.numeric(eligible))
      }
    }
  }
})

test_that("latent RR sectioned spatial scores and AI match dense observation formulas", {
  set.seed(729)
  dataset <- expand.grid(environment=factor(letters[1:3]),
    row=factor(seq_len(3)), range=factor(seq_len(3)))
  dataset$id <- factor(rep(seq_len(5), length.out=nrow(dataset)))
  dataset$response <- rnorm(nrow(dataset))
  dataset <- dataset[sample(setdiff(seq_len(nrow(dataset)), c(2L, 10L, 23L))), ]
  distanceRow <- abs(outer(as.integer(dataset$row), as.integer(dataset$row), "-"))
  distanceRange <- abs(outer(as.integer(dataset$range), as.integer(dataset$range), "-"))
  environment <- as.integer(dataset$environment)
  for(selected in list(letters[1:3], c("c", "a"))){
    residualCovariance <- function(parameters){
      scales <- exp(c(parameters[1L], parameters[1L] + parameters[2:3]))
      nSelected <- length(selected)
      rangeRho <- tanh(parameters[3L + seq_len(nSelected)])
      rowRho <- tanh(parameters[3L + nSelected + seq_len(nSelected)])
      result <- matrix(0, nrow(dataset), nrow(dataset))
      for(section in seq_len(3)){
        rows <- which(environment == section)
        listed <- match(letters[section], selected)
        result[rows, rows] <- scales[section] * if(is.na(listed)){
          diag(length(rows))
        }else{
          rangeRho[listed]^distanceRange[rows, rows] *
            rowRho[listed]^distanceRow[rows, rows]
        }
      }
      result
    }
    for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
      fitModel <- function(returnParam){
        mmes(response~environment, random=~vsm(rrm(environment, 1L), ism(id), sigma2=0.9), rcov=~dsumm(vsm(ar1m(range, rho=0.2), ar1m(row, rho=-0.15),
                          ism(units), sigma2=0.6),
                      by=environment, levels=selected), data=dataset, control=mmes_control(solve=list(method=solver), task=list(objective=if(reml) "reml" else "ml"), convergence=list(maxIters=3L), output=list(inverse=0L, setup=returnParam), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      }
      setup <- fitModel(TRUE)
      fit <- fitModel(FALSE)
      reference <- .rr_ai_dense_reference(setup,
        unlist(fit$covParWorking, use.names=FALSE), reml, residualCovariance)
      expect_equal(reference$score, reference$numericalScore, tolerance=1e-6)
      expect_equal(as.numeric(fit$engineDiagnostics$scoreWorking),
                   reference$score, tolerance=2e-6)
      expect_equal(as.matrix(fit$engineDiagnostics$averageInformationWorking),
                   reference$information, tolerance=2e-6)
      expect_true(all(diff(as.numeric(fit$llik)) >= -1e-7))
      expect_identical(fit$engineDiagnostics$residualEvaluation, "observed-blocks")
      expect_identical(fit$engineDiagnostics$residualEngine, "grid")
      expect_true(fit$engineDiagnostics$sharedConditionalOperators)
      expect_true(all(is.finite(fit$engineDiagnostics$aiElapsedSeconds)))
      expect_lte(fit$engineDiagnostics$largestResidualBlock, 9)
    }
  }
})

test_that("latent RR residual workspace follows observed blocks, not the full grid", {
  set.seed(730)
  dataset <- expand.grid(environment=factor(letters[1:9]), cell=seq_len(30))
  dataset$row <- factor((dataset$cell - 1L) %% 5L + 1L +
                          (as.integer(dataset$environment) - 1L) * 5L)
  dataset$range <- factor((dataset$cell - 1L) %/% 5L + 1L +
                            (as.integer(dataset$environment) - 1L) * 6L)
  dataset$id <- factor((seq_len(nrow(dataset)) - 1L) %% 7L + 1L)
  dataset$response <- rnorm(nrow(dataset))
  fit <- mmes(response~environment, random=~vsm(rrm(environment, 1L), ism(id)), rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=environment), data=dataset, control=mmes_control(convergence=list(maxIters=1L), output=list(inverse=0L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(tail(fit$covStruct, 1L)[[1L]]$dim, 21870L)
  expect_equal(fit$engineDiagnostics$largestResidualBlock, 30)
  expect_identical(fit$engineDiagnostics$residualEvaluation, "observed-blocks")
  expect_true(all(is.finite(fit$engineDiagnostics$scoreWorking)))
  expect_true(all(is.finite(fit$engineDiagnostics$averageInformationWorking)))
})
