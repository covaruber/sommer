.em_direction_data <- function(){
  set.seed(152)
  dataset <- expand.grid(id=factor(seq_len(8L)), trait=factor(seq_len(3L)),
                         replicate=factor(seq_len(2L)))
  effects <- matrix(rnorm(24L), 8L, 3L)
  dataset$y <- as.numeric(dataset$trait) / 3 +
    effects[cbind(as.integer(dataset$id), as.integer(dataset$trait))] +
    rnorm(nrow(dataset), sd=0.7)
  dataset
}

.em_information_at <- function(setup, fit){
  covariance <- setup$covStruct
  attr(covariance, "directionBlend") <- FALSE
  attr(covariance, "emInformationOnly") <- TRUE
  for(term in seq_along(covariance)){
    covariance[[term]]$par <- as.numeric(fit$covParWorking[[term]])
  }
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    0, 1, FALSE, computeCi=0L, solver=setup$solver, reml=setup$REML,
    responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
    preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
}

.em_dense_step <- function(setup, response){
  X <- as.matrix(setup$X)
  blocks <- lapply(setup$Z, as.matrix)
  design <- do.call(cbind, c(list(X), blocks))
  covariance <- lapply(setup$covStruct, function(term){
    sommer:::evaluate_covstruct_cpp(term, term$par)
  })
  residual <- covariance[[length(covariance)]]
  residualCovariance <- matrix(0, length(response), length(response))
  residualRows <- split(seq_along(response), setup$residualBlock)
  for(rows in residualRows){
    index <- setup$residualIndex[rows]
    residualCovariance[rows, rows] <- residual[index, index, drop=FALSE]
  }
  transform <- if(setup$useH) t(chol(as.matrix(setup$W))) else diag(length(response))
  residualPrecision <- crossprod(transform, solve(residualCovariance) %*% transform)
  prior <- matrix(0, ncol(design), ncol(design))
  termIndices <- list()
  first <- ncol(X) + 1L
  for(term in seq_along(setup$Ai)){
    selected <- which(setup$Zind == term)
    levels <- ncol(blocks[[selected[1L]]])
    indices <- first:(first + length(selected) * levels - 1L)
    prior[indices, indices] <- kronecker(solve(covariance[[term]]),
                                        as.matrix(setup$Ai[[term]]))
    termIndices[[term]] <- matrix(indices, levels)
    first <- max(indices) + 1L
  }
  posterior <- solve(crossprod(design, residualPrecision %*% design) + prior)
  mean <- drop(posterior %*% crossprod(design, residualPrecision %*% response))
  if(!setup$REML){
    random <- seq.int(ncol(X) + 1L, ncol(design))
    conditional <- crossprod(design, residualPrecision %*% design) + prior
    posterior[,] <- 0
    posterior[random, random] <- solve(conditional[random, random, drop=FALSE])
  }
  updated <- lapply(seq_along(termIndices), function(term){
    indices <- termIndices[[term]]
    U <- matrix(mean[indices], nrow(indices))
    precision <- as.matrix(setup$Ai[[term]])
    moment <- crossprod(U, precision %*% U)
    for(i in seq_len(ncol(indices))) for(j in seq_len(ncol(indices))){
      moment[i, j] <- moment[i, j] +
        sum(precision * posterior[indices[, i], indices[, j], drop=FALSE])
    }
    moment / nrow(indices)
  })
  error <- drop(transform %*% (response - drop(design %*% mean)))
  transformedDesign <- transform %*% design
  moment <- matrix(0, nrow(residual), ncol(residual))
  for(rows in residualRows){
    index <- setup$residualIndex[rows]
    local <- transformedDesign[rows, , drop=FALSE]
    moment[index, index] <- moment[index, index, drop=FALSE] +
      tcrossprod(error[rows]) + local %*% posterior %*% t(local)
  }
  updated[[length(updated) + 1L]] <- moment / length(residualRows)
  updated
}

.em_native_fit <- function(setup, generic=FALSE){
  covariance <- setup$covStruct
  attr(covariance, "emGenericMoments") <- generic
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    2L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    1, 1, FALSE, computeCi=0L, solver=setup$solver, reml=setup$REML,
    responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
    preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
}

test_that("cached EM and PX moments match generic panels across covariance shapes", {
  dataset <- .em_direction_data()
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")) for(method in c("em", "pxem")){
    for(shape in c("diagonal", "us", "fixed", "weighted", "weighted-panels",
                  "correlated-residual")){
      arguments <- list(fixed=y~trait,
        random=~vsm(dsm(trait), ism(id), sigma2=0.8),
        rcov=~vsm(ism(units), sigma2=0.6),
        data=dataset,
        control=mmes_control(solve=list(method=solver),
          task=list(objective=if(reml) "reml" else "ml"),
          update=list(recovery=method, blendWeight=1),
          proposal=list(stepWeight=1), convergence=list(maxIters=2L),
          output=list(setup=TRUE),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      if(shape == "us") arguments$random <- ~vsm(usm(trait), ism(id), sigma2=0.8)
      if(shape == "fixed"){
        arguments$random <- ~vsm(ism(id), sigma2=0.8, fixedSigma2=TRUE) +
          vsm(ism(replicate), sigma2=0.3)
      }
      if(shape %in% c("weighted", "weighted-panels")){
        n <- nrow(dataset)
        arguments$W <- Matrix::Matrix(outer(seq_len(n), seq_len(n),
          function(i, j) 0.25^abs(i-j)), sparse=TRUE)
      }
      if(shape == "weighted-panels"){
        arguments$control <- do.call(mmes_control, utils::modifyList(
          unclass(arguments$control), list(execution=list(denseMemoryMB=0.08, memoryPolicy="limit")),
          keep.null=TRUE))
      }
      if(shape == "correlated-residual"){
        arguments$rcov <- ~vsm(ar1m(trait, rho=0.2), ism(units), sigma2=0.6)
      }
      setup <- do.call(mmes, arguments)
      cached <- .em_native_fit(setup)
      generic <- .em_native_fit(setup, generic=TRUE)
      for(field in c("theta", "InfMat", "bu", "llik")){
        tolerance <- if(method == "em" && shape != "correlated-residual") 1e-9 else 2e-6
        expect_equal(cached[[field]], generic[[field]], tolerance=tolerance,
          info=paste(solver, method, shape, field))
      }
      diagnostics <- cached$engineDiagnostics$blend
      expect_gt(cached$engineDiagnostics$emInverseTraceCalls, 0)
      if(shape != "correlated-residual"){
        expect_equal(diagnostics$scalarResidualTraceSteps, 1L)
        expect_lt(diagnostics$momentPanelSolves, generic$engineDiagnostics$blend$momentPanelSolves)
      } else {
        expect_equal(diagnostics$scalarResidualTraceSteps, 0L)
      }
      if(shape == "diagonal"){
        expect_equal(diagnostics$randomMomentPairs, 3)
        expect_equal(diagnostics$randomOffDiagonalPairsSkipped, 3)
        expect_equal(cached$engineDiagnostics$emInverseTraceFallbackPanels, 0)
        expect_equal(diagnostics$momentPanelSolves, 0L)
        if(method == "pxem") expect_gt(diagnostics$scalarResidualContractions, 0L)
      }
      if(shape == "weighted-panels" && method == "pxem"){
        expect_equal(diagnostics$scalarResidualContractions, 1L)
        expect_gt(diagnostics$momentPanelSolves, 0L)
      }
    }
  }
})

test_that("one true EM step matches independent dense REML and ML moments", {
  dataset <- .em_direction_data()
  for(reml in c(TRUE, FALSE)) for(shape in c("scalar", "us", "diag")){
    random <- switch(shape,
      scalar=~vsm(ism(id), sigma2=0.8) + vsm(ism(replicate), sigma2=0.3),
      us=~vsm(usm(trait), ism(id), sigma2=0.8),
      diag=~vsm(dsm(trait), ism(id), sigma2=0.8))
    arguments <- list(fixed=y~trait, random=random,
      rcov=~vsm(ism(units), sigma2=0.6), data=dataset,
      control=mmes_control(convergence=list(maxIters=2L),
        update=list(blendWeight=1),
        proposal=list(stepWeight=1),
        task=list(objective=if(reml) "reml" else "ml"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    setupArguments <- arguments
    setupArguments$control <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(output=list(setup=TRUE)),
      keep.null=TRUE))
    setup <- do.call(mmes, setupArguments)
    reference <- .em_dense_step(setup, dataset$y)
    if(shape == "diag") reference[[1L]] <- diag(diag(reference[[1L]]))
    fitted <- do.call(mmes, arguments)
    for(term in seq_along(reference)){
      expect_equal(fitted$theta[[term]], reference[[term]], tolerance=1e-8,
                   ignore_attr=TRUE)
    }
    expect_equal(fitted$engineDiagnostics$blend$closedFormSteps, 1L)
    expect_equal(fitted$engineDiagnostics$blend$generalizedSteps, 0L)
    expect_equal(fitted$engineDiagnostics$blend$aiEvaluations, 0L)
    expect_true(fitted$engineDiagnostics$blend$finalInformationExact)
  }
})

test_that("ML trials reprofile fixed effects and preserve D-first and inverse-mode results", {
  dataset <- .em_direction_data()
  X <- model.matrix(~trait, dataset)
  Z <- model.matrix(~id-1, dataset)
  weights <- Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
    diagonals=list(rep(0.1, nrow(dataset) - 1L), rep(1.5, nrow(dataset)),
      rep(0.1, nrow(dataset) - 1L)))
  transform <- t(chol(as.matrix(weights)))
  residualCovariance <- solve(crossprod(transform))
  for(solver in c("ldlt", "cholmod")) for(method in c("em", "pxem")){
    reference <- NULL
    for(path in c("off", "force")) for(mode in 0:2){
      fit <- mmes(y~trait, random=~vsm(ism(id), sigma2=0.8),
        rcov=~vsm(ism(units), sigma2=0.6), data=dataset, W=weights,
        control=mmes_control(solve=list(method=solver, dFirstML=path),
          task=list(objective="ml"), update=list(recovery=method,
            blendWeight=c(1, 0.5, 0)),
          convergence=list(maxIters=3L), output=list(inverse=mode),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      V <- fit$theta[[1L]][1L] * tcrossprod(Z) +
        fit$theta[[2L]][1L] * residualCovariance
      ViX <- solve(V, X)
      beta <- drop(solve(crossprod(X, ViX), crossprod(ViX, dataset$y)))
      expect_equal(as.numeric(fit$bu[seq_len(ncol(X)), ]), as.numeric(beta), tolerance=1e-9)
      expect_true(all(diff(as.numeric(fit$llik)) >= -1e-10))
      expect_identical(fit$engineDiagnostics$dFirstMLActive, path == "force" && mode != 1L)
      expect_equal(fit$engineDiagnostics$DsymbolicAnalyses, 1)
      if(is.null(reference)) reference <- fit
      for(field in c("theta", "bu", "llik", "InfMat", "dLu", "theta_se")){
        expect_equal(fit[[field]], reference[[field]], tolerance=1e-8,
          info=paste(solver, method, path, mode, field))
      }
      if(mode == 2L){
        expect_equal(as.matrix(fit$C) %*% as.matrix(fit$Ci),
          diag(nrow(fit$C)) * fit$Cscale, tolerance=1e-8)
      }
    }
  }
})

test_that("generic covariance and tied M-steps increase Q and preserve constraints", {
  dataset <- .em_direction_data()
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")) for(model in c("ar", "residual", "fixed", "tie", "rr")){
    arguments <- list(fixed=y~trait, random=~id, data=dataset,
      control=mmes_control(solve=list(method=solver),
        task=list(objective=if(reml) "reml" else "ml"),
        update=list(blendWeight=1), convergence=list(maxIters=4L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    if(model == "ar") arguments$random <- ~vsm(ar1m(trait), ism(id))
    if(model == "residual") arguments$rcov <- ~vsm(ar1m(trait), ism(units))
    if(model == "fixed") arguments$random <- ~vsm(ar1m(trait, rho=0.2, fixed=TRUE), ism(id))
    if(model == "tie"){
      arguments$random <- ~id + replicate
      arguments$vcc <- data.frame(parameter=c(1L, 2L), group=1L)
    }
    if(model == "rr") arguments$random <- ~vsm(rrm(trait, 1L), ism(id))
    fitted <- do.call(mmes, arguments)
    diagnostics <- fitted$engineDiagnostics$blend
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
    expect_true(all(diagnostics$qGains >= -1e-9))
    expect_equal(diagnostics$aiEvaluations, 0L)
    expect_equal(diagnostics$finalAIRefreshes, 1L)
    if(model %in% c("ar", "residual", "tie", "rr")){
      expect_gt(diagnostics$generalizedSteps, 0L)
      expect_length(diagnostics$numericalMstepConverged, ncol(fitted$llik))
      expect_gt(sum(diagnostics$numericalMstepEvaluations), 0)
      expect_equal(sum(diagnostics$numericalMstepConverged, na.rm=TRUE),
                   diagnostics$numericalMstepConvergences)
    }
    if(model == "fixed") expect_equal(fitted$covPar[[1L]][2L], 0.2)
    if(model == "tie") expect_equal(fitted$covParWorking[[1L]], fitted$covParWorking[[2L]])
  }
})

test_that("weighted relationship EM matches dense moments without a full inverse", {
  dataset <- .em_direction_data()
  precision <- diag(8L) + matrix(0.04, 8L, 8L)
  dimnames(precision) <- list(levels(dataset$id), levels(dataset$id))
  attr(precision, "inverse") <- TRUE
  weightMatrices <- list(
    Matrix::Diagonal(nrow(dataset), seq(0.7, 1.3, length.out=nrow(dataset))),
    Matrix::bdiag(replicate(nrow(dataset) / 2L,
      matrix(c(1.4, -0.2, -0.2, 1.4), 2L), simplify=FALSE)),
    Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
      diagonals=list(rep(0.15, nrow(dataset) - 1L),
        seq(1, 2, length.out=nrow(dataset)), rep(0.15, nrow(dataset) - 1L))))
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")) for(weights in weightMatrices){
    arguments <- list(fixed=y~trait,
      random=~vsm(ism(id), Gu=precision, sigma2=0.8),
      rcov=~vsm(ism(units), sigma2=0.6), W=weights, data=dataset,
      control=mmes_control(convergence=list(maxIters=2L),
        update=list(blendWeight=1),
        task=list(objective=if(reml) "reml" else "ml"),
        solve=list(method=solver),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    setupArguments <- arguments
    setupArguments$control <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(output=list(setup=TRUE)),
      keep.null=TRUE))
    setup <- do.call(mmes, setupArguments)
    reference <- .em_dense_step(setup, dataset$y)
    fitted <- do.call(mmes, arguments)
    for(term in seq_along(reference)){
      expect_equal(fitted$theta[[term]], reference[[term]], tolerance=1e-8, ignore_attr=TRUE)
    }
    expect_equal(dim(fitted$llik), c(1L, 2L))
    expect_equal(ncol(fitted$monitor), ncol(fitted$llik))
  }
})

test_that("all blend endpoints evaluate the same correlated-weight model", {
  dataset <- .em_direction_data()
  n <- nrow(dataset)
  weights <- Matrix::bandSparse(n, k=c(-1L, 0L, 1L),
    diagonals=list(rep(0.15, n - 1L), seq(1, 2, length.out=n), rep(0.15, n - 1L)))
  arguments <- list(fixed=y~trait, random=~vsm(ism(id), sigma2=0.8),
    rcov=~vsm(ism(units), sigma2=0.6), data=dataset, W=weights,
    control=mmes_control(convergence=list(maxIters=1L),
      proposal=list(stepWeight=1), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  fitWith <- function(weight){
    fittedArgs <- arguments
    fittedArgs$control <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(update=list(blendWeight=weight)),
      keep.null=TRUE))
    do.call(mmes, fittedArgs)
  }
  reference <- fitWith(0)
  for(weight in c(0.4, 1)){
    fitted <- fitWith(weight)
    for(name in c("llik", "bu", "InfMat", "InfMatInv")){
      expect_equal(fitted[[name]], reference[[name]], tolerance=1e-9)
    }
  }
})

test_that("complete unstructured and diagonal residuals have exact EM M-steps", {
  dataset <- .em_direction_data()
  for(reml in c(TRUE, FALSE)) for(shape in c("us", "diag")){
    arguments <- list(fixed=y~trait, random=~id, data=dataset,
      rcov=if(shape == "us") ~vsm(usm(trait), ism(units)) else
        ~vsm(dsm(trait), ism(units)),
      control=mmes_control(convergence=list(maxIters=2L),
        update=list(blendWeight=1),
        task=list(objective=if(reml) "reml" else "ml"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    setupArguments <- arguments
    setupArguments$control <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(output=list(setup=TRUE)),
      keep.null=TRUE))
    setup <- do.call(mmes, setupArguments)
    reference <- .em_dense_step(setup, dataset$y)
    if(shape == "diag") reference[[2L]] <- diag(diag(reference[[2L]]))
    fitted <- do.call(mmes, arguments)
    expect_equal(as.numeric(fitted$theta[[1L]]), as.numeric(reference[[1L]]), tolerance=1e-8)
    expect_equal(fitted$covParWorking[[2L]][1L], log(reference[[2L]][1L, 1L]), tolerance=1e-8)
    covariance <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[2L]],
      as.numeric(fitted$covParWorking[[2L]]))
    expect_equal(covariance, reference[[2L]], tolerance=1e-8, ignore_attr=TRUE)
    expect_equal(fitted$engineDiagnostics$blend$closedFormSteps, 1L)
    expect_equal(fitted$engineDiagnostics$blend$generalizedSteps, 0L)
  }
})

test_that("EM fitting covers covariance-factor families and ragged sectioned residuals", {
  set.seed(642)
  dataset <- expand.grid(id=factor(1:5), site=factor(letters[1:4]), replicate=1:2)
  dataset$y <- rnorm(nrow(dataset))
  coordinates <- cbind(as.integer(dataset$site) %% 2, as.integer(dataset$site) %/% 2)
  adjacency <- matrix(1, 4L, 4L) - diag(4L)
  dimnames(adjacency) <- list(levels(dataset$site), levels(dataset$site))
  structures <- list(
    ism(dataset$site), dsm(dataset$site), usm(dataset$site),
    csm(dataset$site), csm(dataset$site, variance="heterogeneous"),
    ar1m(dataset$site, variance="heterogeneous"),
    ar2m(dataset$site), ar3m(dataset$site),
    ar2m(dataset$site, variance="heterogeneous"),
    ar3m(dataset$site, variance="heterogeneous"),
    mam(dataset$site), corgm(dataset$site), fam(dataset$site, 1L),
    antem(dataset$site), toeplitzm(dataset$site),
    sar(dataset$site, adjacency), car(dataset$site, adjacency),
    maternm(coordinates), metricm(coordinates),
    maternm(coordinates, anisotropy="geometric"),
    metricm(coordinates, model="power", anisotropy="product"),
    ownm(dataset$site, fun=function(par){
      tanh(par[1L])^abs(outer(1:4, 1:4, "-"))
    }, par=0.2, par_names="correlation"))
  for(shape in structures){
    fitted <- mmes(y~1, random=~vsm(shape, ism(id)), data=dataset, control=mmes_control(convergence=list(maxIters=2L), update=list(blendWeight=1), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_true(all(is.finite(fitted$llik)), info=shape$covFactor$model)
    expect_true(diff(c(fitted$llik)) >= -1e-9, info=shape$covFactor$model)
    expect_equal(fitted$engineDiagnostics$blend$aiEvaluations, 0L,
                 info=shape$covFactor$model)
  }
  dataset <- dataset[-c(1L, 9L), , drop=FALSE]
  dataset$replicate <- factor(dataset$replicate)
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
  fitted <- mmes(y~replicate, random=~id, rcov=~dsumm(vsm(ar1m(site), ism(units)), by=replicate), data=dataset, control=mmes_control(convergence=list(maxIters=3L), update=list(blendWeight=1), task=list(objective=if(reml) "reml" else "ml"), solve=list(method=solver), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(all(diff(c(fitted$llik)) >= -1e-9))
  expect_true(all(fitted$engineDiagnostics$blend$qGains >= -1e-9))
  }
})

test_that("mixed schedules skip unused EM or AI work and return clean inference", {
  dataset <- .em_direction_data()
  for(reml in c(TRUE, FALSE)) for(solver in c("ldlt", "cholmod")){
  fit <- function(weights, iterations=length(weights), computeCi=0L){
    mmes(y~trait, random=~id, data=dataset, control=mmes_control(convergence=list(maxIters=iterations), update=list(blendWeight=weights), output=list(inverse=computeCi), task=list(objective=if(reml) "reml" else "ml"), solve=list(method=solver), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  mixed <- fit(c(1, 0.5, 0, 0))
  diagnostics <- mixed$engineDiagnostics$blend
  expect_equal(diagnostics$emEvaluations, 2L)
  expect_equal(diagnostics$aiEvaluations, 2L)
  expect_true(all(diff(c(mixed$llik)) >= -1e-9))
  for(mode in 1:2){
    uncertainty <- fit(c(1, 0.5, 0, 0), computeCi=mode)
    expect_equal(uncertainty$theta_se, mixed$theta_se, tolerance=1e-9)
    expect_equal(uncertainty$theta, mixed$theta, tolerance=1e-9)
  }
  one <- fit(1)
  expect_equal(one$engineDiagnostics$blend$emEvaluations, 0L)
  expect_equal(one$engineDiagnostics$blend$aiEvaluations, 0L)
  expect_false(one$convergence)
  ai <- fit(0, iterations=4L)
  expect_false(isTRUE(ai$engineDiagnostics$blend$active))
  }
})

test_that("zero-weight iterations reuse native AI and match a pure-AI restart", {
  dataset <- .em_direction_data()
  for(solver in c("ldlt", "cholmod")) for(method in c("em", "pxem")){
    for(shape in c("diagonal", "us", "correlated-residual", "ragged-residual")){
      arguments <- list(fixed=y~trait,
        random=~vsm(dsm(trait), ism(id), sigma2=0.8),
        rcov=~vsm(ism(units), sigma2=0.6), data=dataset,
        control=mmes_control(solve=list(method=solver),
          convergence=list(maxIters=4L),
          update=list(recovery=method, blendWeight=c(1, 0, 0, 0)),
          proposal=list(stepWeight=1),
          diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)))
      if(shape == "us") arguments$random <- ~vsm(usm(trait), ism(id), sigma2=0.8)
      if(shape %in% c("correlated-residual", "ragged-residual")){
        arguments$rcov <- ~vsm(ar1m(trait, rho=0.2), ism(units), sigma2=0.6)
      }
      if(shape == "ragged-residual") arguments$data <- dataset[-1L, , drop=FALSE]
      setupArguments <- arguments
      setupArguments$control <- do.call(mmes_control, utils::modifyList(
        unclass(arguments$control), list(output=list(setup=TRUE)),
        keep.null=TRUE))
      setup <- do.call(mmes, setupArguments)
      prefixArguments <- arguments
      prefixArguments$control <- do.call(mmes_control, utils::modifyList(
        unclass(arguments$control),
        list(convergence=list(maxIters=2L), update=list(blendWeight=c(1, 0))),
        keep.null=TRUE))
      prefix <- do.call(mmes, prefixArguments)
      covariance <- setup$covStruct
      attr(covariance, "directionBlend") <- FALSE
      for(term in seq_along(covariance)){
        covariance[[term]]$par <- as.numeric(prefix$covParWorking[[term]])
      }
      reference <- sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai,
        setup$yvar, setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
        3L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
        rep(0, 3L), rep(1, 3L), FALSE, computeCi=0L, solver=setup$solver,
        reml=TRUE, responsePrepared=setup$responsePrepared,
        preparedMean=setup$preparedMean, preparedSd=setup$preparedSd,
        preparedIntercept=setup$preparedIntercept)
      fitted <- do.call(mmes, arguments)
      for(field in c("theta", "bu", "InfMat", "dLu")){
        expect_equal(fitted[[field]], reference[[field]], tolerance=2e-6,
          ignore_attr=TRUE, info=paste(solver, method, shape, field))
      }
      expect_equal(fitted$llik[, -1L, drop=FALSE], reference$llik, tolerance=2e-6)
      diagnostics <- fitted$engineDiagnostics$blend
      expect_equal(diagnostics$nativeAIEvaluations, 2L)
      expect_equal(diagnostics$latentAIEvaluations, 0L)
      expect_equal(diagnostics$emEvaluations, 1L)
      expect_gte(diagnostics$acceptedAIFitReuses, 2L)
      expect_equal(diagnostics$likelihoodEvaluations, 4L + diagnostics$lineSearchHalvings)
      parameters <- sum(lengths(lapply(setup$covStruct, function(term) term$par)))
      expect_equal(fitted$engineDiagnostics$profile$counters$sensitivityRHSColumns,
        (diagnostics$likelihoodEvaluations - 1L) * parameters)
    }
  }
})

test_that("exact EM moments respect their memory budget", {
  dataset <- .em_direction_data()
  expect_error(mmes(y~trait, random=~id, data=dataset,
    control=mmes_control(execution=list(denseMemoryMB=0.001),
      convergence=list(maxIters=2L), update=list(blendWeight=1),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "EM moments")
  expect_error(mmes(y~trait, random=~id, data=dataset,
    control=mmes_control(execution=list(denseMemoryMB=0.003),
      convergence=list(maxIters=2L), update=list(blendWeight=1),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "EM moments")
})

test_that("full EM skips sensitivity solves and duplicate likelihood fits", {
  dataset <- .em_direction_data()
  for(solver in c("ldlt", "cholmod")){
    fitted <- mmes(y~trait, random=~id, data=dataset,
      control=mmes_control(solve=list(method=solver),
        convergence=list(maxIters=4L), update=list(blendWeight=1),
        diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE)))
    diagnostics <- fitted$engineDiagnostics
    expect_equal(diagnostics$blend$iterationSensitivityColumns, 0)
    expect_equal(diagnostics$blend$aiEvaluations, 0L)
    expect_equal(diagnostics$blend$emEvaluations, 3L)
    expect_equal(diagnostics$blend$likelihoodEvaluations, 5L)
    expect_equal(diagnostics$coefficientEvaluations, 5)
    if(solver == "ldlt") expect_equal(diagnostics$CsymbolicAnalyses, 1)
    expect_gt(diagnostics$profile$counters$sensitivityRHSColumns, 0)
    expect_equal(nrow(diagnostics$profile$secondsByIteration), ncol(fitted$llik))
    expect_true(all(is.finite(fitted$InfMat)))
  }
})

test_that("blending interpolates proposals rather than information matrices", {
  dataset <- .em_direction_data()
  fit <- function(weight, iterations, damping=1){
    mmes(y~trait, random=~vsm(ism(id), sigma2=0.8), rcov=~vsm(ism(units), sigma2=0.6), data=dataset, control=mmes_control(convergence=list(maxIters=iterations), update=list(blendWeight=weight), proposal=list(stepWeight=damping), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  initial <- fit(0, 1L)
  em <- fit(1, 2L)
  mixed <- fit(0.5, 2L, 0.05)
  start <- unlist(initial$covParWorking)
  expected <- start + 0.05 * (0.5 * (unlist(em$covParWorking) - start) -
    0.5 * drop(solve(initial$InfMat, initial$dLu)))
  expect_equal(mixed$engineDiagnostics$blend$lineSearchHalvings, 0L)
  expect_equal(unlist(mixed$covParWorking), expected, tolerance=1e-9)
})

test_that("fixed models skip all proposals and latent RR preserves tied scales", {
  dataset <- .em_direction_data()
  fixed <- mmes(y~trait, random=~vsm(ism(id), sigma2=0.8, fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=0.6, fixedSigma2=TRUE), data=dataset, control=mmes_control(convergence=list(maxIters=3L), update=list(blendWeight=1), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(fixed$convergence)
  expect_equal(ncol(fixed$llik), 1L)
  expect_equal(fixed$engineDiagnostics$blend$emEvaluations, 0L)
  expect_equal(fixed$engineDiagnostics$blend$aiEvaluations, 0L)
  expect_equal(as.numeric(unlist(fixed$theta)), c(0.8, 0.6))
  setup <- mmes(y~trait, random=~vsm(rrm(trait, 1L), ism(id)) + replicate, data=dataset, control=mmes_control(output=list(setup=TRUE), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  scales <- which(setup$vcParams$structure %in% 1:2 & setup$vcParams$kind == "scale")
  for(solver in c("ldlt", "cholmod")){
    fitted <- mmes(y~trait, random=~vsm(rrm(trait, 1L), ism(id)) + replicate, data=dataset, vcc=data.frame(parameter=scales, group=1L), control=mmes_control(convergence=list(maxIters=3L), update=list(blendWeight=1), solve=list(method=solver), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(fitted$covParWorking[[1L]][1L], fitted$covParWorking[[2L]][1L])
    expect_equal(fitted$theta_se[scales[1L], scales[1L]],
                 fitted$theta_se[scales[2L], scales[2L]], tolerance=1e-9)
    expect_equal(fitted$theta_se[scales[1L], scales[1L]],
                 fitted$theta_se[scales[1L], scales[2L]], tolerance=1e-9)
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
    expect_true(all(is.finite(fitted$InfMat)))
  }
})

test_that("large damped EM proposals backtrack without decreasing likelihood", {
  dataset <- .em_direction_data()
  fitted <- mmes(y~trait, random=~vsm(ism(id), sigma2=0.8), rcov=~vsm(ism(units), sigma2=0.6), data=dataset, control=mmes_control(convergence=list(maxIters=3L), update=list(blendWeight=1), proposal=list(stepWeight=20), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_gt(fitted$engineDiagnostics$blend$lineSearchHalvings, 0L)
  expect_true(all(diff(c(fitted$llik)) >= -1e-9))
})

test_that("EM and blended fits reach the standalone AI optimum", {
  dataset <- .em_direction_data()
  set.seed(412)
  dataset$y <- as.integer(dataset$trait) +
    2 * rnorm(nlevels(dataset$id))[as.integer(dataset$id)] + rnorm(nrow(dataset), sd=0.6)
  fit <- function(weight){
    mmes(y~trait, random=~id, data=dataset, control=mmes_control(convergence=list(maxIters=150L, likelihoodTol=1e-8, parameterTol=1e-8), update=list(blendWeight=weight), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  ai <- fit(0)
  for(weight in c(1, 0.5)){
    fitted <- fit(weight)
    expect_true(fitted$convergence)
    expect_equal(fitted$theta, ai$theta, tolerance=1e-5)
    expect_equal(as.numeric(fitted$bu), as.numeric(ai$bu), tolerance=1e-5)
    expect_equal(tail(c(fitted$llik), 1L), tail(c(ai$llik), 1L), tolerance=1e-8)
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
  }
})

test_that("final scores, information, and uncertainty match accepted parameters", {
  dataset <- .em_direction_data()
  for(solver in c("ldlt", "cholmod")) for(shape in c("scalar", "us", "ar", "rr")){
    random <- switch(shape, scalar=~id, us=~vsm(usm(trait), ism(id)),
      ar=~vsm(ar1m(trait), ism(id)), rr=~vsm(rrm(trait, 1L), ism(id)))
    arguments <- list(fixed=y~trait, random=random, data=dataset,
      control=mmes_control(solve=list(method=solver),
        convergence=list(maxIters=4L),
        update=list(blendWeight=c(1, 0.5, 0, 0)),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    if(shape == "scalar"){
      arguments$W <- diag(seq(0.7, 1.3, length.out=nrow(dataset))) +
        matrix(0.01, nrow(dataset), nrow(dataset))
    }
    setupArguments <- arguments
    setupArguments$control <- do.call(mmes_control, utils::modifyList(
      unclass(arguments$control), list(output=list(setup=TRUE)),
      keep.null=TRUE))
    setup <- do.call(mmes, setupArguments)
    fitted <- do.call(mmes, arguments)
    reference <- .em_information_at(setup, fitted)
    for(field in c("InfMat", "theta_se", "dLu")){
      expect_equal(fitted[[field]], reference[[field]], tolerance=1e-9)
    }
    expect_equal(as.numeric(fitted$bu), as.numeric(reference$bu), tolerance=1e-9)
    expect_equal(tail(c(fitted$llik), 1L), tail(c(reference$llik), 1L), tolerance=1e-9)
  }
})

test_that("blended directions respect variance floors, fixed scales, and ratios", {
  dataset <- .em_direction_data()
  for(tied in c(FALSE, TRUE)){
    random <- if(tied) ~vsm(ism(id), sigma2=2e-6) +
      vsm(ism(replicate), sigma2=1e-6) else ~vsm(ism(id), sigma2=1e-6)
    fitted <- mmes(y~trait, random=random, data=dataset, vcc=if(tied) data.frame(parameter=1:2, group=1L, scale=c(2, 1)) else NULL, control=mmes_control(update=list(blendWeight=0.5), convergence=list(maxIters=30L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_gte(fitted$covPar[[1L]][1L], 1e-6 - 1e-14)
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
    if(tied){
      expect_equal(fitted$covPar[[1L]][1L], 2 * fitted$covPar[[2L]][1L],
                   tolerance=1e-12)
      expect_gte(fitted$covPar[[2L]][1L], 1e-6 - 1e-14)
    }
  }
  fixed <- mmes(y~trait, random=~vsm(ism(id), sigma2=5e-7, fixedSigma2=TRUE), data=dataset, control=mmes_control(update=list(blendWeight=0.5), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(fixed$covPar[[1L]][1L], 5e-7, tolerance=1e-15)
})
