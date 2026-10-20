.px_data <- function(){
  set.seed(414)
  dataset <- expand.grid(id=factor(1:7), site=factor(letters[1:4]), replicate=1:2)
  dataset$y <- rep(rnorm(7L), 8L) + rnorm(nrow(dataset), sd=0.7) +
    as.integer(dataset$site) / 5
  dataset
}

.px_setup <- function(arguments){
  setupArguments <- arguments
  setupArguments$control <- do.call(mmes_control, utils::modifyList(
    unclass(arguments$control), list(output=list(setup=TRUE)), keep.null=TRUE))
  do.call(mmes, setupArguments)
}

.px_dense_objective <- function(setup, y){
  X <- as.matrix(setup$X)
  Z <- lapply(setup$Z, as.matrix)
  design <- do.call(cbind, c(list(X), Z))
  structures <- setup$covStruct
  covariance <- lapply(structures, function(term){
    sommer:::evaluate_covstruct_cpp(term, term$par)
  })
  residual <- matrix(0, length(y), length(y))
  rows <- split(seq_along(y), setup$residualBlock)
  for(index in rows){
    local <- setup$residualIndex[index]
    residual[index, index] <- covariance[[length(covariance)]][local, local, drop=FALSE]
  }
  transform <- if(setup$useH) t(chol(as.matrix(setup$W))) else diag(length(y))
  precision <- crossprod(transform, solve(residual, transform))
  prior <- matrix(0, ncol(design), ncol(design))
  indices <- list()
  first <- ncol(X) + 1L
  for(term in seq_along(setup$Ai)){
    selected <- which(setup$Zind == term)
    width <- ncol(Z[[selected[1L]]])
    index <- seq.int(first, length.out=width * length(selected))
    prior[index, index] <- kronecker(solve(covariance[[term]]), as.matrix(setup$Ai[[term]]))
    indices[[term]] <- matrix(index, width)
    first <- max(index) + 1L
  }
  posterior <- solve(crossprod(design, precision %*% design) + prior)
  mean <- drop(posterior %*% crossprod(design, precision %*% y))
  if(setup$control$task$objective != "reml"){
    random <- seq.int(ncol(X) + 1L, ncol(design))
    conditional <- crossprod(design, precision %*% design) + prior
    posterior[,] <- 0
    posterior[random, random] <- solve(conditional[random, random, drop=FALSE])
  }
  moments <- lapply(seq_along(indices), function(term){
    index <- indices[[term]]
    U <- matrix(mean[index], nrow(index))
    A <- as.matrix(setup$Ai[[term]])
    moment <- crossprod(U, A %*% U)
    for(i in seq_len(ncol(index))) for(j in seq_len(ncol(index))){
      moment[i, j] <- moment[i, j] +
        sum(A * posterior[index[, i], index[, j], drop=FALSE])
    }
    moment
  })
  function(point, logExpansion=rep(0, length(indices))){
    candidate <- lapply(seq_along(structures), function(term){
      sommer:::evaluate_covstruct_cpp(structures[[term]], point[[term]])
    })
    expandedDesign <- design
    value <- 0
    for(term in seq_along(indices)){
      G <- candidate[[term]] / exp(2 * logExpansion[term])
      value <- value + nrow(indices[[term]]) * as.numeric(determinant(G, logarithm=TRUE)$modulus) +
        sum(diag(solve(G, moments[[term]])))
      expandedDesign[, indices[[term]]] <- design[, indices[[term]], drop=FALSE] *
        exp(logExpansion[term])
    }
    error <- drop(transform %*% (y - expandedDesign %*% mean))
    D <- transform %*% expandedDesign
    for(index in rows){
      local <- setup$residualIndex[index]
      R <- candidate[[length(candidate)]][local, local, drop=FALSE]
      moment <- tcrossprod(error[index]) + D[index, , drop=FALSE] %*% posterior %*%
        t(D[index, , drop=FALSE])
      value <- value + as.numeric(determinant(R, logarithm=TRUE)$modulus) +
        sum(diag(solve(R, moment)))
    }
    value
  }
}

test_that("scalar PX reduction matches independent expanded REML and ML objectives", {
  dataset <- .px_data()
  for(reml in c(TRUE, FALSE)) for(weighted in c(FALSE, TRUE)){
    arguments <- list(fixed=y~site, random=~vsm(ism(id), sigma2=0.8) +
      vsm(ism(replicate), sigma2=0.4), rcov=~vsm(ism(units), sigma2=0.6),
      data=dataset, control=mmes_control(
        task=list(objective=if(reml) "reml" else "ml"),
        update=list(recovery="pxem", blendWeight=1),
        convergence=list(maxIters=2L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    if(weighted) arguments$W <- Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
      diagonals=list(rep(0.1, nrow(dataset) - 1L), rep(1.5, nrow(dataset)),
        rep(0.1, nrow(dataset) - 1L)))
    setup <- .px_setup(arguments)
    objective <- .px_dense_objective(setup, dataset$y)
    initial <- lapply(setup$covStruct, function(term) term$par)
    native <- do.call(mmes, arguments)
    diagnostics <- native$engineDiagnostics$blend
    fitted <- lapply(native$covParWorking, as.numeric)
    expansion <- diagnostics$pxLogExpansion[, 1L]
    q <- objective(fitted, expansion)
    expect_equal(diagnostics$qGains[1L], (objective(initial) - q) / 2, tolerance=1e-7)
    seedArguments <- arguments
    seedArguments$control$update$recovery <- "em"
    seed <- do.call(mmes, seedArguments)
    seedPoint <- lapply(seed$covParWorking, as.numeric)
    reference <- optim(c(unlist(seedPoint), 0, 0), function(point){
      objective(as.list(point[1:3]), point[4:5])
    }, method="BFGS", control=list(reltol=1e-12, maxit=1000))
    expect_equal(q, reference$value, tolerance=1e-6)
    expect_equal(unlist(fitted), reference$par[1:3], tolerance=1e-4, ignore_attr=TRUE)
    expect_equal(diagnostics$pxQGains[1L], (objective(seedPoint) - q) / 2, tolerance=1e-7)
    expect_gt(diagnostics$pxQGains[1L], 0)
    expect_true(all(diff(c(native$llik)) >= -1e-9))
    expect_equal(diagnostics$aiEvaluations, 0L)
    expect_equal(diagnostics$pxEvaluations, 1L)
    expect_equal(diagnostics$method, "pxem")
    reducedCovariance <- expandedCovariance <- matrix(0, nrow(dataset), nrow(dataset))
    for(term in seq_along(setup$Ai)){
      G <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[term]], fitted[[term]])
      relationship <- solve(as.matrix(setup$Ai[[term]]))
      Z <- as.matrix(setup$Z[[which(setup$Zind == term)[1L]]])
      reducedCovariance <- reducedCovariance + G[1L, 1L] * Z %*% relationship %*% t(Z)
      expandedZ <- exp(expansion[term]) * Z
      expandedCovariance <- expandedCovariance + G[1L, 1L] / exp(2 * expansion[term]) *
        expandedZ %*% relationship %*% t(expandedZ)
    }
    expect_equal(expandedCovariance, reducedCovariance, tolerance=1e-12)
  }
})

test_that("PX uses loading-aware latent moments for reduced-rank covariance", {
  dataset <- .px_data()
  for(reml in c(TRUE, FALSE)) for(weighted in c(FALSE, TRUE)) for(correlated in c(FALSE, TRUE)){
  arguments <- list(fixed=y~site, random=~vsm(rrm(site, 2L), ism(id), sigma2=0.8),
    rcov=if(correlated) ~vsm(ar1m(site), ism(units), sigma2=0.6) else
      ~vsm(ism(units), sigma2=0.6), data=dataset,
    control=mmes_control(task=list(objective=if(reml) "reml" else "ml"),
      update=list(recovery="pxem", blendWeight=1),
      convergence=list(maxIters=2L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  if(weighted) arguments$W <- Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
    diagonals=list(rep(0.1, nrow(dataset) - 1L), rep(1.5, nrow(dataset)),
      rep(0.1, nrow(dataset) - 1L)))
  setup <- .px_setup(arguments)
  factor <- setup$covStruct[[1L]]$factors[[1L]]
  loading <- function(par){
    shape <- par[seq.int(factor$par_start, factor$par_end)]
    L <- matrix(0, factor$dim, factor$order)
    reference <- shape[factor$rr_row == 1L & factor$rr_col == 1L]
    L[cbind(factor$rr_row, factor$rr_col)] <- ifelse(factor$rr_diag,
      exp(shape - reference), shape * exp(-reference))
    L
  }
  design <- function(par){
    L <- loading(par)
    latent <- lapply(seq_len(ncol(L)), function(column){
      Reduce(`+`, lapply(seq_len(nrow(L)), function(row){
        L[row, column] * as.matrix(setup$Z[[row]])
      }))
    })
    do.call(cbind, c(list(as.matrix(setup$X)), latent))
  }
  initial <- lapply(setup$covStruct, function(term) term$par)
  D <- design(initial[[1L]])
  transform <- if(weighted) t(chol(as.matrix(setup$W))) else diag(nrow(dataset))
  residualRows <- split(seq_len(nrow(dataset)), setup$residualBlock)
  residualCovariance <- function(par){
    shape <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[2L]], par)
    covariance <- matrix(0, nrow(dataset), nrow(dataset))
    for(rows in residualRows){
      local <- setup$residualIndex[rows]
      covariance[rows, rows] <- shape[local, local, drop=FALSE]
    }
    covariance
  }
  precision <- crossprod(transform, solve(residualCovariance(initial[[2L]]), transform))
  nFixed <- ncol(setup$X)
  latent <- seq.int(nFixed + 1L, ncol(D))
  prior <- matrix(0, ncol(D), ncol(D))
  prior[latent, latent] <- diag(length(latent)) / exp(initial[[1L]][1L])
  posterior <- solve(crossprod(D, precision %*% D) + prior)
  mean <- drop(posterior %*% crossprod(D, precision %*% dataset$y))
  if(!reml){
    conditional <- crossprod(D, precision %*% D) + prior
    posterior[,] <- 0
    posterior[latent, latent] <- solve(conditional[latent, latent, drop=FALSE])
  }
  randomMoment <- sum(mean[latent]^2) + sum(diag(posterior)[latent])
  objective <- function(point, expansion=0){
    latentVariance <- exp(point[[1L]][1L] - 2 * expansion)
    residual <- residualCovariance(point[[2L]])
    expanded <- design(point[[1L]])
    expanded[, latent] <- expanded[, latent, drop=FALSE] * exp(expansion)
    error <- transform %*% (dataset$y - drop(expanded %*% mean))
    expanded <- transform %*% expanded
    length(latent) * log(latentVariance) + randomMoment / latentVariance +
      as.numeric(determinant(residual, logarithm=TRUE)$modulus) +
      sum(diag(solve(residual, tcrossprod(error) + expanded %*% posterior %*% t(expanded))))
  }
  fitted <- do.call(mmes, arguments)
  point <- lapply(fitted$covParWorking, as.numeric)
  diagnostics <- fitted$engineDiagnostics$blend
  expect_equal(diagnostics$qGains[1L],
    (objective(initial) - objective(point, diagnostics$pxLogExpansion[1L, 1L])) / 2,
    tolerance=1e-7, ignore_attr=TRUE)
  expected <- exp(point[[1L]][1L]) * tcrossprod(loading(point[[1L]]))
  actual <- sommer:::evaluate_covstruct_cpp(setup$covStruct[[1L]], point[[1L]])
  expect_equal(actual, expected, tolerance=1e-10, ignore_attr=TRUE)
  expect_equal(diagnostics$latentMomentPreparations, 1L)
  }
})

test_that("PX reduction preserves fixed parameters and cross-term constraints", {
  dataset <- .px_data()
  for(model in c("fixed-shape", "fixed-scale", "random-tie", "residual-tie", "rr-tie")){
    arguments <- list(fixed=y~site, random=~id + replicate,
      data=dataset, control=mmes_control(update=list(recovery="pxem",
        blendWeight=1), convergence=list(maxIters=3L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    if(model == "fixed-shape"){
      arguments$random <- ~vsm(ar1m(site, rho=0.2, fixed=TRUE), ism(id))
    }
    if(model == "fixed-scale"){
      arguments$random <- ~vsm(ar1m(site), ism(id), sigma2=0.8, fixedSigma2=TRUE)
    }
    if(model == "rr-tie") arguments$random <- ~vsm(rrm(site, 1L), ism(id)) + replicate
    setup <- .px_setup(arguments)
    scales <- which(setup$vcParams$kind == "scale")
    if(model %in% c("random-tie", "residual-tie", "rr-tie")){
      selected <- if(model == "residual-tie") scales[c(1L, length(scales))] else scales[1:2]
      arguments$vcc <- data.frame(parameter=selected, group=1L, scale=c(2, 1))
    }
    fitted <- do.call(mmes, arguments)
    diagnostics <- fitted$engineDiagnostics$blend
    expect_true(all(diff(c(fitted$llik)) >= -1e-9), info=model)
    expect_true(all(diagnostics$qGains >= -1e-9), info=model)
    if(model == "fixed-shape") expect_equal(fitted$covPar[[1L]][2L], 0.2)
    if(model == "fixed-scale") expect_equal(fitted$covPar[[1L]][1L], 0.8)
    if(model %in% c("random-tie", "residual-tie", "rr-tie")){
      reported <- unlist(fitted$covPar)
      expect_equal(reported[selected[1L]], 2 * reported[selected[2L]],
        tolerance=1e-10, ignore_attr=TRUE)
    }
    if(model != "rr-tie"){
      referenceSetup <- .px_setup(arguments)
      objective <- .px_dense_objective(referenceSetup, dataset$y)
      initial <- lapply(referenceSetup$covStruct, function(term) term$par)
      # The first update is checked separately from the final fitted point.
      firstArguments <- arguments
      firstArguments$control$convergence$maxIters <- 2L
      first <- do.call(mmes, firstArguments)
      q <- objective(lapply(first$covParWorking, as.numeric),
        first$engineDiagnostics$blend$pxLogExpansion[, 1L])
      expect_equal(first$engineDiagnostics$blend$qGains[1L],
        (objective(initial) - q) / 2, tolerance=1e-7)
    }
  }
  fixed <- mmes(y~site, random=~vsm(ism(id), sigma2=0.8, fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=0.6, fixedSigma2=TRUE), data=dataset, control=mmes_control(update=list(recovery="pxem", blendWeight=1), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(fixed$convergence)
  expect_equal(fixed$engineDiagnostics$blend$pxEvaluations, 0L)
  expect_equal(as.numeric(unlist(fixed$theta)), c(0.8, 0.6))
})

test_that("PX schedules skip endpoint work and refresh final inference", {
  dataset <- .px_data()
  for(shape in c("scalar", "us", "rr")){
    random <- switch(shape, scalar=~id, us=~vsm(usm(site), ism(id)),
      rr=~vsm(rrm(site, 1L), ism(id)))
    arguments <- list(fixed=y~site, random=random, data=dataset,
      control=mmes_control(update=list(recovery="pxem",
        blendWeight=c(1, 0.4, 0, 0)), convergence=list(maxIters=4L),
        solve=list(method="ldlt"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    setup <- .px_setup(arguments)
    fitted <- do.call(mmes, arguments)
    diagnostics <- fitted$engineDiagnostics$blend
    expect_equal(diagnostics$pxEvaluations, 2L)
    expect_equal(diagnostics$emEvaluations, 2L)
    expect_equal(diagnostics$aiEvaluations, 2L)
    expect_equal(diagnostics$pxMstepEvaluations[3:4], c(0L, 0L))
    expect_equal(diagnostics$pxLogExpansion[, 3:4, drop=FALSE], matrix(0, 1L, 2L))
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
    covariance <- setup$covStruct
    attr(covariance, "directionBlend") <- FALSE
    attr(covariance, "emInformationOnly") <- TRUE
    for(term in seq_along(covariance)){
      covariance[[term]]$par <- as.numeric(fitted$covParWorking[[term]])
    }
    reference <- sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
      setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
      1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
      0, 1, FALSE, computeCi=0L, solver=setup$control$solve$method,
      reml=setup$control$task$objective == "reml",
      responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
      preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
    for(field in c("InfMat", "theta_se", "dLu")){
      expect_equal(fitted[[field]], reference[[field]], tolerance=1e-9)
    }
    expect_equal(as.numeric(fitted$bu), as.numeric(reference$bu), tolerance=1e-9)
    expect_equal(tail(c(fitted$llik), 1L), tail(c(reference$llik), 1L), tolerance=1e-9)
  }
  for(solver in c("ldlt", "cholmod")){
    fitted <- mmes(y~site, random=~id, data=dataset, control=mmes_control(update=list(recovery="pxem", blendWeight=1), convergence=list(maxIters=3L), solve=list(method=solver), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    expect_equal(fitted$engineDiagnostics$blend$iterationSensitivityColumns, 0)
    expect_equal(fitted$engineDiagnostics$blend$finalAIRefreshes, 1L)
    if(solver == "ldlt") expect_equal(fitted$engineDiagnostics$CsymbolicAnalyses, 1)
  }
})

test_that("PX handles residual-only models and enforces its additional memory budget", {
  dataset <- .px_data()
  fit <- function(method){
    mmes(y~site, data=dataset, control=mmes_control(update=list(recovery=method, blendWeight=1), convergence=list(maxIters=3L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  ordinary <- fit("em")
  expanded <- fit("pxem")
  expect_equal(expanded$llik, ordinary$llik)
  expect_equal(expanded$covPar, ordinary$covPar)
  expect_equal(expanded$engineDiagnostics$blend$pxEvaluations, 0L)
  expect_error(mmes(y~site, random=~id, data=dataset,
    control=mmes_control(execution=list(denseMemoryMB=0.001),
      update=list(recovery="pxem", blendWeight=1),
      convergence=list(maxIters=2L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "EM moments")
  # Allow EM's stored moments and solve panels, but not PX's extra contractions.
  budget <- (5 * 2 + 4 * (4 + 7) * 32 + 20) * 8 / 1e6
  limited <- function(method){
    mmes(y~site, random=~id, data=dataset, control=mmes_control(
      execution=list(denseMemoryMB=budget, memoryPolicy="limit"),
      update=list(recovery=method, blendWeight=1),
      convergence=list(maxIters=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  expect_s3_class(limited("em"), "mmes")
  expect_error(limited("pxem"), "EM moments")
})

test_that("PX blends reduced directions rather than expanded parameters or information", {
  dataset <- .px_data()
  fit <- function(weight, iterations, damping=1){
    mmes(y~site, random=~vsm(ism(id), sigma2=0.8), rcov=~vsm(ism(units), sigma2=0.6), data=dataset, control=mmes_control(convergence=list(maxIters=iterations), update=list(recovery="pxem", blendWeight=weight), proposal=list(stepWeight=damping), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  initial <- fit(0, 1L)
  expanded <- fit(1, 2L)
  mixed <- fit(0.4, 2L, 0.01)
  start <- unlist(initial$covParWorking)
  expected <- start + 0.01 * (0.4 * (unlist(expanded$covParWorking) - start) -
    0.6 * drop(solve(initial$InfMat, initial$dLu)))
  expect_equal(mixed$engineDiagnostics$blend$lineSearchHalvings, 0L)
  expect_equal(unlist(mixed$covParWorking), expected, tolerance=1e-9)
  large <- fit(1, 3L, 20)
  expect_gt(large$engineDiagnostics$blend$lineSearchHalvings, 0L)
  expect_true(all(diff(c(large$llik)) >= -1e-9))
})

test_that("PX and PX-AI reach the standalone AI optimum from a small starting variance", {
  set.seed(193)
  dataset <- expand.grid(id=factor(1:30), replicate=1:4)
  dataset$y <- rnorm(30L, sd=1.2)[dataset$id] + rnorm(nrow(dataset), sd=0.7)
  fit <- function(method, weight){
    mmes(y~1, random=~vsm(ism(id), sigma2=0.01), rcov=~vsm(ism(units), sigma2=1), data=dataset, control=mmes_control(convergence=list(maxIters=150L, likelihoodTol=1e-8, parameterTol=1e-8), update=list(recovery=method, blendWeight=weight), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  ai <- fit("em", 0)
  for(weight in c(1, 5)){
    fitted <- fit("pxem", weight)
    expect_true(fitted$convergence)
    expect_equal(fitted$theta, ai$theta, tolerance=1e-5)
    expect_equal(as.numeric(fitted$bu), as.numeric(ai$bu), tolerance=1e-5)
    expect_equal(tail(c(fitted$llik), 1L), tail(c(ai$llik), 1L), tolerance=1e-8)
    expect_true(all(diff(c(fitted$llik)) >= -1e-9))
  }
})
