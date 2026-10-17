test_that("correlation variance modes create consistent covariance descriptors", {
  levels3 <- factor(rep(letters[1:3], each=2))
  levels4 <- factor(rep(1:4, each=2))
  levels5 <- factor(rep(1:5, each=2))

  structures <- list(
    csm(levels3, variance="heterogeneous"),
    ar1m(levels4, variance="heterogeneous"),
    ar2m(levels4, variance="heterogeneous"),
    ar3m(levels5, variance="heterogeneous")
  )

  for(structure in structures){
    descriptor <- structure$covFactor
    expect_true(length(descriptor$par) == length(descriptor$free))
    expect_true(length(descriptor$par) == length(descriptor$par_names))
  }

  for(structure in structures[2:4]){
    descriptor <- structure$covFactor
    covariance <- descriptor$evaluator$fun(descriptor$par)
    expect_true(isSymmetric(covariance))
    expect_true(all(is.finite(covariance)))
    expect_gt(min(eigen(covariance, symmetric=TRUE, only.values=TRUE)$values), 0)
  }
})

test_that("native covariance reporting returns model-scale parameters", {
  DT_example <- get("DT_example", envir=asNamespace("sommer"))
  envFixed <- rep(FALSE, nlevels(DT_example$Env) - 1L)
  envFixed[1L] <- TRUE

  fitDiag <- mmes(
    Yield ~ Env,
    random=~vsm(dsm(Env, fixed=envFixed), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  diagNative <- covparams_mmes(fitDiag, 1L)
  expect_equal(nrow(diagNative), nlevels(DT_example$Env))
  expect_true(all(grepl("^variance\\[", diagNative$parameter)))
  expect_false(any(grepl("ratio", diagNative$parameter)))
  expect_equal(diagNative$estimate, diag(fitDiag$theta[[1]]), tolerance=1e-8)
  expect_equal(fitDiag$covParNative, covparams_mmes_se(fitDiag))

  reported <- as.numeric(fitDiag$covPar[[1]])
  scale <- reported[1L]
  ratios <- c(1, reported[-1L])
  localCovariance <- fitDiag$theta_se[seq_along(reported), seq_along(reported),
                                     drop=FALSE]
  free <- fitDiag$covStruct[[1]]$free
  localCovariance[!free,] <- 0
  localCovariance[,!free] <- 0
  analyticJacobian <- matrix(0, nlevels(DT_example$Env), length(reported))
  analyticJacobian[,1L] <- ratios
  if(length(reported) > 1L){
    analyticJacobian[cbind(2:length(ratios), 2:length(reported))] <- scale
  }
  analyticSE <- sqrt(pmax(diag(analyticJacobian %*% localCovariance %*%
                                 t(analyticJacobian)), 0))

  fitUs <- mmes(
    Yield ~ Env,
    random=~vsm(usm(Env), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  usNative <- covparams_mmes(fitUs, 1L)
  expected <- fitUs$theta[[1]][lower.tri(fitUs$theta[[1]], diag=TRUE)]
  expect_equal(usNative$estimate, expected, tolerance=1e-8)
  expect_false(any(grepl("chol", usNative$parameter)))
  usSE <- covparams_mmes_se(fitUs, 1L)
  expect_true(all(is.finite(usSE$StdError)))
  expect_equal(usSE$parameter, usNative$parameter)

  fitCorrelation <- mmes(
    Yield ~ Env,
    random=~vsm(csm(Env, rho=0.2, fixed=TRUE), ism(Name)),
    rcov=~units,
    data=DT_example,
    verbose=FALSE,
    nIters=3
  )
  correlationSE <- covparams_mmes_se(fitCorrelation, 1L)
  expect_equal(correlationSE$StdError[correlationSE$parameter == "rho"], 0)
  expect_true(is.na(correlationSE$Zratio[correlationSE$parameter == "rho"]))
})

test_that("ownm native reporting callback is dispatched generically", {
  group <- factor(rep(c("a", "b"), each=4))
  id <- factor(rep(seq_len(4), 2))
  y <- seq_along(group)
  covariance <- function(par){
    rho <- tanh(par[1])
    matrix(c(1, rho, rho, 1), 2, 2)
  }
  reporter <- function(scale, par, factor, absorb_scale=TRUE){
    c(variance=scale, correlation=tanh(par[1]))
  }

  fit <- mmes(
    y ~ 1,
    random=~vsm(ownm(group, fun=covariance, par=0.1,
                     par_names="rho_work", native_report=reporter), ism(id)),
    rcov=~units,
    data=data.frame(y, group, id),
    verbose=FALSE,
    nIters=2
  )
  native <- covparams_mmes(fit, 1L)
  expect_equal(native$parameter, c("variance", "correlation"))
  expect_true(all(is.finite(native$estimate)))
  nativeSE <- covparams_mmes_se(fit, 1L)
  expect_true(all(is.finite(nativeSE$StdError)))
})

test_that("all covariance-factor families provide finite native reports", {
  x3 <- factor(rep(letters[1:3], each=2))
  x4 <- factor(rep(letters[1:4], each=2))
  coords <- matrix(c(0,0, 1,0, 0,1), ncol=2, byrow=TRUE)
  W <- matrix(c(0,1,1, 1,0,1, 1,1,0), 3, 3)
  dimnames(W) <- list(letters[1:3], letters[1:3])

  structures <- list(
    ism(x3), dsm(x3), usm(x3),
    ar1m(x3), ar1m(x3, variance="heterogeneous"),
    ar2m(x4), ar2m(x4, variance="heterogeneous"),
    ar3m(x4), ar3m(x4, variance="heterogeneous"),
    csm(x3), csm(x3, variance="heterogeneous"),
    mam(x3), corgm(x3), fam(x3, 1), antem(x3), rrm(x3, 1),
    maternm(coords), toeplitzm(x3), sar(x3, W), car(x3, W),
    metricm(coords), metricm(coords, model="power", anisotropy="product"),
    maternm(coords, anisotropy="geometric")
  )

  for(structure in structures){
    factor <- structure$covFactor
    reporter <- factor$native_report$fun
    natural <- factor$par
    if(length(natural)){
      transform <- factor$report$transform
      for(k in seq_along(natural)){
        natural[k] <- switch(transform[k],
          identity=natural[k], exp=exp(natural[k]), tanh=tanh(natural[k]),
          bounded_logit={
            p <- stats::plogis(natural[k])
            factor$report$lower[k] +
              (factor$report$upper[k] - factor$report$lower[k]) * p
          })
      }
    }
    values <- reporter(scale=2, par=natural, factor=factor, absorb_scale=TRUE)
    expect_true(length(values) > 0L, info=factor$model)
    expect_true(all(is.finite(values)), info=factor$model)
    expect_true(all(nzchar(names(values))), info=factor$model)

    factor$par_start <- 2L
    factor$par_end <- length(natural) + 1L
    fitted <- structure(list(
      covStruct=stats::setNames(list(list(
        factors=list(factor), free=c(TRUE, factor$free)
      )), "term"),
      covPar=stats::setNames(list(c(2, natural)), "term"),
      theta_se=diag(length(natural) + 1L)
    ), class="mmes")
    publicValues <- covparams_mmes(fitted)
    publicSE <- covparams_mmes_se(fitted)
    expect_equal(publicValues$estimate, as.numeric(values), info=factor$model)
    expect_true(all(is.finite(publicSE$StdError)), info=factor$model)
  }
})

test_that("normalized SAR and CAR precision paths match the observation-space reference", {
  locations <- letters[1:4]
  adjacency <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), 4L, byrow=TRUE)
  dimnames(adjacency) <- list(locations, locations)
  dataset <- expand.grid(id=factor(seq_len(7L)), site=factor(locations), replicate=1:2)
  dataset$response <- sin(seq_len(nrow(dataset)) / 3) + as.integer(dataset$id) / 7
  relationship <- diag(7L) + matrix(0.08, 7L, 7L)
  dimnames(relationship) <- list(levels(dataset$id), levels(dataset$id))
  attr(relationship, "inverse") <- TRUE
  for(shape in list(sar(dataset$site, adjacency, rho=0.2, fixed=TRUE),
                    car(dataset$site, adjacency, rho=0.2, fixed=TRUE))){
    factor <- shape$covFactor
    covariance <- factor$evaluator$fun(factor$par)
    if(identical(factor$model, "sar")){
      rho <- -factor$rho_bound + 2 * factor$rho_bound *
        stats::plogis(factor$par[1L])
      system <- diag(nrow(adjacency)) - rho * adjacency
      first <- c(1, rep(0, nrow(adjacency) - 1L))
      normalization <- sum(solve(t(system), first)^2)
      expectedPrecision <- normalization * crossprod(system)
      expectedLogDet <- -2 * as.numeric(determinant(system, logarithm=TRUE)$modulus) -
        nrow(system) * log(normalization)
    }else{
      rho <- factor$rho_lower + (factor$rho_upper - factor$rho_lower) *
        stats::plogis(factor$par[1L])
      system <- diag(factor$degree) - rho * adjacency
      normalization <- solve(system)[1L, 1L]
      expectedPrecision <- normalization * system
      expectedLogDet <- -as.numeric(determinant(system, logarithm=TRUE)$modulus) -
        nrow(system) * log(normalization)
    }
    expect_equal(solve(covariance), expectedPrecision, tolerance=1e-10)
    expect_equal(as.numeric(determinant(covariance, logarithm=TRUE)$modulus),
      expectedLogDet, tolerance=1e-10)
    step <- 1e-6
    eta <- factor$par[1L]
    if(identical(factor$model, "sar")){
      derivativeRho <- 2 * factor$rho_bound * stats::plogis(eta) *
        (1 - stats::plogis(eta))
      systemDerivative <- -derivativeRho * adjacency
      first <- c(1, rep(0, nrow(adjacency) - 1L))
      normalizationVector <- solve(t(system), first)
      normalizationDerivativeVector <- solve(t(system),
        derivativeRho * t(adjacency) %*% normalizationVector)
      normalizationDerivative <- 2 * sum(
        normalizationVector * normalizationDerivativeVector)
      precisionDerivative <- normalizationDerivative * crossprod(system) +
        normalization * (crossprod(systemDerivative, system) +
          crossprod(system, systemDerivative))
      logDetDerivative <- -2 * sum(diag(solve(system, systemDerivative))) -
        nrow(system) * normalizationDerivative / normalization
    }else{
      derivativeRho <- (factor$rho_upper - factor$rho_lower) *
        stats::plogis(eta) * (1 - stats::plogis(eta))
      systemDerivative <- -derivativeRho * adjacency
      normalizationVector <- solve(system, c(1, rep(0, nrow(adjacency) - 1L)))
      normalizationDerivativeVector <- solve(system,
        derivativeRho * adjacency %*% normalizationVector)
      normalizationDerivative <- normalizationDerivativeVector[1L]
      precisionDerivative <- normalizationDerivative * system +
        normalization * systemDerivative
      logDetDerivative <- -sum(diag(solve(system, systemDerivative))) -
        nrow(system) * normalizationDerivative / normalization
    }
    covarianceAt <- function(value){
      factor$evaluator$fun(value)
    }
    numericPrecisionDerivative <-
      (solve(covarianceAt(eta + step)) - solve(covarianceAt(eta - step))) /
      (2 * step)
    numericLogDetDerivative <-
      (as.numeric(determinant(covarianceAt(eta + step), logarithm=TRUE)$modulus) -
       as.numeric(determinant(covarianceAt(eta - step), logarithm=TRUE)$modulus)) /
      (2 * step)
    expect_equal(precisionDerivative, numericPrecisionDerivative, tolerance=1e-7)
    expect_equal(logDetDerivative, numericLogDetDerivative, tolerance=1e-7)

    fitModel <- function(solver){
      mmes(response~site, random=~vsm(shape, ism(id), Gu=relationship),
        rcov=~units, data=dataset, solver=solver,
        nIters=1L, verbose=FALSE, dateWarning=FALSE)
    }
    coefficient <- fitModel("ldlt")
    supernodal <- fitModel("cholmod")
    expect_true(is.finite(coefficient$llik[1L]))
    for(component in c("llik", "theta", "bu", "InfMat", "dLu")){
      expect_equal(supernodal[[component]], coefficient[[component]], tolerance=1e-7)
    }
  }
})

test_that("heterogeneous AR1 reports rho and environment variances", {
  DT_example <- get("DT_example", envir=asNamespace("sommer"))
  fit <- mmes(
    Yield ~ Env,
    random=~vsm(ar1m(Env, variance="heterogeneous"), ism(Name)),
    rcov=~vsm(dsm(Env), ism(units)),
    data=DT_example,
    verbose=FALSE,
    nIters=5
  )
  native <- covparams_mmes(fit, 1L)
  varianceRows <- grepl("^variance\\[", native$parameter)
  expect_equal(native$parameter[1L], "rho")
  expect_equal(native$estimate[varianceRows], diag(fit$theta[[1]]), tolerance=1e-8)
  expect_false(any(grepl("ratio|sigma2", native$parameter)))
  nativeSE <- covparams_mmes_se(fit, 1L)
  expect_true(all(is.finite(nativeSE$StdError)))
})


test_that("covariance representations are shared across random and residual factors", {
  x3 <- factor(letters[1:3])
  x4 <- factor(letters[1:4])
  kinds <- vapply(list(
    ism(x3), dsm(x3), fam(x3, 1), rrm(x3, 1), csm(x3),
    ar1m(x3), ar2m(x4), antem(x4), mam(x4), usm(x3)
  ), function(structure) structure$covFactor$representation$kind,
  character(1))

  expect_identical(unname(kinds), c(
    "diagonal", "diagonal", "lowrank_diagonal", "lowrank_diagonal",
    "compound_symmetry", "banded_precision", "autoregressive_precision",
    "modified_cholesky_precision", "banded_covariance", "dense_cholesky"
  ))

  DT_example <- get("DT_example", envir=asNamespace("sommer"))
  setup <- mmes(
    Yield ~ Env,
    random=~vsm(usm(Env), ism(Name)),
    rcov=~vsm(ar1m(Env), ism(units)),
    data=DT_example,
    returnParam=TRUE,
    verbose=FALSE
  )
  expect_identical(setup$covStruct[[1]]$factors[[1]]$representation$kind,
                   "dense_cholesky")
  expect_identical(setup$covStruct[[2]]$factors[[1]]$representation$kind,
                   "banded_precision")
})

test_that("sparse factor kernels preserve covariance fits and uncertainty", {
  set.seed(303L)
  dataset <- expand.grid(id = factor(seq_len(8L)), level = factor(seq_len(8L)), rep = seq_len(2L))
  dataset$section <- interaction(dataset$id, dataset$rep)
  trueCovariance <- 0.4 ^ abs(outer(seq_len(8L), seq_len(8L), "-"))
  effects <- matrix(rnorm(8L * 8L), 8L, 8L) %*% chol(trueCovariance)
  dataset$response <- effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
    rnorm(nrow(dataset), sd = 0.5) + as.integer(dataset$level) / 8
  adjacency <- as.matrix(Matrix::bandSparse(8L, k = 1L, diagonals = list(rep(1, 7L)), symmetric = TRUE))
  dimnames(adjacency) <- rep(list(levels(dataset$level)), 2L)
  directed <- adjacency
  directed[1L, 2L] <- -0.2
  shapes <- list(
    ar2m(dataset$level, pacf = c(0.2, -0.1)),
    antem(dataset$level, order = 3L),
    sar(dataset$level, directed, rho = 0.2),
    car(dataset$level, adjacency, rho = 0.9, fixed = TRUE))
  cases <- list(
    list(solver = "ldlt", REML = TRUE, computeCi = 0L),
    list(solver = "ldlt", REML = FALSE, computeCi = 1L, missing = TRUE),
    list(solver = "cholmod", REML = TRUE, computeCi = 2L, weighted = TRUE),
    list(solver = "ldlt", REML = TRUE, computeCi = 0L, henderson = FALSE),
    list(solver = "ldlt", REML = TRUE, computeCi = 2L, residual = TRUE, missing = TRUE))
  fields <- c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se", "uPevList")
  for(shapeIndex in seq_along(shapes)){
    shape <- shapes[[shapeIndex]]
    for(configuration in cases[c(1L, 2L + (shapeIndex - 1L) %% 4L)]){
      modelData <- dataset
      if(isTRUE(configuration$missing)) modelData$response[c(3L, 17L, 55L)] <- NA_real_
      fitWith <- function(mode){
        arguments <- list(fixed = response ~ level, random = ~vsm(shape, ism(id)), rcov = ~units,
          data = modelData, REML = configuration$REML, computeCi = configuration$computeCi,
          getPEV = configuration$computeCi > 0L, nIters = 2L, verbose = FALSE, dateWarning = FALSE,
          engine = mmes_control(solver = configuration$solver,
            experimental = list(sparseFactors = mode, blockSchur = "off")))
        if(isTRUE(configuration$weighted)) arguments$W <- Matrix::Diagonal(nrow(dataset), x = seq(0.8, 1.2, length.out = nrow(dataset)))
        if(isFALSE(configuration$henderson)) arguments$henderson <- FALSE
        if(isTRUE(configuration$residual)){
          arguments$random <- ~id
          arguments$rcov <- ~vsm(shape, ism(section))
        }
        do.call(mmes, arguments)
      }
      reference <- fitWith("off")
      candidate <- fitWith("force")
      for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-6, info = shape$covFactor$model) }
      if(configuration$computeCi == 2L){
        expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance = 1e-6)
        if(!isFALSE(configuration$henderson)){
          product <- as.matrix(candidate$C %*% candidate$Ci)
          expect_lt(max(abs(product / mean(diag(product)) - diag(nrow(product)))), 1e-7)
        }
      }
    }
  }
  for(shape in shapes[c(1L, 2L, 3L)]){
    fitWith <- function(mode){
      mmes(response ~ level, random = ~vsm(shape, ism(id)), rcov = ~units, data = dataset,
           nIters = 2L, nThreads = 1L, computeCi = 0L, verbose = FALSE, dateWarning = FALSE,
           engine = mmes_control(solver = "pcg", experimental = list(sparseFactors = mode, matrixFreePcg = "force")))
    }
    reference <- fitWith("off")
    candidate <- fitWith("force")
    for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-6) }
  }
  for(shape in shapes[c(1L, 3L)]){
    fitWith <- function(mode){
      mmes(response ~ level, random = ~vsm(shape, ism(id)), rcov = ~units, data = dataset,
           nIters = 4L, nThreads = 1L, computeCi = 0L, verbose = FALSE, dateWarning = FALSE,
           engine = mmes_control(solver = "ldlt", experimental = list(sparseFactors = mode, blockSchur = "off")))
    }
    reference <- fitWith("off")
    candidate <- fitWith("force")
    for(field in fields){ expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-6) }
    expect_equal(candidate$convergence, reference$convergence)
  }
})

test_that("saved built-in factors acquire missing family capabilities", {
  levels <- factor(letters[1:4])
  adjacency <- as.matrix(Matrix::bandSparse(4L, k = 1L,
    diagonals = list(rep(1, 3L)), symmetric = TRUE))
  dimnames(adjacency) <- rep(list(base::levels(levels)), 2L)
  shapes <- list(antem(levels), fam(levels, 1L), rrm(levels, 1L),
                 sar(levels, adjacency, rho = 0.2), car(levels, adjacency, rho = 0.2))
  for(shape in shapes){
    legacy <- shape$covFactor
    legacy$representation$innovation <- NULL
    legacy$representation$loadings <- NULL
    legacy$representation$spatial <- NULL
    if(legacy$model %in% c("sar", "car")) legacy$representation <- list(kind = "dense")
    compiled <- sommer:::.compile_covfactor(legacy)
    expect_identical(compiled$representation, shape$covFactor$representation)
  }
})

test_that("complete innovation residual blocks preserve sparse derivative outputs", {
  set.seed(303L)
  dataset <- expand.grid(id = factor(seq_len(8L)), level = factor(seq_len(4L)), rep = 1:2)
  dataset$section <- interaction(dataset$id, dataset$rep)
  dataset$response <- rnorm(nrow(dataset))
  shape <- antem(dataset$level, order = 1L)
  fitWith <- function(mode){
    mmes(response ~ 1, random = ~id, rcov = ~vsm(shape, ism(section)),
         data = dataset, nIters = 2L, nThreads = 1L, computeCi = 0L,
         getPEV = FALSE, verbose = FALSE, dateWarning = FALSE,
         engine = mmes_control(solver = "ldlt",
           experimental = list(sparseFactors = mode, blockSchur = "off")))
  }
  reference <- fitWith("off")
  candidate <- fitWith("force")
  for(field in c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")){
    expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7)
  }
})

test_that("factor capabilities are label independent and spatial products have exact scores", {
  set.seed(304L)
  dataset <- expand.grid(id = factor(seq_len(8L)), site = factor(letters[1:4]),
                         environment = factor(seq_len(3L)), rep = 1:2)
  dataset$response <- rnorm(nrow(dataset))
  adjacency <- as.matrix(Matrix::bandSparse(4L, k = 1L,
    diagonals = list(rep(1, 3L)), symmetric = TRUE))
  dimnames(adjacency) <- rep(list(levels(dataset$site)), 2L)
  shapes <- list(sar(dataset$site, adjacency, rho = 0.2),
                 car(dataset$site, adjacency, rho = 0.2),
                 antem(dataset$site), fam(dataset$site, 1L), rrm(dataset$site, 1L))
  for(shape in shapes){
    for(mode in c("off", "force")){
      setup <- mmes(response ~ 1, random = ~vsm(shape, ar1m(environment, rho = 0.3), ism(id)),
        data = dataset, nIters = 2L, nThreads = 1L, computeCi = 0L,
        getPEV = FALSE, verbose = FALSE, dateWarning = FALSE, returnParam = TRUE,
        engine = mmes_control(solver = "ldlt",
          experimental = list(sparseFactors = mode, blockSchur = "off")))
      runEngine <- function(descriptors){
        sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
          setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
          setup$nIters, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv,
          descriptors, setup$blendWeight, setup$stepWeight, FALSE,
          responsePrepared = setup$responsePrepared, preparedMean = setup$preparedMean,
          preparedSd = setup$preparedSd, preparedIntercept = setup$preparedIntercept)
      }
      reference <- runEngine(setup$covStruct)
      anonymous <- setup$covStruct
      anonymous[[1L]]$factors[[1L]]$model <- NULL
      candidate <- runEngine(anonymous)
      for(field in c("llik", "monitor", "theta", "bu", "InfMat", "dLu", "theta_se")){
        expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7)
      }
      if(identical(shape$covFactor$model, "ante")){
        custom <- anonymous
        descriptor <- custom[[1L]]$factors[[1L]]
        evaluator <- local({
          dimension <- descriptor$dim
          innovation <- descriptor$representation$innovation
          function(parameters){
            transition <- diag(dimension)
            transition[cbind(innovation$rows, innovation$columns)] <-
              -parameters[seq_len(innovation$coefficient_count)]
            variances <- c(1, exp(parameters[-seq_len(innovation$coefficient_count)]))
            inverse <- solve(transition)
            inverse %*% diag(variances) %*% t(inverse)
          }
        })
        descriptor$evaluator <- list(backend = "R", fun = evaluator)
        descriptor$derivative <- list(backend = "numeric", rel_step = 1e-6)
        descriptor$ante_ncoef <- descriptor$ante_row <- descriptor$ante_col <- NULL
        sommer:::.validate_covfactor(descriptor)
        custom[[1L]]$factors[[1L]] <- descriptor
        customFit <- runEngine(custom)
        for(field in c("llik", "monitor", "theta", "bu", "InfMat", "dLu", "theta_se")){
          expect_equal(customFit[[field]], candidate[[field]], tolerance = 1e-6)
        }
        invalid <- descriptor
        invalid$representation$innovation$rows[1L] <- 0L
        expect_error(sommer:::.validate_covfactor(invalid), "Invalid modified-Cholesky")
      }
      if(shape$covFactor$model %in% c("sar", "car")){
        generic <- anonymous
        generic[[1L]]$factors[[1L]]$representation <- list(kind = "dense")
        baseline <- runEngine(generic)
        for(field in c("llik", "monitor", "theta", "bu", "InfMat", "dLu", "theta_se")){
          expect_equal(candidate[[field]], baseline[[field]], tolerance = 1e-6)
        }
        invalid <- anonymous[[1L]]$factors[[1L]]
        invalid$representation$spatial$upper <- invalid$representation$spatial$lower
        expect_error(sommer:::.validate_covfactor(invalid), "Invalid normalized spatial")
      }
    }
  }
})

test_that("large ANTE sparse derivative storage preserves fit outputs", {
  set.seed(2401L)
  dataset <- expand.grid(id = factor(seq_len(16L)), level = factor(seq_len(48L)), replicate = 1:2)
  dataset$response <- rnorm(nrow(dataset)) + as.integer(dataset$id) / 5 + as.integer(dataset$level) / 20
  shape <- antem(dataset$level, order = 2L)
  fitWith <- function(mode){
    mmes(response ~ 1, random = ~vsm(shape, ism(id)), rcov = ~units, data = dataset,
         nIters = 2L, computeCi = 0L, verbose = FALSE, dateWarning = FALSE,
         engine = mmes_control(solver = "ldlt", experimental = list(blockSchur = "off", sparseFactors = mode)))
  }
  reference <- fitWith("off")
  candidate <- fitWith("force")
  for(field in c("llik", "monitor", "theta", "covPar", "bu", "InfMat", "dLu", "theta_se")){
    expect_equal(candidate[[field]], reference[[field]], tolerance = 1e-7)
  }
})
