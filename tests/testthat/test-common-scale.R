.common_scale_data <- function(){
  set.seed(913)
  dataset <- expand.grid(id=factor(seq_len(24L)), level=factor(seq_len(3L)),
                         replicate=seq_len(3L))
  effects <- matrix(rnorm(72L), 24L) %*%
    chol(0.4^abs(outer(seq_len(3L), seq_len(3L), "-")))
  dataset$h <- factor(rep(seq_len(12L), length.out=nrow(dataset)))
  dataset$unit <- interaction(dataset$id, dataset$replicate)
  dataset$y <- as.integer(dataset$level) / 3 +
    effects[cbind(as.integer(dataset$id), as.integer(dataset$level))] +
    rnorm(nrow(dataset), sd=0.8)
  dataset
}

.common_scale_native <- function(setup, point=NULL, shape=NULL, profile=FALSE,
                                 information="none", inverse=0L){
  covariance <- setup$covStruct
  if(!is.null(point)){
    offset <- 0L
    for(term in seq_along(covariance)){
      size <- length(covariance[[term]]$par)
      covariance[[term]]$par <- point[offset + seq_len(size)]
      offset <- offset + size
    }
  }
  attr(covariance, "directionBlend") <- FALSE
  attr(covariance, "commonScale") <- "off"
  attr(covariance, "emInformationOnly") <- TRUE
  attr(covariance, "covarianceInformation") <- information
  if(!is.null(shape)) attr(covariance, "aiTangent") <- shape
  attr(covariance, "profileCommonScale") <- profile
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    0, 1, FALSE, computeCi=inverse, solver=setup$solver, reml=setup$REML,
    responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
    preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
}

test_that("common-scale control uses update module and validates canonical requests", {
  expect_identical(mmes_control()$update$commonScale, "off")
  expect_identical(mmes_control(update=list(commonScale="AUTO"))$update$commonScale,
    "auto")
  expect_identical(mmes_control(update=list(commonScale="force"))$update$commonScale,
    "force")
  expect_error(mmes_control(update=list(commonScale="invalid")),
    "Experimental toggle")
  control <- mmes_control()
  control$update$commonScale <- "invalid"
  expect_error(mmes(y~level, random=~id, data=.common_scale_data(), control=control),
    "Experimental toggle")
})

test_that("profile likelihood, score and AI match independent same-point evaluations", {
  dataset <- .common_scale_data()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    setup <- mmes(y~level, random=~vsm(ar1m(level, rho=0.25), ism(id)), data=dataset, control=mmes_control(output=list(setup=TRUE), task=list(objective=if(reml) "reml" else "ml"), solve=list(method=solver), diagnostics=list(verbose=FALSE)))
    # The residual is the reference scale; random log-scale and rho remain shapes.
    point <- c(log(1.5), atanh(0.25), 0)
    shape <- diag(3)[, 1:2, drop=FALSE]
    common <- c(1, 0, 1)
    base <- .common_scale_native(setup, point)
    profiled <- .common_scale_native(setup, point, shape, TRUE)
    degrees <- nrow(dataset) - if(reml) ncol(setup$X) else 0L
    expect_equal(profiled$.commonScaleDF, degrees)
    expect_equal(profiled$.commonScale, profiled$.commonScaleQuadratic / degrees)
    scale <- profiled$.commonScale
    samePoint <- .common_scale_native(setup, point + common * log(scale))
    expect_equal(profiled$llik, samePoint$llik, tolerance=1e-9)
    expect_equal(profiled$dLu, samePoint$dLu, tolerance=1e-8)
    expect_equal(as.numeric(profiled$.commonScaleQuadraticVector),
      as.numeric(base$InfMat %*% common), tolerance=1e-8)
    cross <- crossprod(shape, samePoint$InfMat %*% common)
    schur <- crossprod(shape, samePoint$InfMat %*% shape) -
      tcrossprod(cross) / drop(crossprod(common, samePoint$InfMat %*% common))
    expect_equal(profiled$InfMat, schur, tolerance=1e-8)
    expect_equal(drop(crossprod(common, profiled$dLu)), 0, tolerance=1e-8)
    expect_equal(profiled$engineDiagnostics$nativeAISensitivityColumns, 2L)
    objective <- function(logScale){
      as.numeric(.common_scale_native(setup, point + common * logScale)$llik)
    }
    optimum <- optimize(objective, c(log(scale)-1, log(scale)+1), maximum=TRUE,
      tol=1e-9)
    expect_equal(optimum$maximum, log(scale), tolerance=1e-6)
    expect_equal(optimum$objective, as.numeric(profiled$llik), tolerance=1e-9)
    step <- 1e-5
    derivative <- vapply(seq_len(ncol(shape)), function(k){
      up <- .common_scale_native(setup, point + step * shape[, k], shape, TRUE)
      down <- .common_scale_native(setup, point - step * shape[, k], shape, TRUE)
      (as.numeric(up$llik) - as.numeric(down$llik)) / (2 * step)
    }, numeric(1))
    expect_equal(derivative, drop(-0.5 * crossprod(shape, profiled$dLu)), tolerance=1e-6)
  }
})

test_that("profiled fits preserve converged estimates and exact final uncertainty", {
  dataset <- .common_scale_data()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    model <- function(mode, information="full", inverse=2L, returnParam=FALSE){
      mmes(y~level, random=~vsm(ar1m(level), ism(id)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(output=list(covariance=information), update=list(commonScale=mode), diagnostics=list(profile="force"))), list(convergence=list(maxIters=80L, likelihoodTol=1e-11, parameterTol=1e-9), task=list(objective=if(reml) "reml" else "ml"), output=list(inverse=inverse, setup=returnParam), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))
    }
    ordinary <- model("off")
    profiled <- model("force")
    expect_true(ordinary$convergence)
    expect_true(profiled$convergence)
    expect_equal(tail(c(profiled$llik), 1L), tail(c(ordinary$llik), 1L), tolerance=1e-8)
    for(field in c("theta", "covPar", "bu", "theta_se", "Ci", "uPevList")){
      expect_equal(profiled[[field]], ordinary[[field]], tolerance=2e-5, info=field)
    }
    diagnostics <- profiled$engineDiagnostics$commonScale
    expect_equal(diagnostics$shapeParameters, diagnostics$independentParameters - 1L)
    expect_equal(diagnostics$profileSensitivityColumns,
      diagnostics$profileEvaluations * diagnostics$shapeParameters)
    expect_equal(profiled$engineDiagnostics$coefficientEvaluations,
      diagnostics$profileEvaluations + diagnostics$finalInformationEvaluations)
    expect_lte(diagnostics$finalProjectedScoreNorm, diagnostics$scoreTolerance)
    expect_true(all(diff(as.numeric(profiled$llik)) >= -1e-9))
    point <- unlist(profiled$covParWorking, use.names=FALSE)
    exact <- .common_scale_native(model("off", returnParam=TRUE), point,
      information="full", inverse=2L)
    for(field in c("InfMat", "theta_se", "dLu", "Ci")){
      expect_equal(profiled[[field]], exact[[field]], tolerance=1e-9, info=field)
    }
    omitted <- model("force", "none")
    for(field in c("theta", "bu", "Ci", "uPevList", "llik")){
      expect_equal(omitted[[field]], profiled[[field]], tolerance=1e-10, info=field)
    }
    expect_null(omitted$theta_se)
    expect_null(omitted$InfMat)
    expect_error(vpredict(omitted, ~V1/(V1+V3)),
      "refit with output=list\\(covariance='full'\\)")
    expect_false(any(startsWith(names(profiled), ".commonScale")))
    expect_false(".aiTangent" %in% names(profiled))
  }
})

test_that("scale-only fits profile exactly and retain original observations under compression", {
  dataset <- .common_scale_data()
  for(reml in c(TRUE, FALSE)) for(compression in c("off", "groups", "classes")){
    for(information in c("full", "none")){
      fit <- mmes(y~level, data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(output=list(covariance=information), update=list(commonScale="force"), representation=list(compression=compression))), list(convergence=list(maxIters=1L), task=list(objective=if(reml) "reml" else "ml"), output=list(inverse=2L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))
      residuals <- lm(y~level, dataset)$residuals
      degrees <- nrow(dataset) - if(reml) 3L else 0L
      expect_true(fit$convergence)
      expect_equal(as.numeric(fit$covPar[[1L]]), sum(residuals^2) / degrees,
        tolerance=1e-10)
      expect_equal(fit$engineDiagnostics$commonScale$degreesOfFreedom, degrees)
      expect_equal(fit$engineDiagnostics$commonScale$profileSensitivityColumns, 0)
      expect_equal(length(fitted(fit)), nrow(dataset))
    }
  }
})

test_that("correlated weighted residuals retain exact profile derivatives", {
  dataset <- .common_scale_data()
  for(reml in c(TRUE, FALSE)){
    setup <- mmes(y~level, random=~id, rcov=~vsm(ar1m(level), ism(unit)), data=dataset, W=diag(rep(c(0.8, 1, 1.2), length.out=nrow(dataset))), control=mmes_control(task=list(objective=if(reml) "reml" else "ml"), output=list(setup=TRUE), diagnostics=list(verbose=FALSE)))
    expect_equal(ncol(setup$X), 3L)
    point <- c(log(1.5), 0, atanh(0.2))
    shape <- diag(3)[, c(1L, 3L), drop=FALSE]
    common <- c(1, 1, 0)
    profiled <- .common_scale_native(setup, point, shape, TRUE)
    ordinary <- .common_scale_native(setup,
      point + common * log(profiled$.commonScale))
    expect_equal(profiled$.commonScaleDF, nrow(dataset) - if(reml) 3L else 0L)
    expect_equal(profiled$llik, ordinary$llik, tolerance=1e-9)
    expect_equal(profiled$dLu, ordinary$dLu, tolerance=1e-8)
    cross <- crossprod(shape, ordinary$InfMat %*% common)
    expected <- crossprod(shape, ordinary$InfMat %*% shape) -
      tcrossprod(cross) / drop(crossprod(common, ordinary$InfMat %*% common))
    expect_equal(profiled$InfMat, expected, tolerance=1e-8)
  }
})

test_that("random-residual variance-ratio ties reduce to analytic scale-only fitting", {
  dataset <- .common_scale_data()
  for(reml in c(TRUE, FALSE)) for(information in c("full", "none")){
    model <- function(mode, iterations){
      mmes(y~level, random=~id, data=dataset,
        vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)),
        control=mmes_control(output=list(covariance=information, inverse=2L),
          update=list(commonScale=mode),
          task=list(objective=if(reml) "reml" else "ml"),
          convergence=list(maxIters=iterations, likelihoodTol=1e-11,
            parameterTol=1e-9),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
    }
    profiled <- model("force", 1L)
    ordinary <- model("off", 60L)
    expect_true(profiled$convergence)
    expect_equal(profiled$engineDiagnostics$commonScale$shapeParameters, 0L)
    expect_equal(profiled$engineDiagnostics$commonScale$profileSensitivityColumns, 0)
    expect_equal(profiled$covPar[[1L]], 2 * profiled$covPar[[2L]], ignore_attr=TRUE)
    expect_equal(profiled$covPar, ordinary$covPar, tolerance=2e-5)
    expect_equal(profiled$Ci, ordinary$Ci, tolerance=2e-5)
    expect_equal(tail(c(profiled$llik), 1L), tail(c(ordinary$llik), 1L), tolerance=1e-8)
  }
})

test_that("profile final information supports prediction and covariance-dependent inference", {
  dataset <- .common_scale_data()
  model <- function(mode, inverse=2L){
    mmes(y~level, random=~id, data=dataset, control=mmes_control(convergence=list(maxIters=60L, likelihoodTol=1e-11, parameterTol=1e-9), output=list(inverse=inverse), update=list(commonScale=mode), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  ordinary <- model("off")
  profiled <- model("force")
  expect_true(profiled$convergence)
  expect_equal(predict(profiled, D="id")$pvals, predict(ordinary, D="id")$pvals,
    tolerance=2e-5)
  expect_equal(vpredict(profiled, ~V1/(V1+V2)), vpredict(ordinary, ~V1/(V1+V2)),
    tolerance=2e-5)
  for(method in c("none", "satterthwaite", "kr")){
    expect_equal(as.data.frame(wald_mmes(profiled, denDF=method)),
      as.data.frame(wald_mmes(ordinary, denDF=method)), tolerance=2e-5)
  }
  noInverse <- model("force", 0L)
  expect_equal(noInverse$Ci, model("off", 0L)$Ci)
  expect_equal(noInverse$covPar, profiled$covPar, tolerance=1e-10)
})

test_that("unsupported scale profiling requests fail explicitly or report auto fallback", {
  dataset <- .common_scale_data()
  base <- list(fixed=y~level, random=~id, data=dataset)
  settings <- list(
    list(control=list(solve=list(method="pcg"))),
    list(control=list(formulation="direct")),
    list(control=list(task=list(mode="solve"))),
    list(family=poisson()),
    list(control=list(update=list(blendWeight=0.3))),
    list(control=list(curvature=list(method="hybrid-ai"))),
    list(control=list(representation=list(equations="saddle-point"))))
  for(setting in settings){
    override <- setting$control
    if(is.null(override)) override <- list()
    setting$control <- NULL
    control <- utils::modifyList(unclass(mmes_control(
      update=list(commonScale="force"), convergence=list(maxIters=1L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
      override, keep.null=TRUE)
    expect_error({
      setting$control <- do.call(mmes_control, control)
      do.call(mmes, c(base, setting))
    }, "commonScale")
  }
  expect_error(mmes(y~level, random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE), data=dataset, control=mmes_control(output=list(setup=TRUE), update=list(commonScale="force"))), "every product-level")
  expect_error(mmes(y~level, random=~vsm(rrm(level, 1L), ism(id)), data=dataset, control=mmes_control(output=list(setup=TRUE), update=list(commonScale="force"))), "latent reduced-rank")
  dataset$alias <- as.integer(dataset$level)
  expect_error(mmes(y~level+alias, random=~id, data=dataset, control=mmes_control(output=list(setup=TRUE), update=list(commonScale="force"))), "full-column-rank")
  fallback <- mmes(y~level, random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE), data=dataset, control=mmes_control(convergence=list(maxIters=1L), update=list(commonScale="auto"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(fallback$engineDiagnostics$commonScale$active)
  expect_match(fallback$engineDiagnostics$commonScale$fallbackReason, "every product-level")
  dataset$constant <- rep(c(1, 2, 3), each=24L, length.out=nrow(dataset))
  expect_error(mmes(constant~level, data=dataset, control=mmes_control(update=list(commonScale="force"), diagnostics=list(verbose=FALSE, dateWarning=FALSE))), "positive|variance floor")
})
