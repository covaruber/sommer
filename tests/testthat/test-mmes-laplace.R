laplaceExampleData <- function(){
  set.seed(912)
  group <- factor(rep(seq_len(30L), each=6L))
  x <- rnorm(length(group))
  groupEffect <- rnorm(nlevels(group), sd=0.45)
  probability <- plogis(-0.35 + 0.55 * x + groupEffect[group])
  data.frame(y=rbinom(length(group), 1L, probability), x=x, group=group)
}

laplaceControl <- function(...){
  mmes_control(
    task=list(objective="ml", glmm=list(method="laplace", maxit=2000L)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE),
    ...
  )
}

test_that("native Laplace GLMMs expose likelihood and prediction methods", {
  data <- laplaceExampleData()
  fit <- mmes(y~x, random=~group, data=data, family=binomial(),
    control=laplaceControl())

  expect_s3_class(fit, "mmes.glmm.laplace")
  expect_identical(fit$uPevMode, 0L)
  expect_null(fit$uPevList)
  expect_null(fit$fixedCovariance)
  expect_false("hessian" %in% names(fit$optimizer))
  expect_true(all(is.na(fit$betaSummary[, "Std. Error"])))
  expect_true(all(is.na(fit$covParNative$StdError)))
  expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
    control=laplaceControl(output=list(inverse=1L))), "defer PEV")
  expect_true(is.finite(fit$llik))
  expect_equal(fit$AIC, AIC(fit))
  expect_equal(fit$BIC, BIC(fit))
  expect_equal(as.numeric(coef(fit)), as.numeric(fit$b))
  expect_equal(fitted(fit), fit$fitted.values)
  expect_equal(fitted(fit, type="link"), fit$linear.predictors)
  expect_equal(residuals(fit), fit$y - fit$fitted.values)
  expect_equal(summary(fit)$logo$Method, "Laplace")
  expect_equal(summary(fit)$logo$logLik, fit$llik)
  expect_true(is.logical(fit$convergence))
  expect_identical(fit$REML, FALSE)
  expect_equal(as.numeric(logLik(fit)), fit$llik)
  expect_equal(AIC(fit), fit$AIC)
  expect_identical(fit$executionPlan$glmm$selected, "Native Laplace")
  expect_match(paste(capture.output(print(fit)), collapse="\n"),
    "Generalized Linear Mixed Model fit by\\s+Laplace")
  newdata <- data.frame(x=0, group=factor("1", levels=levels(data$group)))
  expect_length(predict(fit, newdata=newdata, re.form=NA), 1L)
  expect_length(predict(fit, newdata=newdata), 1L)
  predictionSE <- predict(fit, newdata=newdata, se.fit=TRUE)
  expect_named(predictionSE, c("fit", "se.fit", "cov.fit"))
  expect_true(is.finite(predictionSE$se.fit))
  expect_gt(predictionSE$se.fit, 0)
  predictionCovariance <- predict(fit, newdata=newdata, se.fit=TRUE,
    cov.fit=TRUE)
  expect_equal(predictionCovariance$cov.fit,
    matrix(predictionSE$se.fit^2, 1L, 1L,
      dimnames=dimnames(predictionCovariance$cov.fit)), tolerance=1e-8)
  responseSE <- predict(fit, newdata=newdata, type="response",
    re.form=NA, se.fit=TRUE)
  expect_equal(responseSE$se.fit, abs(fit$family$mu.eta(
    predict(fit, newdata=newdata, type="link", re.form=NA))) *
      predict(fit, newdata=newdata, type="link", re.form=NA, se.fit=TRUE)$se.fit)
  expect_error(predict(fit, newdata=newdata, cov.fit=TRUE),
    "requires se.fit=TRUE")
  expect_length(predict(fit, newdata=newdata, re.form="population",
    type="response"), 1L)
  populationSE <- predict(fit, newdata=newdata, re.form="population",
    type="response", se.fit=TRUE, cov.fit=TRUE)
  expect_true(is.finite(populationSE$se.fit))
  expect_equal(diag(populationSE$cov.fit), populationSE$se.fit^2,
    tolerance=1e-8)
  populationReplicates <- simulate(fit, nsim=400L, seed=821L,
    re.form="population")
  populationMeans <- predict(fit, re.form="population", type="response")
  expect_lt(max(abs(rowMeans(populationReplicates) - populationMeans)), 0.1)
  expect_error(predict(fit, newdata=newdata, re.form="population",
    type="link"), "response-scale marginal mean")
  unseenData <- data.frame(x=0,
    group=factor("new", levels=c(levels(data$group), "new")))
  expect_error(predict(fit, newdata=unseenData),
    "unseen random-effect level")
  expect_length(predict(fit, newdata=unseenData, re.form=NA), 1L)
  expect_true(is.finite(predict(fit, newdata=unseenData, re.form=NA,
    se.fit=TRUE)$se.fit))
  expect_length(predict(fit, newdata=unseenData, re.form="population",
    type="response"), 1L)
  expect_true(is.finite(predict(fit, newdata=newdata, re.form="population",
    type="response", se.fit=TRUE)$se.fit))
  expect_equal(predict(fit, type="link"), fit$linear.predictors)
  set.seed(140)
  priorRandomState <- .Random.seed
  simulations <- simulate(fit, nsim=3, seed=22)
  expect_equal(.Random.seed, priorRandomState)
  expect_equal(simulations, simulate(fit, nsim=3, seed=22))
  expect_equal(dim(simulations), c(fit$nobs, 3L))
  expect_true(all(as.matrix(simulations) %in% c(0, 1)))
  populationSimulations <- simulate(fit, nsim=3, seed=22,
    re.form="population")
  expect_equal(populationSimulations,
    simulate(fit, nsim=3, seed=22, re.form="population"))
  expect_false(isTRUE(all.equal(simulations, populationSimulations)))
  expect_error(simulate(fit, re.form="population", maxMemoryMB=0),
    "positive")
  samplingFactors <- sommer:::.laplace_random_effect_sampling_factors(fit, 256)
  set.seed(140)
  randomDraws <- replicate(4000L,
    sommer:::.laplace_draw_random_effects(samplingFactors))
  descriptorCovariance <- sommer:::evaluate_covstruct_cpp(
    fit$covStruct[[1L]], fit$covStruct[[1L]]$par)
  relationshipCovariance <- solve(as.matrix(fit$randomRelationships[[1L]]))
  expect_equal(cov(t(randomDraws)),
    kronecker(descriptorCovariance, relationshipCovariance),
    tolerance=0.05)
  set.seed(29)
  firstQuantiles <- residuals(fit, type="quantile")
  set.seed(29)
  expect_equal(firstQuantiles, residuals(fit, type="quantile"))
  expect_true(all(is.finite(firstQuantiles)))
  expect_length(fit$laplaceDiagnostics$gradient, fit$npar)
  expect_true(all(is.finite(fit$laplaceDiagnostics$gradient)))
  freeCovarianceRows <- fit$covParNative$free
  if(fit$laplaceDiagnostics$uncertaintyValid){
    expect_true(all(is.finite(fit$covParNative$StdError[freeCovarianceRows])))
  }else{
    expect_true(all(is.na(fit$covParNative$StdError[freeCovarianceRows])))
  }
  postFitCovariance <- covparams_mmes_se(fit)
  expect_equal(postFitCovariance$StdError,
    summary(fit)$varcomp$StdError)
  restoredFit <- unserialize(serialize(fit, NULL))
  expect_equal(predict(restoredFit), predict(fit))
  expect_equal(logLik(restoredFit), logLik(fit))
  expect_false(any(vapply(fit, typeof, character(1)) == "externalptr"))
})

test_that("native Laplace inference exposes guarded covariance and Wald tests", {
  data <- laplaceExampleData()
  fixedCovariance <- 0.2
  randomModel <- ~vsm(ism(group), sigma2=fixedCovariance, fixedSigma2=TRUE)
  fit <- mmes(y~x, random=randomModel, data=data, family=binomial(),
    control=laplaceControl())
  smaller <- mmes(y~1, random=randomModel, data=data, family=binomial(),
    control=laplaceControl())
  likelihoodInputs <- unserialize(serialize(fit[c("y", "offset",
    "laplaceFamilyIndex", "priorWeights", "familySpecifications",
    "randomDesign", "randomRelationships", "covStruct")], NULL))
  expect_s3_class(anova(smaller, fit), "anova")

  expect_null(fit$uPevList)
  diagonalPEV <- postPEV(fit, mode=1L)
  fullPEV <- postPEV(fit, mode=2L)
  expect_identical(diagonalPEV$uPevMode, 1L)
  expect_true(all(is.finite(diagonalPEV$uPevList[[1L]])))
  expect_equal(as.numeric(diagonalPEV$uPevList[[1L]]),
    unname(diag(fullPEV$Ci)[ncol(fit$fixedDesign) +
      seq_len(ncol(fit$randomDesign))]),
    tolerance=1e-8)
  expect_identical(fit[c("y", "offset", "laplaceFamilyIndex",
    "priorWeights", "familySpecifications",
    "randomDesign", "randomRelationships", "covStruct")], likelihoodInputs)
  expect_identical(smaller$y, fit$y)
  expect_identical(smaller$offset, fit$offset)
  expect_identical(smaller$laplaceFamilyIndex, fit$laplaceFamilyIndex)
  expect_identical(smaller$priorWeights, fit$priorWeights)
  expect_identical(smaller$familySpecifications, fit$familySpecifications)
  expect_identical(deparse(smaller$randomFormula), deparse(fit$randomFormula))
  expect_identical(smaller$randomDesign, fit$randomDesign)
  expect_identical(smaller$randomRelationships, fit$randomRelationships)
  expect_identical(smaller$covStruct, fit$covStruct)
  expect_equal(dim(fullPEV$Ci),
    c(ncol(fit$fixedDesign) + ncol(fit$randomDesign),
      ncol(fit$fixedDesign) + ncol(fit$randomDesign)))
  expect_identical(postPEV(fullPEV, mode=0L)$uPevMode, 0L)
  expect_false(fit$laplaceDiagnostics$uncertaintyComputed)
  fixedCovariance <- vcov(fit)
  expect_equal(unname(diag(fixedCovariance)),
    unname(summary(fit)$betas$Std.Error^2),
    tolerance=1e-10)
  expect_true(all(is.finite(diag(fixedCovariance))))
  expect_true(summary(fit)$laplaceDiagnostics$uncertaintyComputed)
  expect_true(all(is.finite(summary(fit)$betas$Std.Error)))
  expect_true(all(is.na(fit$covParNative$StdError)))
  waldTable <- anova(fit)
  expect_s3_class(waldTable, "anova")
  expect_true(all(c("Df", "Chisq", "Pr(>Chisq)") %in% names(waldTable)))
  expect_true(all(is.finite(waldTable$Chisq)))
  nativeWald <- wald_mmes(fit, ssType="conditional", terms="x")
  expect_s3_class(nativeWald, "wald.mmes")
  expect_equal(rownames(nativeWald), "x")
  expect_error(wald_mmes(fit, terms="x"),
    "use ssType='conditional'")
  expect_error(vcov(fit, complete=TRUE), "does not accept additional arguments")
  populationVariance <- sommer:::.laplace_population_variance(
    fit$randomDesign, fit$covStruct, fit$randomRelationships,
    fit$randomBlockWidths)
  fixedPredictor <- as.numeric(fit$fixedDesign[1L, , drop=FALSE] %*% fit$b) +
    fit$offset[1L]
  independentPopulationMean <- stats::integrate(function(randomValue)
    plogis(fixedPredictor + sqrt(populationVariance[1L]) * randomValue) *
      stats::dnorm(randomValue), lower=-8, upper=8)$value
  expect_equal(predict(fit, re.form="population", type="response")[1L],
    independentPopulationMean, tolerance=1e-8)

  comparison <- anova(smaller, fit)
  expect_s3_class(comparison, "anova")
  expect_equal(comparison[2L, "Chi Df"], 1)
  expect_true(is.finite(comparison[2L, "Pr(>Chisq)"]))
  expect_error(anova(fit, smaller), "object as the smaller model")
})

test_that("native covariance delta-method derivatives follow descriptor scales", {
  expect_equal(sommer:::.mmes_laplace_vc_derivative(2, "scale", NA, NA), 2)
  expect_equal(sommer:::.mmes_laplace_vc_derivative(2, "exp", NA, NA), 2)
  expect_equal(sommer:::.mmes_laplace_vc_derivative(0.3, "identity", NA, NA), 1)
  expect_equal(sommer:::.mmes_laplace_vc_derivative(0.4, "tanh", NA, NA),
    1 - 0.4^2)
  expect_equal(sommer:::.mmes_laplace_vc_derivative(0.25, "bounded_logit", 0, 1),
    0.25 * 0.75)
})

test_that("native likelihood includes the random-effect Hessian correction", {
  fixedDesign <- Matrix::Matrix(matrix(1, nrow=3L, ncol=1L), sparse=TRUE)
  randomDesign <- Matrix::Matrix(matrix(1, nrow=3L, ncol=1L), sparse=TRUE)
  randomPrecision <- Matrix::Diagonal(1L, x=2)
  response <- c(0, 1, 1)
  fixedEffects <- 0.2
  offset <- c(-0.1, 0, 0.15)
  evaluation <- sommer:::sommer_laplace_glmm_eval(
    fixedDesign, randomDesign, randomPrecision, response, offset,
    fixedEffects, 0, "binomial", log(2)
  )

  mode <- as.numeric(evaluation$mode)
  eta <- as.numeric(offset + fixedEffects + mode)
  probability <- plogis(eta)
  score <- sum(response - probability) - 2 * mode
  hessian <- 2 + sum(probability * (1 - probability))
  conditionalLogLik <- sum(response * eta - log1p(exp(eta)))
  expectedLaplace <- conditionalLogLik - mode^2 +
    0.5 * log(2) - 0.5 * log(hessian)

  expect_true(evaluation$valid)
  expect_equal(score, 0, tolerance=1e-7)
  expect_equal(evaluation$logDetHessian, log(hessian), tolerance=1e-7)
  expect_equal(as.numeric(evaluation$logLik), as.numeric(expectedLaplace),
    tolerance=1e-7)
})

test_that("Laplace backend supports Poisson models with formula offsets", {
  set.seed(57)
  group <- factor(rep(seq_len(24L), each=5L))
  exposure <- runif(length(group), 0.5, 2)
  x <- rnorm(length(group))
  randomEffect <- rnorm(nlevels(group), sd=0.3)
  data <- data.frame(group=group, exposure=exposure, x=x)
  data$count <- rpois(nrow(data),
    exposure * exp(0.2 + 0.25 * x + randomEffect[group]))
  fit <- mmes(count~x+offset(log(exposure)), random=~group,
    data=data, family=poisson(), control=laplaceControl())

  expect_true(is.finite(fit$llik))
  expect_true(all(fitted(fit) > 0))
  expect_equal(fitted(fit, type="link"), fit$linear.predictors)
})

test_that("native Laplace uses sommer structured random covariance descriptors", {
  data <- laplaceExampleData()
  data$x <- factor(rep(letters[1:3], length.out=nrow(data)))
  fit <- mmes(y~x, random=~vsm(usm(x), ism(group)),
    data=data, family=binomial(), control=laplaceControl())

  expect_s3_class(fit, "mmes.glmm.laplace")
  expect_true(is.finite(fit$llik))
  expect_gt(nrow(fit$covParNative), 1L)
  expect_true(all(is.finite(fit$covParNative$estimate)))
  expect_true(all(is.na(fit$betaSummary[, "Std. Error"])))
  expect_length(predict(fit, newdata=data.frame(
    x=factor("a"), group=factor("1", levels=levels(data$group))),
    re.form=NA), 1L)
  expect_error(predict(fit, re.form="population", type="response"),
    "independent random-intercept")
})

test_that("native Laplace converges with supplied relationship precision", {
  set.seed(22)
  individual <- factor(rep(seq_len(10L), each=8L))
  value <- rnorm(length(individual))
  precisionFactor <- matrix(rnorm(100), 10L, 10L)
  precision <- crossprod(precisionFactor) / 10 + diag(1, 10L)
  rownames(precision) <- colnames(precision) <- as.character(seq_len(10L))
  attr(precision, "inverse") <- TRUE
  randomEffect <- rnorm(10L, sd=0.35)
  probability <- plogis(-0.2 + 0.4 * value + randomEffect[individual])
  data <- data.frame(y=rbinom(length(individual), 1L, probability),
    value=value, individual=individual)
  fit <- mmes(y~value, random=~vsm(ism(individual), Gu=precision), data=data,
    family=binomial(), control=laplaceControl())

  expect_true(fit$convergence)
  expect_true(is.finite(fit$llik))
  expect_equal(fit$optimizer$convergence, 0L)
  expect_true(all(is.na(fit$betaSummary[, "Std. Error"])))
  samplingFactors <- sommer:::.laplace_random_effect_sampling_factors(fit, 256)
  set.seed(23)
  randomDraws <- replicate(6000L,
    sommer:::.laplace_draw_random_effects(samplingFactors))
  descriptorCovariance <- sommer:::evaluate_covstruct_cpp(
    fit$covStruct[[1L]], fit$covStruct[[1L]]$par)
  expect_equal(cov(t(randomDraws)),
    kronecker(descriptorCovariance, solve(precision)), tolerance=0.1)
})

test_that("Laplace backend rejects unsupported model and engine combinations", {
  data <- laplaceExampleData()
  expect_error(
    mmes(y~x, random=~group, data=data, family=binomial(),
      control=mmes_control(task=list(glmm=list(method="laplace")),
        diagnostics=list(verbose=FALSE))),
    "objective='ml'"
  )
  expect_error(
    mmes(y~x, random=~group, data=data, family=binomial(),
      control=laplaceControl(formulation="direct")),
    "formulation='direct'"
  )
  expect_error(
    mmes(y~x, random=~group, rcov=~vsm(ism(group)), data=data,
      family=binomial(), control=laplaceControl()),
    "structured rcov terms are not yet supported"
  )
  data$trials <- 4L
  data$successes <- rbinom(nrow(data), data$trials, 0.4)
  expect_error(mmes(y~x, random=~group, data=data,
    family=binomial("cauchit"), control=laplaceControl()),
    "Unsupported native Laplace family/link")
})

test_that("Laplace task controls validate method and optimizer budget", {
  expect_error(mmes_control(task=list(glmm=list(method="invalid"))),
    "task\\$glmm\\$method")
  expect_error(mmes_control(task=list(glmm=list(maxit=0L))),
    "task\\$glmm\\$maxit")
  expect_error(mmes_control(task=list(glmm=list(modeMaxit=0L))),
    "task\\$glmm\\$modeMaxit")
  expect_error(mmes_control(task=list(glmm=list(modeTol=0))),
    "task\\$glmm\\$modeTol")
})

nativeLaplaceFixture <- function(familyName="binomial"){
  list(
    fixedDesign=Matrix::Matrix(cbind(1, c(-1, 0, 1, 0.5)), sparse=TRUE),
    randomDesign=Matrix::Matrix(matrix(c(1, 1, 0, 0, 0, 1, 1, -1),
      nrow=4L), sparse=TRUE),
    randomPrecision=Matrix::Matrix(matrix(c(2, 0.3, 0.3, 1.5), 2L),
      sparse=TRUE),
    response=if(familyName == "binomial") c(0, 1, 1, 0) else c(0, 2, 3, 1),
    offset=c(-0.1, 0, 0.2, -0.2),
    fixedEffects=c(0.2, 0.4), modeStart=c(0, 0), family=familyName,
    logDetPrecision=log(2 * 1.5 - 0.3^2), tolerance=1e-10
  )
}

test_that("native multivariate likelihood agrees with independent dense calculations", {
  for(familyName in c("binomial", "poisson")){
    fixture <- nativeLaplaceFixture(familyName)
    evaluation <- do.call(sommer:::sommer_laplace_glmm_eval, fixture)
    fixedPredictor <- as.numeric(fixture$fixedDesign %*% fixture$fixedEffects) +
      fixture$offset
    randomDesign <- as.matrix(fixture$randomDesign)
    randomPrecision <- as.matrix(fixture$randomPrecision)
    conditionalDensity <- function(randomMode){
      predictor <- fixedPredictor + as.numeric(randomDesign %*% randomMode)
      logDensity <- if(familyName == "binomial"){
        dbinom(fixture$response, 1L, plogis(predictor), log=TRUE)
      }else dpois(fixture$response, exp(predictor), log=TRUE)
      sum(logDensity) - drop(crossprod(randomMode, randomPrecision %*% randomMode)) / 2
    }
    independentMode <- optim(c(0, 0), function(randomMode)
      -conditionalDensity(randomMode), method="BFGS",
      control=list(reltol=1e-12))
    mode <- as.numeric(evaluation$mode)
    predictor <- fixedPredictor + as.numeric(randomDesign %*% mode)
    expectedMean <- if(familyName == "binomial") plogis(predictor) else exp(predictor)
    expectedWeight <- if(familyName == "binomial"){
      expectedMean * (1 - expectedMean)
    }else expectedMean
    expectedHessian <- randomPrecision +
      crossprod(randomDesign, randomDesign * expectedWeight)
    numericalHessian <- optimHess(mode, function(randomMode)
      -conditionalDensity(randomMode))
    expectedLogDet <- as.numeric(determinant(expectedHessian, logarithm=TRUE)$modulus)
    expectedLikelihood <- conditionalDensity(mode) +
      fixture$logDetPrecision / 2 - expectedLogDet / 2

    expect_true(evaluation$valid)
    expect_equal(mode, independentMode$par, tolerance=1e-5)
    expect_equal(numericalHessian, expectedHessian, tolerance=1e-5)
    expect_equal(evaluation$logDetHessian, expectedLogDet, tolerance=1e-9)
    expect_equal(evaluation$logLik, expectedLikelihood, tolerance=1e-9)
    expect_lt(evaluation$modeScore, 1e-8)
  }
})

test_that("native mode starts and evaluation order do not change the likelihood", {
  fixture <- nativeLaplaceFixture()
  baseline <- do.call(sommer:::sommer_laplace_glmm_eval, fixture)
  for(startingMode in list(c(3, -2), c(-4, 2), as.numeric(baseline$mode))){
    fixture$modeStart <- startingMode
    evaluation <- do.call(sommer:::sommer_laplace_glmm_eval, fixture)
    expect_true(evaluation$valid)
    expect_equal(evaluation$logLik, baseline$logLik, tolerance=1e-9)
    expect_equal(as.numeric(evaluation$mode), as.numeric(baseline$mode),
      tolerance=1e-8)
  }
  fixture$modeStart <- c(0, 0)
  repeated <- do.call(sommer:::sommer_laplace_glmm_eval, fixture)
  expect_identical(repeated$logLik, baseline$logLik)
})

test_that("native validation and mode failures have explicit contracts", {
  fixture <- nativeLaplaceFixture("poisson")
  for(invalidCount in c(-1, 0.5, -0.5)){
    fixture$response[1L] <- invalidCount
    expect_error(do.call(sommer:::sommer_laplace_glmm_eval, fixture),
      "non-negative integers")
  }
  fixture <- nativeLaplaceFixture()
  fixture$maxModeIterations <- 1L
  evaluation <- do.call(sommer:::sommer_laplace_glmm_eval, fixture)
  expect_false(evaluation$valid)
  expect_identical(evaluation$reason, "mode iteration limit reached")
  expect_equal(evaluation$iterations, 1L)
  fixture <- nativeLaplaceFixture()
  fixture$randomPrecision[1L, 2L] <- 0.1
  expect_error(do.call(sommer:::sommer_laplace_glmm_eval, fixture), "symmetric")
  fixture <- nativeLaplaceFixture()
  fixture$logDetPrecision <- 0
  expect_error(do.call(sommer:::sommer_laplace_glmm_eval, fixture),
    "log determinant is inconsistent")
  fixture <- nativeLaplaceFixture()
  fixture$response[1L] <- NA_real_
  expect_error(do.call(sommer:::sommer_laplace_glmm_eval, fixture), "finite")
})

test_that("Laplace preparation preserves filtering, offsets and missing-data policies", {
  data <- laplaceExampleData()
  data$y[c(2L, 11L)] <- NA_real_
  data$x[7L] <- NA_real_
  data$group[14L] <- NA
  included <- complete.cases(data)
  fit <- mmes(y~x, random=~group, data=data, family=binomial(),
    naMethodX="omit", naMethodRandom="omit", naMethodR="fail",
    control=laplaceControl())
  completeFit <- mmes(y~x, random=~group, data=data[included, ],
    family=binomial(), control=laplaceControl())
  expect_identical(fit$obsInfo$included, included)
  expect_equal(fit$y, data$y[included])
  expect_equal(fit$offset, numeric(sum(included)))
  expect_equal(fit$llik, completeFit$llik, tolerance=1e-8)
  expect_equal(fit$b, completeFit$b, tolerance=1e-6)
  expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
    naMethodX="fail", control=laplaceControl()), "fixed-effect variables")
  expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
    naMethodRandom="fail", control=laplaceControl()), "random-effect variables")
  data$exposure <- seq(0.5, 2, length.out=nrow(data))
  fitOffset <- mmes(y~x+offset(log(exposure)), random=~group,
    data=data, family=binomial(), control=laplaceControl())
  expect_equal(fitOffset$offset, log(data$exposure[included]))
  data$x[7L] <- Inf
  expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
    control=laplaceControl()), "finite retained")
})

test_that("Laplace predictions preserve training contrasts, levels and transformed terms", {
  data <- laplaceExampleData()
  data$treatment <- factor(rep(letters[1:3], length.out=nrow(data)),
    levels=letters[1:4])
  fit <- mmes(y~treatment+poly(x, 2), random=~group, data=data,
    contrasts=list(treatment="contr.sum"), family=binomial(),
    control=laplaceControl())
  expect_equal(fit$contrasts$treatment, "contr.sum")
  expect_equal(fit$fixedXlevels$treatment, letters[1:3])
  selectedRows <- c(1L, 7L, 12L)
  expect_equal(predict(fit, newdata=data[selectedRows, ], re.form=NA),
    as.numeric(fit$fixedDesign[selectedRows, , drop=FALSE] %*% fit$b),
    tolerance=1e-8)
  newdata <- data.frame(treatment=factor("a"), x=0)
  expect_length(predict(fit, newdata=newdata, re.form=NA), 1L)
  expect_error(predict(fit, newdata=data.frame(treatment="d", x=0),
    re.form=NA), "new level")
})

test_that("Laplace residuals and prediction options have explicit semantics", {
  data <- laplaceExampleData()
  fit <- mmes(y~x, random=~group, data=data, family=binomial(),
    control=laplaceControl())
  expect_equal(residuals(fit, type="working"),
    (fit$y - fitted(fit)) / fit$family$mu.eta(fit$linear.predictors))
  expect_equal(residuals(fit, type="pearson"),
    (fit$y - fitted(fit)) / sqrt(fit$family$variance(fitted(fit))))
  expect_equal(sum(residuals(fit, type="deviance")^2), fit$deviance)
  expect_equal(predict(fit, type="response"), fitted(fit))
  expect_equal(predict(fit, re.form=NA),
    as.numeric(fit$fixedDesign %*% fit$b) + fit$offset)
  expect_equal(predict(fit, PEV=FALSE, VarU=FALSE), predict(fit))
  expect_error(predict(fit, D="group"), "aggregation and prediction uncertainty")
  expect_error(predict(fit, PEV=TRUE), "aggregation and prediction uncertainty")
  expect_true(is.finite(predict(fit, se.fit=TRUE)$se.fit[1L]))
  expect_error(predict(fit, re.form=~group), "re.form must be NULL")
  expect_error(predict(fit, re.form=NA_real_), "re.form must be NULL")
})

test_that("Laplace rejects unused Gaussian controls rather than ignoring them", {
  data <- laplaceExampleData()
  for(extraControl in list(
    list(convergence=list(maxIters=2L)),
    list(proposal=list(stepWeight=0.5)),
    list(output=list(setup=TRUE)),
    list(execution=list(memoryPolicy="limit")))){
    control <- do.call(laplaceControl, extraControl)
    expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
      control=control), "does not use these Gaussian engine controls")
  }
})

test_that("Laplace uses the public dense panel memory budget", {
  fixtureData <- laplaceExampleData()
  defaultFit <- mmes(y~x, random=~group, data=fixtureData,
    family=binomial(), control=laplaceControl())
  budgetFit <- mmes(y~x, random=~group, data=fixtureData,
    family=binomial(), control=laplaceControl(execution=list(denseMemoryMB=0.003)))
  expect_equal(budgetFit$llik, defaultFit$llik, tolerance=1e-8)
  expect_equal(budgetFit$b, defaultFit$b, tolerance=1e-7)
  expect_error(mmes(y~x, random=~group, data=fixtureData,
    family=binomial(), control=laplaceControl(execution=list(denseMemoryMB=1e-8))),
    "dense memory budget")
  expect_error(predict(budgetFit, se.fit=TRUE),
    "solve panel exceeds execution\\$denseMemoryMB")
  expect_error(predict(budgetFit, se.fit=TRUE, cov.fit=TRUE),
    "prediction covariance exceeds execution\\$denseMemoryMB")
})

test_that("Laplace curvature validity is distinct from optimizer convergence", {
  curvature <- sommer:::.mmes_laplace_curvature(diag(c(2, 4)))
  expect_true(curvature$valid)
  expect_equal(curvature$covariance, diag(c(0.5, 0.25)))
  for(invalidHessian in list(diag(c(1, -1)), diag(c(1, 0)),
    diag(c(1, 1e-12)), matrix(NA_real_, 2L, 2L))){
    expect_false(sommer:::.mmes_laplace_curvature(invalidHessian)$valid)
  }
  data <- laplaceExampleData()
  limitedControl <- mmes_control(
    task=list(objective="ml", glmm=list(method="laplace", maxit=1L)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  warningMessages <- character()
  limitedFit <- withCallingHandlers(
    mmes(y~x, random=~group, data=data, family=binomial(),
      control=limitedControl),
    warning=function(condition){
      warningMessages <<- c(warningMessages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("optimizer did not converge", warningMessages)))
  expect_false(limitedFit$convergence)
  expect_false(limitedFit$laplaceDiagnostics$outerConverged)
  expect_true(limitedFit$laplaceDiagnostics$modeConverged)
  expect_true(all(is.na(limitedFit$betaSummary[, "Std. Error"])))
})

test_that("Laplace covariance starts preserve the marginal optimum", {
  data <- laplaceExampleData()
  lowStart <- mmes(y~x, random=~vsm(ism(group), sigma2=0.2), data=data,
    family=binomial(), control=laplaceControl())
  highStart <- mmes(y~x, random=~vsm(ism(group), sigma2=0.8), data=data,
    family=binomial(), control=laplaceControl())
  expect_true(lowStart$convergence)
  expect_true(highStart$convergence)
  expect_equal(lowStart$llik, highStart$llik, tolerance=1e-7)
  expect_equal(lowStart$b, highStart$b, tolerance=1e-4)
})

test_that("Laplace reports a failed initial mode instead of returning a fit", {
  data <- laplaceExampleData()
  control <- mmes_control(task=list(objective="ml",
    glmm=list(method="laplace", modeMaxit=1L)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  expect_error(mmes(y~x, random=~group, data=data, family=binomial(),
    control=control), "starting parameters are invalid: mode iteration limit")
})

test_that("Laplace permits models with no free parameters", {
  data <- laplaceExampleData()
  fit <- mmes(y~0, random=~vsm(ism(group), sigma2=0.5, fixedSigma2=TRUE),
    data=data, family=binomial(), control=laplaceControl())
  expect_true(fit$convergence)
  expect_equal(fit$npar, 0L)
  expect_true(is.finite(fit$llik))
  expect_equal(predict(fit, re.form=NA), numeric(nrow(data)))
})

test_that("Laplace boundary diagnostics respect descriptor transforms", {
  data <- laplaceExampleData()
  reference <- with(data, vsm(ism(group), sigma2=1))$covStruct
  expect_identical(sommer:::.mmes_laplace_boundary(list(reference),
    list(reference)), character())
  boundary <- reference
  boundary$par[1L] <- log(1e-10)
  expect_length(sommer:::.mmes_laplace_boundary(list(boundary),
    list(reference)), 1L)
  largerScale <- reference
  largerScale$par[1L] <- 30
  expect_identical(sommer:::.mmes_laplace_boundary(list(largerScale),
    list(reference)), character())
})

laplacePreparedFixture <- function(familyName="binomial", structured=FALSE){
  fixtureData <- laplaceExampleData()[seq_len(24L), ]
  fixtureData$group <- droplevels(fixtureData$group)
  fixtureData$environment <- factor(rep(c("a", "b"), length.out=nrow(fixtureData)))
  if(familyName == "poisson"){
    fixtureData$y <- rep(c(0, 2, 3, 1), length.out=nrow(fixtureData))
  }
  control <- laplaceControl()
  prepared <- sommer:::.mmes_prepare_model(y~x,
    if(structured) ~vsm(usm(environment), ism(group)) else ~group,
    NULL, fixtureData, if(familyName == "binomial") binomial() else poisson(),
    "exclude", list(naMethodX="exclude", naMethodRandom="exclude",
      naMethodR="exclude", contrasts=NULL, engine=attr(control, "native")))
  descriptor <- prepared$covStruct[[1L]]
  descriptor$par[1L] <- log(0.4)
  if(structured) descriptor$par[-1L] <- 0.15
  list(prepared=prepared, data=fixtureData, descriptor=descriptor)
}

test_that("modular native objective gradients include mode and determinant dependence", {
  for(familyName in c("binomial", "poisson")){
    for(structured in c(FALSE, TRUE)){
      fixture <- laplacePreparedFixture(familyName, structured)
      prepared <- fixture$prepared
      context <- sommer:::sommer_laplace_objective_create(prepared$X,
        do.call(cbind, prepared$Z), fixture$data$y,
        numeric(nrow(fixture$data)), familyName, list(fixture$descriptor),
        prepared$Ai)
      parameters <- c(-0.2, 0.35,
        fixture$descriptor$par[fixture$descriptor$free])
      evaluation <- sommer:::sommer_laplace_objective_eval(context,
        parameters, tolerance=1e-10)
      numericalGradient <- vapply(seq_along(parameters), function(position){
        upperParameters <- lowerParameters <- parameters
        upperParameters[position] <- parameters[position] + 1e-4
        lowerParameters[position] <- parameters[position] - 1e-4
        upperEvaluation <- sommer:::sommer_laplace_objective_eval(context,
          upperParameters, tolerance=1e-10, gradient=FALSE)
        lowerEvaluation <- sommer:::sommer_laplace_objective_eval(context,
          lowerParameters, tolerance=1e-10, gradient=FALSE)
        expect_true(upperEvaluation$valid)
        expect_true(lowerEvaluation$valid)
        (upperEvaluation$logLik - lowerEvaluation$logLik) / 2e-4
      }, numeric(1))
      expect_true(evaluation$valid)
      expect_equal(as.numeric(evaluation$gradient), numericalGradient,
        tolerance=2e-5)
      repeated <- sommer:::sommer_laplace_objective_eval(context,
        parameters, tolerance=1e-10)
      expect_equal(repeated$logLik, evaluation$logLik, tolerance=1e-12)
      expect_equal(as.numeric(repeated$gradient),
        as.numeric(evaluation$gradient), tolerance=1e-10)
      expect_equal(repeated$work$relationshipFactorizations, 1L)
      expect_equal(repeated$work$symbolicAnalyses, 1L)
      expect_gt(repeated$work$descriptorReuses, 0L)
    }
  }
})

test_that("preparation-only bridge preserves the actual family without Gaussian fitting", {
  fixture <- laplacePreparedFixture()
  expect_equal(as.numeric(fixture$prepared$yvar), fixture$data$y)
  expect_s3_class(fixture$prepared$modelFrame, "data.frame")
  expect_null(fixture$prepared$llik)
  expect_null(fixture$prepared$responsePrepared)
  expect_equal(ncol(fixture$prepared$X), 2L)
})

test_that("native mode scores converge below likelihood rounding precision", {
  for(structured in c(FALSE, TRUE)){
    fixture <- laplacePreparedFixture(structured=structured)
    prepared <- fixture$prepared
    context <- sommer:::sommer_laplace_objective_create(prepared$X,
      do.call(cbind, prepared$Z), fixture$data$y, numeric(nrow(fixture$data)),
      "binomial", list(fixture$descriptor), prepared$Ai)
    evaluation <- sommer:::sommer_laplace_objective_eval(context,
      c(-0.2001, 0.35, fixture$descriptor$par[fixture$descriptor$free]),
      tolerance=1e-12)
    expect_true(evaluation$valid, info=evaluation$reason)
    expect_lte(evaluation$modeScore,
      1e-12 * (1 + max(abs(evaluation$mode))))
  }
})

test_that("native objective caches modes only for identical parameters and controls", {
  fixture <- laplacePreparedFixture()
  prepared <- fixture$prepared
  makeContext <- function(memoryMB=1000){
    sommer:::sommer_laplace_objective_create(prepared$X,
      do.call(cbind, prepared$Z), fixture$data$y, numeric(nrow(fixture$data)),
      "binomial", list(fixture$descriptor), prepared$Ai, memoryMB)
  }
  parameters <- c(-0.2, 0.35, fixture$descriptor$par[fixture$descriptor$free])
  context <- makeContext()
  initial <- sommer:::sommer_laplace_objective_eval(context, parameters,
    gradient=FALSE)
  scored <- sommer:::sommer_laplace_objective_eval(context, parameters)
  expect_equal(scored$work$modeEvaluations, initial$work$modeEvaluations)
  expect_equal(scored$work$descriptorUpdates, initial$work$descriptorUpdates)
  expect_length(scored$gradient, length(parameters))
  expect_equal(scored$logLik, initial$logLik)
  expect_error(sommer:::sommer_laplace_objective_eval(context, parameters,
    tolerance=0), "controls must be positive")
  limited <- sommer:::sommer_laplace_objective_eval(context, parameters,
    maxModeIterations=1L)
  expect_false(limited$valid)
  expect_match(limited$reason, "iteration limit")
  restored <- sommer:::sommer_laplace_objective_eval(context, parameters)
  fresh <- sommer:::sommer_laplace_objective_eval(makeContext(), parameters)
  expect_equal(restored$logLik, fresh$logLik, tolerance=1e-12)
  expect_equal(restored$gradient, fresh$gradient, tolerance=1e-12)
  expect_gt(restored$work$modeEvaluations, scored$work$modeEvaluations)
  tighter <- sommer:::sommer_laplace_objective_eval(context, parameters,
    tolerance=1e-10)
  expect_gt(tighter$work$modeEvaluations, restored$work$modeEvaluations)
  expect_lte(tighter$modeScore, 1e-9)
  expect_error(sommer:::sommer_laplace_objective_eval(makeContext(1e-8),
    parameters), "dense memory budget")
  narrow <- sommer:::sommer_laplace_objective_eval(makeContext(0.001),
    parameters)
  expect_equal(narrow$gradient, fresh$gradient, tolerance=1e-12)
})

test_that("native mode factors reanalyze changed covariance sparsity", {
  fixture <- laplacePreparedFixture(structured=TRUE)
  prepared <- fixture$prepared
  fixture$descriptor$par[-1L] <- 0
  context <- sommer:::sommer_laplace_objective_create(prepared$X,
    do.call(cbind, prepared$Z), fixture$data$y, numeric(nrow(fixture$data)),
    "binomial", list(fixture$descriptor), prepared$Ai)
  diagonalParameters <- c(-0.2, 0.35,
    fixture$descriptor$par[fixture$descriptor$free])
  diagonal <- sommer:::sommer_laplace_objective_eval(context, diagonalParameters)
  correlatedParameters <- diagonalParameters
  correlationPosition <- 2L +
    which(grepl("chol\\[", fixture$descriptor$par_names))[1L]
  correlatedParameters[correlationPosition] <- 0.2
  correlated <- sommer:::sommer_laplace_objective_eval(context, correlatedParameters)
  restored <- sommer:::sommer_laplace_objective_eval(context, diagonalParameters)
  expect_true(diagonal$valid)
  expect_true(correlated$valid)
  expect_true(restored$valid)
  expect_gt(correlated$work$symbolicAnalyses, diagonal$work$symbolicAnalyses)
  expect_gt(restored$work$symbolicAnalyses, correlated$work$symbolicAnalyses)
  expect_equal(restored$logLik, diagonal$logLik, tolerance=1e-12)
  expect_equal(restored$gradient, diagonal$gradient, tolerance=1e-12)
})

test_that("descriptor-native precision matches independent crossed and AR objectives", {
  fixture <- laplacePreparedFixture("poisson")
  fixture$data$environment <- factor(rep(letters[1:3],
    length.out=nrow(fixture$data)))
  control <- laplaceControl()
  for(randomFormula in list(
    ~vsm(ar1m(environment, rho=0.3), ism(group), sigma2=0.4) + environment,
    ~vsm(usm(leg(x, 1)), ism(group), sigma2=0.4),
    ~vsm(ar1m(environment, rho=-0.2, fixed=TRUE), ism(group),
      sigma2=0.4, fixedSigma2=TRUE))){
    prepared <- sommer:::.mmes_prepare_model(y~x, randomFormula, NULL,
      fixture$data, poisson(), "exclude", list(naMethodX="exclude",
        naMethodRandom="exclude", naMethodR="exclude", contrasts=NULL,
        engine=attr(control, "native")))
    descriptors <- prepared$covStruct[seq_along(prepared$Ai)]
    suppliedRelationships <- lapply(prepared$Ai, function(relationship){
      dimension <- nrow(relationship)
      Matrix::Matrix(diag(1.5, dimension) +
        0.1 * outer(seq_len(dimension), seq_len(dimension)) / dimension^2,
        sparse=TRUE)
    })
    independentBlocks <- lapply(seq_along(descriptors), function(termPosition){
      descriptor <- descriptors[[termPosition]]
      covariance <- sommer:::evaluate_covstruct_cpp(descriptor, descriptor$par)
      kronecker(solve(covariance), as.matrix(suppliedRelationships[[termPosition]]))
    })
    independentPrecision <- Matrix::bdiag(independentBlocks)
    randomDesign <- do.call(cbind, prepared$Z)
    context <- sommer:::sommer_laplace_objective_create(prepared$X,
      randomDesign, fixture$data$y, numeric(nrow(fixture$data)), "poisson",
      descriptors, suppliedRelationships)
    parameters <- c(-0.2, 0.35, unlist(lapply(descriptors,
      function(descriptor) descriptor$par[descriptor$free]), use.names=FALSE))
    evaluation <- sommer:::sommer_laplace_objective_eval(context,
      parameters, tolerance=1e-10)
    independent <- sommer:::sommer_laplace_glmm_eval(prepared$X, randomDesign,
      independentPrecision, fixture$data$y, numeric(nrow(fixture$data)),
      parameters[1:2], numeric(ncol(randomDesign)), "poisson",
      as.numeric(determinant(independentPrecision, logarithm=TRUE)$modulus),
      tolerance=1e-10)
    expect_true(evaluation$valid)
    expect_true(independent$valid)
    expect_equal(evaluation$logLik, independent$logLik, tolerance=1e-10)
    expect_equal(evaluation$mode, independent$mode, tolerance=1e-9)
    numericalGradient <- vapply(seq_along(parameters), function(parameterPosition){
      upperParameters <- lowerParameters <- parameters
      upperParameters[parameterPosition] <- parameters[parameterPosition] + 1e-4
      lowerParameters[parameterPosition] <- parameters[parameterPosition] - 1e-4
      upper <- sommer:::sommer_laplace_objective_eval(context, upperParameters,
        tolerance=1e-10, gradient=FALSE)
      lower <- sommer:::sommer_laplace_objective_eval(context, lowerParameters,
        tolerance=1e-10, gradient=FALSE)
      expect_true(upper$valid)
      expect_true(lower$valid)
      (upper$logLik - lower$logLik) / 2e-4
    }, numeric(1))
    expect_equal(as.numeric(evaluation$gradient), numericalGradient,
      tolerance=2e-5)
    expect_equal(evaluation$work$relationshipFactorizations, length(descriptors))
  }
})

test_that("modular binomial scores retain small tail probabilities", {
  fixture <- laplacePreparedFixture()
  prepared <- fixture$prepared
  for(tailPredictor in c(-40, 40, -700, 700)){
    observationTotal <- nrow(fixture$data)
    responseValue <- as.numeric(tailPredictor > 0)
    context <- sommer:::sommer_laplace_objective_create(
      Matrix::Matrix(matrix(1, observationTotal, 1L), sparse=TRUE),
      Matrix::Matrix(0, observationTotal, ncol(prepared$Z[[1L]]), sparse=TRUE),
      rep(responseValue, observationTotal), rep(tailPredictor, observationTotal),
      "binomial", list(fixture$descriptor), prepared$Ai)
    evaluation <- sommer:::sommer_laplace_objective_eval(context,
      c(0, fixture$descriptor$par[fixture$descriptor$free]))
    expectedScore <- observationTotal * if(tailPredictor > 0)
      plogis(-tailPredictor) else -plogis(tailPredictor)
    expect_true(evaluation$valid)
    expect_equal(as.numeric(evaluation$gradient)[1L] / expectedScore, 1,
      tolerance=1e-12)
    expect_equal(evaluation$logLik,
      -observationTotal * log1p(exp(-abs(tailPredictor))), tolerance=1e-12)
  }
})
