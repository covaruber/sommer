nativeFamilyCases <- function(){
  list(
    logit=list(family=binomial("logit"), response=c(0, 2, 5), trials=c(5, 5, 5)),
    probit=list(family=binomial("probit"), response=c(0, 2, 5), trials=c(5, 5, 5)),
    cloglog=list(family=binomial("cloglog"), response=c(0, 2, 5), trials=c(5, 5, 5)),
    poisson=list(family=poisson(), response=c(0, 2, 4), trials=c(1, 1, 1)),
    ziPoisson=list(family=zipoisson(), response=c(0, 2, 4), trials=c(1, 1, 1)),
    nb=list(family=nbinom2(2.3), response=c(0, 2, 4), trials=c(1, 1, 1)),
    gaussian=list(family=gaussian(), response=c(-1.2, 0.4, 2.1), trials=c(1, 1, 1)),
    gamma=list(family=Gamma("log"), response=c(0.2, 1.3, 2.7), trials=c(1, 1, 1)))
}

nativeFamilyControl <- function(){
  mmes_control(task=list(objective="ml", glmm=list(method="laplace")),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

independentFamilyDensity <- function(case, predictor, logParameter){
  familyCode <- case$family$family
  meanValue <- case$family$linkinv(predictor)
  if(familyCode == "binomial"){
    return(sum(dbinom(case$response, case$trials, meanValue, log=TRUE)))
  }
  if(familyCode == "poisson") return(sum(dpois(case$response, meanValue, log=TRUE)))
  if(familyCode == "ziPoisson"){
    zeroProbability <- plogis(logParameter)
    countMean <- exp(predictor)
    logDensity <- ifelse(case$response == 0,
      log(zeroProbability + (1 - zeroProbability) * exp(-countMean)),
      log1p(-zeroProbability) +
        dpois(case$response, countMean, log=TRUE))
    return(sum(logDensity))
  }
  if(grepl("^Negative Binomial", familyCode)){
    return(sum(dnbinom(case$response, size=exp(logParameter), mu=meanValue, log=TRUE)))
  }
  if(familyCode == "gaussian"){
    return(sum(dnorm(case$response, meanValue, sd=exp(logParameter / 2), log=TRUE)))
  }
  sum(dgamma(case$response, shape=exp(-logParameter),
    scale=meanValue * exp(logParameter), log=TRUE))
}

test_that("native families match independent densities and mixed derivatives", {
  predictor <- c(-1.2, 0.3, 1.1)
  for(case in nativeFamilyCases()){
    specification <- sommer:::.laplace_family_specs(case$family, list())$specifications
    if(specification[[1L]]$family %in% c("gaussian", "Gamma")){
      specification[[1L]]$logParameter <- log(0.7)
    }
    evaluate <- function(predictorValue, logParameter=specification[[1L]]$logParameter){
      localSpecification <- specification
      localSpecification[[1L]]$logParameter <- logParameter
      sommer:::sommer_glmm_family_eval(case$response, predictorValue,
        localSpecification, rep(1L, 3L), case$trials)
    }
    evaluation <- evaluate(predictor)
    expect_equal(evaluation$logDensity, independentFamilyDensity(case, predictor,
      specification[[1L]]$logParameter), tolerance=1e-12)
    for(recordPosition in seq_along(predictor)){
      upperPredictor <- lowerPredictor <- predictor
      upperPredictor[recordPosition] <- predictor[recordPosition] + 1e-4
      lowerPredictor[recordPosition] <- predictor[recordPosition] - 1e-4
      upper <- evaluate(upperPredictor)
      lower <- evaluate(lowerPredictor)
      expect_equal(as.numeric(evaluation$score)[recordPosition],
        (independentFamilyDensity(case, upperPredictor, specification[[1L]]$logParameter) -
         independentFamilyDensity(case, lowerPredictor, specification[[1L]]$logParameter)) / 2e-4,
        tolerance=1e-6)
      expect_equal(as.numeric(evaluation$curvature)[recordPosition],
        -(as.numeric(upper$score)[recordPosition] - as.numeric(lower$score)[recordPosition]) / 2e-4,
        tolerance=1e-6)
      expect_equal(as.numeric(evaluation$curvatureD1)[recordPosition],
        (as.numeric(upper$curvature)[recordPosition] -
          as.numeric(lower$curvature)[recordPosition]) / 2e-4, tolerance=1e-6)
      expect_equal(as.numeric(evaluation$curvatureD2)[recordPosition],
        (as.numeric(upper$curvatureD1)[recordPosition] -
          as.numeric(lower$curvatureD1)[recordPosition]) / 2e-4, tolerance=1e-6)
    }
    if(specification[[1L]]$free){
      upper <- evaluate(predictor, specification[[1L]]$logParameter + 1e-4)
      lower <- evaluate(predictor, specification[[1L]]$logParameter - 1e-4)
      expect_equal(sum(evaluation$parameterScore),
        (independentFamilyDensity(case, predictor, specification[[1L]]$logParameter + 1e-4) -
         independentFamilyDensity(case, predictor, specification[[1L]]$logParameter - 1e-4)) / 2e-4,
        tolerance=1e-6)
      expect_equal(as.numeric(evaluation$scoreParameter),
        as.numeric(upper$score - lower$score) / 2e-4, tolerance=1e-6)
      expect_equal(as.numeric(evaluation$curvatureParameter),
        as.numeric(upper$curvature - lower$curvature) / 2e-4, tolerance=1e-6)
    }
  }
})

test_that("native expanded family scores differentiate the full marginal objective", {
  for(case in nativeFamilyCases()){
    dataset <- data.frame(response=rep(case$response, 4L),
      x=rep(c(-1, 0, 1), 4L), group=factor(rep(seq_len(4L), each=3L)))
    specifications <- sommer:::.laplace_family_specs(case$family, list())$specifications
    descriptor <- with(dataset, vsm(ism(group), sigma2=0.4))$covStruct
    fixedDesign <- Matrix::Matrix(cbind(1, dataset$x), sparse=TRUE)
    randomDesign <- Matrix::sparse.model.matrix(~group-1, dataset)
    context <- sommer:::sommer_laplace_objective_create(fixedDesign, randomDesign,
      dataset$response, rep(0.1, nrow(dataset)), "", list(descriptor),
      list(Matrix::Diagonal(4L)), familySpecs=specifications,
      familyIndex=rep(1L, nrow(dataset)), trials=rep(case$trials, 4L))
    parameters <- c(0.2, 0.15, descriptor$par[descriptor$free],
      if(specifications[[1L]]$free) specifications[[1L]]$logParameter)
    evaluation <- sommer:::sommer_laplace_objective_eval(context, parameters,
      tolerance=1e-11)
    numericalGradient <- vapply(seq_along(parameters), function(parameterPosition){
      upperParameters <- lowerParameters <- parameters
      upperParameters[parameterPosition] <- parameters[parameterPosition] + 1e-4
      lowerParameters[parameterPosition] <- parameters[parameterPosition] - 1e-4
      upper <- sommer:::sommer_laplace_objective_eval(context, upperParameters,
        tolerance=1e-11, gradient=FALSE)
      lower <- sommer:::sommer_laplace_objective_eval(context, lowerParameters,
        tolerance=1e-11, gradient=FALSE)
      expect_true(upper$valid, info=upper$reason)
      expect_true(lower$valid, info=lower$reason)
      (upper$logLik - lower$logLik) / 2e-4
    }, numeric(1))
    expect_true(evaluation$valid, info=evaluation$reason)
    expect_lte(max(abs(as.numeric(evaluation$gradient) - numericalGradient)), 2e-5)
  }
})

test_that("zero-inflated Poisson fits, predicts, simulates and reports uncertainty", {
  set.seed(154)
  group <- factor(rep(seq_len(24L), each=8L))
  x <- rnorm(length(group))
  randomEffect <- rnorm(nlevels(group), sd=0.3)
  countMean <- exp(-0.15 + 0.35 * x + randomEffect[group])
  zeroInflated <- rbinom(length(group), size=1, prob=0.27) == 1
  response <- rpois(length(group), countMean)
  response[zeroInflated] <- 0L
  data <- data.frame(response=response, x=x, group=group)
  fit <- mmes(response~x,
    random=~vsm(ism(group), sigma2=0.09, fixedSigma2=TRUE), data=data,
    family=zipoisson(pi=0.15),
    control=nativeFamilyControl())

  expect_true(fit$convergence)
  expect_equal(fit$familyParameters$parameter, "zeroInflation")
  expect_gt(fit$familyParameters$estimate, 0)
  expect_lt(fit$familyParameters$estimate, 1)
  expect_equal(fitted(fit), (1 - fit$familyParameters$estimate) *
    exp(fit$linear.predictors), tolerance=1e-10)
  expect_equal(predict(fit, type="response"), fitted(fit), tolerance=1e-10)
  expect_equal(fit$family$linkinv(fit$linear.predictors), fitted(fit),
    tolerance=1e-10)
  expect_true(all(is.finite(residuals(fit, type="quantile"))))
  simulations <- simulate(fit, nsim=1000L, seed=221L)
  expectedZeroRate <- mean(fit$familyParameters$estimate +
    (1 - fit$familyParameters$estimate) * exp(-fit$fitted.values /
      (1 - fit$familyParameters$estimate)))
  expect_lt(abs(mean(simulations == 0) - expectedZeroRate), 0.06)
  expect_true(is.finite(summary(fit)$familyParameters$StdError))
  expect_error(mmes(response~x, random=~group, data=data,
    family=zipoisson(), control=mmes_control(
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "supported only by native Laplace")
})

test_that("zero-inflated Poisson approaches the Poisson special case", {
  set.seed(165)
  group <- factor(rep(seq_len(20L), each=6L))
  x <- rnorm(length(group))
  randomEffect <- rnorm(nlevels(group), sd=0.3)
  data <- data.frame(x=x, group=group)
  data$count <- rpois(nrow(data), exp(0.2 + 0.3*x + randomEffect[group]))
  random <- ~vsm(ism(group), sigma2=0.09, fixedSigma2=TRUE)
  poissonFit <- mmes(count~x, random=random, data=data, family=poisson(),
    control=nativeFamilyControl())
  zeroInflatedFit <- mmes(count~x, random=random, data=data,
    family=zipoisson(pi=1e-8),
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace",
      familyParameters=list(zeroInflation=1e-8, estimate=FALSE))),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))

  expect_true(poissonFit$convergence)
  expect_true(zeroInflatedFit$convergence)
  expect_equal(zeroInflatedFit$llik, poissonFit$llik, tolerance=1e-7)
  expect_equal(zeroInflatedFit$b, poissonFit$b, tolerance=1e-6)
  expect_equal(zeroInflatedFit$fitted.values, poissonFit$fitted.values,
    tolerance=1e-7)
  expect_equal(zeroInflatedFit$npar, poissonFit$npar)
})

test_that("grouped binomial forms preserve constants filtering and residual units", {
  set.seed(323)
  dataset <- data.frame(x=rnorm(90L), group=factor(rep(seq_len(15L), each=6L)))
  dataset$trials <- rep(c(3L, 6L, 10L), length.out=nrow(dataset))
  groupEffect <- rnorm(15L, sd=0.4)
  dataset$successes <- rbinom(nrow(dataset), dataset$trials,
    plogis(-0.3 + 0.4 * dataset$x + groupEffect[dataset$group]))
  dataset$failures <- dataset$trials - dataset$successes
  dataset$x[c(2L, 13L)] <- NA
  for(linkName in c("logit", "probit", "cloglog")){
    fit <- mmes(binm(successes, trials=trials)~x,
      random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
      data=dataset, family=binomial(linkName), control=nativeFamilyControl())
    matrixFit <- mmes(cbind(successes, failures)~x,
      random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
      data=dataset, family=binomial(linkName), control=nativeFamilyControl())
    expect_true(fit$convergence)
    expect_equal(matrixFit$llik, fit$llik, tolerance=1e-8)
    expect_equal(matrixFit$b, fit$b, tolerance=1e-7)
    expect_equal(fit$priorWeights, dataset$trials[fit$obsInfo$included])
    expect_equal(fit$y, dataset$successes[fit$obsInfo$included] / fit$priorWeights)
    expect_equal(residuals(fit, "pearson"),
      (fit$y - fitted(fit)) / sqrt(fitted(fit) * (1-fitted(fit)) / fit$priorWeights))
    expect_equal(sum(residuals(fit, "deviance")^2), fit$deviance)
    specification <- fit$familySpecifications
    evaluation <- sommer:::sommer_glmm_family_eval(round(fit$y * fit$priorWeights),
      fit$linear.predictors, specification, fit$laplaceFamilyIndex, fit$priorWeights)
    expectedLikelihood <- sum(dbinom(round(fit$y * fit$priorWeights),
      fit$priorWeights, fitted(fit), log=TRUE)) -
      sum(fit$randomMode^2) / (2 * 0.16) +
      length(fit$randomMode) * log(1 / 0.16) / 2 - fit$logDetHessian / 2
    expect_equal(fit$llik, expectedLikelihood, tolerance=1e-9)
    expect_equal(evaluation$logDensity, sum(dbinom(round(fit$y * fit$priorWeights),
      fit$priorWeights, fitted(fit), log=TRUE)), tolerance=1e-10)
  }
})

test_that("native Gaussian likelihood agrees with exact marginal GLS", {
  set.seed(114)
  dataset <- data.frame(x=rnorm(40L), group=factor(rep(seq_len(8L), each=5L)))
  dataset$response <- 0.4 + 0.5 * dataset$x + rnorm(8L, sd=sqrt(0.3))[dataset$group] +
    rnorm(40L, sd=sqrt(0.7))
  control <- mmes_control(task=list(objective="ml", glmm=list(method="laplace",
    familyParameters=list(dispersion=0.7, estimate=FALSE))),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  fit <- mmes(response~x, random=~vsm(ism(group), sigma2=0.3, fixedSigma2=TRUE),
    data=dataset, family=gaussian(), control=control)
  design <- cbind(1, dataset$x)
  randomDesign <- model.matrix(~group-1, dataset)
  marginalCovariance <- diag(0.7, 40L) + 0.3 * tcrossprod(randomDesign)
  precision <- solve(marginalCovariance)
  expectedBeta <- solve(crossprod(design, precision %*% design),
    crossprod(design, precision %*% dataset$response))
  discrepancy <- dataset$response - drop(design %*% expectedBeta)
  expectedLikelihood <- -0.5 * (40L * log(2*pi) +
    as.numeric(determinant(marginalCovariance, logarithm=TRUE)$modulus) +
    drop(crossprod(discrepancy, precision %*% discrepancy)))
  expectedBetaCovariance <- solve(crossprod(design, precision %*% design))
  expect_equal(as.numeric(fit$b), as.numeric(expectedBeta), tolerance=1e-7)
  expect_equal(fit$llik, expectedLikelihood, tolerance=1e-9)
  expect_equal(unname(vcov(fit)), expectedBetaCovariance, tolerance=1e-6)
  prediction <- predict(fit, re.form=NA, se.fit=TRUE)
  expect_equal(prediction$se.fit^2,
    diag(design %*% expectedBetaCovariance %*% t(design)),
    tolerance=1e-6)
  expect_equal(fit$dispersion, 0.7)
  expect_equal(fit$npar, 2L)

  set.seed(815)
  simulationCount <- 1000L
  randomEffects <- matrix(rnorm(8L * simulationCount, sd=sqrt(0.3)),
    nrow=8L)
  simulatedResponses <- outer(0.4 + 0.5 * dataset$x,
    rep(1, simulationCount)) +
    randomDesign %*% randomEffects +
    matrix(rnorm(40L * simulationCount, sd=sqrt(0.7)), nrow=40L)
  simulatedEstimates <- expectedBetaCovariance %*%
    crossprod(design, precision %*% simulatedResponses)
  nominalIntervals <- abs(sweep(simulatedEstimates, 1L, c(0.4, 0.5), "-")) <=
    1.96 * sqrt(diag(vcov(fit)))
  expect_true(all(abs(rowMeans(nominalIntervals) - 0.95) < 0.03))
})

test_that("native NB2 and Gamma fits estimate shape and dispersion", {
  set.seed(817)
  dataset <- data.frame(x=rnorm(300L), group=factor(rep(seq_len(30L), each=10L)))
  groupEffect <- rnorm(30L, sd=0.25)
  meanValue <- exp(0.3 + 0.35 * dataset$x + groupEffect[dataset$group])
  dataset$count <- rnbinom(nrow(dataset), size=2.5, mu=meanValue)
  dataset$positive <- rgamma(nrow(dataset), shape=5, scale=meanValue / 5)
  nbFit <- mmes(count~x, random=~vsm(ism(group), sigma2=0.0625, fixedSigma2=TRUE),
    data=dataset, family=nbinom2(1.5), control=nativeFamilyControl())
  gammaFit <- mmes(positive~x, random=~vsm(ism(group), sigma2=0.0625, fixedSigma2=TRUE),
    data=dataset, family=Gamma("log"), control=nativeFamilyControl())
  expect_true(nbFit$convergence)
  expect_true(gammaFit$convergence)
  expect_lt(abs(nbFit$theta.nb - 2.5), 1.2)
  expect_lt(abs(gammaFit$dispersion - 0.2), 0.07)
  expect_lt(max(abs(nbFit$b - c(0.3, 0.35))), 0.3)
  expect_lt(max(abs(gammaFit$b - c(0.3, 0.35))), 0.2)
  expect_equal(nbFit$npar, 3L)
  expect_true(nbFit$familyParameters$estimated)
  expect_true(gammaFit$familyParameters$estimated)
  expect_equal(nbFit$family$variance(2), 2 + 4 / nbFit$theta.nb)
  expect_equal(residuals(gammaFit, "pearson"),
    (gammaFit$y - fitted(gammaFit)) / (fitted(gammaFit) * sqrt(gammaFit$dispersion)))
  expect_equal(sum(residuals(gammaFit, "deviance")^2), gammaFit$deviance)
  expect_equal(AIC(nbFit), -2*nbFit$llik + 2*nbFit$npar)
  fixedFit <- mmes(count~x, random=~vsm(ism(group), sigma2=0.0625, fixedSigma2=TRUE),
    data=dataset, family=nbinom2(2.5),
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace",
      familyParameters=list(theta=2.5, estimate=FALSE))),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(fixedFit$theta.nb, 2.5)
  expect_equal(fixedFit$npar, 2L)
  expect_false(fixedFit$familyParameters$estimated)
  expect_gte(nbFit$llik, fixedFit$llik - 1e-7)
  nbSummaryParameters <- summary(nbFit)$familyParameters
  expect_equal(nbSummaryParameters[, c("trait", "parameter", "estimate",
    "estimated")], nbFit$familyParameters[, c("trait", "parameter", "estimate",
    "estimated")])
  expect_true(all(is.finite(nbSummaryParameters$StdError)))
  expect_match(paste(capture.output(print(nbFit)), collapse="\n"),
    "Family parameters")
})

test_that("mixed native fits estimate separate family parameters with shared covariance", {
  set.seed(207)
  dataset <- expand.grid(record=seq_len(5L), group=seq_len(20L),
    trait=c("count", "continuous"))
  dataset$group <- factor(dataset$group)
  dataset$trait <- factor(dataset$trait)
  dataset$x <- rnorm(nrow(dataset))
  shapeCovariance <- matrix(c(1, 0.25, 0.25, 1), 2L)
  latentEffects <- matrix(rnorm(40L), 20L) %*% chol(0.1 * shapeCovariance)
  traitIndex <- as.integer(dataset$trait)
  eta <- ifelse(dataset$trait == "count", 0.2, 0.5) + 0.3 * dataset$x +
    latentEffects[cbind(as.integer(dataset$group), traitIndex)]
  dataset$response <- numeric(nrow(dataset))
  countRows <- dataset$trait == "count"
  dataset$response[countRows] <- rnbinom(sum(countRows), size=2.5, mu=exp(eta[countRows]))
  dataset$response[!countRows] <- rnorm(sum(!countRows), eta[!countRows], sqrt(0.7))
  fit <- mmes(response~trait-1+x,
    random=~vsm(usm(trait, theta=shapeCovariance, fixed=rep(TRUE, 2L)),
      ism(group), sigma2=0.1, fixedSigma2=TRUE),
    data=dataset, family=familym("trait", count=nbinom2(2), continuous=gaussian()),
    control=nativeFamilyControl())
  expect_true(fit$convergence)
  expect_equal(fit$npar, 5L)
  expect_true(all(fit$familyParameters$estimated))
  expect_equal(fit$familyParameters$parameter, c("theta", "dispersion"))
  expect_lt(abs(fit$dispersion["continuous"] - 0.7), 0.3)
  expect_gt(fit$familyParameters$estimate[1L], 0.8)
  expect_lt(fit$familyParameters$estimate[1L], 6)
  expect_equal(length(fit$laplaceDiagnostics$gradient), fit$npar)
  context <- sommer:::sommer_laplace_objective_create(fit$fixedDesign,
    fit$randomDesign, fit$y, fit$offset, "", fit$covStruct,
    list(Matrix::Diagonal(20L)), familySpecs=fit$familySpecifications,
    familyIndex=fit$laplaceFamilyIndex, trials=fit$priorWeights)
  parameters <- fit$optimizer$par
  numericalGradient <- vapply(seq_along(parameters), function(parameterPosition){
    upperParameters <- lowerParameters <- parameters
    upperParameters[parameterPosition] <- parameters[parameterPosition] + 1e-4
    lowerParameters[parameterPosition] <- parameters[parameterPosition] - 1e-4
    upper <- sommer:::sommer_laplace_objective_eval(context, upperParameters,
      tolerance=1e-11, gradient=FALSE)
    lower <- sommer:::sommer_laplace_objective_eval(context, lowerParameters,
      tolerance=1e-11, gradient=FALSE)
    expect_true(upper$valid, info=upper$reason)
    expect_true(lower$valid, info=lower$reason)
    (upper$logLik - lower$logLik) / 2e-4
  }, numeric(1))
  exact <- sommer:::sommer_laplace_objective_eval(context, parameters, tolerance=1e-11)
  expect_lte(max(abs(as.numeric(exact$gradient) - numericalGradient)), 2e-5)
})

test_that("expanded binomial tails preserve density and finite derivative output", {
  for(linkName in c("logit", "probit", "cloglog")){
    for(predictor in c(-40, 40)){
      specification <- sommer:::.laplace_family_specs(binomial(linkName), list())$specifications
      for(response in c(0, 5)){
        evaluation <- sommer:::sommer_glmm_family_eval(response, predictor,
          specification, 1L, 5)
        logSuccess <- switch(linkName,
          logit=plogis(predictor, log.p=TRUE),
          probit=pnorm(predictor, log.p=TRUE),
          cloglog=log(-expm1(-exp(predictor))))
        logFailure <- switch(linkName,
          logit=plogis(predictor, lower.tail=FALSE, log.p=TRUE),
          probit=pnorm(predictor, lower.tail=FALSE, log.p=TRUE),
          cloglog=-exp(predictor))
        expect_equal(evaluation$logDensity, if(response == 0) 5*logFailure else 5*logSuccess,
          tolerance=1e-10)
        expect_true(all(is.finite(unlist(evaluation))))
      }
    }
  }
})

test_that("mixed native likelihoods index all agreed families and share latent effects", {
  cases <- nativeFamilyCases()
  dataset <- expand.grid(replicate=seq_len(3L), trait=names(cases), group=seq_len(3L))
  dataset$trait <- factor(dataset$trait, levels=names(cases))
  dataset$group <- factor(dataset$group)
  dataset$x <- rep(c(-0.8, 0, 0.8), length.out=nrow(dataset))
  dataset$response <- vapply(seq_len(nrow(dataset)), function(recordPosition)
    cases[[as.character(dataset$trait[recordPosition])]]$response[dataset$replicate[recordPosition]],
    numeric(1))
  dataset$trials <- vapply(as.character(dataset$trait), function(trait) cases[[trait]]$trials[1L], numeric(1))
  for(trait in c("logit", "probit", "cloglog")){
    dataset$response[dataset$trait == trait] <- rep(c(1, 2, 4), 3L)
  }
  dataset$x[5L] <- NA
  dataset$trials[8L] <- NA
  familyObject <- do.call(familym, c(list(by="trait", trials="trials"), lapply(cases, `[[`, "family")))
  settings <- list(nb=list(theta=2.3, estimate=FALSE),
    ziPoisson=list(zeroInflation=0.2, estimate=FALSE),
    gamma=list(dispersion=0.7, estimate=FALSE),
    gaussian=list(dispersion=0.8, estimate=FALSE))
  control <- mmes_control(task=list(objective="ml",
    glmm=list(method="laplace", familyParameters=settings)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  fit <- mmes(response~trait-1+x,
    random=~vsm(usm(trait,
      fixed=rep(TRUE, length(cases) * (length(cases) + 1L) / 2L - 1L)), ism(group),
      sigma2=0.1, fixedSigma2=TRUE),
    data=dataset, family=familyObject, control=control)
  expect_true(fit$convergence)
  expect_equal(length(fitted(fit)), nrow(dataset)-2L)
  expect_equal(fit$laplaceFamilyIndex, match(as.character(dataset$trait[fit$obsInfo$included]), names(cases)))
  expectedMeans <- sommer:::.pql_apply(fit$laplaceFamilies, fit$laplaceFamilyIndex,
    "linkinv", fit$linear.predictors)
  expect_equal(fitted(fit), expectedMeans, tolerance=1e-8)
  expect_equal(sum(residuals(fit, "deviance")^2), fit$deviance)
  newdata <- dataset[seq_len(7L), ]
  newdata$x <- 0.1
  linkPrediction <- predict(fit, newdata=newdata, re.form=NA)
  responsePrediction <- predict(fit, newdata=newdata, re.form=NA, type="response")
  expect_equal(responsePrediction, sommer:::.pql_apply(fit$laplaceFamilies,
    match(as.character(newdata$trait), names(cases)), "linkinv", linkPrediction))
  conditionalPrediction <- predict(fit, newdata=newdata)
  expect_length(conditionalPrediction, nrow(newdata))
  expect_true(all(is.finite(conditionalPrediction)))
  expect_equal(predict(unserialize(serialize(fit, NULL)), type="response"), fitted(fit))
  expect_true(is.finite(logLik(fit)))
  expect_equal(fit$familyParameters$estimate, c(0.2, 2.3, 0.8, 0.7))
  simulated <- as.matrix(simulate(fit, nsim=2, seed=329))
  populationSimulated <- as.matrix(simulate(fit, nsim=2, seed=329,
    re.form="population"))
  binomialRows <- fit$laplaceFamilyIndex %in% 1:3
  countRows <- fit$laplaceFamilyIndex %in% 4:5
  gammaRows <- fit$laplaceFamilyIndex == 8L
  expect_true(all(simulated[binomialRows, ] >= 0 &
    simulated[binomialRows, ] <= 1))
  expect_equal(simulated[binomialRows, ] * fit$priorWeights[binomialRows],
    round(simulated[binomialRows, ] * fit$priorWeights[binomialRows]))
  expect_true(all(simulated[countRows, ] >= 0 &
    simulated[countRows, ] == floor(simulated[countRows, ])))
  expect_true(all(simulated[gammaRows, ] > 0))
  expect_equal(dim(populationSimulated), c(fit$nobs, 2L))
  expect_true(all(populationSimulated[binomialRows, ] >= 0 &
    populationSimulated[binomialRows, ] <= 1))
  quantileResiduals <- residuals(fit, type="quantile")
  expect_length(quantileResiduals, fit$nobs)
  expect_true(all(is.finite(quantileResiduals)))
  unknownData <- newdata
  unknownData$trait <- as.character(unknownData$trait)
  unknownData$trait[1L] <- "unknown"
  expect_error(predict(fit, newdata=unknownData, re.form=NA, type="response"),
    "new levels|Traits without a family")
  expect_error(mmes(response~trait-1+x, random=~group, data=dataset,
    family=familyObject, naMethodX="fail", control=control), "Missing")
  expect_equal(fit$priorWeights, dataset$trials[fit$obsInfo$included])
})

test_that("expanded native family contracts reject unsupported or invalid data", {
  expect_error(binm(c(1.5, 1), trials=c(3, 3)), "finite integers")
  expect_error(familym("trait", a=binomial(), a=poisson()), "Families must be named")
  expect_error(sommer:::.laplace_family_specs(Gamma(), list()), "Unsupported")
  expect_error(sommer:::.laplace_family_specs(poisson("identity"), list()), "Unsupported")
  expect_error(sommer:::.laplace_family_specs(nbinom2(), list(theta=-1)), "positive")
  expect_error(sommer:::.laplace_family_specs(binomial(), list(dispersion=1)), "Unknown")
  expect_error(sommer:::.laplace_family_specs(nbinom2(), list(estimate=NA)), "TRUE/FALSE")
  expect_error(sommer:::.laplace_family_specs(
    familym("trait", a=poisson()), list(b=list())), "Unknown familyParameters")
  expect_error(nbinom2(link="identity"), "supports link='log'")
  expect_error(zipoisson(pi=0), "strictly between 0 and 1")
  expect_error(zipoisson(link="identity"), "supports link='log'")
  expect_error(sommer:::.laplace_family_specs(zipoisson(),
    list(zeroInflation=1)), "strictly between 0 and 1")
  expect_error(familym("trait", a=poisson(), trials=NA_character_), "column name")
  invalidData <- data.frame(response=c(1, 2, 3, 4), group=factor(c("a", "a", "b", "b")))
  expect_error(mmes(response~1, random=~group, data=invalidData, family=poisson(),
    control=mmes_control(task=list(glmm=list(familyParameters=list(theta=1))))),
    "native Laplace only")
  specification <- sommer:::.laplace_family_specs(Gamma("log"), list())$specifications
  expect_error(sommer:::sommer_glmm_family_eval(0, 0, specification, 1L, 1), "strictly positive")
  specification <- sommer:::.laplace_family_specs(binomial(), list())$specifications
  expect_error(sommer:::sommer_glmm_family_eval(1.5, 0, specification, 1L, 3),
    "integer successes")
  expect_error(sommer:::sommer_glmm_family_eval(1, 0, specification, NA_integer_, 3),
    "positive integers")
})
