ordinalKernelSpec <- function(link="logit", thresholds=c(-1, 0.4, 1.8)){
  list(list(family="ordinal", link=link, thresholds=thresholds,
    logParameter=0, free=FALSE))
}

ordinalKernelEval <- function(response, predictor, specification=ordinalKernelSpec(),
                             trials=rep(1, length(response))){
  sommer:::sommer_glmm_family_eval(response, predictor, specification,
    rep(1L, length(response)), trials)
}

test_that("ordinal scaffolding rejects misleading scalar-family calculations", {
  for(link in c("logit", "probit")){
    family <- ordinal(link, c(-1, 1))
    expect_s3_class(family, "family")
    expect_equal(family$thresholds, c(-1, 1))
    expect_error(family$variance(2), "category-aware native Laplace")
    expect_error(family$dev.resids(1, 2, 1), "category-aware native Laplace")
    expect_error(glm(y~1, data=data.frame(y=1:3), family=family),
      "category-aware native Laplace")
    expect_identical(sommer:::.laplace_family_specs(family, list())$
      specifications[[1L]]$family, "ordinal")
  }
  for(thresholds in list(numeric(), c(1, 1), c(2, 1), c(0, NA), Inf)){
    expect_error(ordinal(thresholds=thresholds), "strictly increasing")
  }
  expect_error(ordinal("cloglog"), "supports link")
  data <- data.frame(y=rep(1:3, 3))
  expect_error(mmes(y~1, data=data, family=ordinal(),
    control=mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "not supported by PQL")
  expect_error(mmes(y~1, data=data, family=familym("trait", score=ordinal()),
    control=mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "not supported by PQL")
})

test_that("ordinal category probabilities, expected codes and derivatives agree", {
  thresholds <- c(-1, 0.4, 1.8)
  response <- 1:4
  for(link in c("logit", "probit")){
    specification <- ordinalKernelSpec(link, thresholds)
    cdf <- if(link == "logit") plogis else pnorm
    for(predictorValue in c(-3, -0.5, 0, 1, 3)){
      predictor <- rep(predictorValue, 4)
      probabilities <- diff(c(0, cdf(thresholds - predictorValue), 1))
      evaluation <- ordinalKernelEval(response, predictor, specification)
      expect_equal(evaluation$logDensity, sum(log(probabilities)), tolerance=1e-10)
      expect_equal(as.numeric(evaluation$mean),
        rep(sum(response * probabilities), 4), tolerance=1e-12)
      individual <- vapply(response, function(category)
        exp(ordinalKernelEval(category, predictorValue, specification)$logDensity),
        numeric(1))
      expect_lte(max(abs(individual - probabilities)), 1e-12)
      expect_lte(abs(sum(individual) - 1), 1e-12)
      expect_true(all(evaluation$curvature >= 0))
      step <- 1e-4
      upper <- ordinalKernelEval(response, predictor + step, specification)
      lower <- ordinalKernelEval(response, predictor - step, specification)
      independentLogDensity <- function(eta){
        log(diff(c(0, cdf(thresholds - eta), 1)))
      }
      expect_lte(max(abs(as.numeric(evaluation$score) -
        (independentLogDensity(predictorValue + step) -
          independentLogDensity(predictorValue - step)) / (2 * step))), 1e-6)
      expect_lte(max(abs(evaluation$curvature +
        (upper$score - lower$score) / (2 * step))), 1e-6)
      expect_lte(max(abs(evaluation$curvatureD1 -
        (upper$curvature - lower$curvature) / (2 * step))), 1e-6)
      expect_lte(max(abs(evaluation$curvatureD2 -
        (upper$curvatureD1 - lower$curvatureD1) / (2 * step))), 1e-6)
      expect_equal(as.numeric(evaluation$parameterScore), rep(0, 4))
    }
  }
})

test_that("two-category ordinal kernels reduce to binary likelihoods", {
  predictor <- c(-8, -1, 0, 2, 8)
  response <- c(1, 2, 1, 2, 2)
  for(link in c("logit", "probit")){
    ordinal <- ordinalKernelEval(response, predictor, ordinalKernelSpec(link, 0.7))
    binarySpecification <- sommer:::.laplace_family_specs(binomial(link),
      list())$specifications
    binary <- sommer:::sommer_glmm_family_eval(response - 1, predictor - 0.7,
      binarySpecification, rep(1L, 5), rep(1, 5))
    for(component in c("logDensity", "score", "curvature",
                      "curvatureD1", "curvatureD2")){
      expect_lte(max(abs(ordinal[[component]] - binary[[component]])), 1e-12)
    }
    expect_equal(ordinal$mean, 1 + binary$mean, tolerance=1e-12)
  }
})

test_that("ordinal intervals remain finite in both tails and for narrow cutpoints", {
  for(link in c("logit", "probit")){
    density <- if(link == "logit") dlogis else dnorm
    for(location in c(-12, 0, 12)){
      for(width in c(1e-5, 0.5)){
        thresholds <- c(location, location + width)
        evaluation <- ordinalKernelEval(2, 0, ordinalKernelSpec(link, thresholds))
        # Quadrature is independent of subtracting cumulative probabilities.
        reference <- integrate(density, thresholds[1], thresholds[2],
          rel.tol=1e-10, abs.tol=0)$value
        expect_lte(abs(evaluation$logDensity - log(reference)), 1e-8)
        for(component in c("score", "curvature", "curvatureD1", "curvatureD2")){
          expect_true(all(is.finite(evaluation[[component]])))
        }
        if(link == "probit"){
          moments <- vapply(1:4, function(order){
            integrate(function(x){
              hermite <- switch(order, x, x^2 - 1, x^3 - 3*x,
                x^4 - 6*x^2 + 3)
              hermite * dnorm(x)
            }, thresholds[1], thresholds[2], rel.tol=1e-10,
              abs.tol=0)$value / reference
          }, numeric(1))
          logDerivatives <- c(moments[1], moments[2] - moments[1]^2,
            moments[3] - 3*moments[2]*moments[1] + 2*moments[1]^3,
            moments[4] - 4*moments[3]*moments[1] - 3*moments[2]^2 +
              12*moments[2]*moments[1]^2 - 6*moments[1]^4)
          expect_lte(max(abs(c(as.numeric(evaluation$score),
            -as.numeric(evaluation$curvature),
            -as.numeric(evaluation$curvatureD1),
            -as.numeric(evaluation$curvatureD2)) - logDerivatives)), 1e-6)
        }
      }
    }
    for(predictor in c(-40, 40)){
      evaluation <- ordinalKernelEval(1:4, rep(predictor, 4),
        ordinalKernelSpec(link))
      expect_true(is.finite(evaluation$logDensity))
      expect_true(all(is.finite(evaluation$score)))
      expect_true(all(is.finite(evaluation$curvature)))
    }
    if(link == "probit"){
      for(location in c(-40, 40)){
        evaluation <- ordinalKernelEval(2, 0,
          ordinalKernelSpec(link, c(location, location + 1e-5)))
        expect_equal(as.numeric(evaluation$score), location + 5e-6,
          tolerance=1e-8)
        expect_equal(as.numeric(evaluation$curvature), 1, tolerance=1e-9)
        expect_lt(abs(evaluation$curvatureD1), 1e-9)
        expect_lt(abs(evaluation$curvatureD2), 1e-9)
      }
    }
  }
})

test_that("ordinal rows retain their cutpoints in a mixed native family evaluation", {
  ordinalSpecification <- ordinalKernelSpec("probit")
  poissonSpecification <- sommer:::.laplace_family_specs(poisson(),
    list())$specifications
  response <- c(1, 0, 3, 5)
  predictor <- c(-0.2, 0.3, 0.4, 0.5)
  evaluation <- sommer:::sommer_glmm_family_eval(response, predictor,
    c(ordinalSpecification, poissonSpecification), c(1L, 2L, 1L, 2L),
    rep(1, 4))
  ordinal <- ordinalKernelEval(response[c(1, 3)], predictor[c(1, 3)],
    ordinalSpecification)
  poisson <- sommer:::sommer_glmm_family_eval(response[c(2, 4)],
    predictor[c(2, 4)], poissonSpecification, rep(1L, 2), rep(1, 2))
  expect_equal(evaluation$logDensity, ordinal$logDensity + poisson$logDensity)
  for(component in c("mean", "score", "curvature", "curvatureD1", "curvatureD2")){
    expect_equal(as.numeric(evaluation[[component]])[c(1, 3)],
      as.numeric(ordinal[[component]]))
    expect_equal(as.numeric(evaluation[[component]])[c(2, 4)],
      as.numeric(poisson[[component]]))
  }
})

test_that("native ordinal validation rejects invalid cutpoints and responses", {
  for(thresholds in list(numeric(), c(1, 1), c(1, -1), NA_real_, Inf)){
    expect_error(ordinalKernelEval(1, 0, ordinalKernelSpec(thresholds=thresholds)),
      "finite, strictly increasing thresholds")
  }
  for(response in c(0, -1, 1.5, 5)){
    expect_error(ordinalKernelEval(response, 0), "integer category codes")
  }
  expect_error(ordinalKernelEval(NA_real_, 0), "must be finite")
  expect_error(ordinalKernelEval(1, 0, trials=2), "only supported for binomial")
  expect_error(ordinalKernelEval(1, 0, ordinalKernelSpec("cloglog")),
    "Unsupported native Laplace")
})

ordinalEstimatedContext <- function(data, link="logit", thresholds=c(-1, 0.4, 1.8),
                                    fixed=cbind(data$x), fixedVariance=FALSE,
                                    memoryMB=1000, variance=0.16){
  specification <- ordinalKernelSpec(link, thresholds)
  specification[[1]]$free <- TRUE
  descriptor <- with(data, vsm(ism(group), sigma2=variance,
    fixedSigma2=fixedVariance))$covStruct
  context <- sommer:::sommer_laplace_objective_create(
    Matrix::Matrix(fixed, sparse=TRUE),
    Matrix::sparse.model.matrix(~group-1, data), data$y, rep(0, nrow(data)),
    "", list(descriptor), list(Matrix::Diagonal(nlevels(data$group))),
    memoryMB=memoryMB, familySpecs=specification,
    familyIndex=rep(1L, nrow(data)), trials=rep(1, nrow(data)))
  list(context=context, covariance=descriptor$par[descriptor$free],
    cutpoints=c(thresholds[1], log(diff(thresholds))))
}

test_that("estimated ordinal cutpoints enforce identification and observed categories", {
  data <- data.frame(y=rep(1:4, 5), x=seq(-1, 1, length.out=20),
    group=factor(rep(1:5, each=4)))
  expect_error(ordinalEstimatedContext(data, fixed=cbind(1, data$x)),
    "confounded with the fixed design")
  expect_error(ordinalEstimatedContext(data,
    fixed=cbind(data$x, 1 - data$x)), "confounded with the fixed design")
  expect_error(ordinalEstimatedContext(data,
    fixed=Matrix::sparse.model.matrix(~group-1, data)),
    "confounded with the fixed design")
  expect_error(ordinalEstimatedContext(data, memoryMB=1e-6),
    "identification check exceeds the dense memory budget")
  interceptOnly <- ordinalEstimatedContext(data,
    fixed=matrix(numeric(), nrow(data), 0), fixedVariance=TRUE)
  value <- sommer:::sommer_laplace_objective_eval(interceptOnly$context,
    interceptOnly$cutpoints, tolerance=1e-11)
  expect_true(value$valid, info=value$reason)
  expect_equal(length(value$gradient), 3L)
  data$y[data$y == 2] <- 1
  expect_error(ordinalEstimatedContext(data), "every declared category")
})

test_that("joint ordinal cutpoint scores differentiate the full marginal likelihood", {
  data <- data.frame(y=rep(c(1, 2, 3, 4, 2, 3), 5),
    x=rep(c(-1, 0, 1), 10), group=factor(rep(1:5, each=6)))
  for(link in c("logit", "probit")){
    for(thresholds in list(c(-1, 0.4, 1.8), 0.2, c(-0.2, -0.199, 1),
                          c(-0.2, -0.19999, 1))){
      dataset <- data
      if(length(thresholds) == 1) dataset$y <- ifelse(data$y <= 2, 1, 2)
      setup <- ordinalEstimatedContext(dataset, link, thresholds)
      parameters <- c(0.25, setup$covariance, setup$cutpoints)
      evaluate <- function(parameters, gradient=TRUE){
        sommer:::sommer_laplace_objective_eval(setup$context, parameters,
          tolerance=1e-11, gradient=gradient)
      }
      evaluation <- evaluate(parameters)
      expect_true(evaluation$valid, info=evaluation$reason)
      numerical <- vapply(seq_along(parameters), function(position){
        upper <- lower <- parameters
        upper[position] <- upper[position] + 1e-4
        lower[position] <- lower[position] - 1e-4
        high <- evaluate(upper, FALSE)
        low <- evaluate(lower, FALSE)
        expect_true(high$valid, info=high$reason)
        expect_true(low$valid, info=low$reason)
        (high$logLik - low$logLik) / 2e-4
      }, numeric(1))
      expect_lte(max(abs(as.numeric(evaluation$gradient) - numerical)), 2e-5)
      for(value in c(-1000, 1000)){
        invalid <- parameters
        invalid[length(invalid)] <- value
        if(length(thresholds) > 1){
          rejected <- evaluate(invalid)
          expect_false(rejected$valid)
          expect_match(rejected$reason, "family.*numerical range")
        }
      }
      restored <- evaluate(parameters)
      expect_equal(restored$logLik, evaluation$logLik, tolerance=1e-12)
      expect_equal(restored$gradient, evaluation$gradient, tolerance=1e-12)
      expect_error(evaluate(parameters[-1]), "incompatible length")
    }
  }
})

test_that("joint native ordinal optimization agrees with independent fixed-effect likelihood", {
  set.seed(3191)
  n <- 320
  data <- data.frame(x=rnorm(n), group=factor(rep(1:40, each=8)))
  for(link in c("logit", "probit")){
    latent <- 0.6 * data$x + if(link == "logit") rlogis(n) else rnorm(n)
    data$y <- as.integer(cut(latent, c(-Inf, -1, 0.2, 1.3, Inf)))
    setup <- ordinalEstimatedContext(data, link, fixedVariance=TRUE, variance=1e-8)
    # A vanishing, fixed random variance supplies an independent polr limit.
    parameters <- c(0, setup$cutpoints)
    evaluate <- function(parameters){
      value <- sommer:::sommer_laplace_objective_eval(setup$context, parameters,
        tolerance=1e-11)
      if(!value$valid) stop(value$reason)
      value
    }
    native <- optim(parameters, function(x){
      value <- sommer:::sommer_laplace_objective_eval(setup$context, x,
        tolerance=1e-11, gradient=FALSE)
      if(!value$valid) return(Inf)
      -value$logLik
    },
      function(x) -as.numeric(evaluate(x)$gradient), method="BFGS",
      control=list(reltol=1e-11, maxit=200))
    expect_equal(native$convergence, 0L)
    expect_lt(max(abs(evaluate(native$par)$gradient)), 1e-3)
    reference <- MASS::polr(ordered(y)~x, data=data,
      method=if(link == "logit") "logistic" else "probit",
      control=list(reltol=1e-11, maxit=200))
    cutpoints <- cumsum(c(native$par[2], exp(native$par[3:4])))
    expect_lte(abs(native$value + as.numeric(logLik(reference))), 1e-5)
    expect_lte(max(abs(c(native$par[1], cutpoints) -
      c(coef(reference), reference$zeta))), 1e-4)
  }
})

test_that("joint cutpoint parameters preserve mixed-family gradient ordering", {
  data <- data.frame(y=rep(c(1, 2, 1, 2, 3, 4, -0.5, 0.8), 5),
    x=rep(c(-1, 0, 1, 0.5), 10),
    group=factor(rep(1:5, each=8)))
  index <- rep(c(1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L), 5)
  specifications <- c(ordinalKernelSpec("logit", 0.1),
    ordinalKernelSpec("probit"), sommer:::.laplace_family_specs(gaussian(),
      list(dispersion=0.7))$specifications)
  specifications[[1]]$free <- specifications[[2]]$free <- TRUE
  descriptor <- with(data, vsm(ism(group), sigma2=0.16))$covStruct
  create <- function(fixed){
    sommer:::sommer_laplace_objective_create(
      Matrix::Matrix(fixed, sparse=TRUE),
      Matrix::sparse.model.matrix(~group-1, data), data$y, rep(0, nrow(data)),
      "", list(descriptor), list(Matrix::Diagonal(5)),
      familySpecs=specifications, familyIndex=index, trials=rep(1, nrow(data)))
  }
  context <- create(cbind(data$x))
  parameters <- c(0.25, descriptor$par[descriptor$free], 0.1,
    -1, log(diff(c(-1, 0.4, 1.8))), log(0.7))
  evaluate <- function(parameters){
    sommer:::sommer_laplace_objective_eval(context, parameters, tolerance=1e-11)
  }
  result <- evaluate(parameters)
  expect_true(result$valid, info=result$reason)
  numerical <- vapply(seq_along(parameters), function(position){
    upper <- lower <- parameters
    upper[position] <- upper[position] + 1e-4
    lower[position] <- lower[position] - 1e-4
    high <- evaluate(upper)
    low <- evaluate(lower)
    expect_true(high$valid, info=high$reason)
    expect_true(low$valid, info=low$reason)
    (high$logLik - low$logLik) / 2e-4
  }, numeric(1))
  expect_lte(max(abs(as.numeric(result$gradient) - numerical)), 2e-5)
  expect_error(create(cbind(data$x, as.numeric(index == 1))),
    "confounded with the fixed design")
  expect_error(create(cbind(data$x, as.numeric(index <= 2))),
    "confounded with the fixed design")
  invalid <- parameters
  invalid[3] <- 0.9
  invalid[5] <- -1000
  rejected <- evaluate(invalid)
  expect_false(rejected$valid)
  expect_match(rejected$reason, "family parameters outside numerical range")
  restored <- evaluate(parameters)
  expect_equal(restored$logLik, result$logLik, tolerance=1e-12)
  expect_equal(restored$gradient, result$gradient, tolerance=1e-12)
})

test_that("fixed-cutpoint ordinal derivatives propagate through the Laplace objective", {
  data <- data.frame(y=rep(c(1, 2, 3, 4, 2, 3), 4),
    x=rep(c(-1, 0, 1), 8), group=factor(rep(1:4, each=6)))
  descriptor <- with(data, vsm(ism(group), sigma2=0.4))$covStruct
  for(link in c("logit", "probit")){
    context <- sommer:::sommer_laplace_objective_create(
      Matrix::Matrix(cbind(1, data$x), sparse=TRUE),
      Matrix::sparse.model.matrix(~group-1, data), data$y, rep(0, nrow(data)),
      "", list(descriptor), list(Matrix::Diagonal(4)),
      familySpecs=ordinalKernelSpec(link), familyIndex=rep(1L, nrow(data)),
      trials=rep(1, nrow(data)))
    parameters <- c(0.2, 0.15, descriptor$par[descriptor$free])
    evaluation <- sommer:::sommer_laplace_objective_eval(context, parameters,
      tolerance=1e-11)
    expect_true(evaluation$valid, info=evaluation$reason)
    numerical <- vapply(seq_along(parameters), function(position){
      upper <- lower <- parameters
      upper[position] <- upper[position] + 1e-4
      lower[position] <- lower[position] - 1e-4
      upperFit <- sommer:::sommer_laplace_objective_eval(context, upper,
        tolerance=1e-11, gradient=FALSE)
      lowerFit <- sommer:::sommer_laplace_objective_eval(context, lower,
        tolerance=1e-11, gradient=FALSE)
      expect_true(upperFit$valid, info=upperFit$reason)
      expect_true(lowerFit$valid, info=lowerFit$reason)
      (upperFit$logLik - lowerFit$logLik) / 2e-4
    }, numeric(1))
    expect_lte(max(abs(as.numeric(evaluation$gradient) - numerical)), 2e-5)
  }
})
