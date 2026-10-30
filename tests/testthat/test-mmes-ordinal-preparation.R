ordinalPreparationData <- function(){
  data.frame(y=ordered(rep(c("low", "medium", "high"), 8),
      levels=c("low", "medium", "high")),
    x=seq(-1, 1, length.out=24),
    treatment=factor(rep(c("a", "b"), each=12)),
    group=factor(rep(1:6, each=4)), exposure=seq(0.5, 1, length.out=24))
}

prepareOrdinalFixture <- function(data=ordinalPreparationData(), fixed=y~x,
                                   family=ordinal(), settings=list(),
                                   naMethodY="exclude", naMethodX="exclude",
                                   naMethodRandom="exclude", contrasts=NULL){
  control <- mmes_control(task=list(objective="ml", glmm=list(method="laplace",
    familyParameters=settings)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  random <- as.formula("~group", env=environment(fixed))
  sommer:::.mmes_prepare_model(fixed, random, NULL, data, family, naMethodY,
    list(naMethodX=naMethodX, naMethodRandom=naMethodRandom, naMethodR="exclude",
      contrasts=contrasts, engine=attr(control, "native")),
    familyParameters=control$task$glmm$familyParameters)
}

test_that("ordinal preparation preserves ordered labels and removes only the intercept column", {
  data <- ordinalPreparationData()
  original <- data
  prepared <- prepareOrdinalFixture(data, y~x+treatment+offset(log(exposure)))
  expect_identical(data, original)
  expect_identical(prepared$ordinal$levels, levels(data$y))
  expect_equal(prepared$ordinal$response, as.numeric(data$y))
  expect_identical(which(prepared$obsInfo$included), seq_len(nrow(data)))
  full <- model.matrix(y~x+treatment+offset(log(exposure)), data)
  expect_equal(unname(as.matrix(prepared$X)), unname(full[, -1, drop=FALSE]))
  expect_identical(prepared$ordinal$columns, 2:3)
  expect_identical(prepared$ordinal$fullColumns, c("Intercept", "x", "treatmentb"))
  expect_identical(prepared$ordinal$xlevels$treatment, levels(data$treatment))
  expect_equal(prepared$ordinal$assignment, c(1L, 2L))
  expect_equal(prepared$partitionsX$x, matrix(1L))
  expect_equal(prepared$partitionsX$treatment, matrix(2L))
  expect_false("1" %in% names(prepared$partitionsX))
  expect_equal(model.offset(prepared$modelFrame), log(data$exposure))
  expect_equal(prepared$ordinal$specification$thresholds, qlogis(c(1/3, 2/3)))
  expect_true(prepared$ordinal$specification$free)
  expect_false(prepared$ordinal$specification$supplied)
  expect_equal(prepared$ordinal$counts, rep(8L, 3))
  expect_equal(sommer:::.laplace_ordinal_thresholds(prepared$ordinal$initialParameters),
    qlogis(c(1/3, 2/3)))

  noIntercept <- prepareOrdinalFixture(data, y~0+x)
  expect_equal(as.matrix(noIntercept$X), as.matrix(prepareOrdinalFixture(data)$X))
  expect_error(prepareOrdinalFixture(data, y~0+treatment),
    "confounded with the fixed design")
  data$constant <- 1
  expect_error(prepareOrdinalFixture(data, y~constant+x),
    "confounded with the fixed design")
  data$duplicate <- data$x
  expect_error(prepareOrdinalFixture(data, y~x+duplicate), "full-rank")
  interceptOnly <- prepareOrdinalFixture(data, y~1)
  expect_equal(ncol(interceptOnly$X), 0L)
  expect_length(interceptOnly$partitionsX, 0)
})

test_that("ordinal controls use retained category counts and honor fixed cutpoints", {
  data <- ordinalPreparationData()
  data$x[c(1, 4)] <- NA_real_
  data$y[2] <- NA
  data$group[3] <- NA
  for(link in c("logit", "probit")){
    prepared <- prepareOrdinalFixture(data, family=ordinal(link))
    included <- setdiff(seq_len(nrow(data)), 1:4)
    expect_equal(which(prepared$obsInfo$included), included)
    expect_equal(prepared$ordinal$response, as.numeric(data$y[included]))
    counts <- tabulate(as.integer(data$y[included]), 3)
    expect_equal(prepared$ordinal$counts, counts)
    cumulative <- head(cumsum(counts)/sum(counts), -1)
    expected <- if(link == "logit") qlogis(cumulative) else qnorm(cumulative)
    expect_equal(prepared$ordinal$specification$thresholds, expected)
    expect_equal(nrow(prepared$X), length(included))
    expect_equal(nrow(prepared$Z[[1]]), length(included))
    omit <- prepareOrdinalFixture(data, family=ordinal(link), naMethodY="omit",
      naMethodX="omit", naMethodRandom="omit")
    expect_equal(omit$X, prepared$X)
    expect_equal(omit$ordinal$specification, prepared$ordinal$specification)
  }
  expect_error(prepareOrdinalFixture(data, naMethodX="fail"), "fixed-effect variables")

  supplied <- prepareOrdinalFixture(family=ordinal(thresholds=c(-2, 1)))
  expect_equal(supplied$ordinal$specification$thresholds, c(-2, 1))
  expect_true(supplied$ordinal$specification$supplied)
  controlled <- prepareOrdinalFixture(family=ordinal(thresholds=c(-2, 1)),
    settings=list(thresholds=c(-1, 2), estimate=FALSE))
  expect_equal(controlled$ordinal$specification$thresholds, c(-1, 2))
  expect_false(controlled$ordinal$specification$free)
  expect_length(controlled$ordinal$initialParameters, 0)
  expect_identical(colnames(controlled$X), c("Intercept", "x"))
  expect_true("1" %in% names(controlled$partitionsX))
  expect_error(prepareOrdinalFixture(settings=list(estimate=FALSE)),
    "Fixed ordinal cutpoints require")
  expect_error(prepareOrdinalFixture(settings=list(thresholds=0)), "one fewer")
  expect_error(prepareOrdinalFixture(settings=list(unknown=TRUE)),
    "Unknown familyParameters")
  expect_error(prepareOrdinalFixture(settings=list(estimate=1)), "single TRUE/FALSE")
  for(value in list(c(0, 0), c(1, -1), NA_real_, Inf, 1+1i, matrix(c(-1, 1)))){
    expect_error(prepareOrdinalFixture(settings=list(thresholds=value)),
      "finite and strictly increasing")
  }
})

test_that("ordinal categories cannot be silently reordered or collapsed", {
  data <- ordinalPreparationData()
  data$y <- ordered(as.character(data$y), levels=c("high", "medium", "low"))
  prepared <- prepareOrdinalFixture(data)
  expect_identical(prepared$ordinal$levels, c("high", "medium", "low"))
  expect_equal(prepared$ordinal$response, as.numeric(data$y))

  for(value in list(as.character(data$y), factor(data$y, ordered=FALSE),
                    matrix(as.integer(data$y)))){
    data$y <- value
    expect_error(prepareOrdinalFixture(data), "ordered factors or plain numeric")
  }
  for(value in list(rep(1L, 24), rep(c(1L, 3L), 12),
                    rep(c(0, 1, 2), 8), rep(c(1, 2.5, 3), 8),
                    rep(c(1, Inf, 3), 8), rep(NA_real_, 24))){
    data$y <- value
    expect_error(prepareOrdinalFixture(data), "Ordinal")
  }
  data <- ordinalPreparationData()
  data$y <- as.integer(data$y)
  numericPrepared <- prepareOrdinalFixture(data)
  expect_identical(numericPrepared$ordinal$levels, c("1", "2", "3"))
  expect_equal(numericPrepared$ordinal$response, data$y)
  data$y <- ordered(data$y, levels=1:4)
  expect_error(prepareOrdinalFixture(data), "every declared category")
  fixed <- prepareOrdinalFixture(data, settings=list(thresholds=c(-1, 0, 1),
    estimate=FALSE))
  expect_identical(fixed$ordinal$levels, as.character(1:4))
  expect_equal(fixed$ordinal$counts, c(8L, 8L, 8L, 0L))

  data <- ordinalPreparationData()
  data$x[data$y == "medium"] <- NA
  expect_error(prepareOrdinalFixture(data), "retained observations in every declared category")
  data$y <- as.integer(data$y)
  expect_error(prepareOrdinalFixture(data), "retained observations in every declared category")
  expect_error(prepareOrdinalFixture(family=familym("trait", score=ordinal())),
    "ordinal familym")
})

test_that("ordinal preparation respects environments, transformations and contrasts", {
  data <- ordinalPreparationData()
  prepared <- prepareOrdinalFixture(data, y~I(x^2)+treatment,
    contrasts=list(treatment="contr.sum"))
  full <- model.matrix(y~I(x^2)+treatment, data,
    contrasts.arg=list(treatment="contr.sum"))
  expect_equal(unname(as.matrix(prepared$X)), unname(full[, -1, drop=FALSE]))
  expect_identical(prepared$ordinal$contrasts, list(treatment="contr.sum"))
  fixed <- with(data, y~x+offset(log(exposure)))
  environment(fixed) <- list2env(as.list(data), parent=environment())
  noData <- prepareOrdinalFixture(data=NULL, fixed=fixed)
  expect_equal(as.matrix(noData$X), as.matrix(prepareOrdinalFixture(data)$X))
  expect_identical(noData$ordinal$levels, levels(data$y))
  expect_identical(noData$ordinal$originalFixed, fixed)
})

test_that("ordinal coordinate layout and preparation connect to native joint optimization", {
  prepared <- prepareOrdinalFixture()
  specification <- prepared$ordinal$specification
  expect_identical(sommer:::.laplace_family_count(specification), 2L)
  proposal <- sommer:::.laplace_proposal_specs(list(specification),
    c(0.1, -1, log(3)), list(fixedCount=1, covarianceCount=0))
  expect_equal(proposal[[1]]$thresholds, c(-1, 2))
  expect_equal(proposal[[1]]$logParameter, 0)
  for(coordinates in list(numeric(), NA_real_, c(0, -1000), c(0, 1000), 1+1i,
                          matrix(c(0, 1)))){
    expect_error(sommer:::.laplace_ordinal_thresholds(coordinates), "Ordinal cutpoint")
  }
  expect_error(sommer:::.laplace_ordinal_coordinates(c(-1e308, 1e308)),
    "gaps exceed")
  expect_error(sommer:::.laplace_ordinal_coordinates(NULL), "nonempty")
  expect_error(sommer:::.laplace_family_count(list(family="ordinal", free=TRUE)),
    "response-prepared cutpoints")

  descriptor <- prepared$covStruct[[1]]
  context <- sommer:::sommer_laplace_objective_create(prepared$X,
    do.call(cbind, prepared$Z), prepared$ordinal$response,
    rep(0, nrow(prepared$X)), "", list(descriptor), prepared$Ai,
    familySpecs=list(specification), familyIndex=rep(1L, nrow(prepared$X)),
    trials=rep(1, nrow(prepared$X)))
  parameters <- c(rep(0, ncol(prepared$X)), descriptor$par[descriptor$free],
    prepared$ordinal$initialParameters)
  evaluation <- sommer:::sommer_laplace_objective_eval(context, parameters)
  expect_true(evaluation$valid, info=evaluation$reason)
  expect_true(all(is.finite(evaluation$gradient)))
  expect_equal(length(evaluation$gradient), length(parameters))
  fixed <- specification
  fixed$free <- FALSE
  expect_identical(sommer:::.laplace_family_count(fixed), 0L)
})

test_that("prepared ordinal factor models reproduce independent likelihood fits", {
  set.seed(7325)
  data <- data.frame(x=rnorm(240), treatment=factor(rep(c("a", "b"), 120)),
    group=factor(rep(1:30, each=8)), exposure=seq(0.5, 1, length.out=240))
  for(link in c("logit", "probit")){
    latent <- 0.6 * data$x + 0.2 * (data$treatment == "b") +
      log(data$exposure) + if(link == "logit") rlogis(240) else rnorm(240)
    data$y <- ordered(cut(latent, c(-Inf, -1, 0.8, Inf),
      labels=c("small", "average", "large")),
      levels=c("small", "average", "large"))
    prepared <- prepareOrdinalFixture(data,
      y~x+treatment+offset(log(exposure)), family=ordinal(link))
    descriptor <- prepared$covStruct[[1]]
    descriptor$free[] <- FALSE
    descriptor$par[1] <- log(1e-8)
    context <- sommer:::sommer_laplace_objective_create(prepared$X,
      do.call(cbind, prepared$Z), prepared$ordinal$response,
      model.offset(prepared$modelFrame)[prepared$obsInfo$included],
      "", list(descriptor), prepared$Ai,
      familySpecs=list(prepared$ordinal$specification),
      familyIndex=rep(1L, nrow(prepared$X)), trials=rep(1, nrow(prepared$X)))
    evaluate <- function(parameters, gradient=TRUE){
      sommer:::sommer_laplace_objective_eval(context, parameters,
        tolerance=1e-11, gradient=gradient)
    }
    initial <- c(rep(0, ncol(prepared$X)), prepared$ordinal$initialParameters)
    if(link == "logit"){
      extreme <- evaluate(c(20.24324, 7.418257, -4.781367, 15.95651), FALSE)
      expect_true(extreme$valid, info=extreme$reason)
      expect_true(any(extreme$observationCurvature == 0))
      expect_true(is.finite(extreme$logLik))
    }
    fit <- nlminb(initial, function(parameters){
      value <- evaluate(parameters, FALSE)
      if(!value$valid) return(Inf)
      -value$logLik
    }, function(parameters) -as.numeric(evaluate(parameters)$gradient),
    control=list(rel.tol=1e-10, iter.max=200, eval.max=1000))
    reference <- MASS::polr(y~x+treatment+offset(log(exposure)), data=data,
      method=if(link == "logit") "logistic" else "probit",
      control=list(reltol=1e-11, maxit=200))
    expect_equal(fit$convergence, 0L)
    expect_lt(max(abs(evaluate(fit$par)$gradient)), 1e-3)
    expect_lte(abs(fit$objective + as.numeric(logLik(reference))), 1e-5)
    expect_lte(max(abs(fit$par[1:2] - coef(reference))), 1e-4)
    thresholds <- sommer:::.laplace_ordinal_thresholds(fit$par[-(1:2)])
    expect_lte(max(abs(thresholds - reference$zeta)), 1e-4)
    expect_identical(prepared$ordinal$levels, levels(data$y))
  }
})

test_that("native GLMM sparse scaling safely handles zero-curvature rows", {
  random <- Matrix::sparseMatrix(i=1:12, j=rep(1:3, each=4), x=1)
  fixed <- Matrix::Matrix(matrix(1, 12, 1), sparse=TRUE)
  for(allSaturated in c(TRUE, FALSE)){
    offset <- if(allSaturated) rep(1000, 12) else rep(c(1000, 0), 6)
    response <- as.numeric(offset > 0)
    evaluation <- sommer:::sommer_laplace_glmm_eval(fixed, random,
      Matrix::Diagonal(3), response, offset, 0, rep(0, 3), "binomial", 0)
    expect_true(evaluation$valid, info=evaluation$reason)
    expect_true(is.finite(evaluation$logLik))
    expect_true(any(evaluation$observationCurvature == 0))
    if(allSaturated){
      expect_equal(as.numeric(evaluation$observationCurvature), rep(0, 12))
      expect_equal(evaluation$logLik, 0, tolerance=1e-12)
      expect_equal(as.numeric(evaluation$mode), rep(0, 3))
    }else expect_true(any(evaluation$observationCurvature > 0))
  }
})
