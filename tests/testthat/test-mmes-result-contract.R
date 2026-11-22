test_that("the documented result schema remains stable", {
  expected <- c(
    "call", "method", "control", "executionPlan", "engine",
    "fixedFormula", "randomFormula", "residualFormula", "family", "contrasts",
    "b", "u", "bu", "uList", "randomMode", "groupParameters",
    "fitted.values", "linear.predictors", "residuals", "y", "offset", "priorWeights",
    "data", "dataOriginal", "obsInfo", "fixedDesign", "randomDesign",
    "randomMetadata", "randomRelationships", "predictionSettings",
    "Dtable", "partitions", "partitionsX",
    "nobs", "npar", "llik", "AIC", "BIC", "deviance", "dispersion",
    "REML", "likelihoodMethod", "convergence", "monitor",
    "theta", "covPar", "covParNative", "covParCovariance", "covStruct",
    "vcParams", "theta_se", "InfMat", "Ci", "Ci_11", "Cchol", "VarU",
    "uPevList", "toBoundary", "remlOmissions",
    "nonlinear", "nonlinearParameters", "nonlinearCoordinates",
    "nonlinearRandom", "nonlinearInference", "specification", "familyParameters",
    "familySpecifications", "zeroCoefficients", "thresholds", "category.probabilities",
    "optimizer", "engineDiagnostics", "pqlMonitor", "pqlConverged",
    "laplaceDiagnostics", "nonlinearDiagnostics", "specificationDiagnostics")
  expect_identical(sommer:::.mmes_result_fields, expected)
})

test_that("fitted model routes share the documented result fields", {
  set.seed(741)
  data <- data.frame(group=factor(rep(seq_len(8), each=12)),
    x=rep(seq(0, 3, length.out=12), 8))
  data$y <- 3*exp(-.6*data$x)+rep(rnorm(8, sd=.2), each=12)+
    rnorm(nrow(data), sd=.12)
  data$count <- rpois(nrow(data),
    exp(.5+.25*data$x+rep(rnorm(8, sd=.3), each=12)))
  control <- mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  nonlinear_control <- mmes_control(optimizer=list(maxit=200L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  ml <- mmes_method(objective="ml")
  laplace <- mmes_method(objective="ml", integration="laplace")
  curve <- list(growth=nonlinear_curve("x"))
  start <- list(parameters=c(asymptote=3, rate=.6))

  fits <- list(
    linear=mmes(y~x, data=data, control=control),
    mixed=mmes(y~x, random=~group, data=data, control=control),
    direct=mmes(y~x, random=~group, data=data,
      control=mmes_control(formulation="direct",
        diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    glm=mmes(count~x, data=data, family=poisson(), method=ml, control=control),
    pql=mmes(count~x, random=~group, data=data, family=poisson(),
      control=control),
    glmm=mmes(count~x, random=~group, data=data, family=poisson(),
      method=laplace, control=control),
    nonlinear=mmes(y~0+curve_term("growth"), data=data, curves=curve,
      start=start, method=ml, control=nonlinear_control),
    nonlinear_mixed=mmes(y~0+curve_term("growth"), random=~group, data=data,
      curves=curve, start=start,
      method=mmes_method(objective="ml", integration="exact"), control=nonlinear_control),
    nonlinear_glmm=mmes(count~0+curve_term("growth"), random=~group, data=data,
      family=poisson(), curves=curve, start=start, method=laplace,
      control=nonlinear_control)
  )

  fields <- sommer:::.mmes_result_fields
  expect_false(anyDuplicated(fields) > 0L)
  for(name in names(fits)){
    fit <- fits[[name]]
    expect_true(inherits(fit, "mmes"), info=name)
    expect_true(setequal(intersect(names(fit), fields), fields), info=name)
    expect_true(all(fields %in% names(unserialize(serialize(fit, NULL)))), info=name)
    expect_equal(fit$nobs, nrow(data), info=name)
    expect_identical(fit$fixedFormula[[2L]], as.name(if(name %in%
      c("glm", "pql", "glmm", "nonlinear_glmm")) "count" else "y"), info=name)
    expect_true(is.logical(fit$convergence), info=name)
    expect_equal(as.numeric(fit$b), as.numeric(coef(fit))[seq_along(fit$b)],
      info=name)
    expect_true(length(fitted(fit)) > 0L, info=name)
    expect_true(length(residuals(fit)) > 0L, info=name)
    expect_equal(fit$residuals, residuals(fit, type="response"), info=name)
  }
  for(name in c("linear", "mixed", "direct")){
    expect_equal(fits[[name]]$fitted.values, fitted(fits[[name]]), info=name)
    expect_equal(fits[[name]]$linear.predictors, fitted(fits[[name]]), info=name)
    expect_identical(fits[[name]]$family$family, "gaussian", info=name)
    expect_null(fits[[name]]$nonlinearParameters, info=name)
  }
  for(name in c("glm", "glmm", "nonlinear", "nonlinear_mixed", "nonlinear_glmm")){
    expect_equal(fits[[name]]$AIC, AIC(fits[[name]]), info=name)
    expect_equal(fits[[name]]$BIC, BIC(fits[[name]]), info=name)
  }
  expect_true(all(is.na(fits$pql$llik)))
  expect_true(is.na(fits$pql$AIC))
  expect_true(is.na(fits$pql$BIC))
  expect_null(fits$linear$randomFormula)
  expect_null(fits$glm$u)
  expect_null(fits$nonlinear$Ci)
  expect_null(fits$nonlinear_mixed$Ci)
  expect_true(length(fits$nonlinear_mixed$nonlinearParameters) > 0L)
  expect_identical(fits$linear$likelihoodMethod, "Gaussian REML")
  expect_identical(fits$pql$likelihoodMethod, "PQL (REML working fits)")
  expect_identical(fits$glmm$likelihoodMethod, "Laplace ML")
  expect_false(fits$nonlinear$REML)
  for(name in c("mixed", "direct", "glmm", "nonlinear", "nonlinear_mixed")){
    cleared <- postPEV(fits[[name]], mode=0)
    expect_true(all(fields %in% names(cleared)), info=name)
  }
  for(name in c("mixed", "direct", "glmm")){
    enriched <- postPEV(fits[[name]], mode=1)
    expect_true(all(fields %in% names(enriched)), info=name)
    expect_equal(fitted(enriched), fitted(fits[[name]]), info=name)
  }
  for(mode in 0:1){
    with_variance <- postVarU(fits$mixed, mode=mode)
    expect_true(all(fields %in% names(with_variance)), info=as.character(mode))
  }
  nonlinear_post <- post_nonlinear(fits$nonlinear_mixed, mode=0)
  expect_true(all(fields %in% names(nonlinear_post)))
  expect_null(nonlinear_post$uPevList)
})

test_that("standardization retains existing values and skips non-fit results", {
  existing <- structure(list(b=c(a=1), u=c(g=2), y=c(1, 2),
    fitted.values=c(1, 2), residuals=c(0, 0), linear.predictors=c(3, 4), AIC=NA_real_,
    llik=NA_real_, npar=2L, nobs=2L, nonlinearParameters=c(a=3)),
    class="mmes")
  normalized <- sommer:::.mmes_standardize_result(existing, y~x, ~g, ~units,
    poisson())
  for(name in names(existing))
    expect_identical(normalized[[name]], existing[[name]], info=name)
  expect_true(all(sommer:::.mmes_result_fields %in% names(normalized)))
  expect_null(normalized[["Ci"]])
  expect_identical(normalized$family$family, "poisson")
  setup <- list(X=matrix(1), y=1)
  expect_identical(sommer:::.mmes_standardize_result(setup, y~1, NULL, NULL,
    gaussian()), setup)
  solution <- structure(setup, class="mmesSolve")
  expect_identical(sommer:::.mmes_standardize_result(solution, y~1, NULL, NULL,
    gaussian()), solution)
})

test_that("categorical models keep model-aware fitted values and residuals", {
  set.seed(742)
  data <- data.frame(group=factor(rep(seq_len(8), each=12)),
    x=rep(seq(-1, 1, length.out=12), 8))
  data$ordinal <- ordered(cut(.5*data$x+rlogis(nrow(data)),
    c(-Inf, -.5, .7, Inf)), labels=c("low", "middle", "high"))
  data$nominal <- factor(sample(c("a", "b", "c"), nrow(data), replace=TRUE))
  control <- mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  method <- mmes_method(objective="ml", integration="laplace")
  ordinal_fit <- mmes(ordinal~x,
    random=~vsm(ism(group), sigma2=.3, fixedSigma2=TRUE),
    data=data, family=ordinal(), method=method, control=control,
    familyParameters=list(thresholds=c(-.5, .7), estimate=FALSE))
  nominal_fit <- mmes(nominal~x,
    random=~vsm(ism(group), sigma2=.3, fixedSigma2=TRUE),
    data=data, family=multinomial(reference="c"), method=method, control=control)
  for(fit in list(ordinal_fit, nominal_fit)){
    expect_true(all(sommer:::.mmes_result_fields %in% names(fit)))
    expect_null(fit$residuals)
    expect_equal(nrow(fitted(fit)), nrow(data))
    expect_equal(nrow(as.matrix(residuals(fit))), nrow(data))
    expect_equal(fitted(unserialize(serialize(fit, NULL))), fitted(fit))
    expect_equal(predict(fit, newdata=data[1:3,]), fitted(fit)[1:3,],
      tolerance=1e-8)
  }
})

test_that("parameter random effects and missing observations retain their meaning", {
  set.seed(743)
  data <- expand.grid(x=seq(-1.8, 1.8, length.out=7),
    group=factor(letters[1:3]), environment=factor(c("A", "B")))
  data$y <- 1+.3*(data$environment=="B")+1.7*tanh(.9*data$x)+
    rnorm(nrow(data), sd=.18)
  data$y[1] <- NA
  fit <- mmes(y~environment+a*tanh(k*x), data=data,
    parameters=list(a=parameter(), k=parameter(link="log",
      random=~vsm(ism(group), sigma2=.025, fixedSigma2=TRUE))),
    start=c(a=1.5, k=1), method=mmes_method(objective="ml"),
    control=mmes_control(optimizer=list(maxit=200L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(all(sommer:::.mmes_result_fields %in% names(fit)))
  expect_equal(fit$nobs, nrow(data)-1L)
  expect_equal(fit$u, unlist(fit$randomMode, use.names=TRUE))
  expect_equal(fit$residuals, residuals(fit))
  expect_equal(predict(fit, target="conditional"), fitted(fit), tolerance=1e-8)
})

test_that("mixed-family categorical results keep their existing representations", {
  set.seed(744)
  data <- data.frame(x=rep(seq(-1, 1, length.out=24), 2),
    group=factor(rep(rep(seq_len(6), each=4), 2)),
    trait=factor(rep(c("count", "ord"), each=24), levels=c("count", "ord")))
  data$y <- c(rpois(24, exp(.2+.3*data$x[1:24])),
    as.integer(cut(.5*data$x[25:48]+rlogis(24), c(-Inf, -.5, .7, Inf))))
  fit <- mmes(y~0+trait+trait:x,
    random=~vsm(ism(group), sigma2=.3, fixedSigma2=TRUE), data=data,
    family=familym("trait", count=poisson(), ord=ordinal()),
    familyParameters=list(ord=list(thresholds=c(-.5, .7), estimate=FALSE)),
    method=mmes_method(objective="ml", integration="laplace"),
    control=mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(all(sommer:::.mmes_result_fields %in% names(fit)))
  expect_null(fit$residuals)
  expect_named(fitted(fit), c("count", "ord"))
  expect_named(residuals(fit), c("count", "ord"))
  expect_equal(fitted(unserialize(serialize(fit, NULL))), fitted(fit))
})

test_that("categorical nonlinear fits retain model-aware outputs", {
  set.seed(746)
  data <- data.frame(group=factor(rep(seq_len(8), each=12)),
    x=rep(seq(0, 3, length.out=12), 8))
  data$y <- ordered(cut(1.5*exp(-.6*data$x)+rlogis(nrow(data)),
    c(-Inf, -.5, .7, Inf)), labels=c("low", "middle", "high"))
  fit <- mmes(y~0+curve_term("growth"),
    random=~vsm(ism(group), sigma2=.1, fixedSigma2=TRUE), data=data,
    family=ordinal(), curves=list(growth=nonlinear_curve("x")),
    start=list(parameters=c(asymptote=1.5, rate=.6)),
    familyParameters=list(thresholds=c(-.5, .7), estimate=FALSE),
    method=mmes_method(objective="ml", integration="laplace"),
    control=mmes_control(optimizer=list(maxit=200L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_s3_class(fit, "mmes.nonlinear")
  expect_true(all(sommer:::.mmes_result_fields %in% names(fit)))
  expect_null(fit$residuals)
  expect_equal(nrow(fitted(fit)), nrow(data))
  expect_length(residuals(fit), nrow(data))
  expect_equal(fitted(unserialize(serialize(fit, NULL))), fitted(fit))
})

test_that("normalization preserves lean output and non-fit tasks", {
  set.seed(745)
  data <- data.frame(group=factor(rep(seq_len(6), each=6)),
    x=rep(seq_len(6), 6))
  data$y <- 1+.3*data$x+rep(rnorm(6), each=6)+rnorm(nrow(data), sd=.2)
  data$y[1] <- NA
  diagnostics <- list(verbose=FALSE, dateWarning=FALSE)
  fit <- mmes(y~x, random=~group, data=data,
    control=mmes_control(output=list(covariance="none"), diagnostics=diagnostics))
  expect_equal(fit$nobs, nrow(data)-1L)
  expect_equal(fit$fitted.values, fitted(fit))
  expect_equal(fit$residuals, residuals(fit))
  expect_true(all(is.na(fit$covParNative$StdError)))
  expect_null(fit$InfMat)
  expect_error(covparams_mmes_se(fit), "Covariance information was omitted")
  setup <- mmes(y~x, random=~group, data=data,
    control=mmes_control(output=list(setup=TRUE), diagnostics=diagnostics))
  expect_false(inherits(setup, "mmes"))
  expect_false("fitted.values" %in% names(setup))
  solved <- mmes(y~x, random=~group, data=data, covPar=fit,
    method=mmes_method(task="solve"),
    control=mmes_control(diagnostics=diagnostics))
  expect_s3_class(solved, "mmesSolve")
  expect_false("fitted.values" %in% names(solved))
  expect_equal(as.numeric(solved$fitted), as.numeric(fitted(fit)), tolerance=1e-5)
})
