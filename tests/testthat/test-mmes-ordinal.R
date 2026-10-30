ordinalFitFixture <- function(link="logit", fixedThresholds=FALSE, estimateVariance=FALSE){
  set.seed(781)
  data <- data.frame(x=rnorm(240), group=factor(rep(1:30, each=8)),
    treatment=factor(rep(c("a", "b"), 120)), exposure=seq(0.7, 1, length.out=240))
  eta <- 0.6 * data$x + 0.25 * (data$treatment == "b") +
    log(data$exposure) + rnorm(30, sd=0.4)[data$group]
  latent <- eta + if(link == "logit") rlogis(240) else rnorm(240)
  data$y <- ordered(cut(latent, c(-Inf, -1, 0.8, Inf),
    labels=c("low", "medium", "high")), levels=c("low", "medium", "high"))
  control <- mmes_control(task=list(objective="ml", glmm=list(method="laplace",
    familyParameters=if(fixedThresholds)
      list(thresholds=c(-1, 0.8), estimate=FALSE) else list())),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  random <- if(estimateVariance) ~vsm(ism(group), sigma2=0.16, fixedSigma2=FALSE) else
    ~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE)
  fit <- mmes(y~x+treatment+offset(log(exposure)), random=random,
    data=data, family=ordinal(link), control=control)
  list(fit=fit, data=data, control=control)
}

test_that("ordinal fits expose category distributions and defer uncertainty", {
  for(link in c("logit", "probit")){
    fixture <- ordinalFitFixture(link)
    fit <- fixture$fit
    expect_true(fit$convergence)
    expect_s3_class(fit, "mmes.glmm.laplace")
    expect_null(fit$fixedCovariance)
    expect_null(fit$Ci)
    expect_null(fit$uPevList)
    expect_false("hessian" %in% names(fit$optimizer))
    expect_true(all(is.na(fit$thresholds$StdError)))
    expect_equal(nrow(fit$thresholds), 2)
    expect_true(all(diff(fit$thresholds$estimate) > 0))
    expect_false("Intercept" %in% names(coef(fit)))
    expect_equal(fit$npar, 4)
    expect_equal(AIC(fit), -2 * as.numeric(logLik(fit)) + 8)
    probability <- fitted(fit)
    expect_equal(fitted(fit, type="probabilities"), probability)
    expect_equal(predict(fit, type="probabilities"), probability)
    expect_equal(dim(probability), c(240L, 3L))
    expect_identical(colnames(probability), levels(fixture$data$y))
    expect_lte(max(abs(rowSums(probability) - 1)), 1e-12)
    cdf <- if(link == "logit") plogis else pnorm
    cumulative <- vapply(fit$thresholds$estimate,
      function(threshold) cdf(threshold - fit$linear.predictors), numeric(fit$nobs))
    reference <- cbind(cumulative[, 1], cumulative[, 2]-cumulative[, 1],
      1-cumulative[, 2])
    expect_equal(unname(probability), unname(reference), tolerance=1e-12)
    expect_equal(fitted(fit, type="mean"), as.numeric(probability %*% 1:3))
    expect_equal(fitted(fit, type="link"), fit$linear.predictors)
    expect_equal(predict(fit), probability)
    expect_equal(predict(fit, type="mean"), fitted(fit, type="mean"))
    expect_equal(fit$deviance, -2*sum(log(probability[cbind(1:240, fit$y)])),
      tolerance=1e-10)
    expect_equal(sum(residuals(fit, type="deviance")^2), fit$deviance)
    expect_equal(residuals(fit), fit$y - fitted(fit, type="mean"))
    expect_true(all(is.finite(residuals(fit, type="pearson"))))
    expect_true(all(is.finite(residuals(fit, type="working"))))
    set.seed(1)
    expect_true(all(is.finite(residuals(fit, type="quantile"))))
  }
})

test_that("ordinal fitting never invokes post-fit uncertainty routines", {
  local_mocked_bindings(
    .mmes_laplace_uncertainty=function(...) stop("Fit-time parameter uncertainty"),
    .mmes_laplace_postPEV=function(...) stop("Fit-time PEV"),
    .package="sommer")
  local_mocked_bindings(optimHess=function(...) stop("Fit-time numerical Hessian"),
    .package="stats")
  expect_true(ordinalFitFixture()$fit$convergence)
})

test_that("ordinal prediction preserves fitted contrasts, offsets and category labels", {
  fixture <- ordinalFitFixture()
  fit <- fixture$fit
  data <- fixture$data
  rows <- c(5, 1, 10)
  newdata <- data[rows, ]
  expect_equal(unname(predict(fit, newdata=newdata)),
    unname(fitted(fit)[rows, , drop=FALSE]), tolerance=1e-10)
  expect_equal(predict(fit, newdata=newdata, type="link"),
    fit$linear.predictors[rows], tolerance=1e-10)
  newdata$y <- NULL
  expect_equal(unname(predict(fit, newdata=newdata)),
    unname(fitted(fit)[rows, , drop=FALSE]), tolerance=1e-10)
  eta <- as.numeric(fit$fixedDesign %*% fit$b) + fit$offset
  expect_equal(predict(fit, re.form=NA, type="link"), eta)
  newdata$group <- "unseen"
  expect_error(predict(fit, newdata=newdata), "unseen random-effect")
  expect_equal(dim(predict(fit, newdata=newdata, re.form=NA)), c(3L, 3L))
  expect_equal(dim(predict(fit, newdata=newdata, re.form="population")), c(3L, 3L))
  expect_lte(max(abs(rowSums(predict(fit, re.form="population")) - 1)), 1e-12)
  fixed <- ordinalFitFixture(fixedThresholds=TRUE)$fit
  expect_true("Intercept" %in% names(coef(fixed)))
  expect_equal(fixed$thresholds$estimate, c(-1, 0.8))
  expect_false(any(fixed$thresholds$estimated))
  expect_equal(unname(predict(fixed, newdata=data[rows, ])),
    unname(fitted(fixed)[rows, , drop=FALSE]), tolerance=1e-10)
})

test_that("ordinal simulation preserves levels and category frequencies", {
  fit <- ordinalFitFixture()$fit
  set.seed(42)
  before <- .Random.seed
  simulations <- simulate(fit, nsim=500, seed=82)
  expect_identical(.Random.seed, before)
  expect_identical(simulate(fit, nsim=5, seed=82), simulations[, 1:5])
  expect_true(all(vapply(simulations, is.ordered, logical(1))))
  expect_identical(levels(simulations[[1]]), fit$ordinal$levels)
  frequency <- tabulate(unlist(lapply(simulations, as.integer)), 3) /
    (500 * fit$nobs)
  expect_lte(max(abs(frequency - colMeans(fitted(fit)))), 0.01)
  population <- simulate(fit, nsim=300, seed=74, re.form="population")
  frequency <- tabulate(unlist(lapply(population, as.integer)), 3) /
    (300 * fit$nobs)
  expect_lte(max(abs(frequency - colMeans(predict(fit, re.form="population")))), 0.02)
})

test_that("ordinal post-fit uncertainty includes cutpoints and prediction cross-covariances", {
  fit <- ordinalFitFixture()$fit
  summary <- summary(fit)
  expect_output(print(summary), "Ordinal cutpoints")
  expect_equal(nrow(anova(fit)), 2L)
  expect_true(all(is.finite(summary$thresholds$StdError)))
  expect_true(all(summary$thresholds$StdError > 0))
  expect_true(all(is.na(fit$thresholds$StdError)))
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
  jacobian <- sommer:::.laplace_ordinal_jacobian(fit$familySpecifications[[1]])
  expect_equal(uncertainty$thresholdCovariance,
    jacobian %*% uncertainty$parameterCovariance[3:4, 3:4] %*% t(jacobian),
    ignore_attr=TRUE)
  data <- ordinalFitFixture()$data[1:2, ]
  for(type in c("response", "mean", "link")){
    for(random in list(NA, NULL, "population")){
      if(identical(random, "population") && type == "link") next
      prediction <- predict(fit, newdata=data, type=type, re.form=random,
        se.fit=TRUE, cov.fit=TRUE)
      expect_true(all(is.finite(prediction$se.fit)))
      expect_true(all(prediction$se.fit > 0))
      expect_equal(as.numeric(prediction$se.fit)^2, unname(diag(prediction$cov.fit)),
        tolerance=1e-8)
      if(type == "response"){
        expect_equal(dim(prediction$se.fit), c(2L, 3L))
        expect_equal(dim(prediction$cov.fit), c(6L, 6L))
        # Column-major order: category 1 rows, category 2 rows, category 3 rows.
        sums <- cbind(diag(2), diag(2), diag(2))
        expect_lte(max(abs(sums %*% prediction$cov.fit %*% t(sums))), 1e-8)
      }
    }
  }
  pev <- postPEV(fit, mode=2)
  expect_true(all(is.finite(pev$uPevList[[1]])))
  expect_equal(as.numeric(pev$uPevList[[1]]), unname(diag(pev$Ci)[-(1:2)]))
  expect_equal(pev$Ci[1:2, 1:2], vcov(fit), ignore_attr=TRUE)
  expect_null(fit$Ci)
})

test_that("public ordinal fits agree with independent proportional-odds fits in the zero-variance limit", {
  for(link in c("logit", "probit")){
    fixture <- ordinalFitFixture(link)
    fit <- mmes(y~x+treatment+offset(log(exposure)),
      random=~vsm(ism(group), sigma2=1e-8, fixedSigma2=TRUE),
      data=fixture$data, family=ordinal(link), control=fixture$control)
    reference <- MASS::polr(y~x+treatment+offset(log(exposure)), data=fixture$data,
      method=if(link == "logit") "logistic" else "probit",
      control=list(reltol=1e-12))
    expect_true(fit$convergence)
    expect_equal(unname(coef(fit)), unname(coef(reference)), tolerance=1e-4)
    expect_equal(fit$thresholds$estimate, unname(reference$zeta), tolerance=1e-4)
    expect_lte(abs(as.numeric(logLik(fit))-as.numeric(logLik(reference))), 1e-5)
    expect_equal(unname(fitted(fit)), unname(fitted(reference)), tolerance=2e-5)
  }
})

test_that("ordinal intercept-only, fixed thresholds and binary fits are complete", {
  fixture <- ordinalFitFixture()
  for(link in c("logit", "probit")){
    fit <- mmes(y~1, random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
      data=fixture$data, family=ordinal(link), control=fixture$control)
    expect_true(fit$convergence)
    expect_length(coef(fit), 0L)
    expect_equal(dim(vcov(fit)), c(0L, 0L))
    expect_equal(nrow(summary(fit)$betas), 0L)
    expect_equal(dim(predict(fit, newdata=fixture$data[1:2, ])), c(2L, 3L))
    expect_equal(length(sommer:::.boot_mmes_parameters(fit)), 2L)
    expect_true(all(is.finite(postPEV(fit, mode=1)$uPevList[[1]])))
    data <- fixture$data
    data$y <- ordered(ifelse(data$y == "high", "yes", "no"), levels=c("no", "yes"))
    fit <- mmes(y~x, random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
      data=data, family=ordinal(link), control=fixture$control)
    data$binary <- as.integer(data$y == "yes")
    binary <- mmes(binary~x, random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
      data=data, family=binomial(link), control=fixture$control)
    expect_true(fit$convergence)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(binary)), tolerance=1e-7)
    expect_equal(fit$thresholds$estimate, -unname(binary$b[1]), tolerance=1e-5)
    expect_equal(fitted(fit)[, "yes"], fitted(binary), tolerance=1e-5)
  }
  fixed <- ordinalFitFixture(fixedThresholds=TRUE)$fit
  expect_true(all(is.na(summary(fixed)$thresholds$StdError)))
  expect_true(all(is.finite(vcov(fixed))))
  expect_true(all(is.finite(predict(fixed, re.form=NA, se.fit=TRUE)$se.fit)))
  expect_true(all(is.finite(postPEV(fixed, mode=1)$uPevList[[1]])))
})

test_that("ordinal bootstrap and likelihood comparisons preserve the category experiment", {
  fixture <- ordinalFitFixture()
  fit <- fixture$fit
  reduced <- mmes(y~x+offset(log(exposure)),
    random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
    data=fixture$data, family=ordinal(), control=fixture$control)
  expect_true(reduced$convergence)
  expect_equal(nrow(anova(reduced, fit)), 2)
  bootstrap <- boot_mmes(reduced, fit, nsim=3, seed=18)
  expect_equal(nrow(bootstrap$draws), 3)
  expect_true(all(bootstrap$diagnostics$valid))
  expect_true(all(is.finite(bootstrap$draws)))
  expect_match(colnames(bootstrap$draws)[2], "threshold1")
  expect_match(colnames(bootstrap$draws)[3], "logGap2")
  expect_error(sommer:::.boot_mmes_refit(fit, rep(1, fit$nobs), 100),
    "every declared category")
  incompatible <- fit
  incompatible$ordinal$levels <- c("a", "b", "c")
  expect_error(anova(reduced, incompatible), "famil")
  incompatible <- ordinalFitFixture(fixedThresholds=TRUE)$fit
  changed <- incompatible
  changed$familySpecifications[[1]]$thresholds <- c(-1.1, 0.8)
  expect_error(anova(incompatible, changed), "famil")
})

test_that("ordinal covariance-parameter uncertainty and eligibility guards are explicit", {
  fit <- ordinalFitFixture(estimateVariance=TRUE)$fit
  expect_true(fit$convergence)
  expect_equal(fit$npar, 5)
  expect_true(all(is.finite(summary(fit)$thresholds$StdError)))
  expect_true(all(is.finite(predict(fit, se.fit=TRUE)$se.fit)))
  expect_true(all(is.finite(postPEV(fit, mode=1)$uPevList[[1]])))
  unavailable <- fit
  unavailable$convergence <- FALSE
  expect_warning(expect_error(predict(unavailable, se.fit=TRUE),
    "unavailable"), "did not converge")
  limited <- fit
  limited$laplacePostPEVSettings$denseMemoryMB <- 0.00001
  expect_error(fitted(limited), "denseMemoryMB")
  expect_error(predict(limited, re.form="population"), "denseMemoryMB")
  expect_error(vcov(limited), "denseMemoryMB")
  expect_error(postPEV(limited, mode=1), "denseMemoryMB")
  expect_error(predict(fit, type="class"), "arg")
  expect_error(predict(fit, cov.fit=TRUE), "requires se.fit")
  expect_error(predict(fit, type="link", re.form="population"), "response")
  expect_error(sommer:::.laplace_family_specs(familym("trait", first=ordinal(),
    second=poisson()), list()), "Ordinal familym")
})

test_that("ordinal probability uncertainty matches an independent joint delta method", {
  for(link in c("logit", "probit")){
    fixture <- ordinalFitFixture(link)
    fit <- fixture$fit
    data <- fixture$data[1:2, ]
    eta <- predict(fit, newdata=data, type="link", re.form=NA)
    density <- if(link == "logit") dlogis else dnorm
    cutpoint <- fit$thresholds$estimate
    endpoints <- vapply(cutpoint, function(value) density(value - eta), numeric(2))
    dp <- cbind(-endpoints[, 1], endpoints[, 1] - endpoints[, 2], endpoints[, 2])
    jacobian <- matrix(0, 6, 4)
    for(category in 1:3){
      rows <- (category - 1L) * 2L + 1:2
      jacobian[rows, 1:2] <- dp[, category] * as.matrix(fit$fixedDesign[1:2, ])
      natural <- matrix(0, 2, 2)
      if(category < 3) natural[, category] <- endpoints[, category]
      if(category > 1) natural[, category - 1] <- -endpoints[, category - 1]
      gapMap <- rbind(c(1, 0), c(1, diff(cutpoint)))
      jacobian[rows, 3:4] <- natural %*% gapMap
    }
    covariance <- sommer:::.mmes_laplace_uncertainty(fit)$parameterCovariance
    expected <- jacobian %*% covariance %*% t(jacobian)
    prediction <- predict(fit, newdata=data, re.form=NA, se.fit=TRUE, cov.fit=TRUE)
    expect_equal(unname(prediction$cov.fit), unname(expected), tolerance=1e-7)
  }
})

test_that("ordinal population probabilities match independent integration", {
  for(link in c("logit", "probit")){
    fixture <- ordinalFitFixture(link)
    fit <- fixture$fit
    data <- fixture$data[1:2, ]
    eta <- predict(fit, newdata=data, type="link", re.form=NA)
    cdf <- if(link == "logit") plogis else pnorm
    cumulative <- vapply(fit$thresholds$estimate, function(threshold){
      vapply(eta, function(location) integrate(function(u)
        cdf(threshold - location - u) * dnorm(u, sd=0.4), -Inf, Inf,
        rel.tol=1e-10)$value, numeric(1))
    }, numeric(2))
    expected <- cbind(cumulative[, 1], cumulative[, 2] - cumulative[, 1],
      1 - cumulative[, 2])
    expect_equal(unname(predict(fit, newdata=data, re.form="population")),
      unname(expected), tolerance=1e-8)
  }
})

test_that("ordinal public fits retain category schema after filtering", {
  fixture <- ordinalFitFixture()
  data <- fixture$data
  data$y[data$y == "medium"] <- NA
  expect_error(mmes(y~x, random=~group, data=data, family=ordinal(),
    control=fixture$control), "every declared category")
  fixed <- mmes(y~x, random=~vsm(ism(group), sigma2=0.16, fixedSigma2=TRUE),
    data=data, family=ordinal(thresholds=c(-1, 0.8)),
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace",
      familyParameters=list(estimate=FALSE))),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_true(fixed$convergence)
  expect_identical(colnames(fitted(fixed)), levels(data$y))
  expect_true(all(fitted(fixed)[, "medium"] > 0))
  expect_equal(nrow(predict(fixed, newdata=fixture$data)), 240)
  data$y <- as.integer(fixture$data$y)
  numeric <- mmes(y~x+treatment+offset(log(exposure)),
    random=fixture$fit$randomFormula, data=data, family=ordinal(), control=fixture$control)
  expect_equal(unname(fitted(numeric)), unname(fitted(fixture$fit)))
  expect_identical(colnames(fitted(numeric)), c("1", "2", "3"))
})
