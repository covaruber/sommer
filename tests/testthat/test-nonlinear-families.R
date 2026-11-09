nonlinearFamilyFixture <- function(family, response, randomCurve = TRUE, specifications = NULL, index = NULL,
  trials = NULL, data = NULL, random = ~vsm(ism(group), sigma2 = 0.08), curve = "exponential", coordinate = 2L) {
  n <- length(response)
  if (is.null(data))
    data <- data.frame(x = seq(0.1, 3, length.out = n), group = factor(rep(1:2, length.out = n)))
  data$y <- as.numeric(response)
  control <- mmes_control(task = list(objective = "ml"), diagnostics = list(verbose = FALSE, dateWarning = FALSE))
  prepared <- sommer:::.nonlinear_prepare(y ~ 0, random, data, control, "exclude", "exclude", "exclude",
    "exclude", NULL, data$x)
  if (is.null(specifications))
    specifications <- sommer:::.laplace_family_specs(family, list())$specifications
  if (is.null(index))
    index <- rep(0L, n)
  if (is.null(trials))
    trials <- rep(1, n)
  context <- sommer:::sommer_nonlinear_family_create(prepared$X, do.call(cbind, prepared$Z), data$y,
    numeric(n), data$x, curve, if (randomCurve)
      coordinate else 0L, if (randomCurve)
      0L else sum(vapply(prepared$Z, ncol, integer(1))), prepared$covStruct[seq_along(prepared$rtermss)],
    prepared$Ai, specifications, as.integer(index), trials)
  familyStart <- unlist(lapply(specifications, function(s) {
    if (!s$free)
      return(numeric())
    if (s$family == "ordinal")
      sommer:::.laplace_ordinal_coordinates(s$thresholds) else if (!is.null(s$zeroDesign))
      s$zeroCoefficients else s$logParameter
  }))
  covarianceStart <- unlist(lapply(prepared$covStruct[seq_along(prepared$rtermss)], function(s) s$par[s$free]),
    use.names = FALSE)
  list(context = context, prepared = prepared, parameters = c(log(0.7), log(0.8), if (curve == "logistic") 0.5,
    covarianceStart, familyStart))
}

nonlinearFamilyData <- function(family, curveRandom = FALSE) {
  set.seed(715)
  data <- data.frame(x = rep(seq(0, 5, length.out = 16), 12), group = factor(rep(1:12, each = 16)))
  u <- rep(rnorm(12, sd = 0.2), each = 16)
  eta <- if (curveRandom)
    1.4 * exp(-0.65 * exp(u) * data$x) else 1.4 * exp(-0.65 * data$x) + u
  data$zero <- rep(c(-1, 1), length.out = nrow(data))
  code <- if (grepl("^Negative Binomial", family$family))
    "nbinom2" else family$family
  data$y <- switch(code, binomial = rbinom(nrow(data), 1, family$linkinv(eta)), poisson = rpois(nrow(data),
    exp(eta)), nbinom2 = rnbinom(nrow(data), size = 3, mu = exp(eta)), ziPoisson = {
    y <- rpois(nrow(data), exp(eta))
    y[runif(nrow(data)) < 0.25] <- 0
    y
  }, Gamma = rgamma(nrow(data), shape = 4, scale = exp(eta)/4), gaussian = rnorm(nrow(data), eta, sd = 0.5),
    ordinal = {
      p <- cbind(family$linkinv(-0.6 - eta), family$linkinv(0.7 - eta) - family$linkinv(-0.6 -
        eta), 1 - family$linkinv(0.7 - eta))
      code <- 1L + rowSums(runif(nrow(data)) > t(apply(p, 1, cumsum))[, -3, drop = FALSE])
      ordered(c("low", "middle", "high")[code], levels = c("low", "middle", "high"))
    }, multinomial = {
      p <- cbind(exp(eta), exp(eta - 0.7), 1)
      p <- p/rowSums(p)
      code <- 1L + rowSums(runif(nrow(data)) > t(apply(p, 1, cumsum))[, -3, drop = FALSE])
      factor(c("red", "green", "blue")[code], levels = c("red", "green", "blue"))
    })
  data
}

nonlinearFamilyControl <- function(...) {
  mmes_control(task = list(objective = "ml", nonlinear = list(start = c(asymptote = 1.4, rate = 0.65),
    population = list(order = 8L, maxOrder = 64L, tol = 1e-06))), convergence = list(maxIters = 250L),
    diagnostics = list(verbose = FALSE, dateWarning = FALSE), ...)
}

test_that("categorical labels, references and explicit covariance coordinates survive prediction", {
  data <- nonlinearFamilyData(multinomial())
  fit <- mmes(y ~ 1, data = data, family = multinomial(reference = "red"), nonlinear = nonlinear_curve("x",
    random = list(rate = ~vsm(ism(.nominalCategory), ism(group), sigma2 = 0.04, fixedSigma2 = TRUE))),
    control = nonlinearFamilyControl())
  expect_true(fit$convergence)
  expect_identical(fit$familySpecifications[[1]]$reference, "red")
  expect_equal(colnames(fitted(fit)), levels(data$y))
  expect_equal(unname(predict(fit, newdata = data)), unname(fitted(fit)), tolerance = 1e-08)
  expect_true(all(is.finite(predict(fit, newdata = data[1:4, ], re.form = "population"))))
  expect_identical(levels(simulate(fit, seed = 12)[, 1]), levels(data$y))
  expect_error(mmes(y ~ 0, random = ~group, data = transform(data, y = ordered(y)), family = multinomial(),
    nonlinear = nonlinear_curve("x"), control = nonlinearFamilyControl()), "unordered")
  expect_error(mmes(y ~ 0, random = ~group, data = data, family = ordinal(), nonlinear = nonlinear_curve("x"),
    control = nonlinearFamilyControl()), "ordered factors")
  data <- nonlinearFamilyData(ordinal())
  data$x[c(3, 8)] <- NA
  fit <- mmes(y ~ 1, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data, family = ordinal(),
    nonlinear = nonlinear_curve("x"), control = nonlinearFamilyControl())
  expect_equal(nobs(fit), nrow(data) - 2L)
  expect_equal(colnames(fitted(fit)), levels(data$y))
  data <- nonlinearFamilyData(ordinal())
  data$y[data$y == "middle"] <- "high"
  control <- nonlinearFamilyControl()
  control$task$glmm$familyParameters <- list(thresholds = c(-0.6, 0.7), estimate = FALSE)
  fixedThreshold <- mmes(y ~ 0, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data,
    family = ordinal(), nonlinear = nonlinear_curve("x"), control = control)
  expect_equal(colnames(fitted(fixedThreshold)), levels(data$y))
  expect_equal(fixedThreshold$familySpecifications[[1]]$thresholds, c(-0.6, 0.7))
  data$y <- factor(data$y, levels = levels(data$y), ordered=FALSE)
  expect_error(mmes(y ~ 0, random = ~group, data = data, family = multinomial(), nonlinear = nonlinear_curve("x"),
    control = nonlinearFamilyControl()), "every category")
})

test_that("family covariance estimates, RR priors and relationships persist in public fits", {
  data <- nonlinearFamilyData(poisson(), TRUE)
  data$environment <- factor(rep(1:3, length.out = nrow(data)))
  relationship <- matrix(0.15, 12L, 12L)
  diag(relationship) <- 1
  dimnames(relationship) <- list(levels(data$group), levels(data$group))
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  for (curveRandom in c(FALSE, TRUE)) {
    formula <- ~vsm(rrm(environment, k = 1, loadings = matrix(c(1, 0.7, -0.3), 3L, 1L), fixed = rep(TRUE,
      3L)), ism(group), Gu = precision, sigma2 = 0.04)
    model <- nonlinear_curve("x", random = if (curveRandom)
      list(rate = formula) else list())
    fit <- mmes(y ~ 0, random = if (curveRandom)
      NULL else formula, data = data, family = poisson(), nonlinear = model, control = nonlinearFamilyControl())
    expect_true(fit$convergence)
    expect_equal(unname(fit$covStruct[[1]]$par[fit$covStruct[[1]]$free]), unname(fit$optimizer$par[3]),
      tolerance = 1e-12)
    expect_equal(length(fit$randomMode), 36L)
    expect_equal(length(fit$latentRandomMode), 12L)
    expect_equal(predict(fit, newdata = data), fitted(fit), tolerance = 1e-08)
    expect_true(all(is.finite(predict(fit, newdata = data[1:4, ], re.form = "population"))))
    expect_equal(dim(simulate(fit, 2, seed = 20, re.form = "population")), c(nrow(data), 2L))
    expect_equal(predict(unserialize(serialize(fit, NULL)), newdata = data[1:3, ]), predict(fit,
      newdata = data[1:3, ]))
  }
})

test_that("generative population simulation matches nonlinear prior integration", {
  data <- nonlinearFamilyData(multinomial(), TRUE)
  fit <- mmes(y ~ 0, data = data, family = multinomial(reference = "red"), nonlinear = nonlinear_curve("x",
    random = list(rate = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE))), control = nonlinearFamilyControl())
  population <- predict(fit, re.form = "population")
  draws <- simulate(fit, nsim = 1000, seed = 711, re.form = "population")
  frequencies <- vapply(levels(data$y), function(level) mean(as.matrix(draws) == level), numeric(1))
  expect_equal(frequencies, colMeans(population), tolerance = 0.025)
  local <- fit
  local$control$task$nonlinear$population <- list(order = 2L, maxOrder = 4L, tol = 1e-20)
  expect_error(predict(local, re.form = "population"), "did not meet")
  expect_error(predict(fit, newdata = data.frame()), "nonempty")
  expect_error(predict(fit, newdata = 1:3), "data frame")
})

test_that("grouped trials, nuisance choices and missingness retain existing contracts", {
  data <- nonlinearFamilyData(binomial())
  set.seed(701)
  data$trials <- rep(c(3L, 5L), length.out = nrow(data))
  data$success <- rbinom(nrow(data), data$trials, plogis(1.4 * exp(-0.65 * data$x)))
  training <- data
  fit <- function(fixed = binm(success, trials = trials) ~ 0, data = training, control = nonlinearFamilyControl()) mmes(fixed,
    random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data, family = binomial(),
    nonlinear = nonlinear_curve("x"), control = control)
  grouped <- fit()
  counts <- fit(cbind(success, trials - success) ~ 0)
  expect_equal(coef(grouped), coef(counts), tolerance = 1e-08)
  expect_equal(as.numeric(logLik(grouped)), as.numeric(logLik(counts)), tolerance = 1e-08)
  expect_identical(as.numeric(grouped$priorWeights), as.numeric(data$trials))
  draw <- simulate(grouped, seed = 5)[, 1]
  expect_equal(draw * data$trials, round(draw * data$trials), tolerance = 1e-12)
  expect_true(all(draw >= 0 & draw <= 1))
  missing <- data
  missing$x[c(3, 8)] <- NA
  filtered <- fit(data = missing)
  expect_equal(nobs(filtered), nrow(data) - 2L)
  expect_false(any(filtered$obsInfo$included[c(3, 8)]))
  expect_error(mmes(binm(success, trials = trials) ~ 0, random = ~group, data = missing, family = binomial(),
    nonlinear = nonlinear_curve("x"), control = nonlinearFamilyControl(), naMethodX = "fail"), "Missing nonlinear")
  data <- nonlinearFamilyData(nbinom2(theta = 3))
  control <- nonlinearFamilyControl()
  control$task$glmm$familyParameters <- list(theta = 3, estimate = FALSE)
  fixedShape <- mmes(y ~ 0, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data,
    family = nbinom2(theta = 1), nonlinear = nonlinear_curve("x"), control = control)
  expect_equal(fixedShape$family$theta, 3)
  expect_false(fixedShape$familySpecifications[[1]]$free)
  expect_equal(fixedShape$npar, 2L)
  expect_error(predict(fixedShape, newdata = transform(data[1:3, ], group = "unseen")), "unseen")
  expect_error(simulate(fixedShape, nsim = 2, maxMemoryMB = 1e-04), "exceeds")
  expect_error(simulate(fixedShape, seed = Inf), "seed")
  expect_error(simulate(fixedShape, nsim = 1.5), "integer")
  local <- fixedShape
  local$control$execution$denseMemoryMB <- 1e-04
  expect_error(predict(local, newdata = data), "workspace")
})

test_that("family scores retain structured and reduced-rank loading derivatives", {
  data <- expand.grid(x = seq(0.1, 3, length.out = 5), environment = factor(1:3), group = factor(1:2))
  response <- rep(c(0, 1, 2), length.out = nrow(data))
  for (curve in c("exponential", "logistic", "saturation")) {
    for (randomCurve in c(FALSE, TRUE)) {
      fixture <- nonlinearFamilyFixture(poisson(), response, randomCurve, data = data, random = ~vsm(rrm(environment,
        k = 1), ism(group), sigma2 = 0.04), curve = curve)
      evaluate <- function(par) sommer:::sommer_nonlinear_family_eval(fixture$context, par)
      value <- evaluate(fixture$parameters)
      expect_true(value$valid, info = value$reason)
      score <- vapply(seq_along(fixture$parameters), function(j) {
        up <- down <- fixture$parameters
        up[j] <- up[j] + 1e-05
        down[j] <- down[j] - 1e-05
        a <- evaluate(up)
        b <- evaluate(down)
        expect_true(a$valid, info = a$reason)
        expect_true(b$valid, info = b$reason)
        (a$logLik - b$logLik)/2e-05
      }, numeric(1))
      expect_equal(as.numeric(value$gradient), score, tolerance = 5e-05)
      expect_equal(length(value$mode), 6L)
      expect_equal(length(value$latentMode), 2L)
    }
  }
})

test_that("mixed-family models expose complete categorical and scalar outputs", {
  families <- list(count = poisson(), ordered = ordinal(), nominal = multinomial(), zero = zipoisson(zero = ~zero),
    continuous = Gamma("log"), binary = binomial(), normal = gaussian())
  pieces <- lapply(names(families), function(trait) {
    data <- nonlinearFamilyData(families[[trait]], TRUE)
    if (is.factor(data$y))
      data$y <- as.numeric(data$y)
    data$trait <- trait
    data
  })
  data <- do.call(rbind, pieces)
  family <- do.call(familym, c(list(by = "trait"), families))
  fit <- mmes(y ~ 0, data = data, family = family, nonlinear = nonlinear_curve("x", random = list(rate = ~vsm(ism(.familyCoordinate),
    ism(group), sigma2 = 0.04, fixedSigma2 = TRUE))), control = nonlinearFamilyControl())
  expect_true(fit$convergence)
  expect_equal(nobs(fit), nrow(data))
  for (type in c("response", "link")) {
    fitted <- fitted(fit, type = type)
    expect_named(fitted, names(families))
    expect_equal(unname(predict(fit, newdata = data, type = type)), unname(fitted), tolerance = 1e-09)
    population <- predict(fit, re.form = "population", type = type)
    expect_true(all(is.finite(unlist(population))))
    expect_equal(predict(fit, newdata = data, re.form = "population", type = type), population, tolerance = 1e-08)
  }
  expect_equal(dim(fitted(fit)$nominal), c(192L, 3L))
  expect_equal(dim(fitted(fit)$ordered), c(192L, 3L))
  expect_equal(dim(simulate(fit, 2, seed = 10, re.form = "population")), c(nrow(data), 2L))
  expect_true(all(is.finite(unlist(residuals(fit, type = "pearson")))))
  expect_error(residuals(fit, type = "quantile"), "unordered")
  expect_null(fit$nonlinearInference)
  local <- unserialize(serialize(fit, NULL))
  expect_equal(predict(local, newdata = data[1:12, ]), predict(fit, newdata = data[1:12, ]))
})

test_that("population outputs integrate priors, not conditional uncertainty", {
  data <- nonlinearFamilyData(poisson())
  data$groupCurve <- data$group
  fit <- mmes(y ~ 0, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data, family = poisson(),
    nonlinear = nonlinear_curve("x"), control = nonlinearFamilyControl())
  expect_equal(predict(fit, re.form = "population"), predict(fit, re.form = NA) * exp(0.04/2), tolerance = 1e-09)
  expect_equal(predict(fit, re.form = "population", type = "link"), predict(fit, re.form = NA, type = "link"),
    tolerance = 1e-10)
  fit <- mmes(y ~ 0, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data, family = poisson(),
    nonlinear = nonlinear_curve("x", random = list(rate = ~vsm(ism(groupCurve), sigma2 = 0.04, fixedSigma2 = TRUE))),
    control = nonlinearFamilyControl())
  population <- predict(fit, newdata = data[1:4, ], re.form = "population")
  oracle <- vapply(data$x[1:4], function(x) {
    mean <- function(b) exp(fit$nonlinearParameters["asymptote"] * exp(-fit$nonlinearParameters["rate"] *
      exp(b) * x))
    integrate(function(b) mean(b) * dnorm(b, sd = 0.2), -8, 8, rel.tol = 1e-09)$value * exp(0.04/2)
  }, numeric(1))
  expect_equal(population, oracle, tolerance = 2e-06)
})

test_that("all nonlinear scalar likelihoods have complete analytic scores", {
  cases <- list(list(binomial(), rep(c(0, 1, 1, 0, 1), 4)), list(binomial("probit"), rep(c(0, 1, 1,
    0, 1), 4)), list(binomial("cloglog"), rep(c(0, 1, 1, 0, 1), 4)), list(poisson(), rep(c(0, 1,
    2, 4, 1), 4)), list(nbinom2(theta = 2), rep(c(0, 1, 2, 4, 1), 4)), list(zipoisson(), rep(c(0,
    1, 0, 4, 1), 4)), list(Gamma("log"), rep(c(0.4, 1, 2, 4, 1), 4)), list(gaussian(), rep(c(0.4,
    1, 2, 4, 1), 4)))
  for (case in cases) for (randomCurve in c(FALSE, TRUE)) {
    fixture <- nonlinearFamilyFixture(case[[1]], case[[2]], randomCurve)
    evaluate <- function(par) sommer:::sommer_nonlinear_family_eval(fixture$context, par)
    par <- fixture$parameters
    value <- evaluate(par)
    expect_true(value$valid, info = value$reason)
    if (!value$valid)
      next
    step <- 1e-05
    score <- vapply(seq_along(par), function(j) {
      up <- down <- par
      up[j] <- up[j] + step
      down[j] <- down[j] - step
      a <- evaluate(up)
      b <- evaluate(down)
      expect_true(a$valid, info = a$reason)
      expect_true(b$valid, info = b$reason)
      (a$logLik - b$logLik)/(2 * step)
    }, numeric(1))
    expect_equal(as.numeric(value$gradient), score, tolerance = 3e-05)
    point <- as.numeric(value$latentMode) + c(0.03, -0.05)
    conditional <- function(w) sommer:::sommer_nonlinear_family_conditional(fixture$context, par,
      w)
    hessian <- vapply(seq_along(point), function(j) {
      up <- down <- point
      up[j] <- up[j] + step
      down[j] <- down[j] - step
      -as.numeric(conditional(up)$score - conditional(down)$score)/(2 * step)
    }, numeric(length(point)))
    expect_equal(as.matrix(conditional(point)$hessian), hessian, tolerance = 3e-06)
  }
})

test_that("ordinal and nominal block likelihood derivatives are complete", {
  for (link in c("logit", "probit")) {
    family <- ordinal(link = link)
    specifications <- sommer:::.laplace_family_specs(family, list())$specifications
    specifications[[1]]$thresholds <- c(-0.4, 0.6)
    fixture <- nonlinearFamilyFixture(family, rep(1:3, 8), specifications = specifications)
    par <- fixture$parameters
    evaluate <- function(p) sommer:::sommer_nonlinear_family_eval(fixture$context, p)
    value <- evaluate(par)
    expect_true(value$valid, info = value$reason)
    score <- vapply(seq_along(par), function(j) {
      up <- down <- par
      up[j] <- up[j] + 1e-05
      down[j] <- down[j] - 1e-05
      (evaluate(up)$logLik - evaluate(down)$logLik)/2e-05
    }, numeric(1))
    expect_equal(as.numeric(value$gradient), score, tolerance = 3e-05)
  }
  specifications <- sommer:::.laplace_family_specs(multinomial(), list())$specifications
  specifications[[1]]$nominalGroups <- matrix(1:24, ncol = 2, byrow = TRUE)
  fixture <- nonlinearFamilyFixture(multinomial(), rep(rep(1:3, 4), each = 2), specifications = specifications)
  par <- fixture$parameters
  evaluate <- function(p) sommer:::sommer_nonlinear_family_eval(fixture$context, p)
  value <- evaluate(par)
  expect_true(value$valid, info = value$reason)
  score <- vapply(seq_along(par), function(j) {
    up <- down <- par
    up[j] <- up[j] + 1e-05
    down[j] <- down[j] - 1e-05
    (evaluate(up)$logLik - evaluate(down)$logLik)/2e-05
  }, numeric(1))
  expect_equal(as.numeric(value$gradient), score, tolerance = 3e-05)
  conditional <- function(mode) sommer:::sommer_nonlinear_family_conditional(fixture$context, par,
    mode)
  point <- as.numeric(value$latentMode) + c(0.05, -0.03)
  hessian <- vapply(seq_along(point), function(j) {
    up <- down <- point
    up[j] <- up[j] + 1e-05
    down[j] <- down[j] - 1e-05
    -as.numeric(conditional(up)$score - conditional(down)$score)/2e-05
  }, numeric(length(point)))
  expect_equal(as.matrix(conditional(point)$hessian), hessian, tolerance = 2e-06)
  expect_gt(abs(hessian[1, 2]), 1e-05)
  expect_identical(evaluate(par)$logLik, value$logLik)
  expect_error(conditional(numeric(3)), "length")
  invalid <- par
  invalid[1] <- NA
  expect_false(evaluate(invalid)$valid)
})

test_that("covariance constraints preserve family tails and narrow cutpoint scores", {
  data <- data.frame(x = seq(0.1, 3, length.out = 24), group = factor(rep(1:2, each = 12)), cluster = factor(rep(1:3,
    8)))
  fixture <- nonlinearFamilyFixture(nbinom2(theta = 2), rep(0:3, 6), FALSE, data = data, random = ~vsm(ism(group),
    sigma2 = 0.08) + vsm(ism(cluster), sigma2 = 0.16))
  prepared <- fixture$prepared
  mapping <- sommer:::.laplace_covariance_map(prepared$covStruct[seq_along(prepared$rtermss)], prepared$rtermss,
    data.frame(parameter = 1:2, group = 1L, scale = c(1, 2)))
  layout <- list(fixedCount = 2L, covarianceCount = length(mapping$map$initial), covarianceMap = mapping$map)
  par <- c(log(0.7), log(0.8), mapping$map$initial, log(2))
  evaluate <- function(par) sommer:::.laplace_constrained_evaluate(par, layout, function(expanded) sommer:::sommer_nonlinear_family_eval(fixture$context,
    expanded))
  value <- evaluate(par)
  expect_true(value$valid, info = value$reason)
  score <- vapply(seq_along(par), function(j) {
    up <- down <- par
    up[j] <- up[j] + 1e-05
    down[j] <- down[j] - 1e-05
    (evaluate(up)$logLik - evaluate(down)$logLik)/2e-05
  }, numeric(1))
  expect_equal(as.numeric(value$gradient), score, tolerance = 5e-05)
  for (link in c("logit", "probit")) {
    specifications <- sommer:::.laplace_family_specs(ordinal(link), list())$specifications
    specifications[[1]]$thresholds <- c(0, 1e-07)
    fixture <- nonlinearFamilyFixture(ordinal(link), rep(1:3, 8), specifications = specifications)
    par <- fixture$parameters
    evaluate <- function(par) sommer:::sommer_nonlinear_family_eval(fixture$context, par)
    value <- evaluate(par)
    expect_true(value$valid, info = value$reason)
    score <- vapply(seq_along(par), function(j) {
      up <- down <- par
      up[j] <- up[j] + 1e-04
      down[j] <- down[j] - 1e-04
      (evaluate(up)$logLik - evaluate(down)$logLik)/2e-04
    }, numeric(1))
    expect_equal(as.numeric(value$gradient), score, tolerance = 1e-04)
  }
})

for (family in list(binomial(), binomial("probit"), binomial("cloglog"), poisson(), nbinom2(theta = 3),
  zipoisson(), Gamma("log"), ordinal(), ordinal("probit"), multinomial())) {
  test_that(paste("public nonlinear surfaces work for", family$family, family$link), {
    data <- nonlinearFamilyData(family)
    fit <- mmes(y ~ 0, random = ~vsm(ism(group), sigma2 = 0.04, fixedSigma2 = TRUE), data = data,
      family = family, nonlinear = nonlinear_curve("x"), control = nonlinearFamilyControl())
    expect_true(fit$convergence, info = family$family)
    expect_identical(fit$executionPlan$objective, "ml (Laplace family)")
    expect_true(all(is.finite(coef(fit))))
    expect_true(is.finite(as.numeric(logLik(fit))))
    expect_equal(nobs(fit), nrow(data))
    for (type in c("response", "link")) {
      expect_equal(predict(fit, type = type), fitted(fit, type = type), tolerance = 1e-09)
      expect_equal(unname(predict(fit, newdata = data, type = type)), unname(fitted(fit, type = type)),
        tolerance = 1e-09)
      expect_true(all(is.finite(predict(fit, re.form = NA, type = type))))
      population <- predict(fit, re.form = "population", type = type)
      expect_true(all(is.finite(population)))
      expect_equal(unname(predict(fit, newdata = data, re.form = "population", type = type)), unname(population),
        tolerance = 1e-08)
    }
    if (family$family %in% c("ordinal", "multinomial")) {
      expect_equal(unname(rowSums(fitted(fit))), rep(1, nrow(data)), tolerance = 1e-10)
      expect_equal(colnames(fitted(fit)), levels(data$y))
      expect_equal(rowSums(predict(fit, re.form = "population")), rep(1, nrow(data)), tolerance = 1e-08,
        ignore_attr = TRUE)
    }

    expect_false(fit$nonlinearDiagnostics$uncertaintyComputed)
    expect_null(fit$optimizer$hessian)
    types <- if (family$family == "multinomial")
      c("response", "pearson", "deviance") else c("response", "pearson", "deviance", "working", "quantile")
    for (type in types) expect_true(all(is.finite(unlist(residuals(fit, type = type)))))
    set.seed(8)
    before <- .Random.seed
    simulation <- simulate(fit, nsim = 2, seed = 17)
    expect_identical(.Random.seed, before)
    expect_equal(dim(simulation), c(nrow(data), 2L))
    expect_identical(simulate(fit, nsim = 2, seed = 17), simulation)
    expect_equal(dim(simulate(fit, nsim = 2, seed = 17, re.form = "population")), c(nrow(data), 2L))
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(predict(restored, newdata = data[1:4, ]), predict(fit, newdata = data[1:4, ]))
    expect_null(fit$nonlinearInference)
  })
}

test_that("mixed blocks and zero-regression parameters have complete native scores", {
  specifications <- list(ordinal = list(family = "ordinal", link = "probit", logParameter = 0, free = TRUE,
    thresholds = c(-0.5, 0.8)), nominal = list(family = "multinomial", link = "softmax", logParameter = 0,
    free = FALSE, nominalGroups = matrix(13:24, ncol = 2, byrow = TRUE)), count = list(family = "ziPoisson",
    link = "log", logParameter = qlogis(0.25), free = TRUE, zeroDesign = cbind(1, rep(c(-1, 1), 18)),
    zeroOffset = rep(0.1, 36), zeroCoefficients = c(qlogis(0.25), 0.2)))
  response <- c(rep(1:3, 4), rep(rep(1:3, 2), each = 2), rep(c(0, 1, 2), 4))
  fixture <- nonlinearFamilyFixture(poisson(), response, specifications = specifications, index = rep(0:2,
    c(12, 12, 12)))
  evaluate <- function(par) sommer:::sommer_nonlinear_family_eval(fixture$context, par)
  value <- evaluate(fixture$parameters)
  expect_true(value$valid, info = value$reason)
  score <- vapply(seq_along(fixture$parameters), function(j) {
    up <- down <- fixture$parameters
    up[j] <- up[j] + 1e-05
    down[j] <- down[j] - 1e-05
    a <- evaluate(up)
    b <- evaluate(down)
    expect_true(a$valid, info = a$reason)
    expect_true(b$valid, info = b$reason)
    (a$logLik - b$logLik)/2e-05
  }, numeric(1))
  expect_equal(as.numeric(value$gradient), score, tolerance = 5e-05)
})

test_that("nonlinear family Laplace likelihood has normalized integration and Gaussian limits", {
  fixture <- nonlinearFamilyFixture(poisson(), rep(c(0, 1, 2, 1), 4))
  par <- fixture$parameters
  value <- sommer:::sommer_nonlinear_family_eval(fixture$context, par)
  x <- seq(0.1, 3, length.out = 16)
  y <- rep(c(0, 1, 2, 1), 4)
  exact <- sum(vapply(1:2, function(group) {
    rows <- seq.int(group, 16, 2)
    density <- function(w) vapply(w, function(b) sum(dpois(y[rows], exp(0.7 * exp(-0.8 * exp(b) *
      x[rows])), log = TRUE)) + dnorm(b, 0, sqrt(0.08), log = TRUE), numeric(1))
    shift <- max(density(seq(-2, 2, length.out = 100)))
    shift + log(integrate(function(w) exp(density(w) - shift), -Inf, Inf, rel.tol = 1e-09)$value)
  }, numeric(1)))
  expect_lt(abs(value$logLik - exact), 0.015)
  par[3] <- log(1e-10)
  limit <- sommer:::sommer_nonlinear_family_eval(fixture$context, par)
  expect_true(limit$valid, info = limit$reason)
  expect_equal(limit$logLik, sum(dpois(y, exp(0.7 * exp(-0.8 * x)), log = TRUE)), tolerance = 2e-08)
  fixture <- nonlinearFamilyFixture(gaussian(), rep(c(0.4, 1, 2, 4, 1), 4), FALSE)
  value <- sommer:::sommer_nonlinear_family_eval(fixture$context, fixture$parameters)
  n <- 20L
  z <- model.matrix(~factor(rep(1:2, length.out = n)) - 1)
  mean <- 0.7 * exp(-0.8 * seq(0.1, 3, length.out = n))
  variance <- diag(n) + 0.08 * tcrossprod(z)
  residual <- rep(c(0.4, 1, 2, 4, 1), 4) - mean
  exact <- -0.5 * (n * log(2 * pi) + as.numeric(determinant(variance, logarithm = TRUE)$modulus) +
    sum(residual * solve(variance, residual)))
  expect_equal(value$logLik, exact, tolerance = 1e-10)
})

for (family in list(binomial(), poisson(), nbinom2(theta = 3), zipoisson(), Gamma("log"), ordinal(),
  multinomial())) {
  test_that(paste("random curve family methods work for", family$family), {
    data <- nonlinearFamilyData(family, TRUE)
    fit <- mmes(y ~ 0, data = data, family = family, nonlinear = nonlinear_curve("x", random = list(rate = ~vsm(ism(group),
      sigma2 = 0.04, fixedSigma2 = TRUE))), control = nonlinearFamilyControl())
    expect_true(fit$convergence)
    expect_lt(fit$nonlinearDiagnostics$modeScore, 1e-06)
    expect_equal(unname(predict(fit, newdata = data)), unname(fitted(fit)), tolerance = 1e-09)
    expect_true(all(is.finite(predict(fit, re.form = "population"))))
    population <- predict(fit, re.form = "population")
    subset <- if (is.matrix(population))
      population[1:4, , drop = FALSE] else population[1:4]
    expect_equal(unname(predict(fit, newdata = data[1:4, ], re.form = "population")), unname(subset),
      tolerance = 1e-08)
    expect_equal(dim(simulate(fit, nsim = 2, seed = 12, re.form = "population")), c(nrow(data), 2L))
    expect_false(fit$nonlinearDiagnostics$uncertaintyComputed)
  })
}
