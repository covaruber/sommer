zipFormulaData <- function(seed=109L, groups=60L, perGroup=8L){
  set.seed(seed)
  group <- factor(rep(seq_len(groups), each=perGroup))
  x <- rnorm(length(group))
  z <- rnorm(length(group))
  latent <- rnorm(groups, sd=0.4)[group]
  lambda <- exp(0.7 + 0.3*x + latent)
  probability <- plogis(-1 + 0.8*z)
  y <- rpois(length(group), lambda)
  y[rbinom(length(group), 1, probability) == 1] <- 0L
  data.frame(y=y, x=x, z=z, group=group, latent=latent)
}

zipFormulaControl <- function(){
  mmes_control(task=list(objective="ml", glmm=list(method="laplace")),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

test_that("ZIP zero formula validation is bounded and backward compatible", {
  expect_error(zipoisson(zero=y~x), "one-sided")
  expect_error(zipoisson(zero=~(1|group)), "without random effects")
  data <- zipFormulaData()
  old <- mmes(y~x, random=~group, data=data, family=zipoisson(),
    control=zipFormulaControl())
  explicit <- mmes(y~x, random=~group, data=data, family=zipoisson(zero=~1),
    control=zipFormulaControl())
  expect_equal(logLik(old), logLik(explicit), tolerance=1e-10)
  expect_equal(old$familyParameters, explicit$familyParameters)
  expect_null(explicit$zeroCoefficients)
  expect_error(mmes(y~x, random=~group, data=data,
    family=zipoisson(zero=~z), control=mmes_control(task=list(objective="ml",
      glmm=list(method="laplace",
        familyParameters=list(zeroInflation=0.2, estimate=FALSE))))),
    "requires estimate=TRUE")
  expect_error(mmes(y~x, random=~group, data=data,
    family=zipoisson(zero=~z + I(2*z)), control=zipFormulaControl()),
    "full rank")
})

test_that("native zero-design scores differentiate the full Laplace objective", {
  data <- zipFormulaData(groups=4L, perGroup=5L)
  fit <- mmes(y~x, random=~group, data=data, family=zipoisson(zero=~z),
    control=zipFormulaControl())
  context <- sommer:::sommer_laplace_objective_create(
    fit$fixedDesign, fit$randomDesign, fit$laplaceNativeResponse,
    fit$offset, "", fit$covStruct, fit$randomRelationships,
    familySpecs=fit$familySpecifications, familyIndex=fit$laplaceFamilyIndex,
    trials=fit$priorWeights)
  parameters <- fit$optimizer$par
  parameters[] <- c(0.7, 0.3, log(0.16), -1, 0.5)
  evaluate <- function(value) sommer:::sommer_laplace_objective_eval(
    context, value, tolerance=1e-10)
  result <- evaluate(parameters)
  expect_true(result$valid)
  numerical <- vapply(seq_along(parameters), function(position){
    upper <- lower <- parameters
    upper[position] <- upper[position] + 1e-4
    lower[position] <- lower[position] - 1e-4
    (evaluate(upper)$logLik - evaluate(lower)$logLik)/2e-4
  }, numeric(1))
  expect_equal(as.numeric(result$gradient), numerical, tolerance=2e-5)
})

test_that("zero coefficients precede subsequent estimated family parameters", {
  data <- zipFormulaData(groups=10L)
  other <- data
  other$y <- rnbinom(nrow(other), mu=exp(0.7+0.3*other$x+other$latent), size=2)
  data$trait <- "count"
  other$trait <- "overdispersed"
  data <- rbind(data, other)
  data$trait <- factor(data$trait)
  fit <- mmes(y~trait-1+trait:x, random=~group, data=data,
    family=familym("trait", count=zipoisson(zero=~z), overdispersed=nbinom2()),
    control=zipFormulaControl())
  expect_equal(fit$npar, 8)
  context <- sommer:::sommer_laplace_objective_create(
    fit$fixedDesign, fit$randomDesign, fit$laplaceNativeResponse,
    fit$offset, "", fit$covStruct, fit$randomRelationships,
    familySpecs=fit$familySpecifications, familyIndex=fit$laplaceFamilyIndex,
    trials=fit$priorWeights)
  parameters <- c(0.7, 0.7, 0.3, 0.3, log(0.16), -1, 0.5, log(2))
  evaluate <- function(value) sommer:::sommer_laplace_objective_eval(
    context, value, tolerance=1e-10)
  result <- evaluate(parameters)
  expect_true(result$valid)
  numerical <- vapply(seq_along(parameters), function(position){
    upper <- lower <- parameters
    upper[position] <- upper[position] + 1e-4
    lower[position] <- lower[position] - 1e-4
    (evaluate(upper)$logLik-evaluate(lower)$logLik)/2e-4
  }, numeric(1))
  expect_equal(as.numeric(result$gradient), numerical, tolerance=2e-5)
  expect_equal(fit$familyParameters$trait, "overdispersed")
  expect_equal(fit$zeroCoefficients$trait, rep("count", 2))
})

test_that("a no-intercept zero formula keeps its declared basis", {
  data <- zipFormulaData()
  fit <- mmes(y~x, random=~group, data=data, family=zipoisson(zero=~0+z),
    control=zipFormulaControl())
  expect_equal(fit$zeroCoefficients$term, "z")
  expect_equal(fit$npar, 4)
  probability <- plogis(data$z*fit$zeroCoefficients$estimate)
  expect_equal(fit$fitted.values, (1-probability)*exp(fit$linear.predictors))
  expect_equal(predict(fit, newdata=data[1:4, ], type="response"),
    fit$fitted.values[1:4], tolerance=1e-8)
})

test_that("formula ZIP wires fitted means, inference, prediction and diagnostics", {
  data <- zipFormulaData()
  fit <- mmes(y~x, random=~group, data=data, family=zipoisson(zero=~z),
    control=zipFormulaControl())
  expect_true(fit$convergence)
  expect_false(fit$laplaceDiagnostics$uncertaintyComputed)
  expect_null(fit$fixedCovariance)
  expect_null(fit$Ci)
  expect_null(fit$familyParameters)
  expect_equal(fit$npar, 5)
  expect_equal(fit$zeroCoefficients$term, c("(Intercept)", "z"))
  expect_equal(fit$zeroCoefficients$estimate, c(-1, 0.8), tolerance=0.4)
  probability <- plogis(as.numeric(model.matrix(~z, data) %*%
    fit$zeroCoefficients$estimate))
  mean <- (1 - probability)*exp(fit$linear.predictors)
  expect_equal(fit$fitted.values, mean, tolerance=1e-10)
  expect_equal(predict(fit, type="response"), mean, tolerance=1e-10)
  expect_equal(predict(fit, newdata=data[1:4, ], type="response"),
    mean[1:4], tolerance=1e-8)
  summary <- summary(fit)
  expect_true(all(is.finite(summary$zeroCoefficients$StdError)))
  expect_true(all(is.finite(vcov(fit))))
  expect_equal(residuals(fit, type="pearson"),
    (data$y-mean)/sqrt(mean + probability/(1-probability)*mean^2))
  set.seed(12)
  uniforms <- runif(nrow(data))
  count <- exp(fit$linear.predictors)
  lower <- probability + (1-probability)*ppois(data$y-1, count)
  lower[data$y == 0] <- 0
  upper <- probability + (1-probability)*ppois(data$y, count)
  set.seed(12)
  expect_equal(residuals(fit, type="quantile"),
    qnorm(lower + uniforms*(upper-lower)), tolerance=1e-9)
  draws <- simulate(fit, nsim=1000, seed=3)
  expect_equal(mean(as.matrix(draws)), mean(mean), tolerance=0.05)
  population <- predict(fit, newdata=data[1:4, ], type="response", re.form="population")
  expect_equal(population, (1-probability[1:4]) *
    exp(as.numeric(fit$fixedDesign[1:4, ] %*% fit$b) + fit$covParNative$estimate/2))
  for(form in list(NULL, NA, "population")){
    prediction <- predict(fit, newdata=data[1:4, ], type="response",
      re.form=form, se.fit=TRUE, cov.fit=TRUE)
    expect_true(all(is.finite(prediction$se.fit)))
    expect_equal(prediction$se.fit^2, diag(prediction$cov.fit), tolerance=1e-8)
  }
  linkPrediction <- predict(fit, newdata=data[1:4, ], type="link",
    re.form=NA, se.fit=TRUE, cov.fit=TRUE)
  expect_equal(linkPrediction$se.fit^2,
    diag(as.matrix(fit$fixedDesign[1:4, ]) %*% vcov(fit) %*%
      t(as.matrix(fit$fixedDesign[1:4, ]))), tolerance=1e-6,
    ignore_attr=TRUE)
  zeroPrediction <- predict(fit, newdata=data[1:4, ], type="response",
    re.form=NA, se.fit=TRUE, cov.fit=TRUE)
  uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
  mu <- zeroPrediction$fit
  gradient <- cbind(mu*as.matrix(fit$fixedDesign[1:4, ]), 0,
    -mu*probability[1:4]*model.matrix(~z, data[1:4, ]))
  expected <- gradient %*% uncertainty$parameterCovariance %*% t(gradient)
  expect_equal(zeroPrediction$cov.fit, unname(expected), tolerance=1e-6,
    ignore_attr=TRUE)
  newGroup <- data[1:2, ]
  newGroup$group <- c("new-a", "new-b")
  expect_true(all(is.finite(predict(fit, newdata=newGroup, type="response",
    re.form="population", se.fit=TRUE)$se.fit)))
  pev <- postPEV(fit, mode=1)
  expect_true(all(is.finite(unlist(pev$uPevList))))
})

test_that("zero-only missing covariates use the fixed missing-data policy", {
  data <- zipFormulaData()
  data$z[c(2, 9)] <- NA_real_
  fit <- mmes(y~x, random=~group, data=data, family=zipoisson(zero=~z),
    control=zipFormulaControl())
  expect_equal(which(!fit$obsInfo$included), c(2L, 9L))
  expect_true(all(fit$obsInfo$reason[c(2, 9)] == "missing fixed covariate"))
  expect_error(mmes(y~x, random=~group, data=data,
    family=zipoisson(zero=~z), naMethodX="fail", control=zipFormulaControl()),
    "fixed-effect")
  expect_error(predict(fit, newdata=data[2, ], type="response", re.form=NA),
    "missing")
})

test_that("zero designs retain factors, transformations and offsets for new data", {
  data <- zipFormulaData()
  data$f <- factor(rep(c("a", "b"), length.out=nrow(data)))
  data$off <- rep(0.2, nrow(data))
  fit <- mmes(y~x, random=~group, data=data,
    family=zipoisson(zero=~f+I(z^2)+offset(off)), control=zipFormulaControl())
  new <- data[1:5, ]
  new$f <- as.character(new$f)
  probability <- plogis(as.numeric(model.matrix(~f+I(z^2), new) %*%
    fit$zeroCoefficients$estimate)+new$off)
  expect_equal(predict(fit, newdata=new, type="response", re.form=NA),
    (1-probability)*exp(as.numeric(model.matrix(~x, new) %*% fit$b)))
  new$f[1] <- "unseen"
  expect_error(predict(fit, newdata=new, type="response", re.form=NA),
    "new level")
})

test_that("mixed traits restrict zero-only covariates to their ZIP rows", {
  data <- zipFormulaData()
  other <- data
  other$y <- rpois(nrow(other), exp(0.7+0.3*other$x+other$latent))
  data$trait <- "count"
  other$trait <- "ordinary"
  other$z <- NA_real_
  long <- rbind(data, other)
  long$trait <- factor(long$trait)
  fit <- mmes(y~trait-1+trait:x, random=~group, data=long,
    family=familym("trait", count=zipoisson(zero=~z), ordinary=poisson()),
    control=zipFormulaControl())
  expect_equal(fit$nobs, nrow(long))
  expect_equal(fit$zeroCoefficients$trait, rep("count", 2))
  expect_equal(predict(fit, newdata=long[c(2, nrow(data)+2), ],
    type="response"), fit$fitted.values[c(2, nrow(data)+2)], tolerance=1e-8)
  expect_true(all(is.finite(residuals(fit, type="pearson"))))
  expect_true(all(is.finite(as.matrix(simulate(fit, nsim=2, seed=1)))))
  expect_true(all(is.finite(summary(fit)$zeroCoefficients$StdError)))
  expect_equal(predict(fit, newdata=long[nrow(long), , drop=FALSE],
    type="response"), tail(fit$fitted.values, 1), tolerance=1e-8)
})
