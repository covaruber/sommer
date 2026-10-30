multinomialPublicFixture <- function(){
  set.seed(8)
  group <- factor(rep(seq_len(50L), each=6L))
  x <- stats::rnorm(length(group))
  redEffect <- stats::rnorm(nlevels(group), sd=.7)
  blueEffect <- stats::rnorm(nlevels(group), sd=.7)
  redLogit <- .4 + .5*x + redEffect[group]
  blueLogit <- -.2 - .3*x + blueEffect[group]
  probability <- cbind(exp(redLogit), exp(blueLogit), 1)
  probability <- probability / rowSums(probability)
  data <- data.frame(
    y=factor(apply(probability, 1L, function(value)
      sample(c("low", "mid", "high"), 1L, prob=value)),
      levels=c("low", "mid", "high")),
    x=x, group=group)
  control <- mmes_control(task=list(objective="ml",
    glmm=list(method="laplace", maxit=500L)),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  list(data=data, control=control)
}

test_that("public multinomial fits expose probabilities and post-fit methods", {
  fixture <- multinomialPublicFixture()
  fit <- mmes(y~x, random=~group, data=fixture$data,
    family=multinomial(reference="high"), control=fixture$control)

  expect_s3_class(fit, "mmes.glmm.laplace")
  expect_identical(fit$REML, FALSE)
  expect_true(is.finite(fit$llik))
  expect_equal(as.numeric(logLik(fit)), fit$llik)
  expect_equal(dim(fitted(fit)), c(nrow(fixture$data), 3L))
  expect_equal(fitted(fit), fit$category.probabilities)
  expect_equal(unname(rowSums(fitted(fit))), rep(1, nrow(fixture$data)),
    tolerance=1e-12)
  expect_identical(colnames(fitted(fit)), levels(fixture$data$y))
  expect_equal(fitted(fit, type="link"), fit$linear.predictors)
  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, type="probabilities"), fitted(fit))

  conditional <- predict(fit, newdata=fixture$data[1:4, , drop=FALSE])
  expect_equal(conditional, fitted(fit)[1:4, , drop=FALSE], tolerance=1e-10)
  zeroRandom <- predict(fit, newdata=fixture$data[1:4, , drop=FALSE],
    re.form=NA)
  expect_equal(unname(rowSums(zeroRandom)), rep(1, 4L), tolerance=1e-12)
  expect_error(predict(fit, re.form="population"), "not yet supported")
  expect_error(predict(fit, newdata=transform(fixture$data[1:2, ],
    group=factor(c("unseen", "unseen")))), "unseen random-effect")

  summaryFit <- summary(fit)
  expect_true(all(is.finite(summaryFit$betas$Std.Error)))
  expect_equal(dim(vcov(fit)), c(length(fit$b), length(fit$b)))
  interceptFit <- mmes(y~1, random=~group, data=fixture$data,
    family=multinomial(reference="high"), control=fixture$control)
  expect_s3_class(anova(interceptFit, fit), "anova")
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_equal(dim(residuals(fit)), dim(fitted(fit)))
  expect_equal(dim(residuals(fit, type="pearson")), dim(fitted(fit)))
  expect_error(residuals(fit, type="deviance"), "not defined")

  predictionSE <- predict(fit, newdata=fixture$data[1:3, , drop=FALSE],
    se.fit=TRUE, cov.fit=TRUE)
  expect_equal(dim(predictionSE$se.fit), c(3L, 3L))
  expect_equal(dim(predictionSE$cov.fit), c(9L, 9L))
  expect_true(all(is.finite(predictionSE$se.fit)))

  simulation <- simulate(fit, nsim=2L, seed=10L)
  expect_equal(dim(simulation), c(nrow(fixture$data), 2L))
  expect_true(all(vapply(simulation, is.factor, logical(1))))
  expect_true(all(vapply(simulation, function(value)
    identical(levels(value), levels(fixture$data$y)), logical(1))))
  populationSimulation <- simulate(fit, nsim=1L, seed=10L,
    re.form="population")
  expect_true(is.factor(populationSimulation[[1L]]))
  expect_true(identical(levels(populationSimulation[[1L]]),
    levels(fixture$data$y)))
  expect_null(fit$uPevList)
  postFit <- postPEV(fit, mode=1L)
  expect_identical(postFit$uPevMode, 1L)
  expect_equal(dim(postFit$uPevList[[1L]]), c(nlevels(fixture$data$group), 1L))
  postFull <- postPEV(fit, mode=2L)
  expect_equal(dim(postFull$Ci),
    rep(length(fit$b) + length(fit$randomMode), 2L))
  bootstrap <- boot_mmes(fit, nsim=2L, seed=10L, maxit=300L)
  expect_equal(bootstrap$nsim, 2L)
  expect_equal(nrow(bootstrap$diagnostics), 2L)
  expect_equal(nrow(bootstrap$draws), 2L)
})

test_that("multinomial family requires the native Laplace ML route", {
  fixture <- multinomialPublicFixture()
  expect_s3_class(multinomial(), "family")
  expect_error(mmes(y~x, random=~group, data=fixture$data,
    family=multinomial(), control=mmes_control()),
    "supported only by native Laplace")
  expect_error(mmes(y~x, random=~group, data=fixture$data,
    family=multinomial(), control=mmes_control(task=list(objective="reml",
      glmm=list(method="laplace")))),
    "maximum likelihood")
  expect_error(mmes(y~x, random=~group, data=fixture$data,
    family=multinomial(), control=mmes_control(task=list(objective="ml",
      glmm=list(method="laplace")), output=list(inverse=1L))),
    "defers uncertainty")
})
