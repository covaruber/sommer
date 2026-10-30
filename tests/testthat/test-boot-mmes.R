bootstrapFixture <- function(){
  set.seed(617)
  data <- data.frame(group=factor(rep(seq_len(35), each=8)), x=rnorm(280))
  data$y <- rpois(280, exp(0.8+0.12*data$x+rnorm(35, sd=0.6)[data$group]))
  control <- mmes_control(task=list(objective="ml", glmm=list(method="laplace")),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  list(
    data=data, control=control,
    null=mmes(y~1, random=~group, data=data, family=poisson(), control=control),
    alternative=mmes(y~x, random=~group, data=data, family=poisson(), control=control))
}

test_that("bootstrap native refits reproduce the fitted objective and scores", {
  fits <- bootstrapFixture()
  fit <- fits$alternative
  refit <- sommer:::.boot_mmes_refit(fit, fit$y, 1000L)
  expect_true(refit$converged)
  expect_equal(refit$logLik, fit$llik, tolerance=1e-6)
  expect_equal(unname(refit$parameters), unname(fit$optimizer$par), tolerance=1e-3)
  expect_true(is.finite(refit$scoreMax))
  expect_false(fit$laplaceDiagnostics$uncertaintyComputed)
  data <- fits$data
  data$y <- simulate(fit, nsim=1, seed=18, re.form="population")[[1L]]
  public <- mmes(y~x, random=~group, data=data, family=poisson(),
    control=fits$control)
  direct <- sommer:::.boot_mmes_refit(fit, data$y, 1000L)
  expect_equal(direct$logLik, public$llik, tolerance=1e-6)
  expect_equal(unname(direct$parameters), unname(public$optimizer$par),
    tolerance=1e-3)
})

test_that("population bootstrap is reproducible and leaves fitted models unchanged", {
  fit <- bootstrapFixture()$alternative
  fields <- c("optimizer", "b", "covStruct", "randomRelationships",
    "laplaceDiagnostics", "familySpecifications", "fixedCovariance", "Ci")
  before <- unserialize(serialize(fit[fields], NULL))
  set.seed(421)
  state <- .Random.seed
  first <- boot_mmes(fit, nsim=5, seed=39)
  expect_identical(.Random.seed, state)
  second <- boot_mmes(fit, nsim=5, seed=39)
  expect_identical(first$draws, second$draws)
  expect_identical(first$diagnostics, second$diagnostics)
  expect_identical(fit[fields], before)
  expect_equal(first$successful, 5)
  expect_equal(dim(first$draws), c(5L, length(fit$optimizer$par)))
  expect_equal(first$percentile,
    t(apply(first$draws, 2, quantile, probs=c(0.025, 0.975))))
  expect_true(all(first$diagnostics$rejectedEvaluations >= 0))
  expect_null(first$p.value.bounds)
  expect_true(is.na(first$p.value))
})

test_that("nested fixed-effect bootstrap uses null simulations and plus-one tests", {
  fits <- bootstrapFixture()
  result <- boot_mmes(fits$null, fits$alternative, nsim=5, seed=39)
  expect_equal(result$successful, 5)
  expect_equal(result$statistic, 2*(fits$alternative$llik-fits$null$llik))
  expected <- (1+sum(result$diagnostics$statistic >= result$statistic))/6
  expect_equal(result$p.value, expected)
  expect_equal(unname(result$p.value.bounds), rep(expected, 2))
  expect_equal(result$monteCarloSE, sqrt(expected*(1-expected)/5))
  expect_error(boot_mmes(fits$alternative, fits$null, nsim=2), "strictly nested")
  expect_error(boot_mmes(list(), nsim=2), "native Laplace")
  expect_error(boot_mmes(fits$null, nsim=0), "nsim")
  expect_error(boot_mmes(fits$null, seed=Inf), "seed")
  fits$null$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(boot_mmes(fits$null, nsim=2), "output exceeds")
})

test_that("bootstrap refits preserve grouped binomial trial counts", {
  set.seed(117)
  data <- data.frame(group=factor(rep(seq_len(35), each=8)), trials=rep(8, 280))
  probability <- plogis(-0.2+rnorm(35, sd=0.7)[data$group])
  data$successes <- rbinom(280, 8, probability)
  fit <- mmes(binm(successes, trials=trials)~1, random=~group,
    data=data, family=binomial(),
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace")),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(sommer:::.boot_mmes_refit(fit, fit$y, 1000L)$logLik,
    fit$llik, tolerance=1e-6)
  result <- boot_mmes(fit, nsim=3, seed=9)
  expect_equal(result$successful, 3)
})

test_that("failed bootstrap refits remain explicit and bound rather than erase tests", {
  fits <- bootstrapFixture()
  expect_warning(result <- boot_mmes(fits$null, fits$alternative,
    nsim=3, seed=39, maxit=1), "bootstrap replicates failed")
  expect_equal(nrow(result$diagnostics), 3)
  expect_true(any(nzchar(result$diagnostics$error)))
  expect_true(is.na(result$p.value))
  expect_equal(diff(result$p.value.bounds), (3-result$successful)/4,
    ignore_attr=TRUE)
  expect_false(fits$null$laplaceDiagnostics$uncertaintyComputed)
})

test_that("bootstrap supports ZIP zero formulas and mixed nuisance parameters", {
  set.seed(141)
  data <- data.frame(group=factor(rep(seq_len(40), each=8)), x=rnorm(320))
  latent <- rnorm(40, sd=0.5)[data$group]
  data$y <- rpois(320, exp(1+0.2*data$x+latent))
  data$y[rbinom(320, 1, plogis(-1+0.4*data$x)) == 1] <- 0
  other <- data
  other$y <- rnorm(320, 1+0.2*other$x+latent, sd=0.8)
  data$trait <- "count"
  other$trait <- "continuous"
  data <- rbind(data, other)
  data$trait <- factor(data$trait)
  fit <- mmes(y~trait-1+trait:x, random=~group, data=data,
    family=familym("trait", count=zipoisson(zero=~x), continuous=gaussian()),
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace")),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  refit <- sommer:::.boot_mmes_refit(fit, fit$y, 1000L)
  expect_equal(refit$logLik, fit$llik, tolerance=1e-6)
  expect_equal(unname(refit$parameters), unname(fit$optimizer$par), tolerance=1e-3)
  expect_true(all(c("family:count:(Intercept)", "family:count:x",
    "family:continuous:dispersion") %in% names(refit$parameters)))
  result <- boot_mmes(fit, nsim=3, seed=39)
  expect_equal(result$successful, 3)
  expect_equal(ncol(result$draws), length(fit$optimizer$par))
  expect_false(fit$laplaceDiagnostics$uncertaintyComputed)
})

test_that("population simulation does not alter stored relationship factors", {
  fit <- bootstrapFixture()$null
  relationships <- unserialize(serialize(fit$randomRelationships, NULL))
  simulate(fit, nsim=2, seed=12, re.form="population")
  expect_identical(fit$randomRelationships, relationships)
})
