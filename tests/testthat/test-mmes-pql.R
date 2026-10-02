test_that("Gaussian identity fits retain the linear mixed-model path", {
  data(DT_example)
  DT <- DT_example[complete.cases(DT_example[, c("Yield", "Env", "Name")]), ]

  fitDefault <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                     nIters=10, verbose=FALSE)
  fitGaussian <- mmes(Yield~Env, random=~Name, rcov=~units, data=DT,
                      nIters=10, verbose=FALSE, family=gaussian())

  expect_s3_class(fitGaussian, "mmes")
  expect_false(inherits(fitGaussian, "mmes.glmm"))
  expect_equal(fitDefault$theta, fitGaussian$theta, tolerance=1e-12)
})

test_that("PQL binomial and Poisson fits expose response and link scales", {
  set.seed(10)
  group <- factor(rep(seq_len(24), each=5))
  x <- rep(c(0, 1, 0, 1, 0), 24)
  binomialData <- data.frame(
    y=rbinom(length(x), 1, plogis(-0.4 + 0.8 * x)), x=x, group=group
  )
  poissonData <- data.frame(
    y=rpois(length(x), exp(0.2 + 0.35 * x)), x=x, group=group
  )

  binomialFit <- mmes(y~x, random=~group, rcov=~units, data=binomialData,
                      family=binomial(), nIters=8, verbose=FALSE,
                      pqlControl=list(maxit=4))
  poissonFit <- mmes(y~x, random=~group, rcov=~units, data=poissonData,
                     family=poisson(), nIters=8, verbose=FALSE,
                     pqlControl=list(maxit=4))

  for(fit in list(binomialFit, poissonFit)){
    expect_s3_class(fit, "mmes.glmm")
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(fitted(fit, type="link"))))
    expect_true(all(is.finite(residuals(fit, type="deviance"))))
    expect_true(all(is.finite(residuals(fit, type="working"))))
    expect_lte(nrow(fit$pqlMonitor), 4L)
  }
  expect_true(all(fitted(binomialFit) > 0 & fitted(binomialFit) < 1))
  expect_true(all(fitted(poissonFit) > 0))
})

test_that("PQL Poisson fit without random effects matches glm", {
  set.seed(23)
  data <- data.frame(x=rep(c(0, 1), each=100))
  data$y <- rpois(nrow(data), exp(0.2 + 0.6 * data$x))

  glmFit <- glm(y~x, data=data, family=poisson())
  pqlFit <- mmes(y~x, rcov=~units, data=data, family=poisson(),
                 nIters=15, verbose=FALSE,
                 pqlControl=list(maxit=20, tol=1e-8))

  expect_equal(as.numeric(pqlFit$b), unname(coef(glmFit)), tolerance=1e-4)
})

test_that("PQL honors offsets on the link scale", {
  set.seed(41)
  data <- data.frame(x=rep(c(0, 1), each=100), exposure=runif(200, 0.4, 2))
  data$y <- rpois(nrow(data), data$exposure * exp(0.2 + 0.5 * data$x))

  glmFit <- glm(y~x+offset(log(exposure)), data=data, family=poisson())
  pqlFit <- mmes(y~x+offset(log(exposure)), rcov=~units, data=data,
                 family=poisson(), nIters=15, verbose=FALSE,
                 pqlControl=list(maxit=20, tol=1e-8))

  expect_equal(as.numeric(pqlFit$b), unname(coef(glmFit)), tolerance=1e-4)
  expect_equal(fitted(pqlFit, type="link"), pqlFit$offset +
               pqlFit$linear.predictorsNoOffset, tolerance=1e-12)
})

test_that("PQL covariance warm starts agree with converged cold starts", {
  set.seed(72)
  data <- data.frame(group=factor(rep(seq_len(20L), each=8L)), x=rnorm(160L))
  effects <- rnorm(20L, sd=0.4)
  data$y <- rpois(160L, exp(0.2 + 0.3 * data$x + effects[data$group]))
  fit <- function(warm){
    mmes(y~x, random=~group, rcov=~units, data=data, family=poisson(),
      nIters=80, tolParConvLL=1e-9, tolParConvNorm=1e-9, verbose=FALSE,
      dateWarning=FALSE, pqlControl=list(maxit=30, tol=1e-9, warmStart=warm))
  }
  warm <- fit(TRUE)
  cold <- fit(FALSE)
  expect_equal(unname(warm$pqlMonitor[1L, "warmStarted"]), 0)
  expect_true(all(warm$pqlMonitor[-1L, "warmStarted"] == 1))
  expect_true(all(cold$pqlMonitor[, "warmStarted"] == 0))
  expect_equal(unname(warm$pqlMonitor[1L, "symbolicAnalyses"]), 1)
  expect_true(all(warm$pqlMonitor[-1L, "symbolicAnalyses"] == 0))
  expect_true(all(cold$pqlMonitor[, "symbolicAnalyses"] >= 1))
  expect_null(warm$.ldltCache)
  expect_null(attr(warm$covStruct, "ldltCache"))
  expect_equal(as.numeric(warm$b), as.numeric(cold$b), tolerance=1e-4)
  expect_equal(fitted(warm), fitted(cold), tolerance=1e-4)
  expect_equal(unlist(warm$theta), unlist(cold$theta), tolerance=1e-4)
  expect_error(fit(NA), "warmStart")

  cholmodData <- data.frame(x=rep(c(0, 1), each=40L))
  cholmodData$y <- rpois(80L, exp(0.2 + 0.4 * cholmodData$x))
  fitCholmod <- function(warmStart){
    mmes(y~x, rcov=~units, data=cholmodData, family=poisson(),
      solver="cholmod", nIters=4, verbose=FALSE, dateWarning=FALSE,
      pqlControl=list(maxit=6, tol=1e-8, warmStart=warmStart))
  }
  cholmodWarm <- fitCholmod(TRUE)
  cholmodCold <- fitCholmod(FALSE)
  expect_equal(unname(cholmodWarm$pqlMonitor[1L, "symbolicAnalyses"]), 1)
  expect_true(all(cholmodWarm$pqlMonitor[-1L, "symbolicAnalyses"] == 0))
  expect_true(all(cholmodCold$pqlMonitor[, "symbolicAnalyses"] >= 1))
  expect_equal(as.numeric(cholmodWarm$b), as.numeric(cholmodCold$b), tolerance=1e-4)
  expect_equal(unlist(cholmodWarm$theta), unlist(cholmodCold$theta), tolerance=1e-4)
  expect_null(cholmodWarm$.cholmodCache)
  expect_null(attr(cholmodWarm$covStruct, "cholmodCache"))
})

test_that("PQL accepts sparse non-diagonal W after observation filtering", {
  set.seed(42)
  n <- 80
  precision <- Matrix::bandSparse(
    n, k=c(-1, 0, 1),
    diagonals=list(rep(-0.2, n - 1L), rep(1.5, n), rep(-0.2, n - 1L))
  )
  data <- data.frame(y=rpois(n, exp(0.1)), x=rep(c(0, 1), length.out=n))
  data$y[7] <- NA_real_

  fit <- mmes(y~x, rcov=~units, data=data, W=precision, family=poisson(),
              nIters=8, verbose=FALSE, pqlControl=list(maxit=4))

  expect_s3_class(fit, "mmes.glmm")
  expect_length(fitted(fit), n - 1L)
  expect_true(all(is.finite(fitted(fit))))
})
pqlGlmData <- function(){
  set.seed(1); n <- 300
  d <- data.frame(x=rnorm(n))
  eta <- 0.5 + 0.3 * d$x
  d$yg <- rgamma(n, shape=5, rate=5 / exp(eta))
  d$yig <- abs(rnorm(n, exp(eta), 0.3)) + 0.05
  d$tot <- sample(5:15, n, TRUE)
  d$succ <- rbinom(n, d$tot, plogis(eta - 0.5))
  d$fail <- d$tot - d$succ
  d$yb <- rbinom(n, 1, 0.4)
  d
}

pqlFit <- function(formula, family, data, ...){
  mmes(formula, rcov=~units, data=data, family=family, verbose=FALSE,
       nIters=40, tolParConvLL=1e-8, tolParConvNorm=1e-8,
       pqlControl=utils::modifyList(list(maxit=50, tol=1e-12), list(...)))
}

test_that("PQL without random effects reproduces glm for several families and links", {
  d <- pqlGlmData()
  cases <- list(
    list(yg~x, Gamma()), list(yg~x, Gamma("log")),
    list(yig~x, inverse.gaussian("log")), list(yb~x, binomial("probit")),
    list(yb~x, binomial("cloglog"))
  )
  for(case in cases){
    fit <- pqlFit(case[[1]], case[[2]], d)
    ref <- glm(case[[1]], data=d, family=case[[2]])
    expect_equal(as.numeric(fit$b), unname(coef(ref)), tolerance=1e-6)
    expect_equal(fit$dispersion, summary(ref)$dispersion, tolerance=1e-4)
    expect_equal(fit$deviance, deviance(ref), tolerance=1e-6)
  }
})

test_that("binm() carries binomial trials as prior weights", {
  d <- pqlGlmData()
  ref <- glm(cbind(succ, fail)~x, data=d, family=binomial())
  fit <- pqlFit(binm(succ, fail)~x, binomial(), d)
  expect_equal(as.numeric(fit$b), unname(coef(ref)), tolerance=1e-6)
  expect_equal(fit$deviance, deviance(ref), tolerance=1e-6)
  expect_equal(residuals(fit, type="pearson"),
               unname(residuals(ref, type="pearson")), tolerance=1e-5)

  d$p <- binm(d$succ, trials=d$tot)
  sub <- d[-(1:7), ]
  expect_equal(attr(sub$p, "trials"), sub$tot)
  fitSub <- pqlFit(p~x, binomial(), sub)
  refSub <- glm(cbind(succ, fail)~x, data=sub, family=binomial())
  expect_equal(as.numeric(fitSub$b), unname(coef(refSub)), tolerance=1e-6)

  expect_warning(expect_warning(pqlFit(I(succ/tot)~x, binomial(), d), "binm"),
                 "non-integer")
  expect_error(binm(3, 1, 4), "exactly one")
  expect_error(binm(5, trials=4), "successes")
})

test_that("negative binomial theta estimation matches MASS::glm.nb", {
  set.seed(3); n <- 400
  d <- data.frame(x=rnorm(n))
  d$y <- rnbinom(n, size=1.5, mu=exp(0.4 + 0.3 * d$x))
  ref <- MASS::glm.nb(y~x, data=d)
  fit <- pqlFit(y~x, MASS::negative.binomial(1), d, estimateTheta=TRUE, maxit=100,
                tol=1e-10)
  expect_equal(fit$theta.nb, ref$theta, tolerance=1e-5)
  expect_equal(fit$theta.nb.se, ref$SE.theta, tolerance=1e-4)
  expect_equal(as.numeric(fit$b), unname(coef(ref)), tolerance=1e-6)
  expect_equal(fit$dispersion, 1)
  expect_error(pqlFit(y~x, poisson(), d, estimateTheta=TRUE), "negative.binomial")
})

test_that("PQL fits do not report working-model likelihoods", {
  d <- pqlGlmData()
  d$g <- factor(rep(1:30, each=10))
  fit <- mmes(binm(succ, fail)~x, random=~g, data=d, family=binomial(),
              verbose=FALSE)
  expect_true(is.na(fit$AIC))
  expect_true(all(is.na(fit$llik)))
  expect_error(anova(fit, fit), "not available for PQL")
  expect_identical(summary(fit)$logo$Method, "PQL")
})
