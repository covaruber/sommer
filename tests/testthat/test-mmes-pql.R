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