diagnosticData <- function(){
  set.seed(1)
  n <- 60
  DT <- data.frame(y=rnorm(n), g=factor(rep(1:12, 5)), e=factor(rep(1:3, 20)),
    w=runif(n, .5, 2))
  DT$y <- DT$y + rep(rnorm(12), 5)
  DT
}

fitQuiet <- function(...) suppressMessages(suppressWarnings(mmes(...)))

bruteLeverage <- function(fit, R){
  W <- as.matrix(fit$W)
  Ci <- as.matrix(postPEV(fit, 2L)$Ci)
  Tm <- W %*% Ci %*% t(W)
  list(h=diag(Tm %*% solve(R)), v=diag(R - Tm))
}

test_that("fitted and residual panels no longer plot residuals against themselves", {
  DT <- diagnosticData()
  fit <- fitQuiet(y ~ e, random=~g, rcov=~units, data=DT)
  d <- .mmes_diagnostics(fit, stnd=FALSE)
  expect_equal(d$fitted, as.numeric(fitted(fit)))
  expect_equal(d$residuals, as.numeric(residuals(fit)))
  expect_lt(abs(cor(d$fitted, d$residuals)), 0.99)
  expect_error(plot(fit), NA)
})

test_that("leverage matches brute-force algebra for iid, direct and heterogeneous residuals", {
  DT <- diagnosticData()
  n <- nrow(DT)
  fit <- fitQuiet(y ~ e, random=~g, rcov=~units, data=DT)
  b <- bruteLeverage(fit, diag(as.numeric(fit$covPar[[2]]), n))
  d <- .mmes_diagnostics(fit)
  expect_equal(d$leverage, b$h, tolerance=1e-10)
  expect_equal(d$residuals, as.numeric(residuals(fit)) / sqrt(b$v), tolerance=1e-10)

  direct <- fitQuiet(y ~ e, random=~g, rcov=~units, data=DT,
    control=mmes_control(formulation="direct", output=list(coefficientMatrix=TRUE)))
  expect_equal(.mmes_diagnostics(direct)$leverage, d$leverage, tolerance=1e-8)

  het <- fitQuiet(y ~ e, random=~g, rcov=~vsm(dsm(e), ism(units)), data=DT)
  v <- as.numeric(het$covPar[[2]])
  Rh <- diag(c(v[1], v[1] * v[2:3])[as.integer(DT$e)], n)
  expect_equal(.mmes_diagnostics(het)$leverage, bruteLeverage(het, Rh)$h, tolerance=1e-10)
  expect_error(plot(het), NA)
})

test_that("leverage honours the weight matrix W", {
  DT <- diagnosticData()
  fit <- fitQuiet(y ~ e, random=~g, rcov=~units, data=DT, W=diag(1/DT$w))
  R <- as.numeric(fit$covPar[[2]]) * diag(DT$w)
  expect_equal(.mmes_diagnostics(fit)$leverage, bruteLeverage(fit, R)$h, tolerance=1e-10)
})

test_that("plot rejects unsupported fits clearly", {
  skip_if_not(exists(".mmes_leverage_linear", asNamespace("sommer")))
  fit <- fitQuiet(y ~ e, random=~g, rcov=~units, data=diagnosticData())
  fit$residualLayout <- NULL
  expect_error(.mmes_leverage_residual(fit), "residual layout")
})
