.direct_pql_control <- function(formulation="direct", warmStart=TRUE,
                                objective="reml", inverse=0L, estimateTheta=FALSE){
  sommer:::.mmes_control_internal(formulation=formulation, task=list(objective=objective,
    pql=list(maxit=40L, tol=1e-8, warmStart=warmStart, estimateTheta=estimateTheta)),
    convergence=list(maxIters=80L, likelihoodTol=1e-7, parameterTol=1e-7),
    output=list(inverse=inverse),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

.direct_pql_data <- function(){
  set.seed(204)
  data <- data.frame(group=factor(rep(seq_len(12L), each=8L)),
    visit=factor(rep(seq_len(8L), 12L)), x=rnorm(96L),
    exposure=runif(96L, 0.8, 1.3))
  eta <- 0.3 + 0.15*data$x + rep(rnorm(12L, sd=0.3), each=8L)
  data$count <- rpois(96L, data$exposure*exp(eta))
  data$trials <- rep(8:11, 24L)
  data$success <- rbinom(96L, data$trials, plogis(eta))
  data$positive <- rgamma(96L, shape=5, scale=exp(eta)/5)
  data
}

test_that("direct PQL retains the sparse predictor without a Henderson fit", {
  original <- sommer:::.mmes_fit
  formulations <- logical()
  local_mocked_bindings(.mmes_fit=function(...){
    args <- list(...)
    if(isTRUE(args$.internal$pqlInner)){
      formulations <<- c(formulations, args$henderson)
    }
    original(...)
  }, .package="sommer")
  data <- .direct_pql_data()
  fit <- sommer:::.mmes_dispatch(count~x+offset(log(exposure)),
    random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
    data=data, family=poisson(), control=.direct_pql_control())
  expect_s3_class(fit, "mmes.glmm")
  expect_identical(fit$engine, "direct")
  expect_true(length(formulations) > 1L)
  expect_false(any(formulations))
  expect_s4_class(fit$W, "sparseMatrix")
  expect_equal(dim(fit$W), c(96L, length(fit$bu)))
  expect_equal(fitted(fit, type="link"),
    fit$offset + as.numeric(fit$W %*% fit$bu), tolerance=1e-12)
  expect_equal(fitted(fit), exp(fitted(fit, type="link")), tolerance=1e-12)
  expect_equal(fit$linear.predictorsNoOffset,
    as.numeric(fit$W %*% fit$bu), tolerance=1e-12)
  expect_true(fit$pqlConverged)
  expect_true(fit$convergence)
  expect_identical(fit$executionPlan$solver$selected, "direct")
  expect_identical(fit$executionPlan$glmm$selected, "PQL")
  expect_identical(fit$executionPlan$glmm$selectedPath, "pql")
})

test_that("direct PQL without random effects reproduces independent GLMs", {
  data <- .direct_pql_data()
  cases <- list(
    list(count~x+offset(log(exposure)), poisson()),
    list(binm(success, trials=trials)~x, binomial()),
    list(binm(success, trials=trials)~x, binomial("probit")),
    list(binm(success, trials=trials)~x, binomial("cloglog")),
    list(positive~x, Gamma()),
    list(positive~x, Gamma("log")),
    list(positive~x, inverse.gaussian("log")),
    list(positive~x, gaussian("log")),
    list(count~x, quasipoisson()),
    list(count~x, MASS::negative.binomial(2)))
  for(case in cases){
    formula <- case[[1L]]
    family <- case[[2L]]
    fit <- sommer:::.mmes_dispatch(formula, data=data, family=family,
      control=.direct_pql_control())
    frame <- model.frame(formula, data)
    response <- model.response(frame)
    weights <- attr(response, "trials")
    if(is.null(weights)) weights <- rep(1, nrow(frame))
    reference <- glm(as.numeric(response)~x+offset(
      if(is.null(model.offset(frame))) rep(0, nrow(frame)) else model.offset(frame)),
      data=data, weights=weights, family=family,
      control=glm.control(epsilon=1e-12, maxit=100L))
    expect_true(fit$pqlConverged, info=family$family)
    expect_equal(as.numeric(fit$b), unname(coef(reference)), tolerance=2e-6)
    expect_equal(fitted(fit), unname(fitted(reference)), tolerance=2e-6)
    expect_equal(fit$deviance, unname(deviance(reference)), tolerance=2e-6)
    dispersion <- if(sommer:::.pql_fixed_dispersion(family)) 1 else
      summary(reference)$dispersion
    expect_equal(fit$dispersion, dispersion, tolerance=2e-6)
    expect_equal(residuals(fit, type="pearson"),
      unname(residuals(reference, type="pearson")), tolerance=2e-6)
    expect_equal(residuals(fit, type="working"),
      unname(residuals(reference, type="working")), tolerance=2e-6)
  }
})

test_that("direct and Henderson PQL agree with structured covariance and missing rows", {
  data <- .direct_pql_data()
  data$count[c(7L, 51L)] <- NA_real_
  data$x[29L] <- NA_real_
  data$group[33L] <- NA
  data$exposure[75L] <- NA_real_
  data$visit[88L] <- NA
  precision <- Matrix::bdiag(rep(list(Matrix::bandSparse(8L,
    k=c(-1L, 0L, 1L), diagonals=list(rep(-0.15, 7L),
      rep(1.5, 8L), rep(-0.15, 7L)))), 12L))
  relationship <- solve(0.3^abs(outer(seq_len(12L), seq_len(12L), "-")))
  dimnames(relationship) <- list(levels(data$group), levels(data$group))
  attr(relationship, "inverse") <- TRUE
  for(objective in c("ml", "reml")){
    fits <- lapply(c("direct", "henderson"), function(formulation)
      sommer:::.mmes_dispatch(count~x+offset(log(exposure)),
        random=~vsm(ism(group), Gu=relationship, sigma2=0.1, fixedSigma2=TRUE),
        rcov=~vsm(ar1m(visit, rho=0.2, fixed=TRUE), ism(group)),
        W=precision, weights=~group, data=data, family=poisson(),
        control=.direct_pql_control(formulation, objective=objective)))
    direct <- fits[[1L]]
    henderson <- fits[[2L]]
    expect_true(direct$pqlConverged)
    expect_identical(which(!direct$obsInfo$included), c(7L, 29L, 33L, 51L, 75L, 88L))
    expect_identical(direct$obsInfo$included, henderson$obsInfo$included)
    expect_equal(as.numeric(direct$b), as.numeric(henderson$b), tolerance=2e-6)
    expect_equal(as.numeric(direct$u), as.numeric(henderson$u), tolerance=2e-6)
    expect_equal(fitted(direct), fitted(henderson), tolerance=2e-6)
    expect_equal(direct$workingWeights, henderson$workingWeights, tolerance=2e-6)
    expect_equal(direct$workingResponse, henderson$workingResponse, tolerance=2e-6)
    expect_equal(direct$VarBeta,
      sommer:::.wald_fixed_vcov(henderson), tolerance=2e-6, ignore_attr=TRUE)
  }
})

test_that("direct PQL estimates tied covariance parameters and warm starts consistently", {
  data <- .direct_pql_data()
  data$other <- data$visit
  set.seed(809)
  data$count <- rpois(nrow(data), exp(0.5 + 0.1*data$x +
    rnorm(12L, sd=0.6)[data$group] + rnorm(8L, sd=0.4)[data$other]))
  fit <- function(formulation="direct", warmStart=TRUE)
    sommer:::.mmes_dispatch(count~x, random=~group+other, data=data, family=poisson(),
      vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)),
      control=.direct_pql_control(formulation, warmStart))
  warm <- fit()
  cold <- fit(warmStart=FALSE)
  reference <- fit("henderson")
  expect_true(warm$pqlConverged)
  expect_true(cold$pqlConverged)
  expect_equal(as.numeric(warm$covPar[[1L]]),
    2*as.numeric(warm$covPar[[2L]]), tolerance=1e-10)
  expect_equal(unname(warm$monitor[1L, ]), 2*unname(warm$monitor[2L, ]),
    tolerance=1e-10)
  expect_equal(fitted(warm), fitted(cold), tolerance=2e-5)
  expect_equal(fitted(warm), fitted(reference), tolerance=2e-5)
  expect_lt(max(abs(unlist(warm$covPar) - unlist(cold$covPar))), 2e-6)
  expect_equal(unlist(warm$covPar), unlist(reference$covPar), tolerance=2e-5)
  expect_true(all(warm$pqlMonitor[-1L, "warmStarted"] == 1L))
  expect_true(all(cold$pqlMonitor[, "warmStarted"] == 0L))
  expect_null(warm$.ldltCache)
  expect_null(warm$.cholmodCache)
})

test_that("direct PQL supports mixed scalar traits and trait-specific dispersion", {
  data <- .direct_pql_data()
  long <- rbind(transform(data, trait="count", y=count),
    transform(data, trait="positive", y=positive))
  long$trait <- factor(long$trait, levels=c("count", "positive"))
  fits <- lapply(c("direct", "henderson"), function(formulation)
    sommer:::.mmes_dispatch(y~trait+trait:x-1,
      random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
      rcov=~vsm(dsm(trait), ism(units)), data=long,
      family=familym("trait", count=poisson(), positive=Gamma("log")),
      control=.direct_pql_control(formulation)))
  expect_true(fits[[1L]]$pqlConverged)
  expect_equal(fitted(fits[[1L]]), fitted(fits[[2L]]), tolerance=2e-5)
  expect_equal(fits[[1L]]$dispersion, fits[[2L]]$dispersion, tolerance=2e-5)
  expect_equal(unname(fits[[1L]]$dispersion["count"]), 1)
  expect_gt(fits[[1L]]$dispersion["positive"], 0)
})

test_that("direct PQL estimates negative-binomial shape as well as GLM coefficients", {
  set.seed(3)
  data <- data.frame(x=rnorm(200L))
  data$y <- rnbinom(nrow(data), size=1.5, mu=exp(0.4+0.3*data$x))
  reference <- MASS::glm.nb(y~x, data=data)
  fit <- sommer:::.mmes_dispatch(y~x, data=data, family=MASS::negative.binomial(1),
    control=.direct_pql_control(estimateTheta=TRUE))
  expect_true(fit$pqlConverged)
  expect_equal(as.numeric(fit$b), unname(coef(reference)), tolerance=2e-6)
  expect_equal(fit$theta.nb, reference$theta, tolerance=2e-5)
  expect_equal(fit$theta.nb.se, reference$SE.theta, tolerance=2e-5)
  expect_equal(fitted(fit), unname(fitted(reference)), tolerance=2e-6)
})

test_that("direct PQL exposes working-model inference without likelihood claims", {
  data <- .direct_pql_data()
  fit <- sommer:::.mmes_dispatch(count~x, random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
    data=data, family=poisson(), control=.direct_pql_control(inverse=2L))
  standardErrors <- summary(fit)$betas$Std.Error
  expect_length(standardErrors, nrow(fit$b))
  expect_true(all(is.finite(standardErrors)))
  expect_length(unlist(fit$uPevList), nrow(fit$u))
  expect_true(all(is.finite(unlist(fit$uPevList))))
  expect_true(all(unlist(fit$uPevList) >= 0))
  expect_s3_class(wald_mmes(fit), "wald.mmes")
  for(type in c("response", "working", "pearson", "deviance"))
    expect_true(all(is.finite(residuals(fit, type=type))))
  expect_true(all(is.na(fit$llik)))
  expect_true(is.na(fit$AIC))
  expect_true(is.na(fit$BIC))
  expect_error(anova(fit, fit), "not available for PQL")
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(fitted(restored), fitted(fit))
  expect_equal(summary(restored), summary(fit))
})

test_that("direct PQL coefficients and PEV agree with independent working-model algebra", {
  data <- .direct_pql_data()
  fit <- sommer:::.mmes_dispatch(count~x+offset(log(exposure)),
    random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
    data=data, family=poisson(), control=.direct_pql_control(inverse=2L))
  X <- as.matrix(fit$W[, seq_len(nrow(fit$b)), drop=FALSE])
  Z <- as.matrix(fit$W[, nrow(fit$b) + seq_len(nrow(fit$u)), drop=FALSE])
  G <- diag(0.1, ncol(Z))
  V <- Z %*% G %*% t(Z) + diag(1/fit$workingWeights)
  Vi <- solve(V)
  VarBeta <- solve(crossprod(X, Vi %*% X))
  beta <- VarBeta %*% crossprod(X, Vi %*% as.numeric(fit$workingResponse))
  u <- G %*% crossprod(Z, Vi %*% (as.numeric(fit$workingResponse) - X %*% beta))
  P <- Vi - Vi %*% X %*% VarBeta %*% t(X) %*% Vi
  pev <- G - G %*% crossprod(Z, P %*% Z) %*% G
  expect_equal(as.numeric(fit$b), as.numeric(beta), tolerance=1e-9)
  expect_equal(as.numeric(fit$u), as.numeric(u), tolerance=1e-9)
  expect_equal(unname(fit$VarBeta), unname(VarBeta), tolerance=1e-9)
  expect_equal(as.numeric(fit$uPevList[[1L]]), diag(pev), tolerance=1e-9)
})

test_that("direct PQL does not certify failed inner or unfinished outer iterations", {
  data <- .direct_pql_data()
  innerControl <- .direct_pql_control()
  innerControl$convergence$maxIters <- 1L
  innerControl$task$pql$maxit <- 3L
  expect_warning(inner <- sommer:::.mmes_dispatch(count~x, random=~group, data=data,
    family=poisson(), control=innerControl), "Direct PQL convergence was not certified")
  expect_false(inner$convergence)
  expect_false(inner$pqlConverged)
  expect_true(all(inner$pqlMonitor[, "innerConverged"] == 0))
  expect_equal(nrow(inner$pqlMonitor), 3L)
  expect_identical(inner$engineDiagnostics$direct$stopReason, "iteration-limit")
  expect_false(summary(inner)$logo$Converge)

  outerControl <- .direct_pql_control()
  outerControl$task$pql$maxit <- 1L
  expect_warning(outer <- sommer:::.mmes_dispatch(count~x,
    random=~vsm(ism(group), sigma2=0.1, fixedSigma2=TRUE),
    data=data, family=poisson(), control=outerControl),
    "Direct PQL convergence was not certified")
  expect_true(all(outer$pqlMonitor[, "innerConverged"] == 1))
  expect_false(outer$pqlConverged)
  expect_false(outer$convergence)
  expect_equal(fitted(outer, type="link"),
    outer$offset + as.numeric(outer$W %*% outer$bu), tolerance=1e-12)
})

test_that("direct PQL validates weights and supports deferred inference", {
  data <- .direct_pql_data()
  control <- .direct_pql_control()
  expect_error(sommer:::.mmes_dispatch(count~x, data=data, family=poisson(),
    weights=~group, W=Matrix::Matrix(1, nrow(data), nrow(data), sparse=TRUE) +
      Matrix::Diagonal(nrow(data)), control=control),
    "nonzero entries across groups")
  expect_warning(expect_error(sommer:::.mmes_dispatch(count~x, data=data, family=poisson(),
    W=Matrix::Diagonal(nrow(data), x=-1), control=control),
    "positive definite"), "not positive definite")
  expect_error(sommer:::.mmes_dispatch(count~x, data=data, family=poisson(),
    control=sommer:::.mmes_control_internal(formulation="direct", output=list(inverse=1L))),
    "inverse")
  fit <- sommer:::.mmes_dispatch(count~x, data=data, family=poisson(), control=control)
  prediction <- predict(fit, D=matrix(c(1, 0), 1L), VarU=FALSE)
  expect_equal(as.numeric(prediction$vcov), fit$VarBeta[1L, 1L], tolerance=1e-8)
  expect_equal(as.matrix(postPEV(fit, mode=2L)$Ci), fit$VarBeta,
    tolerance=1e-8, ignore_attr=TRUE)
  expect_null(fit$C)
  expect_equal(dim(fit$Ci), c(0L, 0L))
})
