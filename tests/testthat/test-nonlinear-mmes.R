nonlinearExample <- function(curve="exponential"){
  set.seed(738)
  group <- factor(rep(seq_len(20L), each=12L))
  x <- rep(seq(0, 6, length.out=12L), 20L)
  mean <- switch(curve, exponential=4*exp(-0.55*x),
    logistic=4*plogis(1.3*(x-3)), saturation=4*x/(1.5+x))
  data.frame(y=mean+rnorm(20L, sd=0.3)[group]+rnorm(length(x), sd=0.12),
    x=x, group=group)
}

nonlinearControl <- function(start=NULL, ...){
  sommer:::.mmes_control_internal(task=list(objective="ml", nonlinear=list(start=start)),
    convergence=list(maxIters=500L, likelihoodTol=1e-10, parameterTol=1e-8),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE), ...)
}

nonlinearDense <- function(parameters, fixed, random, y, offset, input, curve){
  p <- ncol(fixed)
  k <- if(curve == "logistic") 3L else 2L
  theta <- parameters[p+seq_len(k)]
  a <- exp(theta[1L])
  b <- exp(theta[2L])
  mean <- switch(curve, exponential=a*exp(-b*input),
    logistic=a*plogis(b*(input-theta[3L])), saturation=a*input/(b+input))
  residual <- y-offset-as.numeric(fixed %*% parameters[seq_len(p)])-mean
  covariance <- exp(parameters[p+k+1L])*tcrossprod(as.matrix(random)) +
    diag(exp(parameters[p+k+2L]), length(y))
  factor <- chol(covariance)
  alpha <- backsolve(factor, forwardsolve(t(factor), residual))
  -0.5*(length(y)*log(2*pi)+2*sum(log(diag(factor)))+sum(residual*alpha))
}

test_that("traditional curve contracts and analytic coordinate derivatives agree", {
  inputs <- c(0, 0.4, 1, 3, 6)
  starts <- list(exponential=c(asymptote=4, rate=0.55),
    logistic=c(asymptote=4, rate=1.3, midpoint=3),
    saturation=c(asymptote=4, half=1.5))
  for(curve in names(starts)){
    model <- nonlinear_curve("x", curve)
    theta <- sommer:::.nonlinear_coordinates(model, starts[[curve]])
    value <- sommer:::sommer_nonlinear_curve_eval(curve, inputs, theta)
    expected <- switch(curve, exponential=4*exp(-0.55*inputs),
      logistic=4*plogis(1.3*(inputs-3)), saturation=4*inputs/(1.5+inputs))
    expect_equal(as.numeric(value$mean), expected, tolerance=1e-12)
    for(step in c(1e-4, 1e-5)){
      numerical <- vapply(seq_along(theta), function(index){
        upper <- lower <- theta
        upper[index] <- upper[index]+step
        lower[index] <- lower[index]-step
        as.numeric((sommer:::sommer_nonlinear_curve_eval(curve, inputs, upper)$mean -
          sommer:::sommer_nonlinear_curve_eval(curve, inputs, lower)$mean)/(2*step))
      }, numeric(length(inputs)))
      expect_equal(value$jacobian, numerical, tolerance=1e-7)
    }
    expect_identical(unserialize(serialize(model, NULL)), model)
  }
  expect_error(nonlinear_curve("x", "relu"), "arg")
  expect_error(nonlinear_curve(c("x", "z")), "one numeric")
  model <- nonlinear_curve("x")
  model$parameters <- rev(model$parameters)
  expect_error(sommer:::.nonlinear_validate(model), "modified")
  expect_error(sommer:::sommer_nonlinear_curve_eval("exponential", -1, c(0, 0)),
    "non-negative")
  expect_error(sommer:::sommer_nonlinear_curve_eval("logistic", 1, c(0, 0)),
    "incompatible")
  extreme <- sommer:::sommer_nonlinear_curve_eval("logistic", c(-1000, 1000),
    c(log(4), 0, 0))
  expect_equal(as.numeric(extreme$mean), c(0, 4))
  expect_true(all(is.finite(extreme$jacobian)))
})

test_that("native Gaussian likelihood and all scores match an independent dense oracle", {
  data <- nonlinearExample()[seq_len(36L), ]
  data$z <- rep(c(-1, 0, 1), 12L)
  control <- nonlinearControl()
  fixed <- y~0+z
  for(curve in c("exponential", "logistic", "saturation")){
    prepared <- sommer:::.nonlinear_prepare(fixed, ~group, data, control,
      "exclude", "exclude", "exclude", "exclude", NULL, data$x)
    x <- prepared$X
    z <- do.call(cbind, prepared$Z)
    offset <- seq(-0.1, 0.1, length.out=nrow(data))
    theta <- if(curve == "logistic") c(log(4), log(1.3), 3) else
      c(log(4), log(0.55))
    parameters <- c(0.2, theta, log(0.3^2), log(0.12^2))
    context <- sommer:::sommer_nonlinear_gaussian_create(x, z, data$y,
      offset, data$x, curve, prepared$covStruct[1L], prepared$Ai)
    value <- sommer:::sommer_nonlinear_gaussian_eval(context, parameters)
    oracle <- function(par) nonlinearDense(par, x, z, data$y, offset, data$x, curve)
    expect_true(value$valid)
    expect_equal(value$logLik, oracle(parameters), tolerance=1e-10)
    expect_lt(value$modeScore, 1e-10)
    for(step in c(1e-4, 1e-5)){
      numerical <- vapply(seq_along(parameters), function(index){
        upper <- lower <- parameters
        upper[index] <- upper[index]+step
        lower[index] <- lower[index]-step
        (oracle(upper)-oracle(lower))/(2*step)
      }, numeric(1))
      expect_equal(as.numeric(value$gradient), numerical, tolerance=1e-6)
    }
    expect_error(sommer:::sommer_nonlinear_gaussian_eval(context, parameters[-1L]),
      "incompatible length")
    invalid <- parameters
    invalid[length(invalid)] <- 1000
    rejected <- sommer:::sommer_nonlinear_gaussian_eval(context, invalid)
    expect_false(rejected$valid)
    expect_match(rejected$reason, "residual variance")
  }
})

test_that("Gaussian curve fitting recovers identifiable parameters without inference", {
  targets <- list(exponential=c(asymptote=4, rate=0.55),
    logistic=c(asymptote=4, rate=1.3, midpoint=3),
    saturation=c(asymptote=4, half=1.5))
  for(curve in names(targets)){
    data <- nonlinearExample(curve)
    seed <- .Random.seed
    fit <- sommer:::.mmes_dispatch(y~0, random=~group, data=data,
      nonlinear=nonlinear_curve("x", curve), control=nonlinearControl())
    expect_identical(.Random.seed, seed)
    expect_s3_class(fit, "mmes.nonlinear")
    expect_true(fit$convergence)
    expect_true(fit$nonlinearDiagnostics$identified)
    expect_equal(unname(fit$nonlinearParameters), unname(targets[[curve]]),
      tolerance=0.08)
    expect_lt(max(abs(fit$nonlinearDiagnostics$gradient)), 0.02)
    expect_equal(unname(fit$dispersion), 0.12^2, tolerance=0.15)
    expect_equal(fitted(fit), predict(fit))
    expect_equal(predict(fit, newdata=data), fitted(fit), tolerance=1e-10)
    expect_equal(predict(fit, newdata=data, re.form=NA), fit$fixedMean, tolerance=1e-10)
    expect_equal(predict(fit, re.form="population"), predict(fit, re.form=NA))
    expect_equal(residuals(fit), data$y-fitted(fit))
    expect_equal(as.numeric(logLik(fit)), fit$llik)
    expect_equal(nobs(fit), nrow(data))
    expect_equal(unname(coef(fit)), unname(fit$nonlinearParameters))
    expect_output(print(fit), "exact Gaussian ML")
    expect_null(fit$Ci)
    expect_null(fit$uPevList)
    expect_null(fit$fixedCovariance)
    expect_null(fit$optimizer$hessian)
    expect_false(summary(fit)$diagnostics$uncertaintyComputed)
    expect_identical(fit$executionPlan$nonlinear$selected, "exact Gaussian ML")
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(predict(restored, newdata=data[1:3, ]), predict(fit)[1:3])
    expect_false(any(vapply(fit, typeof, character(1)) == "externalptr"))
    refitted <- eval(fit$call)
    expect_equal(coef(refitted), coef(fit), tolerance=1e-10)
    expect_null(fit$nonlinearInference)
    expect_error(postVarU(fit), "post-fit")
    expect_equal(AIC(fit), -2*fit$llik+2*fit$npar)
    expect_equal(BIC(fit), -2*fit$llik+log(fit$nobs)*fit$npar)
    expect_error(anova(fit), "explicit parm and null")
    expect_error(wald_mmes(fit), "post-fit")
    expect_error(r2(fit), "post-fit")
    expect_error(plot(fit), "post-fit")
    expect_equal(unname(covmatrix_mmes(fit, term=1, se=FALSE)$covariance),
      matrix(fit$covPar[[1L]], 1, 1))
    expect_equal(nrow(simulate(fit,seed=1)),nobs(fit))
  }
})

test_that("traditional mean agrees with Gaussian ML for equivalent centered responses", {
  data <- nonlinearExample()
  fit <- sommer:::.mmes_dispatch(y~0, random=~group, data=data,
    nonlinear=nonlinear_curve("x"), control=nonlinearControl())
  data$mean <- predict(fit, re.form=NA)
  data$centered <- data$y-data$mean
  reference <- sommer:::.mmes_dispatch(centered~0, random=~group, data=data,
    control=sommer:::.mmes_control_internal(task=list(objective="ml"),
      convergence=list(maxIters=100L, likelihoodTol=1e-9, parameterTol=1e-8),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(fit$covPar[[1L]], as.numeric(reference$covPar[[1L]]), tolerance=1e-4)
  expect_equal(unname(fit$dispersion), as.numeric(reference$covPar[[2L]]),
    tolerance=1e-4)
  normalized <- sommer:::.mmes_dispatch(centered~0, random=~group, data=data,
    control=sommer:::.mmes_control_internal(task=list(objective="ml",
      glmm=list(method="laplace", maxit=500L, modeTol=1e-10)),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_equal(fit$llik, normalized$llik, tolerance=1e-7)
  expect_equal(fitted(fit)-data$mean, fitted(normalized), tolerance=1e-5)
})

test_that("curve ML agrees with an independent optimizer", {
  data <- droplevels(nonlinearExample()[1:120, ])
  fit <- sommer:::.mmes_dispatch(y~0, random=~group, data=data,
    nonlinear=nonlinear_curve("x"), control=nonlinearControl())
  oracle <- function(parameters) -nonlinearDense(parameters, fit$fixedDesign,
    fit$randomDesign, data$y, numeric(nrow(data)), data$x, "exponential")
  reference <- nlminb(c(log(3.5), log(0.7), log(0.08), log(0.02)), oracle,
    control=list(iter.max=500L, rel.tol=1e-10, x.tol=1e-8))
  expect_identical(reference$convergence, 0L)
  expect_equal(-reference$objective, fit$llik, tolerance=1e-8)
  expect_equal(unname(reference$par), unname(fit$optimizer$par), tolerance=1e-5)
})

test_that("relationships, covariance shapes and constraints retain their semantics", {
  data <- nonlinearExample()
  data$occasion <- factor(rep(rep(1:3, each=4), 20))
  relationship <- 0.7*diag(20)+0.3
  dimnames(relationship) <- list(levels(data$group), levels(data$group))
  precision <- solve(relationship)
  attr(precision, "inverse") <- TRUE
  # Use an interior AR(1) covariance rather than a random-intercept boundary.
  occasionCovariance <- 0.09*0.4^abs(outer(1:3, 1:3, "-"))
  effects <- as.numeric(t(chol(kronecker(occasionCovariance, relationship))) %*%
    rnorm(60L))
  data$y <- 4*exp(-0.55*data$x)+
    effects[as.integer(data$group)+20L*(as.integer(data$occasion)-1L)]+
    rnorm(nrow(data), sd=0.12)
  fit <- sommer:::.mmes_dispatch(y~0,
    random=~vsm(ar1m(occasion), ism(group), Gu=precision),
    data=data, nonlinear=nonlinear_curve("x"), control=nonlinearControl())
  expect_true(fit$convergence)
  expect_lt(abs(fit$covPar[[1L]][2L]), 0.95)
  expect_equal(predict(fit, newdata=data[1:3, ]), fitted(fit)[1:3], tolerance=1e-10)
  proposal <- c(log(4), log(0.55), log(0.09), atanh(0.4), log(0.0144))
  context <- sommer:::sommer_nonlinear_gaussian_create(fit$fixedDesign,
    fit$randomDesign, fit$y, fit$offset, fit$input, "exponential",
    fit$covStruct, fit$randomRelationships)
  value <- sommer:::sommer_nonlinear_gaussian_eval(context, proposal)
  independent <- function(parameters){
    descriptor <- sommer:::evaluate_covstruct_cpp(fit$covStruct[[1L]], parameters[3:4])
    covariance <- fit$randomDesign %*% kronecker(descriptor, relationship) %*%
      t(fit$randomDesign)+diag(exp(parameters[5L]), nrow(data))
    residual <- data$y-exp(parameters[1L])*exp(-exp(parameters[2L])*data$x)
    factor <- chol(as.matrix(covariance))
    -0.5*(nrow(data)*log(2*pi)+2*sum(log(diag(factor)))+
      sum(forwardsolve(t(factor), residual)^2))
  }
  expect_equal(value$logLik, independent(proposal), tolerance=1e-9)
  expect_equal(fit$llik, independent(fit$optimizer$par), tolerance=1e-8)
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-5
    lower[index] <- lower[index]-1e-5
    expect_equal(as.numeric(value$gradient)[index],
      (independent(upper)-independent(lower))/2e-5, tolerance=1e-5)
  }
  fixed <- sommer:::.mmes_dispatch(y~0, random=~vsm(ism(group), sigma2=0.09, fixedSigma2=TRUE),
    data=data, nonlinear=nonlinear_curve("x"), control=nonlinearControl())
  expect_length(fixed$optimizer$par, 3L)
  expect_equal(fixed$covPar[[1L]], 0.09)
  savedPrediction <- predict(fit, newdata=data[1:3, ])
  precision <- NULL
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(predict(restored, newdata=data[1:3, ]), savedPrediction)
  data$block <- factor(rep(rep(1:4, 3), 20))
  tied <- sommer:::.mmes_dispatch(y~0, random=~group+block, data=data,
    nonlinear=nonlinear_curve("x"),
    vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)),
    control=nonlinearControl())
  expect_true(tied$convergence)
  expect_equal(tied$covPar[[1L]], 2*tied$covPar[[2L]], tolerance=1e-12)
  expect_length(tied$optimizer$par, 4L)
  expect_equal(predict(tied, newdata=data[1:3, ]), fitted(tied)[1:3], tolerance=1e-10)
})

test_that("preparation filters nonlinear inputs and preserves offsets and fixed terms", {
  data <- nonlinearExample()
  data$z <- rep(c(-1, 1), length.out=nrow(data))
  data$o <- seq(0, 0.4, length.out=nrow(data))
  data$y <- data$y+0.7*data$z+data$o
  data$x[2L] <- NA
  data$y[3L] <- NA
  fit <- sommer:::.mmes_dispatch(y~0+z+offset(o), random=~group, data=data,
    nonlinear=nonlinear_curve("x"), control=nonlinearControl())
  expect_equal(which(!fit$obsInfo$included), 2:3)
  expect_identical(fit$obsInfo$reason[2L], "missing fixed covariate")
  expect_equal(unname(fit$b), 0.7, tolerance=0.04)
  expect_equal(predict(fit, newdata=data[-c(2, 3), ]), fitted(fit), tolerance=1e-10)
  expect_error(sommer:::.mmes_dispatch(y~0, random=~group, data=data, nonlinear=nonlinear_curve("x"),
    control=nonlinearControl(), naMethodX="fail"), "Missing values")
  expect_error(sommer:::.mmes_dispatch(y~0, random=~group, data=data, nonlinear=nonlinear_curve("x"),
    control=nonlinearControl(), naMethodX="pass"), "finite retained")
  new <- data[1L, ]
  new$group <- "unseen"
  expect_error(predict(fit, newdata=new), "unseen random-effect level")
  expect_length(predict(fit, newdata=new, re.form=NA), 1L)
  expect_error(predict(fit, newdata=data.frame(z=1)), "input")
  expect_error(predict(fit, newdata=data.frame(x=1, group="1", o=0)),
    "fixed-effect columns")
  expect_error(predict(fit, newdata=data[1L, ], nonsense=1), "Unused")
})

test_that("ordinary intercepts and factor preprocessing survive new-data prediction", {
  data <- nonlinearExample()
  data$treatment <- factor(rep(c("a", "b"), length.out=nrow(data)))
  data$y <- data$y+0.9+0.4*(data$treatment == "b")
  fit <- sommer:::.mmes_dispatch(y~treatment, random=~group, data=data,
    nonlinear=nonlinear_curve("x"), control=nonlinearControl(
      c(asymptote=4, rate=0.55)))
  expect_true(fit$convergence)
  expect_true("Intercept" %in% names(fit$b))
  expect_identical(fit$executionPlan$nonlinear$initialization, "supplied")
  expect_true(all(is.finite(fit$nonlinearDiagnostics$initialParameters)))
  expect_equal(predict(fit, newdata=data), fitted(fit), tolerance=1e-10)
  newdata <- data[c(1, 3), ]
  newdata$treatment <- as.character(newdata$treatment)
  expect_equal(predict(fit, newdata=newdata), fitted(fit)[c(1, 3)], tolerance=1e-10)
  newdata$treatment <- "unseen"
  expect_error(predict(fit, newdata=newdata), "new level")
})

test_that("nonlinear routes reject unsupported models and controls before fitting", {
  data <- nonlinearExample()
  model <- nonlinear_curve("x")
  fit <- function(control=nonlinearControl(), ...)
    sommer:::.mmes_dispatch(y~0, random=~group, data=data, nonlinear=model, control=control, ...)
  expect_error(fit(sommer:::.mmes_control_internal()), "objective='ml'")
  expect_error(fit(nonlinearControl(output=list(inverse=1L))), "post-fit")
  expect_error(fit(nonlinearControl(output=list(setup=TRUE))), "output\\$setup")
  expect_error(fit(family=poisson()), "integer")
  expect_error(fit(rcov=~vsm(ism(units))), "homogeneous")
  expect_error(fit(W=diag(nrow(data))), "W, weights")
  expect_error(fit(covPar=list()), "covPar")
  expect_error(fit(nonlinearControl(execution=list(seed=1))), "execution\\$seed")
  expect_error(fit(nonlinearControl(execution=list(denseMemoryMB=0.001))),
    "workspace")
  expect_error(fit(nonlinearControl(list(asymptote=4, rate=0.5))), "numeric vector")
  expect_error(fit(nonlinearControl(c(asymptote=-4, rate=0.5))), "positive")
  expect_error(fit(nonlinearControl(c(asymptote=4, typo=0.5))), "containing exactly")
  expect_error(sommer:::.mmes_dispatch(y~0, data=data, nonlinear=model, control=nonlinearControl()),
    "random term")
  expect_error(sommer:::.mmes_control_internal(task=list(nonlinear=list(gradient="numeric"))), "analytic")
  expect_error(sommer:::.mmes_dispatch(y~1, data=data, control=nonlinearControl(
    c(asymptote=4, rate=0.5))), "require a nonlinear")
  constant <- data
  constant$x <- 1
  expect_error(sommer:::.mmes_dispatch(y~0, random=~group, data=constant,
    nonlinear=model, control=nonlinearControl()), "distinct inputs")
  exhausted <- sommer:::.mmes_control_internal(task=list(objective="ml"),
    convergence=list(maxIters=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  expect_warning(notConverged <- fit(exhausted), "did not converge")
  expect_false(notConverged$convergence)
  expect_identical(notConverged$executionPlan$status, "not converged")
  descriptor <- with(data, vsm(ar1m(group), ism(group)))$covStruct
  descriptor$factors[[1L]]$derivative$backend <- "numeric"
  expect_error(sommer:::.nonlinear_covariance_eligibility(list(descriptor), "group"),
    "native analytic covariance")
})
