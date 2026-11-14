.laplace_vcc_control <- function(setup=FALSE, parameters=list()){
  sommer:::.mmes_control_internal(task=list(objective="ml", glmm=list(method="laplace",
    maxit=1000L, modeTol=1e-10, familyParameters=parameters)),
    output=list(setup=setup), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

.laplace_vcc_data <- function(){
  set.seed(193)
  data <- expand.grid(group=factor(1:10), block=factor(1:6), rep=1:2)
  data$occasion <- factor(data$rep)
  data$x <- rnorm(nrow(data))
  eta <- 0.6 + 0.2*data$x + rnorm(10, sd=0.65)[data$group] +
    rnorm(6, sd=0.45)[data$block]
  data$y <- rpois(nrow(data), exp(eta))
  data$continuous <- eta + rnorm(nrow(data), sd=0.5)
  data$ordered <- ordered(cut(eta+rlogis(nrow(data)),
    breaks=c(-Inf, -0.5, 1, Inf), labels=c("low", "mid", "high")))
  data
}

.laplace_vcc_evaluator <- function(fit){
  context <- sommer:::sommer_laplace_objective_create(fit$fixedDesign,
    fit$randomDesign, fit$laplaceNativeResponse, fit$offset, "",
    fit$covStruct, fit$randomRelationships, 1000,
    familySpecs=fit$familySpecifications, familyIndex=fit$laplaceFamilyIndex,
    trials=fit$priorWeights)
  function(parameters, gradient=TRUE)
    sommer:::.laplace_constrained_evaluate(parameters, fit$laplaceParameterLayout,
      function(expanded) sommer:::sommer_laplace_objective_eval(context, expanded,
        maxModeIterations=100L, tolerance=1e-10, gradient=gradient))
}

test_that("the Laplace tangent handles signed loadings, ratios and correlation ties", {
  structure <- function(variance, values, kinds, bounds=rep(NA_real_,length(values))){
    list(par=c(log(variance),values), par_names=c("sigma2",names(values)),
      free=rep(TRUE,length(values)+1L), factors=list(list(
        par_start=2L, par_end=length(values)+1L,
        report=list(backend="builtin",transform=kinds,lower=bounds,
          upper=ifelse(is.na(bounds),NA_real_,1)))))
  }
  structures <- list(
    structure(2,c(loading=.3,rho=atanh(.2),ratio=log(1.4)),
      c("identity","tanh","exp")),
    structure(3,c(loading=-.1,rho=atanh(.4),ratio=log(.7)),
      c("identity","tanh","exp")),
    list(par=log(.5),par_names="sigma2",free=TRUE,factors=list()))
  constraints <- data.frame(parameter=c(1,5,2,6,3,7,4,8),
    group=rep(1:4,each=2),scale=c(2,1,-2,1,1,1,1,3))
  mapping <- sommer:::.laplace_covariance_map(structures,c("a","b","c"),constraints)
  expect_equal(mapping$table$start[1],2*mapping$table$start[5])
  expect_equal(mapping$table$start[2],-2*mapping$table$start[6])
  expect_equal(mapping$table$start[3],mapping$table$start[7])
  expect_equal(mapping$table$start[8],3*mapping$table$start[4])
  expect_equal(which(mapping$map$T[,1]!=0),9L)
  layout <- list(fixedCount=2L,covarianceCount=5L,covarianceMap=mapping$map,
    freePositions=lapply(mapping$structures,function(x) which(x$free)),
    parameterStarts=c(0L,4L,8L))
  parameters <- c(.1,.2,mapping$map$initial,.3,.4)
  objective <- function(value) sommer:::.laplace_constrained_evaluate(value,layout,
    function(native) list(valid=TRUE,logLik=-sum(native^2),gradient=-2*native))
  evaluation <- objective(parameters)
  expect_length(evaluation$gradient,length(parameters))
  for(index in seq_along(parameters)){
    upper <- lower <- parameters
    upper[index] <- upper[index]+1e-5
    lower[index] <- lower[index]-1e-5
    expect_equal(unname(evaluation$gradient[index]),
      (objective(upper)$logLik-objective(lower)$logLik)/2e-5,tolerance=1e-8)
  }
  proposal <- sommer:::.laplace_proposal_structures(mapping$structures,parameters,layout)
  expect_equal(unlist(lapply(proposal,`[[`,"par"),use.names=FALSE),
    unname(sommer:::.laplace_expand_parameters(parameters,layout)[3:11]))
  expect_identical(sommer:::.laplace_expand_parameters(parameters,NULL),parameters)
  expect_error(sommer:::.laplace_covariance_map(structures,c("a","b","c"),
    data.frame(parameter=c(3,7),group=1,scale=c(1,2))),"only be constrained to be equal")
  bounded <- list(structure(1,c(rho=0),"bounded_logit",-.5),
    structure(1,c(rho=0),"bounded_logit",-.25))
  expect_error(sommer:::.laplace_covariance_map(bounded,c("a","b"),
    data.frame(parameter=c(2,4),group=1)),"identical bounds")
})

test_that("native Laplace vcc equals independent combined-design reformulations", {
  skip_on_cran()
  data <- .laplace_vcc_data()
  for(ratio in c(1, 2)){
    constraint <- data.frame(parameter=1:2, group=1L, scale=c(ratio, 1))
    fit <- sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
      vcc=constraint, control=.laplace_vcc_control())
    combined <- cbind(sqrt(ratio)*model.matrix(~group-1, data),
      model.matrix(~block-1, data))
    reference <- sommer:::.mmes_dispatch(y~x, random=~vsm(ism(combined)), family=poisson(),
      data=data, control=.laplace_vcc_control())
    expect_true(fit$convergence)
    expect_equal(fit$covPar[[1L]], ratio*fit$covPar[[2L]], tolerance=1e-12)
    expect_equal(fit$llik, reference$llik, tolerance=1e-8)
    expect_equal(fitted(fit), fitted(reference), tolerance=2e-6)
    expect_equal(as.numeric(fit$b), as.numeric(reference$b), tolerance=2e-6)
    expect_equal(fit$npar, 3L)
    expect_equal(attr(logLik(fit), "df"), 3L)
    expect_equal(AIC(fit), -2*fit$llik+6)
    expect_equal(fit$vcParams$group, c(1L, 1L))
    expect_null(fit$fixedCovariance)
    expect_true(all(is.na(fit$covParNative$StdError)))
    expect_false("hessian" %in% names(fit$optimizer))
    uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
    referenceUncertainty <- sommer:::.mmes_laplace_uncertainty(reference)
    expect_true(uncertainty$valid)
    expect_equal(uncertainty$fixedCovariance, referenceUncertainty$fixedCovariance,
      tolerance=3e-4)
    expect_equal(uncertainty$covParNative$StdError[1L],
      ratio*uncertainty$covParNative$StdError[2L], tolerance=1e-10)
    expect_equal(uncertainty$covParCovariance[1L,2L]^2,
      prod(diag(uncertainty$covParCovariance)), tolerance=1e-10)
    expect_equal(uncertainty$covParNative$StdError[2L],
      referenceUncertainty$covParNative$StdError[1L], tolerance=3e-4)
    public <- covparams_mmes_se(fit)
    expect_equal(public$StdError,uncertainty$covParNative$StdError)
    expect_equal(attr(public,"vcov"),uncertainty$covParCovariance)
    expect_equal(summary(fit)$betas$Std.Error,
      unname(sqrt(diag(uncertainty$fixedCovariance))))
    evaluation <- .laplace_vcc_evaluator(fit)
    parameters <- fit$optimizer$par + c(0.1, -0.05, 0.15)
    value <- evaluation(parameters)
    expect_true(value$valid)
    for(index in seq_along(parameters)){
      upper <- lower <- parameters
      upper[index] <- upper[index]+1e-4
      lower[index] <- lower[index]-1e-4
      difference <- (evaluation(upper, FALSE)$logLik -
        evaluation(lower, FALSE)$logLik)/2e-4
      expect_equal(as.numeric(value$gradient)[index], difference, tolerance=2e-5)
    }
    covariance <- as.matrix(postPEV(fit, 2L)$Ci)
    referenceCovariance <- as.matrix(postPEV(reference, 2L)$Ci)
    map <- diag(nrow(covariance))
    map[3:12, 3:12] <- diag(sqrt(ratio), 10L)
    expect_equal(unname(covariance), map %*% referenceCovariance %*% t(map),
      tolerance=3e-5)
    expect_equal(postPEV(fit)$uPevList, postPEV(fit, 2L)$uPevList,
      tolerance=1e-8)
    expect_true(all(is.finite(predict(fit, newdata=data[1:3,],
      se.fit=TRUE, cov.fit=TRUE)$se.fit)))
    expect_true(all(is.finite(predict(fit, newdata=data[1:3,],
      re.form="population", type="response", se.fit=TRUE)$se.fit)))
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(vcov(restored), vcov(fit))
    expect_equal(fitted(restored), fitted(fit))
  }
})

test_that("Laplace constraints expose only random covariance parameters in setup", {
  data <- .laplace_vcc_data()
  local_mocked_bindings(sommer_laplace_objective_create=function(...)
    stop("Native likelihood reached during setup"), .package="sommer")
  setup <- sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
    control=.laplace_vcc_control(setup=TRUE))
  expect_equal(nrow(setup$vcParams), 2L)
  expect_false(any(grepl("units", setup$vcParams$term)))
  names <- paste(setup$vcParams$term, setup$vcParams$parameter, sep=":")
  constrained <- sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
    vcc=data.frame(parameter=names, group="shared", scale=c(2,1)),
    control=.laplace_vcc_control(setup=TRUE))
  expect_equal(constrained$vcParams$start[1], 2*constrained$vcParams$start[2])
  expect_equal(ncol(constrained$covarianceMap$T), 1L)
})

test_that("structured variance and correlation ties match exact Gaussian marginal likelihood", {
  set.seed(978)
  data <- expand.grid(occasion=factor(1:3),group=factor(1:12),
    block=factor(1:6),rep=1:2)
  data$x <- rnorm(nrow(data))
  correlation <- .4^abs(outer(1:3,1:3,`-`))
  group <- matrix(rnorm(36),12,3) %*% chol(.8*correlation)
  block <- matrix(rnorm(18),6,3) %*% chol(.4*correlation)
  data$y <- .5+.2*data$x +
    group[cbind(as.integer(data$group),as.integer(data$occasion))] +
    block[cbind(as.integer(data$block),as.integer(data$occasion))] +
    rnorm(nrow(data),sd=sqrt(.3))
  fit <- sommer:::.mmes_dispatch(y~x, random=~vsm(ar1m(occasion,rho=.2),ism(group))+
      vsm(ar1m(occasion,rho=.2),ism(block)),
    family=gaussian(),data=data,
    vcc=data.frame(parameter=c(1,3,2,4),group=c(1,1,2,2),scale=c(2,1,1,1)),
    control=.laplace_vcc_control(parameters=list(dispersion=.3,estimate=FALSE)))
  expect_true(fit$convergence)
  expect_equal(fit$npar,4L)
  expect_equal(fit$covParNative$estimate[1],2*fit$covParNative$estimate[3],
    tolerance=1e-12)
  expect_equal(fit$covParNative$estimate[2],fit$covParNative$estimate[4],
    tolerance=1e-12)
  evaluation <- .laplace_vcc_evaluator(fit)
  independent <- function(parameters){
    native <- sommer:::.laplace_expand_parameters(parameters,fit$laplaceParameterLayout)
    covariance <- .3*diag(nrow(data))
    for(term in 1:2){
      variance <- exp(native[2+2*term-1])
      rho <- tanh(native[2+2*term])
      # Descriptor Kronecker order is occasion first, then independent group.
      count <- if(term==1) 12L else 6L
      G <- variance*kronecker(rho^abs(outer(1:3,1:3,`-`)),diag(count))
      columns <- if(term==1) seq_len(36) else 36+seq_len(18)
      design <- as.matrix(fit$randomDesign[,columns])
      covariance <- covariance + design %*% G %*% t(design)
    }
    factor <- chol(covariance)
    residual <- data$y-as.numeric(fit$fixedDesign %*% parameters[1:2])
    standardized <- forwardsolve(t(factor),residual)
    -.5*(length(residual)*log(2*pi)+2*sum(log(diag(factor)))+sum(standardized^2))
  }
  proposal <- fit$optimizer$par+c(.1,-.05,.08,-.1)
  expect_equal(evaluation(proposal)$logLik,independent(proposal),tolerance=1e-9)
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    expect_equal(as.numeric(evaluation(proposal)$gradient)[index],
      (independent(upper)-independent(lower))/2e-4,tolerance=2e-5)
  }
  uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
  expect_true(uncertainty$valid)
  expect_equal(uncertainty$covParNative$StdError[1],
    2*uncertainty$covParNative$StdError[3],tolerance=1e-10)
  expect_equal(uncertainty$covParNative$StdError[2],
    uncertainty$covParNative$StdError[4],tolerance=1e-10)
})

test_that("fixed anchors remove constrained coordinates from fitting and uncertainty", {
  data <- .laplace_vcc_data()
  fit <- sommer:::.mmes_dispatch(continuous~x,
    random=~vsm(ism(group), sigma2=0.6, fixedSigma2=TRUE)+block,
    family=gaussian(), data=data,
    vcc=data.frame(parameter=1:2, group=1, scale=c(2,1)),
    control=.laplace_vcc_control(parameters=list(dispersion=0.25, estimate=FALSE)))
  reference <- sommer:::.mmes_dispatch(continuous~x,
    random=~vsm(ism(group), sigma2=0.6, fixedSigma2=TRUE)+
      vsm(ism(block), sigma2=0.3, fixedSigma2=TRUE),
    family=gaussian(), data=data,
    control=.laplace_vcc_control(parameters=list(dispersion=0.25, estimate=FALSE)))
  expect_true(fit$convergence)
  expect_equal(unname(unlist(fit$covPar)), c(0.6,0.3))
  expect_false(any(fit$vcParams$free))
  expect_equal(fit$laplaceParameterLayout$covarianceCount, 0L)
  expect_equal(fit$npar, 2L)
  expect_equal(fit$llik, reference$llik, tolerance=1e-10)
  expect_equal(vcov(fit), vcov(reference), tolerance=1e-6)
  expect_true(all(is.na(sommer:::.mmes_laplace_uncertainty(fit)$covParNative$StdError)))
  expect_equal(postPEV(fit, 2L)$Ci, postPEV(reference, 2L)$Ci, tolerance=1e-8)
})

test_that("covariance constraints preserve grouped trials and estimated family parameters", {
  skip_on_cran()
  data <- .laplace_vcc_data()
  set.seed(881)
  data$success <- rbinom(nrow(data),5,plogis(data$continuous))
  data$failure <- 5-data$success
  data$nb <- rnbinom(nrow(data),size=2,mu=exp(data$continuous))
  data$gamma <- rgamma(nrow(data),shape=3,scale=exp(data$continuous)/3)
  data$zip <- ifelse(runif(nrow(data))<.25,0,rpois(nrow(data),exp(data$continuous)))
  combined <- cbind(sqrt(2)*model.matrix(~group-1,data),model.matrix(~block-1,data))
  cases <- list(
    list(fixed=cbind(success,failure)~x,family=binomial(),parameters=list()),
    list(fixed=nb~x,family=nbinom2(2),parameters=list()),
    list(fixed=gamma~x,family=Gamma("log"),parameters=list()),
    list(fixed=zip~x,family=zipoisson(),parameters=list(zeroInflation=.25,estimate=FALSE)))
  for(case in cases){
    control <- .laplace_vcc_control(parameters=case$parameters)
    fit <- sommer:::.mmes_dispatch(case$fixed,random=~group+block,data=data,family=case$family,
      vcc=data.frame(parameter=1:2,group=1,scale=c(2,1)),control=control)
    reference <- sommer:::.mmes_dispatch(case$fixed,random=~vsm(ism(combined)),data=data,
      family=case$family,control=control)
    expect_true(fit$convergence)
    expect_equal(fit$llik,reference$llik,tolerance=1e-7)
    expect_equal(fitted(fit),fitted(reference),tolerance=3e-5)
    expect_equal(fit$familyParameters,reference$familyParameters,tolerance=3e-5)
    expect_equal(fit$npar,reference$npar)
    expect_equal(length(fit$optimizer$par),fit$npar)
    if(isTRUE(fit$familySpecifications[[1]]$free)){
      uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
      expect_true(uncertainty$valid)
      expect_true(is.finite(uncertainty$familyParameters$StdError))
      expect_equal(uncertainty$covParNative$StdError[1],
        2*uncertainty$covParNative$StdError[2],tolerance=1e-10)
    }
  }
})

test_that("Laplace vcc validates groups, types, fixed anchors and likelihood-only boundaries", {
  data <- .laplace_vcc_data()
  fit <- function(vcc, random=~group+block)
    sommer:::.mmes_dispatch(y~x, random=random, family=poisson(), data=data,
      vcc=vcc, control=.laplace_vcc_control(setup=TRUE))
  expect_error(fit(data.frame(parameter=c(1,3), group=1)), "Unknown parameters")
  expect_error(fit(data.frame(parameter=c(1,1), group=1)), "more than once")
  expect_error(fit(data.frame(parameter=1, group=1)), "single parameter")
  expect_error(fit(data.frame(parameter=1:2, group=NA)), "nonmissing")
  expect_error(fit(data.frame(parameter=c(1,1.5), group=1)), "whole numbers")
  expect_error(fit(data.frame(parameter=1:2, group=1, scale=c(1,-1))), "positive")
  expect_error(fit(data.frame(parameter=1:2, group=1, scale=c(1,0))), "non-zero")
  expect_error(fit(data.frame(parameter=1:2, group=1),
    ~vsm(ism(group), sigma2=0.6, fixedSigma2=TRUE)+
      vsm(ism(block), sigma2=0.3, fixedSigma2=TRUE)), "inconsistent")
  expect_error(fit(data.frame(parameter=c(1,2), group=1),
    ~vsm(ar1m(occasion), ism(group))+block), "different kinds")
  expect_error(sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
    W=Matrix::Diagonal(nrow(data)), vcc=data.frame(parameter=1:2,group=1),
    control=.laplace_vcc_control()), "W or weights")
})

test_that("ordinal and nominal constraints retain full post-fit and refit support", {
  skip_on_cran()
  data <- .laplace_vcc_data()
  set.seed(805)
  logits <- cbind(0.3+2.5*data$continuous, -0.5+2.5*data$continuous, 0)
  probabilities <- exp(logits)/rowSums(exp(logits))
  data$category <- factor(vapply(seq_len(nrow(data)), function(row)
    sample.int(3L, 1L, prob=probabilities[row,]), integer(1L)))
  constraint <- data.frame(parameter=1:2, group=1L, scale=c(2,1))
  combined <- cbind(sqrt(2)*model.matrix(~group-1, data),
    model.matrix(~block-1, data))
  for(kind in c("ordinal", "nominal")){
    family <- if(kind=="ordinal") ordinal() else multinomial(reference="3")
    fixed <- if(kind=="ordinal") ordered~x else category~x
    fit <- sommer:::.mmes_dispatch(fixed, random=~group+block, family=family, data=data,
      vcc=constraint, control=.laplace_vcc_control())
    referenceDesign <- if(kind=="nominal") combined[rep(seq_len(nrow(data)),2L),] else combined
    reference <- sommer:::.mmes_dispatch(fixed, random=~vsm(ism(referenceDesign)), family=family, data=data,
      control=.laplace_vcc_control())
    expect_true(fit$convergence)
    expect_equal(fit$covPar[[1]], 2*fit$covPar[[2]], tolerance=1e-12)
    expect_equal(fit$llik, reference$llik, tolerance=2e-7)
    expect_equal(fitted(fit), fitted(reference), tolerance=3e-5)
    uncertainty <- sommer:::.mmes_laplace_uncertainty(fit)
    expect_true(uncertainty$valid)
    expect_equal(uncertainty$covParNative$StdError[1],
      2*uncertainty$covParNative$StdError[2], tolerance=1e-10)
    for(population in c(FALSE, TRUE)){
      prediction <- predict(fit, newdata=data[1:3,], type="response",
        re.form=if(population) NA else NULL, se.fit=TRUE, cov.fit=TRUE)
      expect_true(all(is.finite(prediction$se.fit)))
      expect_equal(rowSums(prediction$fit), rep(1,3), ignore_attr=TRUE, tolerance=1e-10)
    }
    expect_equal(postPEV(fit)$uPevList, postPEV(fit, 2L)$uPevList,
      tolerance=1e-7)
    refit <- sommer:::.boot_mmes_refit(fit, fit$y, 1000L)
    expect_true(refit$converged)
    expect_equal(refit$logLik, fit$llik, tolerance=1e-6)
    covariance <- refit$parameters[grepl("^covariance:", names(refit$parameters))]
    expect_equal(exp(unname(covariance[1])), 2*exp(unname(covariance[2])), tolerance=1e-10)
    if(kind=="nominal"){
      preparation <- fit$nominalPreparation
      objective <- sommer:::.laplace_nominal_objective(preparation)
      proposal <- fit$optimizer$par + 0.1
      evaluation <- objective$evaluate(proposal)
      for(index in seq_along(proposal)){
        upper <- lower <- proposal
        upper[index] <- upper[index]+1e-4
        lower[index] <- lower[index]-1e-4
        numerical <- (objective$evaluate(upper, FALSE)$logLik-
          objective$evaluate(lower, FALSE)$logLik)/2e-4
        expect_equal(as.numeric(evaluation$gradient)[index], numerical, tolerance=2e-5)
      }
      setup <- sommer:::.mmes_dispatch(fixed, random=~group+block, family=family, data=data,
        vcc=constraint, control=.laplace_vcc_control(setup=TRUE))
      expect_equal(setup$vcParams$group, c(1L,1L))
      expect_equal(setup$parameterLayout$covarianceCount, 1L)
    }
  }
})

test_that("mixed categorical familym respects vcc in prediction, PEV and bootstrap", {
  skip_on_cran()
  data <- .laplace_vcc_data()
  set.seed(195)
  original <- data
  data <- do.call(rbind, lapply(c("count","ord","nom"), function(trait){
    result <- original
    result$trait <- trait
    result$y <- switch(trait, count=original$y, ord=as.integer(original$ordered),
      nom=sample.int(3L,nrow(original),replace=TRUE))
    result
  }))
  data$trait <- factor(data$trait, levels=c("count","ord","nom"))
  family <- familym("trait", count=poisson(), ord=ordinal(),
    nom=multinomial(reference="3"))
  constraint <- data.frame(parameter=1:2, group=1L)
  fit <- sommer:::.mmes_dispatch(y~trait+trait:x-1, random=~group+block, data=data, family=family,
    vcc=constraint, control=.laplace_vcc_control())
  expect_true(fit$convergence)
  expect_equal(fit$covPar[[1]], fit$covPar[[2]], tolerance=1e-12)
  expect_equal(fit$laplaceParameterLayout$covarianceCount, 1L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_equal(postPEV(fit)$uPevList, postPEV(fit,2L)$uPevList, tolerance=1e-7)
  rows <- c(1,2,121,122,241,242)
  for(uncertainty in c("conditional", "full")){
    prediction <- predict(fit, newdata=data[rows,], type="response", se.fit=TRUE,
      cov.fit=TRUE, uncertainty=uncertainty)
    expect_true(all(is.finite(unlist(prediction$se.fit))))
  }
  population <- predict(fit,newdata=data[rows,],type="response",
    re.form="population",se.fit=TRUE,cov.fit=TRUE)
  expect_true(all(is.finite(unlist(population$se.fit))))
  expect_equal(as.numeric(diag(population$cov.fit)),
    as.numeric(unlist(population$se.fit)^2),tolerance=1e-8)
  refit <- sommer:::.boot_mmes_refit(fit, data$y, 1000L)
  expect_true(refit$converged)
  expect_equal(refit$logLik, fit$llik, tolerance=1e-5)
  covariance <- refit$parameters[grepl("^covariance:", names(refit$parameters))]
  expect_equal(covariance[1], covariance[2], ignore_attr=TRUE, tolerance=1e-12)
})

test_that("Laplace constraints persist through bootstrap refits and comparisons", {
  data <- .laplace_vcc_data()
  constraint <- data.frame(parameter=1:2,group=1,scale=c(2,1))
  fit <- sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
    vcc=constraint, control=.laplace_vcc_control())
  smaller <- sommer:::.mmes_dispatch(y~1, random=~group+block, family=poisson(), data=data,
    vcc=constraint, control=.laplace_vcc_control())
  expect_equal(anova(smaller, fit)$`Chi Df`[2], 1L)
  unrestricted <- sommer:::.mmes_dispatch(y~x, random=~group+block, family=poisson(), data=data,
    control=.laplace_vcc_control())
  expect_error(anova(smaller, unrestricted), "vcc constraints")
  parameters <- sommer:::.boot_mmes_parameters(fit)
  covariance <- parameters[grepl("^covariance:", names(parameters))]
  expect_equal(exp(unname(covariance[1])), 2*exp(unname(covariance[2])), tolerance=1e-10)
  refit <- sommer:::.boot_mmes_refit(fit, fit$y, 1000L)
  expect_true(refit$converged)
  covariance <- refit$parameters[grepl("^covariance:", names(refit$parameters))]
  expect_equal(exp(unname(covariance[1])), 2*exp(unname(covariance[2])), tolerance=1e-10)
  expect_equal(refit$logLik, fit$llik, tolerance=1e-7)
  fields <- c("optimizer","b","covStruct","randomRelationships","laplaceDiagnostics",
    "familySpecifications","laplaceParameterLayout","vcParams","vcc","fixedCovariance","Ci")
  before <- unserialize(serialize(fit[fields],NULL))
  set.seed(319)
  state <- .Random.seed
  first <- boot_mmes(fit,nsim=2L,seed=39)
  second <- boot_mmes(fit,nsim=2L,seed=39)
  expect_identical(.Random.seed,state)
  expect_identical(fit[fields],before)
  expect_identical(first$draws,second$draws)
  expect_equal(first$successful,2L)
  covarianceColumns <- grepl("^covariance:",colnames(first$draws))
  draws <- first$draws[,covarianceColumns,drop=FALSE]
  expect_equal(exp(draws[,1]),2*exp(draws[,2]),ignore_attr=TRUE,tolerance=1e-10)
})
