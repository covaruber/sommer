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
  populationSE <- predict(fit,newdata=fixture$data[1:3,],
    re.form="population",se.fit=TRUE,cov.fit=TRUE)
  expect_true(all(is.finite(populationSE$se.fit)))
  expect_equal(as.numeric(diag(populationSE$cov.fit)),as.numeric(populationSE$se.fit^2))
  expect_error(predict(fit, re.form="population", type="link"), "not averaged logits")
  rngBefore <- .Random.seed
  fittedBefore <- fit[c("b","covStruct","randomMode","laplaceDiagnostics")]
  population <- predict(fit, newdata=fixture$data[1:3, , drop=FALSE],
    re.form="population")
  expect_identical(.Random.seed,rngBefore)
  expect_identical(fit[c("b","covStruct","randomMode","laplaceDiagnostics")],
    fittedBefore)
  expect_equal(population,predict(fit,re.form="population")[1:3,,drop=FALSE])
  expect_gt(max(abs(population-zeroRandom[1:3,,drop=FALSE])),1e-5)
  expect_equal(rowSums(population), rep(1,3), ignore_attr=TRUE,
    tolerance=1e-12)
  zeroEta <- predict(fit, newdata=fixture$data[1:3, , drop=FALSE],
    re.form=NA, type="link")
  variance <- sommer:::evaluate_covstruct_cpp(fit$covStruct[[1]],
    fit$covStruct[[1]]$par)[1,1]
  reference <- vapply(seq_len(3L), function(row){
    vapply(seq_len(3L), function(category){
      integrate(function(u){
        eta <- matrix(rep(zeroEta[row,],each=length(u)),ncol=2L)+u
        sommer:::.laplace_nominal_probabilities(eta,fit$nominal)[,category] *
          dnorm(u,sd=sqrt(variance))
      },-Inf,Inf,rel.tol=1e-9)$value
    },numeric(1))
  },numeric(3))
  expect_equal(unname(population), t(reference), tolerance=1e-6)
  unseen <- fixture$data[1:3, , drop=FALSE]
  unseen$group <- factor(rep("new-group",3L))
  expect_equal(predict(fit,newdata=unseen,re.form="population"),population)
  expect_error(predict(fit, newdata=transform(fixture$data[1:2, ],
    group=factor(c("unseen", "unseen")))), "unseen random-effect")

  summaryFit <- summary(fit)
  expect_true(all(is.finite(summaryFit$betas$Std.Error)))
  expect_match(summaryFit$laplaceDiagnostics$uncertaintyStatus,
    "fixed-effect curvature available")
  expect_output(print(summaryFit), "Uncertainty:")
  expect_equal(dim(vcov(fit)), c(length(fit$b), length(fit$b)))
  interceptFit <- mmes(y~1, random=~group, data=fixture$data,
    family=multinomial(reference="high"), control=fixture$control)
  expect_s3_class(anova(interceptFit, fit), "anova")
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_equal(dim(residuals(fit)), dim(fitted(fit)))
  expect_equal(dim(residuals(fit, type="pearson")), dim(fitted(fit)))
  devianceResidual <- residuals(fit, type="deviance")
  observedProbability <- fitted(fit)[cbind(seq_along(fit$y),
    match(as.character(fit$y),colnames(fitted(fit))))]
  expect_equal(unname(devianceResidual^2),-2*log(observedProbability),
    ignore_attr=TRUE,tolerance=1e-12)
  expect_equal(sum(devianceResidual^2),fit$deviance,tolerance=1e-12)

  predictionSE <- predict(fit, newdata=fixture$data[1:3, , drop=FALSE],
    se.fit=TRUE, cov.fit=TRUE)
  expect_equal(dim(predictionSE$se.fit), c(3L, 3L))
  expect_equal(dim(predictionSE$cov.fit), c(9L, 9L))
  expect_true(all(is.finite(predictionSE$se.fit)))
  expect_equal(predictionSE, predict(fit,
    newdata=fixture$data[1:3, , drop=FALSE], se.fit=TRUE, cov.fit=TRUE,
    uncertainty="conditional"))
  expect_error(predict(fit, uncertainty="full"), "requires se.fit")
  expect_error(predict(fit, se.fit=TRUE, uncertainty="invalid"), "arg")

  newdata <- fixture$data[1:2, , drop=FALSE]
  fields <- c("b", "optimizer", "covStruct", "randomMode", "randomRelationships",
    "laplaceDiagnostics", "laplaceModeHessian", "fixedDesign", "randomDesign")
  before <- unserialize(serialize(fit[fields], NULL))
  parameterCovariance <- sommer:::.mmes_laplace_uncertainty(fit)$parameterCovariance
  zeroDesign <- sommer:::.laplace_nominal_prediction_design(fit, newdata,
    includeRandom=FALSE)
  fullZero <- predict(fit, newdata=newdata, re.form=NA, type="link",
    se.fit=TRUE, cov.fit=TRUE, uncertainty="full")
  expectedZero <- as.matrix(zeroDesign$fixed %*%
    parameterCovariance[seq_along(fit$b),seq_along(fit$b)] %*%
    Matrix::t(zeroDesign$fixed))
  expect_equal(unname(fullZero$cov.fit), unname(expectedZero), tolerance=1e-9)
  expect_equal(as.numeric(fullZero$se.fit), sqrt(diag(expectedZero)),
    tolerance=1e-9)
  predictionDesign <- sommer:::.laplace_nominal_prediction_design(fit, newdata)
  objective <- sommer:::.laplace_nominal_objective(fit$nominalPreparation)
  probabilityAt <- function(parameters){
    value <- objective$evaluate(parameters, gradient=FALSE)
    expect_true(value$valid, info=value$reason)
    eta <- matrix(as.numeric(predictionDesign$fixed %*%
      parameters[seq_along(fit$b)]) +
      as.numeric(predictionDesign$random %*% value$mode) +
      predictionDesign$offset, nrow=2L)
    as.numeric(sommer:::.laplace_nominal_probabilities(eta, fit$nominal))
  }
  parameters <- fit$optimizer$par
  jacobian <- vapply(seq_along(parameters), function(position){
    step <- 2e-5 * max(1,abs(parameters[position]))
    plus <- minus <- parameters
    plus[position] <- plus[position]+step
    minus[position] <- minus[position]-step
    (probabilityAt(plus)-probabilityAt(minus))/(2*step)
  }, numeric(6L))
  conditionalProbability <- predict(fit, newdata=newdata)
  softmaxJacobian <- matrix(0,6L,4L)
  for(row in 1:2){
    p <- conditionalProbability[row,]
    softmaxJacobian[row+2*(0:2),row+2*(0:1)] <-
      (diag(p)-tcrossprod(p))[,1:2]
  }
  posteriorEta <- as.matrix(predictionDesign$random %*%
    solve(as.matrix(fit$laplaceModeHessian),
      as.matrix(Matrix::t(predictionDesign$random))))
  expectedProbability <- jacobian %*% parameterCovariance %*% t(jacobian) +
    softmaxJacobian %*% posteriorEta %*% t(softmaxJacobian)
  fullProbability <- predict(fit, newdata=newdata, se.fit=TRUE, cov.fit=TRUE,
    uncertainty="full")
  expect_equal(unname(fullProbability$cov.fit), expectedProbability,
    tolerance=1e-6)
  expect_equal(dim(fullProbability$se.fit), c(2L,3L))
  expect_equal(rowSums(fullProbability$fit), rep(1,2), ignore_attr=TRUE)
  expect_equal(fullProbability$fit, conditionalProbability)
  fullSE <- predict(fit, newdata=newdata, se.fit=TRUE, uncertainty="full")
  expect_equal(fullSE$se.fit, fullProbability$se.fit)
  expect_null(fullSE$cov.fit)
  expect_equal(as.numeric(fullProbability$se.fit),
    sqrt(diag(expectedProbability)), tolerance=1e-6)
  fullLink <- predict(fit, newdata=newdata, type="link", se.fit=TRUE,
    cov.fit=TRUE, uncertainty="full")
  linkCovariance <- sommer:::.laplace_nominal_full_prediction_uncertainty(
    fit, predictionDesign)$cov.fit
  expect_equal(unname(fullLink$cov.fit), unname(linkCovariance))
  expect_equal(as.numeric(fullLink$se.fit), sqrt(diag(linkCovariance)))
  expect_equal(dim(fullLink$se.fit), c(2L,2L))
  expect_identical(fit[fields], before)
  lowMemory <- fit
  lowMemory$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(predict(lowMemory, newdata=newdata, se.fit=TRUE,
    uncertainty="full"), "exceeds")

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

test_that("nominal population quadrature includes cross-category covariance and guards", {
  fixture <- multinomialPublicFixture()
  prepared <- sommer:::.laplace_prepare_nominal(y~x,
    ~vsm(usm(.nominalCategory,theta=matrix(c(1,.35,.35,.8),2)),ism(group),sigma2=.4),
    fixture$data,"high",fixture$control)
  object <- list(nominal=prepared$nominal, b=c(.3,.5,-.2,-.4),
    randomFormula=prepared$nominal$random$formula,
    randomFactorLevels=prepared$randomFactorLevels,
    randomRelationships=prepared$relationships, covStruct=prepared$covStruct,
    randomBlockWidths=prepared$randomBlockWidths,
    laplacePostPEVSettings=list(denseMemoryMB=256))
  rows <- c(1L,301L)
  prediction <- list(fixed=prepared$fixedDesign[rows,,drop=FALSE],
    random=prepared$randomDesign[rows,,drop=FALSE], offset=c(.2,-.15),
    n=1L,rowNames="1")
  result <- sommer:::.laplace_nominal_population_prediction(object,prediction)
  sigma <- sommer:::evaluate_covstruct_cpp(object$covStruct[[1]],
    object$covStruct[[1]]$par)
  root <- t(chol(sigma))
  eta <- as.numeric(prediction$fixed %*% object$b)+prediction$offset
  reference <- vapply(1:3,function(category){
    integrate(function(z1) vapply(z1,function(a){
      integrate(function(z2){
        logits <- cbind(eta[1]+root[1,1]*a+root[1,2]*z2,
          eta[2]+root[2,1]*a+root[2,2]*z2)
        sommer:::.laplace_nominal_probabilities(logits,object$nominal)[,category] *
          dnorm(z2)
      },-8,8,rel.tol=1e-8)$value * dnorm(a)
    },numeric(1)),-8,8,rel.tol=1e-8)$value
  },numeric(1))
  expect_lt(max(abs(as.numeric(result)-reference)),1e-6)
  expect_identical(colnames(result),c("low","mid","high"))
  invalid <- object
  invalid$randomRelationships[[1]][1,2] <- .1
  expect_error(sommer:::.laplace_nominal_population_prediction(invalid,prediction),
    "relationship precision")
  invalid <- object
  invalid$randomFormula <- ~group+x
  expect_equal(sommer:::.laplace_nominal_population_prediction(invalid,prediction),result)
  invalid <- object
  invalid$nominal$predictorLevels <- letters[1:4]
  expect_error(sommer:::.laplace_nominal_population_prediction(invalid,prediction),
    "three non-reference")
  invalid <- object
  invalid$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(sommer:::.laplace_nominal_population_prediction(invalid,prediction),
    "exceeds")
  badPrediction <- prediction
  badPrediction$random <- badPrediction$random*2
  scaled <- sommer:::.laplace_nominal_population_prediction(object,badPrediction)
  expect_equal(as.numeric(rowSums(scaled)),1,tolerance=1e-12)
  expect_gt(max(abs(scaled-result)),1e-4)
  badPrediction <- prediction
  groups <- nrow(object$randomRelationships[[1]])
  badPrediction$random[2L,groups+1L] <- 0
  badPrediction$random[2L,groups+2L] <- 1
  independent <- sommer:::.laplace_nominal_population_prediction(object,badPrediction)
  expect_equal(as.numeric(rowSums(independent)),1,tolerance=1e-12)
  expect_gt(max(abs(independent-result)),1e-4)
  rows <- c(1L,2L,301L,302L)
  varying <- list(fixed=prepared$fixedDesign[rows,,drop=FALSE],
    random=prepared$randomDesign[rows,,drop=FALSE],offset=rep(0,4),
    n=2L,rowNames=c("1","2"))
  varying$random[2L,] <- 0
  varied <- sommer:::.laplace_nominal_population_prediction(object,varying)
  expect_equal(as.numeric(rowSums(varied)),rep(1,2),tolerance=1e-12)
  nearZero <- object
  nearZero$covStruct[[1]]$par[1L] <- -100
  expect_equal(sommer:::.laplace_nominal_population_prediction(nearZero,prediction),
    sommer:::.laplace_nominal_probabilities(matrix(eta,nrow=1L,
      dimnames=list("1",NULL)),object$nominal),tolerance=1e-12)
  highVariance <- object
  highVariance$covStruct[[1]]$par[1] <- log(100)
  expect_error(sommer:::.laplace_nominal_population_prediction(highVariance,prediction),
    "accuracy check")
})

test_that("nominal population quadrature supports independent terms and four categories", {
  fixture <- multinomialPublicFixture()
  fixture$data$other <- factor(rep(1:10,length.out=nrow(fixture$data)))
  build <- function(data,random){
    prepared <- sommer:::.laplace_prepare_nominal(y~x,random,data,"high",
      fixture$control)
    d <- length(prepared$nominal$predictorLevels)
    rows <- (seq_len(d)-1L)*nrow(data)+1L
    list(object=list(nominal=prepared$nominal,b=rep(c(.3,.5),d),
      randomFormula=random,randomFactorLevels=prepared$randomFactorLevels,
      randomRelationships=prepared$relationships,covStruct=prepared$covStruct,
      randomBlockWidths=prepared$randomBlockWidths,
      laplacePostPEVSettings=list(denseMemoryMB=256)),
      prediction=list(fixed=prepared$fixedDesign[rows,,drop=FALSE],
        random=prepared$randomDesign[rows,,drop=FALSE],offset=rep(0,d),
        n=1L,rowNames="1"))
  }
  multiple <- build(fixture$data,
    ~vsm(ism(group),sigma2=.2)+vsm(ism(other),sigma2=.3))
  result <- sommer:::.laplace_nominal_population_prediction(multiple$object,
    multiple$prediction)
  eta <- as.numeric(multiple$prediction$fixed %*% multiple$object$b)
  reference <- vapply(1:3,function(category){
    integrate(function(u){
      logits <- matrix(rep(eta,each=length(u)),ncol=2L)+u
      sommer:::.laplace_nominal_probabilities(logits,multiple$object$nominal)[,category] *
        dnorm(u,sd=sqrt(.5))
    },-Inf,Inf,rel.tol=1e-9)$value
  },numeric(1))
  expect_lt(max(abs(as.numeric(result)-reference)),1e-6)

  fixture$data$y <- factor(rep(c("low","mid","extra","high"),
    length.out=nrow(fixture$data)),levels=c("low","mid","extra","high"))
  four <- build(fixture$data,
    ~vsm(ism(.nominalCategory),ism(group),sigma2=.1))
  result <- sommer:::.laplace_nominal_population_prediction(four$object,
    four$prediction)
  expect_equal(unname(rowSums(result)),1,tolerance=1e-12)
  expect_identical(colnames(result),levels(fixture$data$y))
  expect_true(all(result>0 & result<1))
  expect_lt(max(abs(result[1L,1:3]-result[1L,1L])),1e-10)
  expect_equal(result,sommer:::.laplace_nominal_population_prediction(
    four$object,four$prediction))
})
