nonlinearGraphInferenceControl <- function(start){
  sommer:::.mmes_control_internal(task=list(objective="ml",
    nonlinear=list(start=start)),
    convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

nonlinearGraphInferenceFixture <- function(n=120L){
  set.seed(914L)
  data <- data.frame(x=seq(-1.5,1.5,length.out=n),
    group=factor(rep(seq_len(12L),length.out=n)))
  data$y <- 0.8*data$x+1.1+rnorm(nlevels(data$group),sd=0.25)[data$group]+
    rnorm(n,sd=0.15)
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")),output="output",
    parameters=c("slope","intercept"),
    parameterMap=list(slope=list(type="fixed",transform="log"),
      intercept=list(type="fixed")))
  list(data=data,graph=graph,start=c(slope=log(0.75),intercept=1))
}

test_that("native graph post-fit uncertainty matches an affine Gaussian oracle",{
  fixture <- nonlinearGraphInferenceFixture()
  data <- fixture$data
  fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,
    nonlinear=fixture$graph,
    control=nonlinearGraphInferenceControl(fixture$start))
  expect_true(fit$convergence)
  expect_true(fit$nonlinearDiagnostics$identified)
  post <- post_nonlinear(fit,mode=2L)
  expect_match(post$nonlinearInference$method,"exact Gaussian likelihood")
  expect_equal(post$nonlinearInference$outerScoreMethod,"analytic native score")

  at <- fit$optimizer$par
  z <- as.matrix(fit$randomDesign)
  logLikelihood <- function(value){
    mean <- exp(value[1L])*data$x+value[2L]
    covariance <- exp(value[3L])*tcrossprod(z)+
      diag(exp(value[4L]),nrow(data))
    factor <- chol(covariance)
    residual <- data$y-mean
    alpha <- backsolve(factor,forwardsolve(t(factor),residual))
    -0.5*(nrow(data)*log(2*pi)+2*sum(log(diag(factor)))+
      sum(residual*alpha))
  }
  expected <- solve(-stats::optimHess(at,logLikelihood,
    control=list(ndeps=rep(1e-4,length(at)))))
  expect_equal(unname(vcov(post,type="working")),unname(expected),
    tolerance=2e-3,ignore_attr=TRUE)
  variance <- exp(at[3L])
  conditional <- solve(crossprod(z)/fit$dispersion+
    diag(1/variance,ncol(z)))
  expect_equal(vcov(post,type="conditional"),conditional,
    tolerance=1e-9,ignore_attr=TRUE)
  motion <- post$nonlinearInference$modeSensitivity
  expectedPEV <- conditional+motion%*%vcov(post,type="working")%*%t(motion)
  expect_equal(vcov(post,type="random"),expectedPEV,
    tolerance=1e-9,ignore_attr=TRUE)
  expect_match(post$nonlinearInference$randomMethod,
    "not an exact frequentist PEV")
  expect_equal(unname(coef(post)[["graph:slope"]]),unname(exp(at[1L])),
    tolerance=1e-10)

  newdata <- data.frame(x=c(-0.4,0.2,0.9),
    group=factor(c("1","4","8"),levels=levels(data$group)))
  prediction <- predict(post,newdata,se.fit=TRUE,cov.fit=TRUE)
  expect_equal(prediction$fit,predict(post,newdata),tolerance=1e-10)
  expect_equal(prediction$se.fit^2,diag(prediction$cov.fit),
    tolerance=1e-10,ignore_attr=TRUE)
  fixed <- predict(post,newdata,re.form=NA,se.fit=TRUE)
  expect_equal(fixed$fit,
    exp(at[1L])*newdata$x+at[2L],tolerance=1e-8)

  restored <- unserialize(serialize(fit,NULL))
  restoredPost <- post_nonlinear(restored,mode=2L)
  expect_equal(vcov(restoredPost,type="working"),vcov(post,type="working"),
    tolerance=1e-9,ignore_attr=TRUE)
  expect_equal(predict(restored,newdata,re.form=NA),
    predict(fit,newdata,re.form=NA),tolerance=1e-10)

  boundary <- fit
  boundary$nonlinearDiagnostics$boundaryParameters <- "covariance boundary"
  expect_error(post_nonlinear(boundary),"flagged covariance boundaries")
  settings <- sommer:::.nonlinear_inference_settings(fit)$numerical
  asymmetric <- matrix(c(1,0.2,0.3,1),2L)
  expect_error(sommer:::.nonlinear_positive_inverse(asymmetric,settings,
    "curvature"),"asymmetric")
})

test_that("native graph family predictions use graph link outputs and response integration",{
  fixture <- nonlinearGraphInferenceFixture(96L)
  data <- fixture$data
  eta <- 0.65*data$x-0.15+rnorm(nlevels(data$group),sd=0.35)[data$group]
  set.seed(122L)
  data$y <- rpois(nrow(data),lambda=exp(eta))
  start <- c(slope=log(0.6),intercept=-0.1)
  fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,
    nonlinear=fixture$graph,family=poisson(),
    control=nonlinearGraphInferenceControl(start))
  expect_true(fit$convergence)
  expect_true(fit$nonlinearGraph)
  expect_true(fit$nonlinearFamily)
  post <- post_nonlinear(fit,mode=1L)
  expect_match(post$nonlinearInference$outerScoreMethod,"finite differences")
  expect_true(all(diag(vcov(post,type="all"))>0))
  newdata <- data[seq_len(6L),,drop=FALSE]
  fixedLink <- predict(fit,newdata,re.form=NA,type="link")
  expectedLink <- fit$nonlinearParameters[["slope"]]*newdata$x+
    fit$nonlinearParameters[["intercept"]]
  expect_equal(fixedLink,expectedLink,tolerance=1e-8)
  expect_equal(predict(fit,newdata,re.form=NA,type="response"),
    exp(expectedLink),tolerance=1e-8)
  population <- predict(fit,newdata,re.form="population",type="response")
  expect_true(all(is.finite(population)))
  expect_true(all(population>0))
  simulation <- simulate(fit,nsim=2L,seed=994L,re.form="population")
  expect_equal(dim(simulation),c(nrow(data),2L))
  expect_true(all(as.matrix(simulation)>=0))
})

test_that("native graph inference refines conditional precision without changing the fit",{
  set.seed(75110L)
  data <- expand.grid(x=seq(-1.5,1.5,length.out=12L),
    group=factor(seq_len(4L)))
  eta <- .25+.7*data$x+rnorm(4L,sd=.3)[data$group]
  data$y <- rbinom(nrow(data),1L,plogis(eta))
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")),output="output",
    parameters=c("slope","intercept"),
    parameterMap=list(slope=list(type="fixed"),
      intercept=list(type="random",group="group")),
    random=list(intercept=~vsm(ism(group),sigma2=.09,fixedSigma2=TRUE)))
  fit <- mmes(y~0+graph_term("net"),data=data,family=binomial(),
    predictors=list(net=graph),start=c(slope=.6,intercept=.2),
    method=mmes_method(objective="ml"),
    control=mmes_control(optimizer=list(maxit=200L),
      convergence=list(likelihoodTol=1e-9,parameterTol=1e-8),
      execution=list(denseMemoryMB=256),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_true(fit$convergence)
  original <- list(control=fit$control,numericalControl=fit$.numericalControl,
    parameters=fit$optimizer$par,logLik=fit$llik,mode=fit$randomMode)
  context <- sommer:::.nonlinear_postfit_context(fit,256)
  value <- context$evaluate(fit$optimizer$par)
  expect_equal(value$logLik,fit$llik,tolerance=1e-10)
  post <- post_nonlinear(fit,mode=2L)
  expect_equal(post$nonlinearInference$modeSettings$scoreTol,1e-12)
  expect_equal(post$nonlinearInference$modeSettings$stepTol,1e-12)
  expect_true(all(is.finite(vcov(post,type="working"))))
  expect_true(all(diag(vcov(post,type="working"))>0))
  scalarLaplace <- function(parameters){
    predictor <- parameters[["graph:slope"]]*data$x+
      parameters[["graph:intercept"]]
    sum(vapply(levels(data$group),function(group){
      rows <- data$group==group
      score <- function(mode)
        sum(data$y[rows]-plogis(predictor[rows]+mode))-mode/.09
      mode <- stats::uniroot(score,c(-10,10),tol=1e-12)$root
      probability <- plogis(predictor[rows]+mode)
      hessian <- sum(probability*(1-probability))+1/.09
      sum(stats::dbinom(data$y[rows],1L,probability,log=TRUE))-
        .5*mode^2/.09-.5*log(.09)-.5*log(hessian)
    },numeric(1)))
  }
  oracle <- -stats::optimHess(fit$optimizer$par,scalarLaplace,
    control=list(ndeps=rep(1e-3,length(fit$optimizer$par))))
  expect_equal(post$nonlinearInference$information,oracle,
    tolerance=2e-5,ignore_attr=TRUE)
  prediction <- predict(fit,newdata=data.frame(x=0,
    group=factor("1",levels=levels(data$group))),
    target="zero",type="link",se.fit=TRUE,interval="confidence")
  expect_true(is.finite(prediction$se.fit))
  expect_lt(prediction$lower,prediction$upper)
  strict <- fit
  strict$.numericalControl$task$nonlinear$mode$scoreTol <- 1e-13
  strict$.numericalControl$task$nonlinear$mode$stepTol <- 1e-13
  expect_equal(vcov(post,type="working"),vcov(post_nonlinear(strict,mode=0L),
    type="working"),tolerance=1e-4)
  expect_identical(list(control=fit$control,numericalControl=fit$.numericalControl,
    parameters=fit$optimizer$par,logLik=fit$llik,mode=fit$randomMode),original)
})
