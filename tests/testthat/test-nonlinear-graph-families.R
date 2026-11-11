nonlinearGraphFamilyFixture <- function(){
  set.seed(2047)
  data <- expand.grid(x=seq(-1.2,1.2,length.out=5L),
    group=factor(seq_len(4L)))
  data$group <- factor(rep(seq_len(4L),each=5L))
  data$environment <- factor(rep(seq_len(3L),length.out=nrow(data)))
  eta <- 0.18*data$x+0.35
  data$y <- stats::rpois(nrow(data),exp(eta))
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")
  ),output="output",parameters=c("slope","intercept"))
  control <- sommer:::.mmes_control_internal(
    task=list(objective="ml",nonlinear=list(start=c(slope=0.15,intercept=0.3))),
    convergence=list(maxIters=100L,likelihoodTol=1e-8,parameterTol=1e-7),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  list(data=data,graph=graph,control=control,start=c(slope=0.15,intercept=0.3))
}

test_that("graph family modes use exact conditional curvature and central outer scores",{
  fixture <- nonlinearGraphFamilyFixture()
  data <- fixture$data
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=0.12,
    fixedSigma2=TRUE),data=data,family=poisson(),nonlinear=fixture$graph,
    control=fixture$control)
  expect_true(fit$convergence)
  expect_true(fit$nonlinearGraph)
  expect_null(fit$nonlinearRandom)
  expect_match(fit$nonlinearDiagnostics$outerScore,"central finite differences")
  expect_equal(length(fit$nonlinearDiagnostics$gradient),
    length(fit$optimizer$par))
  expect_error(sommer::mmes(y~0+graph_term("net"),
    random=~vsm(ism(group),sigma2=0.12,fixedSigma2=TRUE),data=data,
    family=poisson(),predictors=list(net=fixture$graph),
    start=list(parameters=fixture$start),method=mmes_method(objective="ml"),
    control=mmes_control(optimizer=list(gradient="analytic"))),
    "optimizer gradient is unavailable")

  context <- sommer:::sommer_nonlinear_graph_family_create(
    fit$fixedDesign,fit$randomDesign,fit$nativeResponse,fit$offset,
    fit$graphData,fit$graphNative$nodes,fit$graphNative$output,
    fit$graphNative$transformed,matrix(0,ncol(fit$randomDesign),
      length(fit$nonlinearCoordinates)),fit$covStruct,
    fit$randomRelationships,fit$familySpecifications,
    as.integer(fit$laplaceFamilyIndex-1L),fit$priorWeights)
  outer <- sommer:::sommer_nonlinear_graph_family_eval(context,
    fit$optimizer$par,TRUE)
  reference <- vapply(seq_along(fit$optimizer$par),function(j){
    upper <- lower <- fit$optimizer$par
    step <- 1e-5*(1+abs(fit$optimizer$par[j]))
    upper[j] <- upper[j]+step
    lower[j] <- lower[j]-step
    (sommer:::sommer_nonlinear_graph_family_eval(context,upper,FALSE)$logLik-
      sommer:::sommer_nonlinear_graph_family_eval(context,lower,FALSE)$logLik)/
      (2*step)
  },numeric(1))
  expect_equal(as.numeric(outer$gradient),reference,tolerance=1e-10)
  expect_match(outer$work$outerScore,"central finite differences")
  conditional <- sommer:::sommer_nonlinear_graph_family_conditional(context,
    fit$optimizer$par,fit$latentRandomMode)
  expect_true(conditional$valid,info=conditional$reason)
  mode <- fit$latentRandomMode
  numerical <- matrix(0,length(mode),length(mode))
  for(j in seq_along(mode)){
    upper <- lower <- mode
    step <- 1e-5*(1+abs(mode[j]))
    upper[j] <- upper[j]+step
    lower[j] <- lower[j]-step
    upperScore <- sommer:::sommer_nonlinear_graph_family_conditional(context,
      fit$optimizer$par,upper)$score
    lowerScore <- sommer:::sommer_nonlinear_graph_family_conditional(context,
      fit$optimizer$par,lower)$score
    numerical[,j] <- -(upperScore-lowerScore)/(2*step)
  }
  expect_equal(as.matrix(conditional$hessian),numerical,tolerance=2e-4)
})

test_that("graph family random coordinates support independent, joint, and reduced-rank priors",{
  fixture <- nonlinearGraphFamilyFixture()
  data <- fixture$data
  parameterMap <- list(slope=list(type="fixed"),
    intercept=list(type="random",group="group"))
  independent <- nonlinear_graph(fixture$graph$nodes,fixture$graph$output,
    fixture$graph$parameters,parameterMap,
    random=list(intercept=~vsm(ism(group),sigma2=0.12,fixedSigma2=TRUE)))
  independentFit <- sommer:::.mmes_dispatch(y~0,data=data,family=poisson(),
    nonlinear=independent,control=fixture$control)
  expect_true(independentFit$convergence)
  expect_equal(length(independentFit$latentRandomMode),4L)

  jointMap <- list(slope=list(type="random",group="group"),
    intercept=list(type="random",group="group"))
  joint <- nonlinear_graph(fixture$graph$nodes,fixture$graph$output,
    fixture$graph$parameters,jointMap,
    joint=list(parameters=c("slope","intercept"),
      random=~vsm(usm(.parameter),ism(group),sigma2=0.12,fixedSigma2=TRUE)))
  jointFit <- sommer:::.mmes_dispatch(y~0,data=data,family=poisson(),
    nonlinear=joint,control=fixture$control)
  expect_true(jointFit$convergence)
  expect_equal(length(jointFit$latentRandomMode),8L)

  reduced <- nonlinear_graph(fixture$graph$nodes,fixture$graph$output,
    fixture$graph$parameters,parameterMap,
    random=list(intercept=~vsm(rrm(environment,k=1,
      loadings=matrix(c(1,0.7,-0.3),3L,1L),fixed=rep(TRUE,3L)),
      ism(group),sigma2=0.12,fixedSigma2=TRUE)))
  reducedFit <- sommer:::.mmes_dispatch(y~0,data=data,family=poisson(),
    nonlinear=reduced,control=fixture$control)
  expect_true(reducedFit$convergence)
  expect_lt(length(reducedFit$latentRandomMode),
    length(reducedFit$randomMode))
})

test_that("graph family preparation covers native links and expanded familym rows",{
  set.seed(612)
  families <- list(logit=binomial("logit"),probit=binomial("probit"),
    cloglog=binomial("cloglog"),counts=poisson(),
    normal=gaussian(),gamma=Gamma("log"),nb=nbinom2(theta=2),
    zip=zipoisson(pi=0.15),ordered=ordinal(thresholds=c(-0.5,0.5)),
    nominal=multinomial())
  pieces <- lapply(names(families),function(trait){
    x <- seq(-1,1,length.out=18L)
    group <- factor(rep(seq_len(6L),each=3L))
    eta <- 0.2*x+0.25
    y <- switch(trait,
      logit=rbinom(length(x),1,plogis(eta)),
      probit=rbinom(length(x),1,pnorm(eta)),
      cloglog=rbinom(length(x),1,-expm1(-exp(eta))),
      counts=rpois(length(x),exp(eta)),
      normal=rnorm(length(x),eta,sd=0.8),
      gamma=rgamma(length(x),shape=3,scale=exp(eta)/3),
      nb=rnbinom(length(x),size=2,mu=exp(eta)),
      zip={value<-rpois(length(x),exp(eta));value[runif(length(x))<0.2]<-0;value},
      ordered=sample(rep(1:3,length.out=length(x))),
      nominal=sample(rep(1:3,length.out=length(x))))
    data.frame(x=x,group=group,trait=trait,y=y)
  })
  data <- do.call(rbind,pieces)
  family <- do.call(familym,c(list(by="trait"),families))
  control <- sommer:::.mmes_control_internal(
    task=list(objective="ml",nonlinear=list(start=c(slope=0.2,intercept=0.25)),
      glmm=list(familyParameters=list(
        gamma=list(dispersion=1,estimate=FALSE),
        nb=list(theta=2,estimate=FALSE),
        zip=list(zeroInflation=0.15,estimate=FALSE)))),
    convergence=list(maxIters=100L,likelihoodTol=1e-7,parameterTol=1e-6),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")
  ),output="output",parameters=c("slope","intercept"))
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=0.1,
    fixedSigma2=TRUE),data=data,family=family,nonlinear=graph,control=control)
  expect_true(fit$convergence)
  expect_identical(fit$nonlinearGraph,TRUE)
  expect_true(nrow(fit$graphData)>nrow(data))
  expect_equal(nrow(fit$graphData),
    nrow(data)+sum(data$trait=="nominal"))
  expect_setequal(unique(names(fit$laplaceFamilies)),names(families))
  expect_true(fit$familySpecifications$normal$free)
  expect_true(is.finite(fit$dispersion[["normal"]]))
  expect_match(fit$nonlinearDiagnostics$outerScore,"central finite differences")
  context <- sommer:::sommer_nonlinear_graph_family_create(
    fit$fixedDesign,fit$randomDesign,fit$nativeResponse,fit$offset,
    fit$graphData,fit$graphNative$nodes,fit$graphNative$output,
    fit$graphNative$transformed,matrix(0,ncol(fit$randomDesign),
      length(fit$nonlinearCoordinates)),fit$covStruct,
    fit$randomRelationships,fit$familySpecifications,
    as.integer(fit$laplaceFamilyIndex-1L),fit$priorWeights)
  conditional <- sommer:::sommer_nonlinear_graph_family_conditional(context,
    fit$optimizer$par,fit$latentRandomMode)
  numerical <- matrix(0,length(fit$latentRandomMode),
    length(fit$latentRandomMode))
  for(j in seq_along(fit$latentRandomMode)){
    upper <- lower <- fit$latentRandomMode
    step <- 1e-5*(1+abs(fit$latentRandomMode[j]))
    upper[j] <- upper[j]+step
    lower[j] <- lower[j]-step
    upperScore <- sommer:::sommer_nonlinear_graph_family_conditional(context,
      fit$optimizer$par,upper)$score
    lowerScore <- sommer:::sommer_nonlinear_graph_family_conditional(context,
      fit$optimizer$par,lower)$score
    numerical[,j] <- -(upperScore-lowerScore)/(2*step)
  }
  expect_equal(as.matrix(conditional$hessian),numerical,tolerance=3e-4)
})
