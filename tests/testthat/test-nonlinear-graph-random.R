nonlinearGraphRandomFixture <- function(){
  group <- factor(rep(seq_len(6L),each=8L))
  x <- rep(seq(-1.5,1.5,length.out=8L),6L)
  data <- data.frame(x=x,group=group)
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x",center=0.1,scale=1.2),
    hidden=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    activation=list(op="tanh",input="hidden"),
    output=list(op="affine",inputs="activation",weights=matrix(3L),bias=4L),
    mean=list(op="identity",input="output")
  ),output="mean",parameters=c("hiddenWeight","hiddenBias",
    "outputWeight","outputBias"),parameterMap=list(
      hiddenWeight=list(type="fixed"),
      hiddenBias=list(type="random",group="group"),
      outputWeight=list(type="fixed"),
      outputBias=list(type="fixed")
    ),random=list(hiddenBias=~vsm(ism(group),sigma2=0.02,
      fixedSigma2=TRUE)))
  coordinates <- c(hiddenWeight=1,hiddenBias=0,outputWeight=1.4,
    outputBias=0.3)
  groupEffect <- c(-0.15,0.08,0.11,-0.05,0.03,-0.02)
  values <- matrix(rep(coordinates,each=nrow(data)),nrow(data),
    dimnames=list(NULL,names(coordinates)))
  values[,"hiddenBias"] <- values[,"hiddenBias"]+groupEffect[group]
  mean <- sommer:::.nonlinear_graph_values(graph,data,values)
  set.seed(732L)
  data$y <- mean+rnorm(nrow(data),sd=0.12)
  control <- sommer:::.mmes_control_internal(task=list(objective="ml",
    nonlinear=list(start=coordinates,
      population=list(order=8L,maxOrder=64L,tol=1e-5))),
    convergence=list(maxIters=100L,likelihoodTol=1e-8,parameterTol=1e-7),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  list(data=data,graph=graph,coordinates=coordinates,truth=mean,control=control)
}

test_that("graph random covariance declarations follow independent and joint contracts",{
  nodes <- list(x=list(op="input",column="x"),
    hidden=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    activation=list(op="tanh",input="hidden"),
    output=list(op="affine",inputs="activation",weights=matrix(3L),bias=4L),
    mean=list(op="identity",input="output"))
  parameters <- c("weight","hiddenBias","outputWeight","outputBias")
  parameterMap <- stats::setNames(lapply(parameters,function(parameter)
    if(parameter %in% c("hiddenBias","outputBias"))
      list(type="random",group="group") else list(type="fixed")),parameters)
  joint <- nonlinear_graph(nodes,"mean",parameters,parameterMap,
    joint=list(parameters=c("hiddenBias","outputBias"),
      random=~vsm(usm(.parameter),ism(group))))
  expect_identical(joint$version,2L)
  expect_identical(sommer:::.nonlinear_graph_validate(joint),joint)
  expect_s3_class(nonlinear_graph(nodes,"mean",parameters,parameterMap),
    "sommer_nonlinear_graph")
  singleRandomMap <- parameterMap
  singleRandomMap$outputBias <- list(type="fixed")
  expect_error(nonlinear_graph(nodes,"mean",parameters,parameterMap,
    random=list(hiddenBias=~vsm(ism(otherGroup)))),"declared exactly once")
  expect_error(nonlinear_graph(nodes,"mean",parameters,singleRandomMap,
    random=list(hiddenBias=~vsm(ism(otherGroup)))),"declared group column")
})

test_that("random graph coordinates use exact conditional curvature and Laplace ML",{
  fixture <- nonlinearGraphRandomFixture()
  data <- fixture$data
  graph <- fixture$graph
  start <- fixture$coordinates
  random <- sommer:::.nonlinear_joint_formula(NULL,graph,data)
  prepared <- sommer:::.nonlinear_prepare(y~0,random,data,fixture$control,
    "exclude","exclude","exclude","exclude",NULL,
    sommer:::.nonlinear_graph_native_layout(graph,data)$inputs)
  layout <- sommer:::.nonlinear_random_layout(graph,prepared)
  expect_true(is.matrix(layout$coordinateMap))
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=graph,control=fixture$control)
  expect_s3_class(fit,"mmes.nonlinear")
  expect_true(fit$convergence)
  expect_lt(fit$nonlinearDiagnostics$modeScore,1e-4)
  expect_match(fit$nonlinearDiagnostics$outerScore,"finite differences")
  expect_true(is.finite(logLik(fit)))
  expect_lt(sqrt(mean((fit$fitted.values-fixture$truth)^2)),0.25)
  residualVariance <- unname(exp(tail(fit$optimizer$par,1L)))
  groupIndex <- as.integer(data$group)
  conditionalObjective <- function(mode){
    parameterRows <- matrix(rep(fit$nonlinearCoordinates,each=nrow(data)),
      nrow(data),dimnames=list(NULL,names(fit$nonlinearCoordinates)))
    parameterRows[,"hiddenBias"] <- parameterRows[,"hiddenBias"]+mode[groupIndex]
    residual <- data$y-sommer:::.nonlinear_graph_values(graph,data,parameterRows,
      fit$control)
    0.5*(sum(residual^2)/residualVariance+sum(mode^2)/0.02)
  }
  referenceMode <- stats::optim(fit$latentRandomMode,conditionalObjective,
    method="BFGS",control=list(reltol=1e-12,maxit=1000L))
  referenceHessian <- stats::optimHess(referenceMode$par,conditionalObjective)
  referenceLogLik <- -referenceMode$value-0.5*nrow(data)*
    (log(2*pi)+log(residualVariance))-0.5*length(referenceMode$par)*log(0.02)-
    0.5*as.numeric(determinant(referenceHessian,logarithm=TRUE)$modulus)
  expect_equal(fit$latentRandomMode,referenceMode$par,tolerance=2e-4)
  expect_equal(as.numeric(logLik(fit)),unname(referenceLogLik),tolerance=2e-5)
  context <- sommer:::sommer_nonlinear_graph_laplace_create(
    fit$fixedDesign,fit$randomDesign,fit$y,fit$offset,fit$graphData,
    fit$graphNative$nodes,fit$graphNative$output,fit$graphNative$transformed,
    fit$nonlinearRandom$coordinateMap*1.0,fit$covStruct,
    fit$randomRelationships,fit$control$execution$denseMemoryMB,
    fit$control$task$nonlinear$mode$maxIters,
    fit$control$task$nonlinear$mode$scoreTol,
    fit$control$task$nonlinear$mode$stepTol)
  at <- fit$optimizer$par
  mode <- fit$latentRandomMode
  score <- sommer:::sommer_nonlinear_graph_laplace_eval(context,at,TRUE)
  referenceScore <- numeric(length(at))
  for(j in seq_along(at)){
    upper <- lower <- at
    step <- 5e-5*(1+abs(at[j]))
    upper[j] <- upper[j]+step
    lower[j] <- lower[j]-step
    referenceScore[j] <- (
      sommer:::sommer_nonlinear_graph_laplace_eval(context,upper,FALSE)$logLik-
      sommer:::sommer_nonlinear_graph_laplace_eval(context,lower,FALSE)$logLik
    )/(2*step)
  }
  expect_equal(as.numeric(score$gradient),referenceScore,tolerance=2e-3)
  conditional <- sommer:::sommer_nonlinear_graph_laplace_conditional(
    context,at,mode)
  expect_true(conditional$valid,info=conditional$reason)
  numericalHessian <- matrix(0,length(mode),length(mode))
  for(j in seq_along(mode)){
    upper <- lower <- mode
    upper[j] <- upper[j]+1e-5
    lower[j] <- lower[j]-1e-5
    numericalHessian[,j] <- -(
      sommer:::sommer_nonlinear_graph_laplace_conditional(context,at,upper)$score-
      sommer:::sommer_nonlinear_graph_laplace_conditional(context,at,lower)$score
    )/(2e-5)
  }
  expect_equal(as.matrix(conditional$hessian),numericalHessian,tolerance=2e-4)
  expect_true(any(abs(predict(fit)-predict(fit,re.form=NA))>1e-6))
  population <- predict(fit,re.form="population")
  expect_equal(length(population),nrow(data))
  expect_true(all(is.finite(population)))
  expect_gt(max(abs(population-predict(fit,re.form=NA))),1e-6)
  quadrature <- sommer:::.laplace_normal_quadrature(64L)
  populationRows <- rep(seq_len(nrow(data)),each=length(quadrature$nodes))
  populationData <- data[populationRows,,drop=FALSE]
  populationParameters <- matrix(rep(fit$nonlinearCoordinates,
    each=length(populationRows)),length(populationRows),
    dimnames=list(NULL,names(fit$nonlinearCoordinates)))
  populationParameters[,"hiddenBias"] <- populationParameters[,"hiddenBias"]+
    rep(sqrt(0.02)*quadrature$nodes,times=nrow(data))
  integratedValues <- sommer:::.nonlinear_graph_values(graph,populationData,
    populationParameters,fit$control)
  referencePopulation <- as.numeric(
    matrix(integratedValues,nrow(data),byrow=TRUE)%*%quadrature$weights)
  expect_equal(as.numeric(population),referencePopulation,tolerance=2e-5)
  saved <- tempfile()
  on.exit(unlink(saved),add=TRUE)
  saveRDS(fit,saved)
  restored <- readRDS(saved)
  expect_equal(predict(restored,newdata=data[1:5,]),
    predict(fit,newdata=data[1:5,]),tolerance=1e-10)
  expect_error(summary(fit,inference=TRUE),"not supported")
  conditional <- simulate(fit,nsim=2L,seed=918L,re.form="conditional")
  populationSimulation <- simulate(fit,nsim=2L,seed=918L,re.form="population")
  expect_identical(populationSimulation,
    simulate(fit,nsim=2L,seed=918L,re.form="population"))
  expect_equal(dim(populationSimulation),c(nrow(data),2L))
  expect_false(isTRUE(all.equal(conditional,populationSimulation)))
})

test_that("joint and reduced-rank graph-coordinate blocks preserve covariance mapping",{
  fixture <- nonlinearGraphRandomFixture()
  graph <- nonlinear_graph(fixture$graph$nodes,fixture$graph$output,
    fixture$graph$parameters,parameterMap=list(
      hiddenWeight=list(type="fixed"),
      hiddenBias=list(type="random",group="group"),
      outputWeight=list(type="fixed"),
      outputBias=list(type="random",group="group")
    ),joint=list(parameters=c("hiddenBias","outputBias"),
      random=~vsm(rrm(.parameter,k=1,loadings=matrix(c(1,.5),2,1),
        fixed=c(TRUE,TRUE)),ism(group),sigma2=0.02,fixedSigma2=TRUE)))
  fit <- sommer:::.mmes_dispatch(y~0,data=fixture$data,nonlinear=graph,
    control=fixture$control)
  expect_true(fit$convergence)
  expect_true(sommer:::.nonlinear_has_reduced_rank(fit$covStruct))
  expect_lt(length(fit$latentRandomMode),length(fit$randomMode))
  expect_true(all(is.finite(predict(fit,re.form="population"))))
  expect_true(all(is.finite(predict(fit,newdata=fixture$data[1:4,]))))
})
