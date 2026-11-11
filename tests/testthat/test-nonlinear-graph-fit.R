nonlinearGraphFitControl <- function(start=NULL){
  sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(start=start)),
    convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

nonlinearGraphFitFixture <- function(hiddenLayers=1L,n=240L){
  if(hiddenLayers==1L){
    parameters <- c(inputWeight=1.2,inputBias=-0.4,
      outputWeight=1.7,outputBias=0.25)
    graph <- nonlinear_graph(list(
      x=list(op="input",column="x",center=0.2,scale=1.1),
      hidden=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
      activation=list(op="tanh",input="hidden"),
      output=list(op="affine",inputs="activation",weights=matrix(3L),bias=4L),
      mean=list(op="identity",input="output")),output="mean",
      parameters=names(parameters))
  }else{
    parameters <- c(inputWeight=1.2,inputBias=-0.4,
      middleWeight=1.1,middleBias=0.2,outputWeight=1.7,outputBias=0.25)
    graph <- nonlinear_graph(list(
      x=list(op="input",column="x",center=0.2,scale=1.1),
      hidden1=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
      activation1=list(op="tanh",input="hidden1"),
      hidden2=list(op="affine",inputs="activation1",weights=matrix(3L),bias=4L),
      activation2=list(op="tanh",input="hidden2"),
      output=list(op="affine",inputs="activation2",weights=matrix(5L),bias=6L),
      mean=list(op="identity",input="output")),output="mean",
      parameters=names(parameters))
  }
  x <- seq(-2,2,length.out=n)
  group <- factor(rep(seq_len(24L),length.out=n))
  data <- data.frame(x=x,z=rep(c(-1,0,1),length.out=n),group=group)
  mean <- nonlinear_graph_eval(graph,data,parameters)$output
  set.seed(193)
  data$y <- mean+rnorm(nlevels(group),sd=0.3)[group]+rnorm(n,sd=0.12)
  list(graph=graph,parameters=parameters,data=data,mean=mean)
}

test_that("graph Gaussian objective and score match an independent dense oracle",{
  fixture <- nonlinearGraphFitFixture(2L,72L)
  nodes <- fixture$graph$nodes
  nodes$output$bias <- 5L
  parameterNames <- fixture$graph$parameters[-6L]
  parameterMap <- stats::setNames(lapply(parameterNames,function(name)
    list(type="fixed")),parameterNames)
  parameterMap$inputWeight <- list(type="fixed",transform="log")
  graph <- nonlinear_graph(nodes,fixture$graph$output,parameterNames,parameterMap)
  graphStart <- fixture$parameters[parameterNames]
  graphStart[["inputWeight"]] <- log(graphStart[["inputWeight"]])
  data <- fixture$data
  control <- nonlinearGraphFitControl(fixture$parameters)
  prepared <- sommer:::.nonlinear_prepare(y~0+z,~group,data,control,
    "exclude","exclude","exclude","exclude",NULL,
    sommer:::.nonlinear_graph_native_layout(graph,data)$inputs)
  native <- sommer:::.nonlinear_graph_native_layout(graph,
    data[prepared$obsInfo$included,,drop=FALSE])
  mapping <- sommer:::.laplace_covariance_map(
    prepared$covStruct[seq_along(prepared$rtermss)],prepared$rtermss,NULL)
  random <- Matrix::Matrix(do.call(cbind,prepared$Z),sparse=TRUE)
  fixed <- prepared$X
  response <- as.numeric(stats::model.response(prepared$modelFrame))
  parameters <- c(fixed=0.13,graphStart,log(0.3^2),log(0.12^2))
  context <- sommer:::sommer_nonlinear_graph_gaussian_create(fixed,random,
    response,rep(0,length(response)),native$inputs,native$nodes,native$output,
    native$transformed,mapping$structures,prepared$Ai)
  denseLogLik <- function(at){
    mean <- nonlinear_graph_eval(graph,
      data[prepared$obsInfo$included,,drop=FALSE],
      at[2L+seq_along(graphStart)-1L])$output+
      as.numeric(fixed%*%at[1L])
    z <- as.matrix(random)
    k <- length(graphStart)
    covariance <- exp(at[k+2L])*tcrossprod(z)+
      diag(exp(at[k+3L]),length(response))
    factor <- chol(covariance)
    residual <- response-mean
    alpha <- backsolve(factor,forwardsolve(t(factor),residual))
    -0.5*(length(response)*log(2*pi)+2*sum(log(diag(factor)))+
      sum(residual*alpha))
  }
  value <- sommer:::sommer_nonlinear_gaussian_eval(context,parameters)
  expect_true(value$valid)
  expect_equal(value$logLik,denseLogLik(parameters),tolerance=1e-10)
  for(step in c(1e-4,1e-5)){
    numerical <- vapply(seq_along(parameters),function(index){
      upper <- lower <- parameters
      upper[index] <- upper[index]+step
      lower[index] <- lower[index]-step
      (denseLogLik(upper)-denseLogLik(lower))/(2*step)
    },numeric(1))
    expect_equal(as.numeric(value$gradient),numerical,tolerance=2e-6)
  }
})

test_that("fixed graph fits preserve linear parity and expose fitted-model methods",{
  fixture <- nonlinearGraphFitFixture(1L,240L)
  data <- fixture$data
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")),output="output",
    parameters=c("slope","intercept"))
  linearData <- data
  set.seed(73)
  linearData$y <- 0.7*data$x+1.1+rnorm(nlevels(data$group),sd=0.3)[data$group]+
    rnorm(nrow(data),sd=0.12)
  starts <- c(slope=0.65,intercept=1)
  set.seed(1001)
  rng <- .Random.seed
  graphFit <- sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=graph,
    control=nonlinearGraphFitControl(starts))
  expect_identical(.Random.seed,rng)
  expect_true(graphFit$convergence)
  expect_true(is.finite(as.numeric(logLik(graphFit))))
  parity <- nonlinear_graph_eval(graph,data.frame(x=c(-2,0,3)),starts)
  expect_equal(parity$output,starts[["slope"]]*c(-2,0,3)+starts[["intercept"]])
  expectedGradient <- cbind(slope=c(-2,0,3),intercept=1)
  expect_equal(parity$gradient,expectedGradient,tolerance=1e-12)
  expect_equal(predict(graphFit,newdata=linearData),fitted(graphFit),
    tolerance=1e-9)
  expect_equal(predict(graphFit,newdata=linearData[1:5,],re.form=NA),
    graphFit$fixedMean[1:5],tolerance=1e-9)
  expect_named(coef(graphFit),c("graph:slope","graph:intercept"))
  expect_named(coef(graphFit,type="all"),
    c("graph:slope","graph:intercept",
      paste0("covariance:",graphFit$randomMetadata,":",graphFit$covStruct[[1]]$par_names),
      "residual:variance"))
  expect_true(inherits(summary(graphFit),"summary.mmes.nonlinear"))
  expect_equal(as.numeric(residuals(graphFit)),
    graphFit$y-graphFit$fitted.values,tolerance=1e-12)
  expect_equal(predict(graphFit,newdata=linearData,type="link"),
    predict(graphFit,newdata=linearData),tolerance=1e-12)
  beforeSimulation <- .Random.seed
  firstSimulation <- simulate(graphFit,nsim=2L,seed=510L)
  expect_identical(.Random.seed,beforeSimulation)
  expect_equal(firstSimulation,simulate(graphFit,nsim=2L,seed=510L))
  populationSimulation <- simulate(graphFit,nsim=2L,seed=511L,
    re.form="population")
  expect_equal(dim(populationSimulation),c(nrow(linearData),2L))
  repeatedFit <- sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=graph,
    control=nonlinearGraphFitControl(starts))
  expect_equal(repeatedFit$optimizer$par,graphFit$optimizer$par,tolerance=1e-10)
  profiled <- profile_nonlinear(graphFit,replications=1L)
  expect_true(all(c("graphEvaluation","marginalObjective","marginalScore") %in%
    profiled$table$phase))
  inferenceReady <- isTRUE(graphFit$nonlinearDiagnostics$identified) &&
    !length(graphFit$nonlinearDiagnostics$boundaryParameters)
  if(inferenceReady){
    inferred <- post_nonlinear(graphFit)
    expect_true(all(is.finite(diag(vcov(inferred,type="working")))))
    profiledInference <- profile_nonlinear(graphFit,replications=1L,inference=TRUE)
    expect_true("postFitInference" %in% profiledInference$table$phase)
  }else{
    expect_error(post_nonlinear(graphFit),
      "locally identifiable|covariance boundaries")
  }
  expect_error(predict(graphFit,newdata=data.frame(z=1)),"missing columns")
  missingPrediction <- linearData[1L,,drop=FALSE]
  missingPrediction$x <- NA_real_
  expect_error(predict(graphFit,newdata=missingPrediction),"finite numeric values")
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=graph,
    control=nonlinearGraphFitControl()),"requires a complete named")
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=graph,
    family=poisson(),control=nonlinearGraphFitControl(starts)),
    "count responses must be non-negative integer counts")
  boundedOutput <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    hidden=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="tanh",input="hidden")),output="output",
    parameters=c("slope","intercept"))
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=boundedOutput,
    control=nonlinearGraphFitControl(c(slope=0.1,intercept=0.1))),
    "requires an identity output")
  randomGraph <- nonlinear_graph(graph$nodes,graph$output,graph$parameters,
    parameterMap=list(slope=list(type="random",group="group"),
      intercept=list(type="fixed")))
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=randomGraph,
    control=nonlinearGraphFitControl(starts)),
    "explicit random or joint covariance formulas")
  smallMemory <- nonlinearGraphFitControl(starts)
  smallMemory$execution$denseMemoryMB <- 0.0001
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=linearData,nonlinear=graph,
    control=smallMemory),"workspace")
  expect_error(sommer:::.mmes_dispatch(y~1,random=~group,data=linearData,nonlinear=graph,
    control=nonlinearGraphFitControl(starts)),"rank deficient")
  incomplete <- linearData
  incomplete$x[1L] <- NA_real_
  filteredFit <- sommer:::.mmes_dispatch(y~0,random=~group,data=incomplete,nonlinear=graph,
    control=nonlinearGraphFitControl(starts))
  expect_equal(nobs(filteredFit),nrow(incomplete)-1L)
  restored <- unserialize(serialize(graphFit,NULL))
  expect_equal(predict(restored,newdata=linearData[1:5,]),
    predict(graphFit,newdata=linearData[1:5,]))
})

test_that("graph composition matches a classical logistic fit",{
  set.seed(917)
  data <- data.frame(x=seq(-4,5,length.out=240L),
    zero=0,group=factor(rep(seq_len(24L),each=10L)))
  data$y <- 4*plogis(1.2*(data$x-0.7))+
    rnorm(nlevels(data$group),sd=0.25)[data$group]+rnorm(nrow(data),sd=0.1)
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    one=list(op="input",column="zero",center=-1),
    linear=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    activation=list(op="tanh",input="linear"),
    mean=list(op="affine",inputs=c("activation","one"),
      weights=matrix(c(3L,3L),ncol=1L),bias=0L),
    output=list(op="identity",input="mean")),output="output",
    parameters=c("halfRate","shift","halfAsymptote"),
    parameterMap=list(halfRate=list(type="fixed"),
      shift=list(type="fixed"),
      halfAsymptote=list(type="fixed",transform="log")))
  graphStart <- c(halfRate=0.6,shift=-0.42,halfAsymptote=log(2))
  graphFit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,nonlinear=graph,
    control=nonlinearGraphFitControl(graphStart))
  curveFit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,
    nonlinear=nonlinear_curve("x","logistic"),
    control=nonlinearGraphFitControl(c(asymptote=4,rate=1.2,midpoint=0.7)))
  expect_true(graphFit$convergence)
  expect_true(curveFit$convergence)
  expect_equal(as.numeric(logLik(graphFit)),as.numeric(logLik(curveFit)),
    tolerance=2e-5)
  expect_equal(predict(graphFit,newdata=data),predict(curveFit,newdata=data),
    tolerance=2e-4)
})

test_that("one- and two-hidden-layer graph fits recover predictions",{
  for(depth in 1:2){
    fixture <- nonlinearGraphFitFixture(depth,240L)
    start <- fixture$parameters+c(rep(0.01,length(fixture$parameters)))
    fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=fixture$data,nonlinear=fixture$graph,
      control=nonlinearGraphFitControl(start))
    expect_true(fit$convergence)
    expect_true(fit$nonlinearDiagnostics$identified)
    expect_equal(predict(fit,newdata=fixture$data),fitted(fit),tolerance=1e-8)
    fixedPrediction <- predict(fit,newdata=fixture$data,re.form=NA)
    expect_lt(sqrt(mean((fixedPrediction-fixture$mean)^2)),0.15)
  }
})
