accelerationControl <- function(warmStart="off",gradient="analytic",trace="panels",
    blocks="off",backend="random",parameters=list()){
  sommer:::.mmes_control_internal(task=list(objective="ml",
    glmm=list(method="laplace",backend=backend,warmStart=warmStart,gradient=gradient,
      trace=trace,blocks=blocks,maxit=1000L,modeTol=1e-11,
      familyParameters=parameters)),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

accelerationFixture <- function(q=12L,n=36L,dense=FALSE){
  set.seed(972)
  data <- data.frame(y=1,id=factor(seq_len(q)))
  setup <- sommer:::.mmes_dispatch(y~1,random=~vsm(ism(id),sigma2=.4),data=data,
    family=poisson(),control=sommer:::.mmes_control_internal(task=list(objective="ml",
      glmm=list(method="laplace")),output=list(setup=TRUE),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  list(X=Matrix::Matrix(cbind(1,rnorm(n)),sparse=TRUE),
    Z=if(dense) Matrix::Matrix(matrix(rnorm(n*q)/sqrt(q),n,q),sparse=TRUE) else
      Matrix::sparseMatrix(i=seq_len(n),j=rep(seq_len(q),length.out=n),x=1,dims=c(n,q)),
    y=rpois(n,2),offset=rep(0,n),descriptor=setup$covStruct[1],
    Q=Matrix::Diagonal(q))
}

accelerationContext <- function(f,warmStart="off",trace="panels",blocks="off",
    backend="random",specs=NULL){
  sommer:::sommer_laplace_objective_create(f$X,f$Z,f$y,f$offset,
    if(is.null(specs)) "poisson" else "",f$descriptor,list(f$Q),1000,
    backend=backend,warmStart=warmStart,trace=trace,blocks=blocks,
    familySpecs=specs,familyIndex=if(!is.null(specs)) rep(1L,length(f$y)) else NULL)
}

accelerationEval <- function(context,parameters=c(.3,-.1,log(.4)),gradient=TRUE,maxit=100L){
  sommer:::sommer_laplace_objective_eval(context,parameters,
    tolerance=1e-11,maxModeIterations=maxit,gradient=gradient)
}

test_that("selected contractions and independent blocks preserve exact sparse algebra", {
  for(dense in c(FALSE,TRUE)){
    f <- accelerationFixture(dense=dense)
    if(!dense) f$Q <- Matrix::bandSparse(12L,k=c(-1,0,1),
      diagonals=list(rep(-.15,11),rep(1.3,12),rep(-.15,11)))
    reference <- accelerationEval(accelerationContext(f))
    h <- as.matrix(reference$modeHessian)
    for(trace in c("selected","panels")) for(blocks in c("on","off")){
      context <- accelerationContext(f,trace=trace,blocks=blocks)
      value <- accelerationEval(context)
      expect_true(value$valid)
      expect_equal(value$logLik,reference$logLik,tolerance=1e-10)
      expect_equal(value$gradient,reference$gradient,tolerance=1e-9)
      expect_equal(sommer:::sommer_laplace_objective_solve(context,diag(12)),
        solve(h),tolerance=1e-10)
      if(trace=="selected") expect_gt(value$work$backend$selectedBuilds,0L)
    }
  }
  f <- accelerationFixture()
  blocked <- accelerationEval(accelerationContext(f,trace="selected",blocks="on"))
  expect_identical(blocked$work$backend$blocks,12L)
  expect_identical(blocked$work$backend$contractionSolves,0L)
  expect_equal(blocked$gradient,accelerationEval(accelerationContext(f))$gradient,
    tolerance=1e-10)
})

test_that("missing selected entries are solved exactly and never treated as zeros", {
  f <- accelerationFixture(q=24L,n=72L)
  f$Q <- Matrix::bandSparse(24,k=c(-1,0,1),
    diagonals=list(rep(-.2,23),rep(1.4,24),rep(-.2,23)))
  f$Z[1,] <- c(1,rep(0,22),1)
  f$y[] <- 0
  f$offset[] <- -1000
  selected <- accelerationEval(accelerationContext(f,trace="selected"))
  reference <- accelerationEval(accelerationContext(f))
  expect_true(selected$valid)
  expect_gt(selected$work$backend$contractionSolves,0L)
  expect_equal(selected$gradient,reference$gradient,tolerance=1e-12)
  expect_equal(selected$logLik,reference$logLik,tolerance=1e-12)
  f <- accelerationFixture(q=4L,n=12L)
  f$Z[1,] <- c(1,0,0,1)
  f$y[] <- 0
  f$offset[] <- 0
  blocks <- accelerationContext(f,trace="selected",blocks="on")
  zero <- accelerationEval(blocks,c(-1000,0,log(.4)))
  coupled <- accelerationEval(blocks,c(.1,0,log(.8)))
  expect_true(coupled$valid)
  expect_equal(coupled$gradient,
    accelerationEval(accelerationContext(f),c(.1,0,log(.8)))$gradient,tolerance=1e-10)
  expect_identical(zero$work$backend$blocks,4L)
  expect_identical(coupled$work$backend$blocks,3L)
})

test_that("warm mode starts preserve likelihood and scores while reducing iterations", {
  f <- accelerationFixture(q=30L,n=180L)
  for(backend in c("random","observation")){
    cold <- accelerationContext(f,backend=backend)
    warm <- accelerationContext(f,warmStart="auto",backend=backend)
    for(step in seq_len(12L)){
      parameters <- c(.3+step*.003,-.1,log(.4)+step*.004)
      baseline <- accelerationEval(cold,parameters)
      value <- accelerationEval(warm,parameters)
      expect_true(value$valid)
      expect_equal(value$logLik,baseline$logLik,tolerance=1e-9)
      expect_equal(value$mode,baseline$mode,tolerance=1e-8)
      expect_equal(value$gradient,baseline$gradient,tolerance=1e-8)
    }
    expect_identical(value$work$warmStart$attempts,11L)
    expect_lt(value$work$warmStart$modeIterations,baseline$work$warmStart$modeIterations)
    repeated <- accelerationEval(warm,parameters)
    expect_identical(repeated$work$warmStart,value$work$warmStart)
    failed <- accelerationEval(warm,parameters+c(2,0,0),maxit=1L)
    expect_false(failed$valid)
    expect_identical(failed$work$warmStart$coldRetries,1L)
    recovered <- accelerationEval(warm,parameters)
    expect_true(recovered$valid)
    expect_equal(recovered$gradient,value$gradient,tolerance=1e-8)
    expect_error(sommer:::sommer_laplace_objective_solve(warm,matrix(0,1,1)),
      "conformable")
  }
})

test_that("ZIP warm starts use the authoritative cold root and preserve signed curvature", {
  f <- accelerationFixture()
  f$y[] <- 0
  f$offset[] <- log(3)
  specs <- list(list(family="ziPoisson",link="log",free=FALSE,
    logParameter=qlogis(.5)))
  cold <- accelerationContext(f,specs=specs)
  auto <- accelerationContext(f,warmStart="auto",specs=specs)
  warm <- accelerationContext(f,warmStart="on",trace="selected",blocks="on",specs=specs)
  for(step in 1:3){
    parameters <- c(step*.01,0,log(.02))
    reference <- accelerationEval(cold,parameters)
    value <- accelerationEval(warm,parameters)
    automatic <- accelerationEval(auto,parameters)
    expect_true(value$valid)
    expect_true(any(value$observationCurvature<0))
    expect_equal(value$logLik,reference$logLik,tolerance=1e-12)
    expect_equal(value$gradient,reference$gradient,tolerance=1e-10)
    expect_equal(automatic$logLik,reference$logLik,tolerance=1e-12)
  }
  expect_false(automatic$work$warmStart$enabled)
  expect_identical(value$work$warmStart$coldChecks,2L)
})

accelerationMixedData <- function(){
  set.seed(628)
  base <- expand.grid(id=factor(1:10),replicate=factor(1:3))
  base$x <- rnorm(nrow(base))
  eta <- .2+.3*base$x+rnorm(10,sd=.4)[base$id]
  out <- do.call(rbind,lapply(c("count","binary","ord","nom","normal"),function(trait){
    data <- base
    data$trait <- trait
    data$y <- switch(trait,count=rpois(nrow(base),exp(eta)),
      binary=rbinom(nrow(base),1,plogis(eta)),
      ord=as.integer(cut(eta+rlogis(nrow(base)),c(-Inf,-.4,.8,Inf))),
      nom=sample.int(3,nrow(base),replace=TRUE),normal=eta+rnorm(nrow(base),sd=.6))
    data
  }))
  out$trait <- factor(out$trait,levels=c("count","binary","ord","nom","normal"))
  out
}

test_that("mixed nominal analytic scores include nuisance, residual and constrained covariance derivatives", {
  data <- accelerationMixedData()
  model <- sommer:::.mmes_dispatch(y~trait+trait:x-1,
    random=~vsm(ism(id),sigma2=.3)+vsm(ism(replicate),sigma2=.15),
    vcc=data.frame(parameter=1:2,group=1L,scale=c(2,1)),data=data,
    family=familym("trait",count=poisson(),binary=binomial(),
      ord=ordinal(),nom=multinomial(),normal=gaussian()),
    control=accelerationControl(
      parameters=list(ord=list(thresholds=c(-.4,.8)))))
  parameters <- model$optimizer$par+rep(c(.1,-.05),length.out=length(model$optimizer$par))
  for(backend in c("random","observation")) for(trace in
      if(backend=="random") c("panels","selected") else "panels"){
    model$laplaceDiagnostics$nativeWork$backend$selected <- backend
    model$laplacePostPEVSettings$trace <- trace
    context <- sommer:::.laplace_objective_context(model)
    evaluate <- function(parameters,gradient=TRUE)
      sommer:::.laplace_constrained_evaluate(parameters,model$laplaceParameterLayout,
        function(expanded) accelerationEval(context,expanded,gradient))
    value <- evaluate(parameters)
    numerical <- vapply(seq_along(parameters),function(position){
      upper <- lower <- parameters
      upper[position] <- upper[position]+1e-4
      lower[position] <- lower[position]-1e-4
      (evaluate(upper,FALSE)$logLik-evaluate(lower,FALSE)$logLik)/2e-4
    },numeric(1))
    expect_true(value$valid)
    expect_lt(max(abs(as.numeric(value$gradient)-numerical)),3e-6)
  }
})

test_that("standalone nominal acceleration preserves analytic scores and mode systems", {
  set.seed(37)
  data <- expand.grid(id=factor(1:12),replicate=1:8)
  data$x <- rnorm(nrow(data))
  data$y <- factor(sample(c("a","b","c"),nrow(data),replace=TRUE))
  reference <- sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),
    data=data,family=multinomial(),control=accelerationControl())
  for(backend in c("random","observation")) for(trace in
      if(backend=="random") c("panels","selected") else "panels"){
    fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),data=data,
      family=multinomial(),control=accelerationControl(warmStart="auto",trace=trace,
        backend=backend,blocks=if(backend=="random") "on" else "off"))
    expect_true(fit$convergence)
    expect_equal(fit$llik,reference$llik,tolerance=1e-8)
    expect_equal(fit$b,reference$b,tolerance=2e-6)
    expect_equal(postPEV(fit)$uPevList,postPEV(reference)$uPevList,tolerance=1e-6)
    expect_gt(fit$laplaceDiagnostics$nativeWork$warmStart$attempts,0L)
    expect_equal(sommer:::.boot_mmes_refit(fit,fit$nominalPreparation$nominal$response,1000)$logLik,
      fit$llik,tolerance=1e-7)
  }
  numeric <- sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),data=data,
    family=multinomial(),control=accelerationControl(gradient="numeric"))
  expect_equal(numeric$llik,reference$llik,tolerance=1e-8)
  expect_equal(numeric$b,reference$b,tolerance=2e-6)
})

test_that("public scalar controls persist across uncertainty, bootstrap and serialization", {
  f <- accelerationFixture(q=10,n=100)
  data <- data.frame(y=f$y,x=f$X[,2],id=factor(rep(1:10,10)))
  fit <- function(control) sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.4,fixedSigma2=TRUE),
    data=data,family=poisson(),control=control)
  reference <- fit(accelerationControl())
  accelerated <- fit(accelerationControl(warmStart="auto",trace="selected",blocks="on"))
  expect_null(accelerated$control$execution$blasThreads)
  expect_null(accelerated$threadingDiagnostics$blasThreadsRequested)
  expect_true(accelerated$convergence)
  expect_equal(accelerated$llik,reference$llik,tolerance=1e-9)
  expect_equal(vcov(accelerated),vcov(reference),tolerance=1e-5)
  expect_equal(postPEV(accelerated,2)$Ci,postPEV(reference,2)$Ci,tolerance=1e-7)
  expect_equal(predict(accelerated,newdata=data[1:2,],se.fit=TRUE),
    predict(reference,newdata=data[1:2,],se.fit=TRUE),tolerance=1e-5)
  expect_equal(sommer:::.boot_mmes_refit(accelerated,data$y,1000)$logLik,
    reference$llik,tolerance=1e-8)
  expect_equal(vcov(unserialize(serialize(accelerated,NULL))),vcov(accelerated))
  expect_null(accelerated$Ci)
  expect_null(accelerated$fixedCovariance)
  expect_equal(fit(accelerationControl(gradient="numeric"))$llik,
    reference$llik,tolerance=1e-8)
  expect_identical(accelerated$executionPlan$glmm$trace,"selected")
  expect_identical(accelerated$executionPlan$glmm$blocks,"on")
  expect_identical(accelerated$executionPlan$glmm$gradient,"analytic")
})

test_that("mixed analytic and numerical full fits retain inference and refit parity", {
  data <- droplevels(subset(accelerationMixedData(),trait %in% c("count","ord","nom")))
  fit <- function(control) sommer:::.mmes_dispatch(y~trait+trait:x-1,
    random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),data=data,
    family=familym("trait",count=poisson(),ord=ordinal(),nom=multinomial()),
    control=control)
  reference <- fit(accelerationControl(gradient="numeric"))
  analytic <- fit(accelerationControl(warmStart="auto",trace="selected",blocks="on"))
  expect_true(analytic$convergence)
  expect_equal(analytic$llik,reference$llik,tolerance=1e-7)
  expect_equal(analytic$b,reference$b,tolerance=1e-5)
  expect_equal(vcov(analytic),vcov(reference),tolerance=1e-4)
  expect_equal(postPEV(analytic)$uPevList,postPEV(reference)$uPevList,tolerance=1e-4)
  expect_equal(sommer:::.boot_mmes_refit(analytic,data$y,1000)$logLik,
    reference$llik,tolerance=1e-6)
  expect_equal(predict(unserialize(serialize(analytic,NULL)),newdata=data[1:2,]),
    predict(analytic,newdata=data[1:2,]))
})

test_that("low-level nominal inverse-subset traces solve unavailable arbitrary derivative entries", {
  q <- 24L
  n <- q/2L
  response <- rep(1:3,length.out=n)
  predictor <- cbind(seq(-.2,.3,length.out=n),seq(.1,-.4,length.out=n))
  design <- Matrix::Diagonal(q)
  precision <- Matrix::bandSparse(q,k=c(-1,0,1),
    diagonals=list(rep(-.2,q-1),rep(1.4,q),rep(-.2,q-1)))
  fixed <- Matrix::Matrix(cbind(rep(c(1,0),each=n),rep(c(0,1),each=n)),sparse=TRUE)
  derivative <- Matrix::Matrix(matrix(1/q,q,q),sparse=TRUE)
  evaluate <- function(trace){
    context <- sommer:::sommer_multinomial_mode_create(response,design,precision,2L,
      trace=trace)
    sommer:::sommer_multinomial_laplace_eval(context,predictor,precision,
      fixed,list(derivative),rep(0,q),100L,1e-11,TRUE)
  }
  reference <- evaluate("panels")
  selected <- evaluate("selected")
  expect_true(selected$valid)
  expect_gt(selected$work$selectedBuilds,0L)
  expect_gt(selected$work$contractionSolves,0L)
  expect_equal(selected$gradient,reference$gradient,tolerance=1e-10)
  expect_equal(selected$logLik,reference$logLik,tolerance=1e-12)
})

test_that("unsupported choices are rejected and exposed by mmes_options", {
  for(name in c("warmStart","gradient","trace","blocks"))
    expect_error(sommer:::.mmes_control_internal(task=list(glmm=setNames(list("invalid"),name))),name)
  f <- accelerationFixture()
  expect_error(accelerationContext(f,backend="observation",trace="selected"),"random-space")
  expect_error(accelerationContext(f,backend="observation",blocks="on"),"random-space")
  data <- data.frame(y=f$y,id=factor(rep(1:12,3)))
  for(choice in list(list(warmStart="on"),list(gradient="numeric"),
      list(trace="selected"),list(blocks="on")))
    expect_error(sommer:::.mmes_dispatch(y~1,random=~id,data=data,family=poisson(),
      control=sommer:::.mmes_control_internal(task=list(glmm=choice))),"Laplace only")
  paths <- sommer:::.mmes_capabilities()$controls$path
  expect_true(all(paste0("control.task.glmm.",c("warmStart","gradient","trace","blocks")) %in% paths))
})
