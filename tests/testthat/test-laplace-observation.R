observationControl <- function(backend, parameters=list(), setup=FALSE){
  mmes_control(task=list(objective="ml",glmm=list(method="laplace",
    backend=backend,modeTol=1e-10,familyParameters=parameters)),
    output=list(setup=setup),diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

observationData <- function(){
  set.seed(187)
  data <- expand.grid(time=factor(1:3),id=factor(1:12),replicate=1:3)
  data$x <- rnorm(nrow(data))
  data$eta <- .4+.2*data$x+rnorm(12,sd=.5)[data$id]
  data$y <- rpois(nrow(data),exp(data$eta))
  data
}

observationNative <- function(n=10L,q=24L){
  set.seed(409)
  data <- data.frame(y=rep(1,q),id=factor(seq_len(q)))
  prepared <- mmes(y~1,random=~vsm(ism(id),sigma2=.4),data=data,
    family=poisson(),control=observationControl("random",setup=TRUE))
  z <- matrix(rnorm(n*q)/sqrt(q),n,q)
  z[n,] <- z[1,]
  list(X=Matrix::Matrix(cbind(1,rnorm(n)),sparse=TRUE),
    Z=Matrix::Matrix(z,sparse=TRUE),y=rpois(n,2),offset=rep(0,n),
    descriptor=prepared$covStruct[1],Q=Matrix::Diagonal(q))
}

observationContext <- function(fixture,backend="observation",memory=1000,specs=NULL){
  sommer:::sommer_laplace_objective_create(fixture$X,fixture$Z,fixture$y,
    fixture$offset,if(is.null(specs)) "poisson" else "",fixture$descriptor,
    list(fixture$Q),memory,familySpecs=specs,
    familyIndex=if(!is.null(specs)) rep(1L,length(fixture$y)) else NULL,
    backend=backend)
}

test_that("observation algebra matches dense mode, determinants, solves and gradients", {
  fixture <- observationNative()
  fixture$Q <- Matrix::bandSparse(24,k=c(-1,0,1),
    diagonals=list(rep(-.1,23),rep(1.2,24),rep(-.1,23)))
  parameters <- c(.3,-.2,log(.4))
  contexts <- lapply(c("random","observation"),function(backend)
    observationContext(fixture,backend))
  values <- lapply(contexts,function(context)
    sommer:::sommer_laplace_objective_eval(context,parameters,tolerance=1e-11))
  expect_true(all(vapply(values,`[[`,logical(1),"valid")))
  expect_equal(values[[2]]$mode,values[[1]]$mode,tolerance=1e-9)
  expect_equal(values[[2]]$logLik,values[[1]]$logLik,tolerance=1e-10)
  expect_equal(values[[2]]$gradient,values[[1]]$gradient,tolerance=1e-8)
  expect_null(values[[2]]$modeHessian)
  expect_identical(values[[2]]$work$backend$selected,"observation")
  prior <- as.matrix(fixture$Q)/.4
  eta <- as.numeric(fixture$X %*% parameters[1:2]+fixture$Z %*% values[[2]]$mode)
  hessian <- prior+crossprod(as.matrix(fixture$Z),exp(eta)*as.matrix(fixture$Z))
  rhs <- cbind(seq_len(24)/24,rep(1,24))
  expect_equal(sommer:::sommer_laplace_objective_solve(contexts[[2]],rhs),
    solve(hessian,rhs),tolerance=1e-9)
  expect_equal(values[[2]]$logDetHessian,
    as.numeric(determinant(hessian,logarithm=TRUE)$modulus),tolerance=1e-10)
  reference <- sum(dpois(fixture$y,exp(eta),log=TRUE))-
    drop(crossprod(values[[2]]$mode,prior %*% values[[2]]$mode))/2+
    as.numeric(determinant(prior,logarithm=TRUE)$modulus)/2-
    as.numeric(determinant(hessian,logarithm=TRUE)$modulus)/2
  expect_equal(values[[2]]$logLik,reference,tolerance=1e-10)
  numerical <- vapply(seq_along(parameters),function(position){
    upper <- lower <- parameters
    upper[position] <- upper[position]+1e-5
    lower[position] <- lower[position]-1e-5
    (sommer:::sommer_laplace_objective_eval(contexts[[2]],upper,
      tolerance=1e-11,gradient=FALSE)$logLik-
     sommer:::sommer_laplace_objective_eval(contexts[[2]],lower,
       tolerance=1e-11,gradient=FALSE)$logLik)/2e-5
  },numeric(1))
  expect_equal(as.numeric(values[[2]]$gradient),numerical,tolerance=1e-7)
  repeated <- sommer:::sommer_laplace_objective_eval(contexts[[2]],parameters)
  cached <- sommer:::sommer_laplace_objective_eval(contexts[[2]],parameters)
  expect_equal(cached$work$modeEvaluations,repeated$work$modeEvaluations)
  expect_equal(cached$work$backend$covarianceBuilds,repeated$work$backend$covarianceBuilds)
  expect_error(observationContext(fixture,memory=1e-5),"workspace")
  expect_error(observationContext(fixture,backend="invalid"),"backend")
  expect_error(sommer:::sommer_laplace_objective_solve(
    observationContext(fixture),matrix(0,24,1)),"accepted mode")
})

test_that("signed and zero family curvature does not require artificial positive weights", {
  fixture <- observationNative()
  fixture$y[] <- 0
  fixture$offset[] <- log(3)
  fixture$descriptor[[1]]$par[] <- log(.02)
  specs <- list(list(family="ziPoisson",link="log",free=FALSE,
    logParameter=qlogis(.5)))
  contexts <- lapply(c("random","observation"),function(backend)
    observationContext(fixture,backend,specs=specs))
  values <- lapply(contexts,function(context)
    sommer:::sommer_laplace_objective_eval(context,c(0,0,log(.02)),tolerance=1e-11))
  expect_true(all(vapply(values,`[[`,logical(1),"valid")))
  expect_true(any(values[[2]]$observationCurvature<0))
  expect_equal(values[[2]]$logLik,values[[1]]$logLik,tolerance=1e-10)
  expect_equal(values[[2]]$gradient,values[[1]]$gradient,tolerance=1e-8)
  expect_equal(values[[2]]$mode,values[[1]]$mode,tolerance=1e-9)
  fixture$y[] <- 0
  fixture$offset[] <- -1000
  zero <- observationContext(fixture)
  value <- sommer:::sommer_laplace_objective_eval(zero,c(0,0,log(.02)))
  expect_true(value$valid)
  expect_equal(as.numeric(value$observationCurvature),rep(0,length(fixture$y)))
  expect_equal(sommer:::sommer_laplace_objective_solve(zero,diag(24)),
    diag(.02,24),tolerance=1e-12)
})

test_that("automatic selection accounts for size, coupling and memory", {
  fixture <- observationNative()
  evaluate <- function(fixture,memory=1000)
    sommer:::sommer_laplace_objective_eval(observationContext(fixture,"auto",memory),
      c(.2,.1,log(.4)),gradient=FALSE)$work$backend
  expect_identical(evaluate(fixture)$selected,"observation")
  fixture$Z <- Matrix::sparseMatrix(i=seq_len(10),j=seq_len(10),x=1,dims=c(10,24))
  expect_identical(evaluate(fixture)$selected,"random")
  fixture <- observationNative(n=20L,q=24L)
  expect_identical(evaluate(fixture)$selected,"random")
  fixture <- observationNative()
  expect_identical(evaluate(fixture,memory=.01)$selected,"random")
  expect_error(mmes_control(task=list(glmm=list(backend="wrong"))),"backend")
  expect_error(mmes(y~x,random=~id,data=observationData(),family=poisson(),
    control=mmes_control(task=list(glmm=list(backend="observation")))),"Laplace only")
})

test_that("public observation fits preserve prediction, covariance, PEV and serialization", {
  data <- observationData()
  fit <- function(backend) mmes(y~x,
    random=~vsm(ar1m(time,rho=.4,fixed=TRUE),ism(id),sigma2=.3,fixedSigma2=TRUE),
    data=data,family=poisson(),control=observationControl(backend))
  random <- fit("random")
  observation <- fit("observation")
  expect_true(observation$convergence)
  expect_equal(observation$llik,random$llik,tolerance=1e-8)
  expect_equal(observation$b,random$b,tolerance=1e-7)
  expect_null(observation$laplaceModeHessian)
  expect_null(observation$Ci)
  expect_null(observation$fixedCovariance)
  expect_identical(observation$executionPlan$glmm$backend$selected,"observation")
  expect_equal(vcov(observation),vcov(random),tolerance=1e-5)
  expect_equal(postPEV(observation,2L)$Ci,postPEV(random,2L)$Ci,tolerance=1e-7)
  expect_equal(postPEV(observation)$uPevList,postPEV(random)$uPevList,tolerance=1e-7)
  for(target in list(NULL,NA,"population")){
    args <- list(newdata=data[1:3,],type="response",se.fit=TRUE,cov.fit=TRUE)
    if(!is.null(target)) args$re.form <- target
    expect_equal(do.call(predict,c(list(object=observation),args)),
      do.call(predict,c(list(object=random),args)),tolerance=1e-5)
  }
  restored <- unserialize(serialize(observation,NULL))
  expect_equal(postPEV(restored)$uPevList,postPEV(observation)$uPevList)
  expect_equal(sommer:::.boot_mmes_refit(observation,data$y,1000)$logLik,
    observation$llik,tolerance=1e-7)
  expect_equal(simulate(observation,nsim=2,seed=8,re.form="population"),
    simulate(random,nsim=2,seed=8,re.form="population"),tolerance=1e-6)
})

test_that("Gaussian residual blocks and family nuisance gradients retain backend parity", {
  data <- observationData()
  data$unit <- interaction(data$id,data$replicate,drop=TRUE)
  set.seed(94)
  responses <- list(gaussian=data$eta+rnorm(nrow(data),sd=.6),
    binomial=rbinom(nrow(data),1,plogis(data$eta)),
    nb=rnbinom(nrow(data),size=3,mu=exp(data$eta)),
    gamma=rgamma(nrow(data),shape=3,scale=exp(data$eta)/3))
  families <- list(gaussian=gaussian(),binomial=binomial(),nb=nbinom2(3),gamma=Gamma("log"))
  for(name in names(responses)){
    data$y <- responses[[name]]
    fit <- function(backend) mmes(y~x,
      random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),
      rcov=if(name=="gaussian") ~vsm(ar1m(time,rho=.25,fixed=TRUE),ism(unit),
        sigma2=.4,fixedSigma2=TRUE) else ~units,
      data=data,family=families[[name]],control=observationControl(backend))
    random <- fit("random")
    observation <- fit("observation")
    expect_true(observation$convergence)
    expect_equal(observation$llik,random$llik,tolerance=1e-7)
    expect_equal(observation$laplaceDiagnostics$gradient,
      random$laplaceDiagnostics$gradient,tolerance=1e-5)
    expect_equal(vcov(observation),vcov(random),tolerance=2e-4)
    expect_equal(postPEV(observation)$uPevList,postPEV(random)$uPevList,tolerance=1e-6)
  }
})

test_that("ordinal and nominal engines share observation algebra and post-fit operators", {
  data <- observationData()
  set.seed(119)
  ordinalResponse <- ordered(cut(data$eta+rlogis(nrow(data)),
    c(-Inf,-.5,.7,Inf)),labels=c("low","middle","high"))
  nominalResponse <- factor(sample(c("a","b","c"),nrow(data),replace=TRUE))
  for(kind in c("ordinal","nominal")){
    data$y <- if(kind=="ordinal") ordinalResponse else nominalResponse
    family <- if(kind=="ordinal") ordinal() else multinomial(reference="c")
    parameters <- if(kind=="ordinal") list(thresholds=c(-.5,.7),estimate=FALSE) else list()
    fit <- function(backend) mmes(y~x,random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),
      data=data,family=family,control=observationControl(backend,parameters))
    random <- fit("random")
    observation <- fit("observation")
    expect_true(observation$convergence)
    expect_equal(observation$llik,random$llik,tolerance=1e-7)
    expect_null(observation$laplaceModeHessian)
    expect_identical(observation$executionPlan$glmm$backend$selected,"observation")
    expect_equal(vcov(observation),vcov(random),tolerance=1e-4)
    expect_equal(postPEV(observation,2L)$Ci,postPEV(random,2L)$Ci,tolerance=1e-4)
    for(target in list(NULL,"population")){
      args <- list(newdata=data[1:2,],type="response",se.fit=TRUE,cov.fit=TRUE)
      if(!is.null(target)) args$re.form <- target
      expect_equal(do.call(predict,c(list(object=observation),args)),
        do.call(predict,c(list(object=random),args)),tolerance=1e-4)
    }
    expect_equal(predict(unserialize(serialize(observation,NULL)),newdata=data[1:2,]),
      predict(observation,newdata=data[1:2,]))
    expect_equal(sommer:::.boot_mmes_refit(observation,data$y,1000)$logLik,
      observation$llik,tolerance=1e-6)
  }
})

test_that("mixed category blocks, vcc and weighted Gaussian blocks preserve the target", {
  original <- observationData()
  original <- droplevels(original[original$replicate==1,])
  set.seed(518)
  data <- do.call(rbind,lapply(c("count","ord","nom"),function(trait){
    result <- original
    result$trait <- trait
    result$y <- switch(trait,count=original$y,
      ord=as.integer(cut(original$eta+rlogis(nrow(original)),c(-Inf,-.5,.7,Inf))),
      nom=sample.int(3,nrow(original),replace=TRUE))
    result
  }))
  data$trait <- factor(data$trait,levels=c("count","ord","nom"))
  fit <- function(backend) mmes(y~trait+trait:x-1,
    random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),
    data=data,family=familym("trait",count=poisson(),ord=ordinal(),nom=multinomial()),
    control=observationControl(backend,list(ord=list(thresholds=c(-.5,.7),estimate=FALSE))))
  random <- fit("random")
  observation <- fit("observation")
  expect_true(observation$convergence)
  expect_equal(observation$llik,random$llik,tolerance=1e-7)
  expect_equal(vcov(observation),vcov(random),tolerance=1e-4)
  expect_equal(postPEV(observation,2)$Ci,postPEV(random,2)$Ci,tolerance=1e-4)
  for(target in list(NULL,"population")){
    args <- list(newdata=data[c(1,37,73),],type="response",se.fit=TRUE,cov.fit=TRUE)
    if(!is.null(target)) args$re.form <- target
    expect_equal(do.call(predict,c(list(object=observation),args)),
      do.call(predict,c(list(object=random),args)),tolerance=1e-4)
  }
  data <- observationData()
  fit <- function(backend) mmes(y~x,
    random=~vsm(ism(id),sigma2=.3)+vsm(ism(replicate),sigma2=.15),
    vcc=data.frame(parameter=1:2,group=1,scale=c(2,1)),
    data=data,family=poisson(),control=observationControl(backend))
  random <- fit("random")
  observation <- fit("observation")
  expect_true(observation$convergence)
  expect_equal(observation$llik,random$llik,tolerance=1e-6)
  expect_equal(observation$covPar,random$covPar,tolerance=1e-5)
  expect_equal(vcov(observation),vcov(random),tolerance=1e-4)
  expect_equal(postPEV(observation)$uPevList,postPEV(random)$uPevList,tolerance=1e-5)
  set.seed(727)
  data$y <- data$eta+rnorm(nrow(data),sd=.5)
  data$w <- seq(.8,1.2,length.out=nrow(data))
  fit <- function(backend) mmes(y~x,
    random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),data=data,
    W=Matrix::Diagonal(x=seq(.9,1.1,length.out=nrow(data))),weights=~w,
    family=gaussian(),control=observationControl(backend))
  random <- fit("random")
  observation <- fit("observation")
  expect_true(observation$convergence)
  expect_equal(observation$llik,random$llik,tolerance=1e-6)
  expect_equal(postPEV(observation,2)$Ci,postPEV(random,2)$Ci,tolerance=1e-6)
})

test_that("public auto-selected observation fits retain represented relationship levels", {
  data <- observationData()
  data <- data[data$time=="1" & data$replicate %in% 1:2,]
  data$id <- factor(data$id,levels=as.character(1:48))
  precision <- .4^abs(outer(1:48,1:48,`-`))
  dimnames(precision) <- list(levels(data$id),levels(data$id))
  attr(precision,"inverse") <- TRUE
  fit <- function(backend) mmes(y~x,random=~vsm(ism(id),Gu=precision,
    sigma2=.3,fixedSigma2=TRUE),data=data,family=poisson(),
    control=observationControl(backend))
  random <- fit("random")
  auto <- fit("auto")
  expect_true(auto$convergence)
  expect_identical(auto$executionPlan$glmm$backend$requested,"auto")
  expect_identical(auto$executionPlan$glmm$backend$selected,"observation")
  expect_null(auto$laplaceModeHessian)
  expect_equal(auto$llik,random$llik,tolerance=1e-8)
  expect_equal(postPEV(auto)$uPevList,postPEV(random)$uPevList,tolerance=1e-7)
  new <- data[1:2,]
  new$id <- factor(c("17","48"),levels=levels(data$id))
  expect_equal(predict(auto,newdata=new,se.fit=TRUE,cov.fit=TRUE),
    predict(random,newdata=new,se.fit=TRUE,cov.fit=TRUE),tolerance=1e-6)
})

test_that("observation errors and post-fit work restore saved BLAS controls", {
  configured <- 4L
  updates <- integer()
  controller <- list(provider="mock",get=function() configured,
    set=function(value){ configured <<- value; updates <<- c(updates,value) })
  local_mocked_bindings(.mmes_blas_controller=function() controller,.package="sommer")
  control <- observationControl("observation")
  control$execution$blasThreads <- 1L
  data <- observationData()
  fit <- mmes(y~x,random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),
    data=data,family=poisson(),control=control)
  expect_identical(configured,4L)
  for(operation in list(vcov,postPEV,function(object)
    predict(object,newdata=data[1:2,],se.fit=TRUE))){
    updates <- integer()
    operation(fit)
    expect_identical(configured,4L)
    expect_true(1L %in% updates)
    expect_identical(tail(updates,1L),4L)
  }
  fit$laplacePostPEVSettings$denseMemoryMB <- 1e-5
  expect_error(postPEV(fit),"workspace|memory")
  expect_identical(configured,4L)
  expect_error(predict(fit,newdata=data[1:2,],se.fit=TRUE),"workspace|memory")
  expect_identical(configured,4L)
})

test_that("structured covariance shape scores and alternative links share the algebra", {
  data <- observationData()
  data <- droplevels(data[data$replicate==1,])
  prepared <- mmes(y~x,random=~vsm(ar1m(time,rho=.4),ism(id),sigma2=.3),
    data=data,family=poisson(),control=observationControl("random",setup=TRUE))
  create <- function(backend,specs=NULL,response=data$y)
    sommer:::sommer_laplace_objective_create(prepared$X,do.call(cbind,prepared$Z),
      response,rep(0,nrow(data)),if(is.null(specs)) "poisson" else "",
      prepared$covStruct[1],prepared$Ai,backend=backend,familySpecs=specs,
      familyIndex=if(!is.null(specs)) rep(1L,nrow(data)) else NULL)
  parameters <- c(.3,.2,prepared$covStruct[[1]]$par)
  contexts <- lapply(c("random","observation"),create)
  values <- lapply(contexts,function(context)
    sommer:::sommer_laplace_objective_eval(context,parameters,tolerance=1e-11))
  expect_true(all(vapply(values,`[[`,logical(1),"valid")))
  expect_equal(values[[2]]$gradient,values[[1]]$gradient,tolerance=1e-8)
  upper <- lower <- parameters
  upper[4] <- upper[4]+1e-5
  lower[4] <- lower[4]-1e-5
  numerical <- (sommer:::sommer_laplace_objective_eval(contexts[[2]],upper,
    gradient=FALSE,tolerance=1e-11)$logLik-
    sommer:::sommer_laplace_objective_eval(contexts[[2]],lower,
      gradient=FALSE,tolerance=1e-11)$logLik)/2e-5
  expect_equal(as.numeric(values[[2]]$gradient)[4],numerical,tolerance=1e-7)
  for(link in c("probit","cloglog")){
    specs <- list(list(family="binomial",link=link,free=FALSE,logParameter=0))
    contexts <- lapply(c("random","observation"),function(backend)
      create(backend,specs,as.numeric(data$y>1)))
    values <- lapply(contexts,function(context)
      sommer:::sommer_laplace_objective_eval(context,parameters,tolerance=1e-11))
    expect_true(all(vapply(values,`[[`,logical(1),"valid")))
    expect_equal(values[[2]]$logLik,values[[1]]$logLik,tolerance=1e-10)
    expect_equal(values[[2]]$gradient,values[[1]]$gradient,tolerance=1e-8)
  }
})
