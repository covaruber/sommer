nonlinearStructureData <- function(groups=3L, axes=3L){
  set.seed(642)
  data <- expand.grid(x=seq(0,5,length.out=8L),environment=seq_len(axes),
    subject=seq_len(groups))
  data$environment <- factor(data$environment)
  data$subject <- factor(data$subject)
  data$occasion <- factor(rep(1:2,length.out=nrow(data)))
  data$y <- 4*exp(-0.6*data$x)+0.08*as.numeric(data$environment)+
    rnorm(nrow(data),sd=0.15)
  data
}

nonlinearStructureContext <- function(random, data=nonlinearStructureData(),
    curveRandom=FALSE, additive=NULL, memoryMB=1000,
    curve="exponential",parameter="rate"){
  model <- nonlinear_curve("x",curve,random=if(curveRandom)
    stats::setNames(list(random),parameter) else list())
  joint <- if(curveRandom) sommer:::.nonlinear_joint_formula(additive,model,data) else random
  control <- mmes_control(task=list(objective="ml"),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  prepared <- sommer:::.nonlinear_prepare(y~0,joint,data,control,
    "exclude","exclude","exclude","exclude",NULL,data$x)
  count <- length(prepared$rtermss)
  structures <- prepared$covStruct[seq_len(count)]
  z <- do.call(cbind,prepared$Z)
  if(curveRandom){
    layout <- sommer:::.nonlinear_random_layout(model,prepared)
    context <- sommer:::sommer_nonlinear_laplace_create(prepared$X,z,data$y,
      numeric(nrow(data)),data$x,curve,layout$coordinate,layout$start,
      structures,prepared$Ai,memoryMB=memoryMB)
    evaluate <- function(par,gradient=TRUE)
      sommer:::sommer_nonlinear_laplace_eval(context,par,gradient)
  }else{
    context <- sommer:::sommer_nonlinear_gaussian_create(prepared$X,z,data$y,
      numeric(nrow(data)),data$x,curve,structures,prepared$Ai,memoryMB)
    evaluate <- function(par,gradient=TRUE)
      sommer:::sommer_nonlinear_gaussian_eval(context,par,gradient)
  }
  theta <- if(curve=="logistic") c(log(4),log(1.2),2.5) else c(log(4),log(0.6))
  parameters <- c(theta,unlist(lapply(structures,
    function(s) s$par[s$free]),use.names=FALSE),log(0.15^2))
  list(context=context,evaluate=evaluate,prepared=prepared,structures=structures,
    parameters=parameters,data=data,design=z,model=model)
}

nonlinearStructureDifference <- function(evaluate,par,step=1e-5){
  vapply(seq_along(par),function(j){
    upper <- lower <- par
    upper[j] <- upper[j]+step; lower[j] <- lower[j]-step
    up <- evaluate(upper); down <- evaluate(lower)
    expect_true(up$valid,info=up$reason)
    expect_true(down$valid,info=down$reason)
    (up$logLik-down$logLik)/(2*step)
  },numeric(1))
}

nonlinearStructureDense <- function(fixture,par){
  structures <- fixture$structures
  position <- 2L
  covariance <- diag(exp(tail(par,1)),nrow(fixture$data))
  for(term in seq_along(structures)){
    descriptor <- structures[[term]]
    count <- sum(descriptor$free)
    descriptor$par[descriptor$free] <- par[position+seq_len(count)]
    position <- position+count
    sigma <- sommer:::evaluate_covstruct_cpp(descriptor,descriptor$par)
    kernel <- kronecker(sigma,solve(as.matrix(fixture$prepared$Ai[[term]])))
    blocks <- which(fixture$prepared$Zind==term)
    design <- as.matrix(do.call(cbind,fixture$prepared$Z[blocks]))
    covariance <- covariance+design %*% kernel %*% t(design)
  }
  residual <- fixture$data$y-exp(par[1])*exp(-exp(par[2])*fixture$data$x)
  factor <- chol(covariance)
  -0.5*(length(residual)*log(2*pi)+2*sum(log(diag(factor)))+
    sum(forwardsolve(t(factor),residual)^2))
}

test_that("native positive-definite covariance structures pass nonlinear acceptance checks",{
  formulas <- list(
    identity=~vsm(ism(subject),sigma2=0.04),
    diagonal=~vsm(dsm(environment),ism(subject),sigma2=0.04),
    unstructured=~vsm(usm(environment),ism(subject),sigma2=0.04),
    autoregressive=~vsm(ar1m(environment),ism(subject),sigma2=0.04),
    compound=~vsm(csm(environment),ism(subject),sigma2=0.04),
    fa=~vsm(fam(environment,k=1),ism(subject),sigma2=0.04),
    homogeneousFA=~vsm(fam(environment,k=1,variance="homogeneous"),
      ism(subject),sigma2=0.04),
    kronecker=~vsm(ar1m(occasion),fam(environment,k=1),
      ism(subject),sigma2=0.04))
  for(formula in formulas){
    for(curveRandom in c(FALSE,TRUE)){
      fixture <- nonlinearStructureContext(formula,curveRandom=curveRandom)
      value <- fixture$evaluate(fixture$parameters)
      expect_true(value$valid,info=value$reason)
      if(!value$valid) next
      expect_equal(as.numeric(value$gradient),
        nonlinearStructureDifference(fixture$evaluate,fixture$parameters),
        tolerance=6e-5)
      if(!curveRandom){
        expect_equal(value$logLik,nonlinearStructureDense(fixture,fixture$parameters),
          tolerance=1e-9)
      }
    }
  }
})

test_that("additive reduced-rank integration matches the full singular covariance model",{
  formulas <- list(
    rank1=~vsm(rrm(environment,k=1),ism(subject),sigma2=0.04),
    rank2=~vsm(rrm(environment,k=2),ism(subject),sigma2=0.04),
    kronecker=~vsm(ar1m(occasion),rrm(environment,k=1),ism(subject),sigma2=0.04),
    twoRR=~vsm(rrm(occasion,k=1),rrm(environment,k=1),ism(subject),sigma2=0.04),
    mixed=~vsm(rrm(environment,k=1),ism(subject),sigma2=0.04)+
      vsm(ism(subject),sigma2=0.02))
  for(formula in formulas){
    fixture <- nonlinearStructureContext(formula)
    value <- fixture$evaluate(fixture$parameters)
    expect_true(value$valid,info=value$reason)
    if(!value$valid) next
    expect_equal(value$logLik,nonlinearStructureDense(fixture,fixture$parameters),
      tolerance=1e-9)
    expect_equal(as.numeric(value$gradient),
      nonlinearStructureDifference(fixture$evaluate,fixture$parameters),
      tolerance=5e-5)
    expect_lt(length(value$latentMode),length(value$mode))
    expect_equal(as.numeric(value$fitted),4*exp(-0.6*fixture$data$x)+
      as.numeric(fixture$design %*% value$mode),tolerance=1e-11)
    expect_equal(fixture$evaluate(fixture$parameters,FALSE)$logLik,value$logLik)
  }
})

test_that("random curve reduced-rank scores include all moving-loading terms",{
  formulas <- list(
    rank1=~vsm(rrm(environment,k=1),ism(subject),sigma2=0.04),
    rank2=~vsm(rrm(environment,k=2),ism(subject),sigma2=0.04),
    kronecker=~vsm(ar1m(occasion),rrm(environment,k=1),ism(subject),sigma2=0.04),
    twoRR=~vsm(rrm(occasion,k=1),rrm(environment,k=1),ism(subject),sigma2=0.04))
  for(formula in formulas){
    fixture <- nonlinearStructureContext(formula,curveRandom=TRUE)
    value <- fixture$evaluate(fixture$parameters)
    expect_true(value$valid,info=value$reason)
    if(!value$valid) next
    for(step in c(1e-4,1e-5)){
      expect_equal(as.numeric(value$gradient),
        nonlinearStructureDifference(fixture$evaluate,fixture$parameters,step),
        tolerance=7e-5)
    }
    expect_lt(length(value$latentMode),length(value$mode))
    expect_equal(as.numeric(value$fitted),4*exp(-0.6*exp(as.numeric(
      fixture$design %*% value$mode))*fixture$data$x),tolerance=1e-11)
    expect_equal(fixture$evaluate(fixture$parameters)$logLik,value$logLik,tolerance=0)
  }
})

test_that("reduced-rank conditional density and Laplace likelihood match latent-space oracles",{
  data <- nonlinearStructureData(groups=1L)
  fixture <- nonlinearStructureContext(
    ~vsm(rrm(environment,k=1,loadings=matrix(c(1,0.65,-0.4),3,1)),
      ism(subject),sigma2=0.04),data=data,curveRandom=TRUE)
  par <- fixture$parameters
  value <- fixture$evaluate(par)
  expect_true(value$valid,info=value$reason)
  loading <- c(1,0.65,-0.4)[data$environment]
  conditional <- function(v){
    mean <- 4*exp(-0.6*exp(loading*v)*data$x)
    sum((data$y-mean)^2)/(2*0.15^2)+v^2/(2*0.04)
  }
  reference <- optimize(conditional,c(-1,1),tol=1e-12)
  expect_equal(as.numeric(value$latentMode),reference$minimum,tolerance=1e-6)
  point <- 0.08
  probe <- sommer:::sommer_nonlinear_conditional_eval(fixture$context,par,point)
  expect_equal(probe$objective,conditional(point),tolerance=1e-11)
  for(step in c(1e-4,1e-5)){
    expect_equal(as.numeric(probe$score),
      -(conditional(point+step)-conditional(point-step))/(2*step),tolerance=2e-6)
    upper <- sommer:::sommer_nonlinear_conditional_eval(fixture$context,par,point+step)
    lower <- sommer:::sommer_nonlinear_conditional_eval(fixture$context,par,point-step)
    expect_equal(as.numeric(probe$hessian),-as.numeric(upper$score-lower$score)/(2*step),
      tolerance=2e-6)
  }
  curvature <- as.numeric(stats::optimHess(reference$minimum,conditional))
  oracle <- -reference$objective-nrow(data)/2*log(2*pi*0.15^2)-
    log(0.04)/2-log(curvature)/2
  expect_equal(value$logLik,oracle,tolerance=2e-5)
  logDensity <- function(v) vapply(v,function(v)
    -conditional(v)-nrow(data)/2*log(2*pi*0.15^2)-log(2*pi*0.04)/2,numeric(1))
  shift <- logDensity(reference$minimum)
  exact <- shift+log(integrate(function(v) exp(logDensity(v)-shift),-2,2,
    rel.tol=1e-10)$value)
  expect_lt(abs(value$logLik-exact),0.015)
})

test_that("FA and reduced-rank fits retain full-coordinate and population predictions",{
  data <- nonlinearStructureData(groups=12L)
  set.seed(738)
  loading <- c(1,0.65,-0.4)
  deviation <- loading[data$environment]*rnorm(12L,sd=0.2)[data$subject]
  data$y <- 4*exp(-0.6*exp(deviation)*data$x)+rnorm(nrow(data),sd=0.12)
  relationship <- 0.9*diag(12)+0.1
  dimnames(relationship) <- list(levels(data$subject),levels(data$subject))
  precision <- solve(relationship)
  attr(precision,"inverse") <- TRUE
  formulas <- list(
    rr=~vsm(rrm(environment,k=1,loadings=matrix(c(1,0.65,-0.4),3,1)),
      ism(subject),Gu=precision,sigma2=0.04,fixedSigma2=TRUE),
    fa=~vsm(fam(environment,k=1),ism(subject),Gu=precision,sigma2=0.04,
      fixedSigma2=TRUE))
  control <- mmes_control(task=list(objective="ml",nonlinear=list(
    start=c(asymptote=4,rate=0.6))),
    convergence=list(maxIters=400L,likelihoodTol=1e-8,parameterTol=1e-7),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  for(name in names(formulas)){
    if(name=="fa"){
      loadings <- matrix(c(0.6,0.2,-0.15),3,1)
      sigma <- tcrossprod(loadings)+diag(0.5,3)
      sigma <- 0.04*sigma/sigma[1,1]
      deviations <- as.numeric(t(chol(kronecker(sigma,relationship))) %*% rnorm(36))
      local <- deviations[(as.integer(data$environment)-1L)*12L+
        as.integer(data$subject)]
      data$y <- 4*exp(-0.6*exp(local)*data$x)+rnorm(nrow(data),sd=0.12)
    }
    fit <- mmes(y~0,data=data,nonlinear=nonlinear_curve("x",
      random=list(rate=formulas[[name]])),control=control)
    expect_true(fit$convergence,info=name)
    expect_equal(unname(coef(fit)),c(4,0.6),tolerance=0.06)
    expect_equal(predict(fit,newdata=data),fitted(fit),tolerance=1e-10)
    expect_equal(predict(fit,newdata=data,re.form=NA),predict(fit,re.form=NA))
    expect_equal(predict(fit,newdata=data,re.form="population"),
      predict(fit,re.form="population"),tolerance=1e-10)
    expect_length(fit$groupParameters[,1],36L)
    expect_length(fit$randomMode,36L)
    expect_length(fit$latentRandomMode,if(name=="rr") 12L else 36L)
    expect_false(fit$nonlinearDiagnostics$uncertaintyComputed)
    expect_null(fit$Ci)
    rows <- c(2,13,22)
    conditional <- predict(fit,newdata=data[rows,])
    population <- predict(fit,newdata=data[rows,],re.form="population")
    sigma <- sommer:::evaluate_covstruct_cpp(fit$covStruct[[1]],fit$covStruct[[1]]$par)
    for(index in seq_along(rows)){
      row <- rows[index]
      variance <- sigma[as.integer(data$environment[row]),as.integer(data$environment[row])]
      oracle <- integrate(function(b) fit$nonlinearParameters[1]*
        exp(-fit$nonlinearParameters[2]*exp(b)*data$x[row])*
        dnorm(b,sd=sqrt(variance)),-5,5,rel.tol=1e-9)$value
      expect_equal(population[index],oracle,tolerance=3e-7)
    }
    restored <- unserialize(serialize(fit,NULL))
    captured <- precision
    precision <- NULL
    expect_equal(predict(restored,newdata=data[rows,]),conditional,tolerance=1e-11)
    expect_equal(predict(restored,newdata=data[rows,],re.form="population"),
      population,tolerance=1e-11)
    precision <- captured
  }
})

test_that("additive reduced-rank modes and covariance constraints preserve existing semantics",{
  data <- nonlinearStructureData(groups=8L)
  data$cluster <- factor(rep(1:2,length.out=nrow(data)))
  control <- mmes_control(task=list(objective="ml",nonlinear=list(
    start=c(asymptote=4,rate=0.6))),convergence=list(maxIters=400L),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  fit <- mmes(y~0,random=~vsm(rrm(environment,k=1),ism(subject),sigma2=0.04),
    nonlinear=nonlinear_curve("x"),data=data,control=control)
  expect_true(fit$convergence)
  expect_identical(fit$likelihoodMethod,"exact Gaussian ML")
  expect_length(fit$randomMode,24L)
  expect_length(fit$latentRandomMode,8L)
  expect_equal(predict(fit,newdata=data),fitted(fit),tolerance=1e-10)
  expect_equal(predict(fit,re.form="population"),predict(fit,re.form=NA))
  saved <- unserialize(serialize(fit,NULL))
  expect_equal(predict(saved,newdata=data[1:3,]),fitted(fit)[1:3])
  fixture <- nonlinearStructureContext(
    ~vsm(rrm(environment,k=1),ism(subject),sigma2=0.04),
    data=data,curveRandom=TRUE,additive=~vsm(ism(cluster),sigma2=0.02))
  mapping <- sommer:::.laplace_covariance_map(fixture$structures,
    fixture$prepared$rtermss,
    data.frame(parameter=c(1L,2L),group=1L,scale=c(1,2)))
  layout <- list(fixedCount=2L,covarianceCount=length(mapping$map$initial),
    covarianceMap=mapping$map)
  par <- c(log(4),log(0.6),mapping$map$initial,log(0.15^2))
  evaluate <- function(par) sommer:::.laplace_constrained_evaluate(par,layout,
    fixture$evaluate)
  value <- evaluate(par)
  expect_true(value$valid,info=value$reason)
  expect_equal(as.numeric(value$gradient),nonlinearStructureDifference(evaluate,par),
    tolerance=8e-5)
  additiveRR <- nonlinearStructureContext(~vsm(ism(subject),sigma2=0.04),
    data=data,curveRandom=TRUE,
    additive=~vsm(rrm(environment,k=1),ism(subject),sigma2=0.02))
  value <- additiveRR$evaluate(additiveRR$parameters)
  expect_true(value$valid,info=value$reason)
  expect_equal(as.numeric(value$gradient),
    nonlinearStructureDifference(additiveRR$evaluate,additiveRR$parameters),
    tolerance=8e-5)
})

test_that("fixed loading maps, invalid latent ranks and resource limits are explicit",{
  fixed <- nonlinearStructureContext(
    ~vsm(rrm(environment,k=1,fixed=rep(TRUE,3)),ism(subject),
      sigma2=0.04,fixedSigma2=TRUE),curveRandom=TRUE)
  expect_length(fixed$parameters,3L)
  value <- fixed$evaluate(fixed$parameters)
  expect_true(value$valid,info=value$reason)
  expect_equal(as.numeric(value$gradient),
    nonlinearStructureDifference(fixed$evaluate,fixed$parameters),tolerance=6e-5)
  fixture <- nonlinearStructureContext(
    ~vsm(rrm(environment,k=2),ism(subject),sigma2=0.04),curveRandom=TRUE)
  par <- fixture$parameters
  # Rank-two leading pivot is the fourth loading after the scale coordinate.
  position <- grep("loading_diag[2,F2]",fixture$structures[[1]]$par_names,fixed=TRUE)
  expect_length(position,1L)
  free <- which(fixture$structures[[1]]$free)
  par[2L+match(position,free)] <- -1000
  value <- fixture$evaluate(par)
  expect_false(value$valid)
  expect_match(value$reason,"loading pivot")
  expect_error(nonlinearStructureContext(
    ~vsm(rrm(environment,k=1),ism(subject)),curveRandom=TRUE,memoryMB=0.001),
    "workspace")
})

test_that("latent loadings compose with other curve parameters and relationship priors",{
  data <- nonlinearStructureData()
  relationship <- 0.8*diag(3)+0.2
  dimnames(relationship) <- list(levels(data$subject),levels(data$subject))
  precision <- solve(relationship)
  attr(precision,"inverse") <- TRUE
  formula <- ~vsm(rrm(environment,k=2),ism(subject),Gu=precision,sigma2=0.04)
  for(case in list(c("exponential","asymptote"),c("logistic","midpoint"),
      c("logistic","rate"),c("saturation","half"))){
    data$y <- switch(case[1],exponential=4*exp(-0.6*data$x),
      logistic=4*plogis(1.2*(data$x-2.5)),saturation=4*data$x/(0.6+data$x))+
      rep(c(-0.02,0.01,0.04),length.out=nrow(data))
    fixture <- nonlinearStructureContext(formula,data=data,curveRandom=TRUE,
      curve=case[1],parameter=case[2])
    value <- fixture$evaluate(fixture$parameters)
    expect_true(value$valid,info=value$reason)
    expect_equal(as.numeric(value$gradient),
      nonlinearStructureDifference(fixture$evaluate,fixture$parameters),
      tolerance=7e-5)
  }
})
