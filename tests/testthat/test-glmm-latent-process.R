latentProcessControl <- function(setup=FALSE,parameters=list()){
  sommer:::.mmes_control_internal(task=list(objective="ml",glmm=list(method="laplace",
    maxit=1000L,modeTol=1e-10,familyParameters=parameters)),
    output=list(setup=setup),diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

latentProcessData <- function(){
  set.seed(917)
  data <- expand.grid(time=factor(1:3),id=factor(1:16),replicate=1:4)
  data$x <- rnorm(nrow(data))
  covariance <- .4*.55^abs(outer(1:3,1:3,`-`))
  latent <- as.vector(t(chol(covariance)) %*% matrix(rnorm(48),3))
  data$eta <- .3+.25*data$x+latent[as.integer(data$time)+3*(as.integer(data$id)-1)]
  data$y <- rpois(nrow(data),exp(data$eta))
  data
}

latentProcessEvaluate <- function(fit,parameters){
  context <- sommer:::sommer_laplace_objective_create(fit$fixedDesign,
    fit$randomDesign,fit$laplaceNativeResponse,fit$offset,"",
    fit$covStruct,fit$randomRelationships,1000,
    familySpecs=fit$familySpecifications,familyIndex=fit$laplaceFamilyIndex,
    trials=fit$priorWeights)
  sommer:::.laplace_constrained_evaluate(parameters,fit$laplaceParameterLayout,
    function(native) sommer:::sommer_laplace_objective_eval(context,native,
      tolerance=1e-10,gradient=TRUE))
}

latentProcessPrior <- function(fit){
  pieces <- lapply(seq_along(fit$covStruct),function(term)
    kronecker(sommer:::evaluate_covstruct_cpp(fit$covStruct[[term]],
      fit$covStruct[[term]]$par),solve(as.matrix(fit$randomRelationships[[term]]))))
  as.matrix(Matrix::bdiag(pieces))
}

test_that("temporal Poisson processes use the existing likelihood and covariance descriptors", {
  data <- latentProcessData()
  missing <- which((data$time=="2" & data$id=="1") | seq_len(nrow(data))==70L)
  data$y[missing] <- NA
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55),ism(id),sigma2=.4),
    data=data,family=poisson(),control=latentProcessControl())
  expect_true(fit$convergence)
  expect_equal(which(!fit$obsInfo$included),missing)
  expect_null(fit$Ci)
  expect_null(fit$fixedCovariance)
  expect_null(fit$uPevList)
  expect_false("hessian" %in% names(fit$optimizer))
  expect_lt(max(abs(latentProcessEvaluate(fit,fit$optimizer$par)$gradient)),.01)
  X <- as.matrix(fit$fixedDesign)
  Z <- as.matrix(fit$randomDesign)
  Q <- solve(latentProcessPrior(fit))
  beta <- as.numeric(fit$b)
  objective <- function(u){
    eta <- as.numeric(X %*% beta+Z %*% u)
    -sum(dpois(fit$y,exp(eta),log=TRUE))+.5*sum(u*(Q %*% u))
  }
  gradient <- function(u) as.numeric(crossprod(Z,exp(X %*% beta+Z %*% u)-fit$y)+Q %*% u)
  reference <- nlminb(rep(0,ncol(Z)),objective,gradient,
    control=list(rel.tol=1e-12,iter.max=1000))
  expect_equal(reference$convergence,0L)
  expect_equal(as.numeric(fit$randomMode),reference$par,tolerance=1e-5)
  rho <- tanh(fit$covStruct[[1]]$par[2])
  expect_equal(fit$randomMode[17],unname(rho/(1+rho^2)*sum(fit$randomMode[c(1,33)])),
    tolerance=1e-8)
  expect_equal(as.numeric(predict(fit,newdata=data[2,],type="link")),
    as.numeric(fit$b[1]+fit$b[2]*data$x[2]+fit$randomMode[17]),tolerance=1e-8)
  eta <- as.numeric(X %*% beta+Z %*% reference$par)
  H <- Q+crossprod(Z,exp(eta)*Z)
  expected <- -objective(reference$par)+.5*as.numeric(determinant(Q,logarithm=TRUE)$modulus)-
    .5*as.numeric(determinant(H,logarithm=TRUE)$modulus)
  expect_equal(fit$llik,expected,tolerance=1e-8)
  proposal <- fit$optimizer$par+.04
  analytic <- latentProcessEvaluate(fit,proposal)$gradient
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    expect_equal(as.numeric(analytic)[index],
      (latentProcessEvaluate(fit,upper)$logLik-latentProcessEvaluate(fit,lower)$logLik)/2e-4,
      tolerance=2e-5)
  }
  prior <- Z %*% latentProcessPrior(fit) %*% t(Z)
  fixed <- predict(fit,re.form=NA,type="link")
  expect_equal(predict(fit,re.form="population",type="response"),
    exp(fixed+diag(prior)/2),tolerance=1e-10)
  retained <- data[fit$obsInfo$included,,drop=FALSE]
  prediction <- predict(fit,newdata=retained[1:4,],re.form="population",
    type="response",se.fit=TRUE,cov.fit=TRUE)
  expect_true(all(is.finite(prediction$se.fit)))
  expect_equal(as.numeric(diag(prediction$cov.fit)),as.numeric(prediction$se.fit^2))
  expect_true(all(is.finite(covparams_mmes_se(fit)$StdError)))
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-8)
  design <- cbind(X,Z)
  information <- crossprod(design,as.numeric(fitted(fit))*design)
  positions <- ncol(X)+seq_len(ncol(Z))
  information[positions,positions] <- information[positions,positions]+Q
  expect_equal(unname(as.matrix(postPEV(fit,2L)$Ci)),
    unname(solve(information)),tolerance=1e-8)
  expect_identical(fit$executionPlan$glmm$latent$scale,"link")
  expect_identical(fit$executionPlan$glmm$latent$uncertainty,"post-fit")
  expect_equal(predict(unserialize(serialize(fit,NULL))),predict(fit))
  expect_equal(sommer:::.boot_mmes_refit(fit,fit$y,1000)$logLik,fit$llik,tolerance=1e-7)
  population <- as.matrix(simulate(fit,nsim=4000,seed=11,re.form="population"))
  expect_lt(max(abs(rowMeans(population)-exp(fixed+diag(prior)/2))/
    (1+exp(fixed+diag(prior)/2))),.1)
  expect_lt(abs(cov(t(population[1:2,]))[1,2]-
    exp(fixed[1]+fixed[2]+sum(diag(prior)[1:2])/2)*(exp(prior[1,2])-1)),.2)
})

test_that("latent covariance panels handle arbitrary loadings, correlations and row groups", {
  data <- latentProcessData()
  setup <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    sigma2=.4,fixedSigma2=TRUE),data=data,family=poisson(),
    control=latentProcessControl(TRUE))
  descriptor <- setup$covStruct[1]
  precision <- .7^abs(outer(1:16,1:16,`-`))
  precision <- Matrix::Matrix(precision,sparse=TRUE)
  set.seed(18)
  design <- Matrix::Matrix(matrix(rnorm(8*48),8,48),sparse=TRUE)
  G <- kronecker(sommer:::evaluate_covstruct_cpp(descriptor[[1]],descriptor[[1]]$par),
    solve(as.matrix(precision)))
  expected <- as.matrix(design %*% G %*% t(design))
  groups <- matrix(1:8,ncol=2,byrow=TRUE)
  result <- sommer:::.laplace_latent_covariance(design,descriptor,list(precision),
    48L,groups,memoryMB=.03)
  for(index in 1:4)
    expect_equal(result[,,index],expected[groups[index,],groups[index,]],tolerance=1e-10)
  expect_equal(sommer:::.laplace_population_variance(design,descriptor,
    list(precision),48L),diag(expected),tolerance=1e-10)
  expect_error(sommer:::.laplace_latent_covariance(design,descriptor,list(precision),
    48L,memoryMB=1e-6),"workspace|denseMemoryMB")
  expect_error(sommer:::.laplace_latent_covariance(design,descriptor,list(precision),
    47L),"align")
})

test_that("binomial and spatial count processes support conditional and population targets", {
  data <- latentProcessData()
  set.seed(192)
  data$y <- rbinom(nrow(data),1,plogis(data$eta))
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    sigma2=.4,fixedSigma2=TRUE),data=data,family=binomial(),
    control=latentProcessControl())
  expect_true(fit$convergence)
  reference <- vapply(predict(fit,re.form=NA,type="link"),function(eta)
    integrate(function(z) plogis(eta+sqrt(.4)*z)*dnorm(z),-9,9,
      rel.tol=1e-10)$value,numeric(1))
  expect_equal(predict(fit,re.form="population",type="response"),reference,tolerance=1e-8)
  known <- data[1:3,]
  expect_equal(predict(fit,newdata=known),predict(fit)[1:3],tolerance=1e-8)
  known$id <- factor("new")
  expect_error(predict(fit,newdata=known),"unseen")
  expect_equal(predict(fit,newdata=known,re.form="population",type="response"),
    reference[1:3],tolerance=1e-8)
  known$time <- factor("4")
  expect_error(predict(fit,newdata=known,re.form="population",type="response"),
    "unseen latent coordinates")
  spatial <- expand.grid(row=factor(1:2),column=factor(1:3),
    site=factor(1:8),replicate=1:3)
  spatial$x <- rnorm(nrow(spatial))
  spatial$y <- rpois(nrow(spatial),exp(.2+.3*spatial$x))
  model <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(row,rho=.3,fixed=TRUE),
    ar1m(column,rho=.6,fixed=TRUE),ism(site),sigma2=.25,fixedSigma2=TRUE),
    data=spatial,family=poisson(),control=latentProcessControl())
  expect_true(model$convergence)
  expect_equal(predict(model,re.form="population",type="response"),
    exp(predict(model,re.form=NA,type="link")+.25/2),tolerance=1e-10)
  expect_true(all(is.finite(predict(model,newdata=spatial[1:3,],
    type="response",se.fit=TRUE)$se.fit)))
})

test_that("ordinal latent processes integrate the fitted correlated prior", {
  data <- latentProcessData()
  set.seed(812)
  data$y <- ordered(cut(data$eta+rnorm(nrow(data)),c(-Inf,-.5,.7,Inf)),
    labels=c("low","medium","high"))
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    sigma2=.4,fixedSigma2=TRUE),data=data,family=ordinal("probit"),
    control=latentProcessControl(parameters=list(thresholds=c(-.5,.7),estimate=FALSE)))
  expect_true(fit$convergence)
  eta <- predict(fit,re.form=NA,type="link")
  cumulative <- cbind(pnorm((-.5-eta)/sqrt(1.4)),pnorm((.7-eta)/sqrt(1.4)))
  expected <- cbind(cumulative[,1],cumulative[,2]-cumulative[,1],1-cumulative[,2])
  expect_equal(unname(predict(fit,re.form="population")),unname(expected),tolerance=1e-10)
  prediction <- predict(fit,newdata=data[1:2,],re.form="population",se.fit=TRUE,
    cov.fit=TRUE)
  expect_true(all(is.finite(prediction$se.fit)))
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-8)
})

test_that("nominal processes integrate row-specific correlated logit loadings", {
  skip_on_cran()
  data <- latentProcessData()
  set.seed(87)
  probability <- cbind(exp(data$eta),exp(-.4*data$eta),1)
  data$y <- factor(apply(probability,1,function(p) sample.int(3,1,prob=p)))
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),
    ism(.nominalCategory),ism(id),sigma2=.2,fixedSigma2=TRUE),
    data=data,family=multinomial(),control=latentProcessControl())
  expect_true(fit$convergence)
  known <- data[1:3,]
  prediction <- predict(fit,newdata=known,re.form="population",se.fit=TRUE,cov.fit=TRUE)
  expect_true(all(is.finite(prediction$se.fit)))
  expect_equal(as.numeric(rowSums(prediction$fit)),rep(1,3),tolerance=1e-12)
  expect_equal(as.numeric(diag(prediction$cov.fit)),as.numeric(prediction$se.fit^2))
  zero <- as.numeric(predict(fit,newdata=data[1,],re.form=NA,type="link"))
  rule <- sommer:::.laplace_normal_quadrature(40L)
  grid <- expand.grid(1:40,1:40)
  logits <- cbind(zero[1]+sqrt(.2)*rule$nodes[grid[,1]],
    zero[2]+sqrt(.2)*rule$nodes[grid[,2]])
  probabilities <- sommer:::.laplace_nominal_probabilities(logits,fit$nominal)
  weights <- rule$weights[grid[,1]]*rule$weights[grid[,2]]
  derivative <- matrix(0,3,2)
  for(category in 1:3) for(coordinate in 1:2)
    derivative[category,coordinate] <- sum(weights*probabilities[,category]*
      ((category==coordinate)-probabilities[,coordinate]))
  jacobian <- derivative %*% as.matrix(fit$fixedDesign[c(1,nrow(data)+1),])
  expected <- jacobian %*% vcov(fit) %*% t(jacobian)
  expect_equal(unname(prediction$cov.fit[c(1,4,7),c(1,4,7)]),
    unname(expected),tolerance=1e-6)
  expect_equal(predict(fit,newdata=known),fitted(fit)[1:3,,drop=FALSE],tolerance=1e-8)
  expect_equal(predict(unserialize(serialize(fit,NULL)),newdata=known,
    re.form="population"),prediction$fit)
  known$time <- factor("4")
  expect_error(predict(fit,newdata=known,re.form="population"),"unseen latent coordinates")
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-8)
})

test_that("mixed count, binary, ordinal and nominal responses share a temporal process", {
  skip_on_cran()
  original <- latentProcessData()
  original <- droplevels(original[original$id %in% as.character(1:8),])
  set.seed(881)
  data <- do.call(rbind,lapply(c("count","binary","ord","nom"),function(trait){
    result <- original
    result$trait <- trait
    result$y <- switch(trait,count=original$y,
      binary=rbinom(nrow(original),1,plogis(original$eta)),
      ord=as.integer(cut(original$eta+rlogis(nrow(original)),c(-Inf,-.5,.7,Inf))),
      nom=sample.int(3,nrow(original),replace=TRUE))
    result
  }))
  data$trait <- factor(data$trait,levels=c("count","binary","ord","nom"))
  family <- familym("trait",count=poisson(),binary=binomial(),ord=ordinal(),
    nom=multinomial())
  fit <- sommer:::.mmes_dispatch(y~trait+trait:x-1,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),
    ism(id),sigma2=.4,fixedSigma2=TRUE),data=data,family=family,
    control=latentProcessControl(parameters=list(ord=list(thresholds=c(-.5,.7),
      estimate=FALSE))))
  expect_true(fit$convergence)
  rows <- c(1,2,97,98,193,194,289,290)
  prediction <- predict(fit,newdata=data[rows,],re.form="population",
    se.fit=TRUE,cov.fit=TRUE)
  expect_true(all(is.finite(unlist(prediction$se.fit))))
  expect_identical(prediction$uncertainty,"full")
  expect_equal(as.numeric(diag(prediction$cov.fit)),
    as.numeric(unlist(prediction$se.fit)^2),tolerance=1e-8)
  zero <- predict(fit,newdata=data[rows,],re.form=NA,type="link")
  expect_equal(as.numeric(prediction$fit$count),
    exp(as.numeric(zero$count)+.4/2),tolerance=1e-10)
  for(trait in c("ord","nom"))
    expect_equal(as.numeric(rowSums(prediction$fit[[trait]])),rep(1,2),tolerance=1e-12)
  nominal <- vapply(1:3,function(category){
    integrate(function(u){
      eta <- matrix(rep(as.numeric(zero$nom[1,]),each=length(u)),ncol=2)+u
      sommer:::.laplace_nominal_probabilities(eta,fit$familySpecifications$nom)[,category]*
        dnorm(u,sd=sqrt(.4))
    },-8,8,rel.tol=1e-9)$value
  },numeric(1))
  expect_lt(max(abs(as.numeric(prediction$fit$nom[1,])-nominal)),1e-6)
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-8)
  expect_equal(sommer:::.boot_mmes_refit(fit,data$y,1000)$logLik,fit$llik,tolerance=1e-6)
  expect_identical(simulate(fit,nsim=2,seed=5,re.form="population"),
    simulate(fit,nsim=2,seed=5,re.form="population"))
})

test_that("relationship precision and random slopes contribute to population means and SE", {
  data <- latentProcessData()
  data$id <- factor(data$id,levels=as.character(1:17))
  precision <- .3^abs(outer(1:17,1:17,`-`))
  dimnames(precision) <- list(levels(data$id),levels(data$id))
  attr(precision,"inverse") <- TRUE
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    Gu=precision,sigma2=.2,fixedSigma2=TRUE),
    data=data,family=poisson(),control=latentProcessControl())
  expect_true(fit$convergence)
  variance <- diag(as.matrix(fit$randomDesign %*% latentProcessPrior(fit) %*%
    t(fit$randomDesign)))
  expect_equal(predict(fit,re.form="population",type="response"),
    exp(predict(fit,re.form=NA,type="link")+variance/2),tolerance=1e-10)
  expect_true(all(is.finite(predict(fit,newdata=data[1:2,],
    re.form="population",type="response",se.fit=TRUE)$se.fit)))
  represented <- data[1,]
  represented$id <- factor("17",levels=levels(data$id))
  expect_equal(as.numeric(predict(fit,newdata=represented,type="link")),
    as.numeric(fit$b[1]+fit$b[2]*represented$x+fit$randomMode[17]),tolerance=1e-8)
  expect_equal(fit$randomMode[17],
    -sum(precision[17,1:16]*fit$randomMode[1:16])/precision[17,17],tolerance=1e-8)
  expect_equal(as.numeric(predict(fit,newdata=represented,re.form="population",
    type="response")),
    exp(as.numeric(fit$b[1]+fit$b[2]*represented$x)+.2*solve(precision)[17,17]/2),
    tolerance=1e-10)
  unknown <- data[1,]
  unknown$id <- factor("new")
  expect_error(predict(fit,newdata=unknown,re.form="population",type="response"),
    "relationship levels")
  slope <- sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.2,fixedSigma2=TRUE)+
    vsm(ism(x),ism(id),sigma2=.1,fixedSigma2=TRUE),data=data,family=poisson(),
    control=latentProcessControl())
  expect_true(slope$convergence)
  variance <- diag(as.matrix(slope$randomDesign %*% latentProcessPrior(slope) %*%
    t(slope$randomDesign)))
  expect_equal(predict(slope,re.form="population",type="response"),
    exp(predict(slope,re.form=NA,type="link")+variance/2),tolerance=1e-10)
  unknown <- data[1:2,]
  unknown$id <- factor("new")
  expect_true(all(is.finite(predict(slope,newdata=unknown,
    re.form="population",type="response",se.fit=TRUE)$se.fit)))
  data$replicate <- factor(data$replicate)
  combined <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    Gu=precision,sigma2=.2,fixedSigma2=TRUE)+
    vsm(ism(replicate),sigma2=.1,fixedSigma2=TRUE),
    data=data,family=poisson(),control=latentProcessControl())
  unknown <- data[1,]
  unknown$replicate <- factor("new")
  expect_equal(as.numeric(predict(combined,newdata=unknown,re.form="population",
    type="response",se.fit=TRUE)$fit),
    exp(as.numeric(combined$b[1]+combined$b[2]*unknown$x)+
      (.2*solve(precision)[1,1]+.1)/2),tolerance=1e-10)
})

test_that("PQL retains its existing structured latent-process workflow", {
  data <- latentProcessData()
  fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
    sigma2=.4,fixedSigma2=TRUE),data=data,family=poisson(),
    control=sommer:::.mmes_control_internal(task=list(objective="ml",pql=list(maxit=100L,tol=1e-7)),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_true(fit$convergence)
  expect_true(all(is.finite(fitted(fit))))
  expect_true(all(is.finite(unlist(postPEV(fit)$uPevList))))
  expect_true(all(is.finite(predict(fit,D="id",VarU=FALSE)$pvals$predicted.value)))
})

test_that("other log-linked families preserve their marginal mean and post-fit contracts", {
  skip_on_cran()
  data <- latentProcessData()
  set.seed(561)
  responses <- list(nb=rnbinom(nrow(data),size=2,mu=exp(data$eta)),
    gamma=rgamma(nrow(data),shape=2,scale=exp(data$eta)/2),
    zip=ifelse(rbinom(nrow(data),1,.25),0,rpois(nrow(data),exp(data$eta))))
  families <- list(nb=nbinom2(2),gamma=Gamma("log"),zip=zipoisson(.25))
  settings <- list(nb=list(theta=2,estimate=FALSE),
    gamma=list(dispersion=.5,estimate=FALSE),
    zip=list(zeroInflation=.25,estimate=FALSE))
  for(name in names(responses)){
    data$y <- responses[[name]]
    fit <- sommer:::.mmes_dispatch(y~x,random=~vsm(ar1m(time,rho=.55,fixed=TRUE),ism(id),
      sigma2=.2,fixedSigma2=TRUE),data=data,family=families[[name]],
      control=latentProcessControl(parameters=settings[[name]]))
    expect_true(fit$convergence)
    expected <- exp(predict(fit,re.form=NA,type="link")+.2/2)
    if(name=="zip") expected <- expected*.75
    expect_equal(predict(fit,re.form="population",type="response"),
      expected,tolerance=1e-10)
    prediction <- predict(fit,newdata=data[1:3,],re.form="population",
      type="response",se.fit=TRUE,cov.fit=TRUE)
    expect_true(all(is.finite(prediction$se.fit)))
    expect_equal(as.numeric(diag(prediction$cov.fit)),as.numeric(prediction$se.fit^2))
  }
})
