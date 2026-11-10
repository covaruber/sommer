gaussianBlockControl <- function(setup=FALSE){
  sommer:::.mmes_control_internal(task=list(objective="ml",glmm=list(method="laplace",
    modeTol=1e-10,maxit=1000L)),output=list(setup=setup),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

gaussianBlockData <- function(){
  set.seed(738)
  data <- expand.grid(time=factor(1:3),id=factor(1:24))
  data$x <- rnorm(nrow(data))
  data$trial <- factor(rep(1:2,each=36))
  R <- .7*.4^abs(outer(1:3,1:3,`-`))
  error <- as.vector(t(matrix(rnorm(nrow(data)),ncol=3) %*% chol(R)))
  data$y <- .4+.2*data$x+rnorm(24,sd=.8)[data$id]+error
  data
}

gaussianBlockObjective <- function(fit,memoryMB=1000){
  context <- sommer:::sommer_laplace_objective_create(fit$fixedDesign,
    fit$randomDesign,fit$laplaceNativeResponse,fit$offset,"",
    fit$covStruct,fit$randomRelationships,memoryMB,
    familySpecs=fit$familySpecifications,familyIndex=fit$laplaceFamilyIndex,
    trials=fit$priorWeights,residualModel=fit$laplaceResidual)
  function(parameters,gradient=TRUE){
    scope <- sommer:::.mmes_blas_scope(fit$control$execution$blasThreads)
    on.exit(scope$restore(),add=TRUE)
    sommer:::.laplace_constrained_evaluate(parameters,fit$laplaceParameterLayout,
      function(native) sommer:::sommer_laplace_objective_eval(context,native,
        tolerance=1e-10,gradient=gradient))
  }
}

test_that("structured Gaussian Laplace matches exact ML and joint coefficient covariance", {
  data <- gaussianBlockData()
  fit <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),
    data=data,family=gaussian(),control=gaussianBlockControl())
  expect_true(fit$convergence)
  expect_equal(fit$npar,5L)
  expect_equal(nrow(fit$vcParams),3L)
  expect_equal(fit$vcParams$term[2:3],rep("residual",2))
  expect_equal(fit$executionPlan$glmm$residual$observations,nrow(data))
  expect_equal(fit$executionPlan$glmm$residual$covarianceParameters,2L)
  expect_identical(fit$executionPlan$glmm$residual$likelihood,
    "joint Gaussian observation block")
  expect_null(fit$familyParameters)
  expect_null(fit$fixedCovariance)
  expect_null(fit$Ci)
  expect_false("hessian" %in% names(fit$optimizer))
  objective <- gaussianBlockObjective(fit)
  independent <- function(parameters){
    native <- sommer:::.laplace_expand_parameters(parameters,fit$laplaceParameterLayout)
    variance <- exp(native[3])
    sigma <- exp(native[4])
    rho <- tanh(native[5])
    R <- kronecker(diag(24),sigma*rho^abs(outer(1:3,1:3,`-`)))
    Z <- as.matrix(fit$randomDesign)
    V <- variance*tcrossprod(Z)+R
    U <- chol(V)
    residual <- data$y-as.numeric(fit$fixedDesign %*% parameters[1:2])
    standardized <- forwardsolve(t(U),residual)
    -.5*(nrow(data)*log(2*pi)+2*sum(log(diag(U)))+sum(standardized^2))
  }
  proposal <- fit$optimizer$par+c(.1,-.05,.08,.04,-.1)
  evaluation <- objective(proposal)
  expect_true(evaluation$valid)
  expect_equal(evaluation$logLik,independent(proposal),tolerance=1e-10)
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    expect_equal(as.numeric(evaluation$gradient)[index],
      (independent(upper)-independent(lower))/2e-4,tolerance=2e-5)
  }
  expect_equal(fit$llik,independent(fit$optimizer$par),tolerance=1e-10)
  expect_lt(max(abs(objective(fit$optimizer$par)$gradient)),.01)
  reference <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),data=data,
    control=sommer:::.mmes_control_internal(task=list(objective="ml"),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
  expect_true(reference$convergence)
  expect_equal(as.numeric(fit$b),as.numeric(reference$b),tolerance=1e-4)
  expect_equal(unname(unlist(fit$covPar)),unname(unlist(reference$covPar)),
    tolerance=2e-4)
  R <- sommer:::.laplace_residual_matrix(fit)
  P <- solve(R)
  conditional <- sommer:::sommer_glmm_family_eval(fit$y,fit$linear.predictors,
    fit$familySpecifications,fit$laplaceFamilyIndex,fit$priorWeights,
    residualModel=fit$laplaceResidual)
  expect_equal(as.numeric(conditional$score),
    as.numeric(P %*% (fit$y-fit$linear.predictors)),tolerance=1e-10)
  expect_equal(as.matrix(conditional$blockCurvature),P,tolerance=1e-10)
  expect_equal(sum(conditional$logDensityRows),conditional$logDensity)
  X <- as.matrix(fit$fixedDesign)
  Z <- as.matrix(fit$randomDesign)
  Q <- diag(1/fit$covPar[[1]][1],ncol(Z))
  design <- cbind(X,Z)
  information <- crossprod(design,P %*% design)
  information[3:ncol(design),3:ncol(design)] <-
    information[3:ncol(design),3:ncol(design)]+Q
  expected <- solve(information)
  expect_equal(unname(as.matrix(postPEV(fit,2L)$Ci)),unname(expected),tolerance=1e-8)
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-10)
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(covparams_mmes_se(fit)$StdError)))
  for(random in list(NULL,NA,"population")){
    prediction <- predict(fit,newdata=data[1:3,],se.fit=TRUE,cov.fit=TRUE,
      re.form=random,type="response")
    expect_true(all(is.finite(prediction$se.fit)))
  }
  expect_equal(predict(unserialize(serialize(fit,NULL))),predict(fit))
  expect_equal(sommer:::.boot_mmes_refit(fit,fit$y,1000L)$logLik,fit$llik,
    tolerance=1e-6)
})

test_that("Gaussian W and block metadata preserve filtering and covariance semantics", {
  data <- gaussianBlockData()
  W <- Matrix::bdiag(replicate(24,matrix(c(2,.1,0,.1,1,.1,0,.1,1.5),3),simplify=FALSE))
  data$y[2] <- NA
  data$x[8] <- NA
  fit <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time,rho=.3),ism(id)),
    W=W,weights=~id,data=data,family=gaussian(),control=gaussianBlockControl())
  expect_true(fit$convergence)
  expect_identical(which(!fit$obsInfo$included),c(2L,8L))
  expect_equal(fit$laplaceResidual$originalRows,which(fit$obsInfo$included))
  expect_equal(nrow(fit$laplaceResidual$W),nrow(data)-2L)
  retained <- data[fit$obsInfo$included,,drop=FALSE]
  explicit <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time,rho=.3),ism(id)),
    W=W[fit$obsInfo$included,fit$obsInfo$included],weights=~id,
    data=retained,family=gaussian(),control=gaussianBlockControl())
  expect_equal(fit$llik,explicit$llik,tolerance=1e-8)
  expect_equal(fitted(fit),fitted(explicit),tolerance=1e-7)
  base <- fit$laplaceResidual
  base$weighted <- FALSE
  base$W <- Matrix::Diagonal(nrow(base$W))
  covariance <- sommer:::sommer_glmm_residual_covariance(base)
  inverse <- solve(chol(as.matrix(fit$laplaceResidual$W)))
  expect_equal(sommer:::.laplace_residual_matrix(fit),
    t(inverse) %*% covariance %*% inverse,tolerance=1e-10)
  evaluate <- gaussianBlockObjective(fit)
  proposal <- fit$optimizer$par+.04
  gradient <- evaluate(proposal)$gradient
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    expect_equal(as.numeric(gradient)[index],
      (evaluate(upper,FALSE)$logLik-evaluate(lower,FALSE)$logLik)/2e-4,
      tolerance=2e-5)
  }
  bad <- W
  bad[1,4] <- bad[4,1] <- .1
  expect_error(sommer:::.mmes_dispatch(y~x,random=~id,W=bad,weights=~id,data=data,
    family=gaussian(),control=gaussianBlockControl()),"across groups")
})

test_that("W alone activates a Gaussian residual descriptor without changing the model scale", {
  data <- gaussianBlockData()
  independent <- sommer:::.mmes_dispatch(y~x,random=~id,data=data,family=gaussian(),
    control=gaussianBlockControl())
  weighted <- sommer:::.mmes_dispatch(y~x,random=~id,data=data,family=gaussian(),
    W=2*Matrix::Diagonal(nrow(data)),weights=~id,control=gaussianBlockControl())
  expect_true(weighted$convergence)
  expect_null(independent$laplaceResidual)
  expect_null(weighted$familyParameters)
  expect_equal(weighted$npar,independent$npar)
  expect_equal(weighted$llik,independent$llik,tolerance=1e-8)
  expect_lt(max(abs(as.numeric(weighted$b)-as.numeric(independent$b))),1e-6)
  expect_equal(diag(sommer:::.laplace_residual_matrix(weighted)),
    rep(independent$dispersion,nrow(data)),tolerance=1e-5)
  expect_equal(postPEV(weighted)$uPevList,postPEV(independent)$uPevList,
    tolerance=1e-6)
})

test_that("Gaussian block simulation and residual diagnostics use the fitted joint covariance", {
  data <- gaussianBlockData()
  fit <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),
    data=data,family=gaussian(),control=gaussianBlockControl())
  R <- sommer:::.laplace_residual_matrix(fit)
  expect_equal(residuals(fit,"pearson"),(fit$y-fitted(fit))/sqrt(diag(R)))
  whitened <- as.numeric(forwardsolve(t(chol(R)),fit$y-fitted(fit)))
  expect_equal(residuals(fit,"deviance"),whitened)
  expect_equal(residuals(fit,"quantile"),residuals(fit,"pearson"),tolerance=1e-8)
  expect_equal(residuals(fit,"working"),residuals(fit,"response"))
  simulations <- as.matrix(simulate(fit,nsim=4000L,seed=19))
  expect_lt(max(abs(rowMeans(simulations)[1:3]-fitted(fit)[1:3])),.07)
  expect_lt(max(abs(cov(t(simulations[1:3,]))-R[1:3,1:3])),.07)
  expect_identical(simulate(fit,nsim=2,seed=19),simulate(fit,nsim=2,seed=19))
  simulations <- as.matrix(simulate(fit,nsim=4000L,seed=19,re.form="population"))
  expect_lt(max(abs(cov(t(simulations[1:3,]))-
    (R[1:3,1:3]+fit$covPar[[1]][1]))),.1)
  lowMemory <- fit
  lowMemory$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(vcov(lowMemory),"memory|Memory")
  expect_error(gaussianBlockObjective(fit,1e-6),"Gaussian observation residual workspace")
  evaluate <- gaussianBlockObjective(fit,.65)
  expect_true(evaluate(fit$optimizer$par,FALSE)$valid)
  expect_error(evaluate(fit$optimizer$par),"Gaussian residual derivative workspace")
  invalid <- fit$laplaceResidual
  invalid$rows[2] <- invalid$rows[1]
  expect_error(sommer:::sommer_glmm_family_eval(fit$y,fit$linear.predictors,
    fit$familySpecifications,fit$laplaceFamilyIndex,fit$priorWeights,
    residualModel=invalid),"Gaussian residual rows are invalid")
})

test_that("Gaussian trait blocks estimate cross-trait residual covariance without duplicate dispersions", {
  set.seed(140)
  records <- 30L
  R <- matrix(c(.8,.3,.3,.5),2)
  error <- matrix(rnorm(2*records),records,2) %*% chol(R)
  effect <- rnorm(10,sd=.7)[rep(1:10,each=3)]
  data <- data.frame(y=as.vector(error+effect),
    trait=factor(rep(c("first","second"),each=records)),
    record=factor(rep(seq_len(records),2)),
    group=factor(rep(rep(1:10,each=3),2)))
  fit <- sommer:::.mmes_dispatch(y~trait-1,random=~group,
    rcov=~vsm(usm(trait),ism(record)),data=data,
    family=familym("trait",first=gaussian(),second=gaussian()),
    control=gaussianBlockControl())
  expect_true(fit$convergence)
  expect_null(fit$familyParameters)
  expect_equal(length(fit$optimizer$par),fit$npar)
  expect_equal(fit$npar,6L)
  covariance <- sommer:::.laplace_residual_matrix(fit)
  expect_gt(abs(covariance[1,records+1]),1e-5)
  expect_equal(covariance[1,2],0)
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(covparams_mmes_se(fit)$StdError)))
  evaluation <- gaussianBlockObjective(fit)
  proposal <- fit$optimizer$par+.04
  native <- evaluation(proposal)
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    difference <- (evaluation(upper,FALSE)$logLik-evaluation(lower,FALSE)$logLik)/2e-4
    expect_equal(as.numeric(native$gradient)[index],difference,tolerance=2e-5)
  }
  simulations <- as.matrix(simulate(fit,nsim=4000L,seed=19))
  selected <- c(1,records+1)
  expect_lt(max(abs(cov(t(simulations[selected,]))-covariance[selected,selected])),.07)
  expect_equal(sommer:::.boot_mmes_refit(fit,fit$y,1000L)$logLik,fit$llik,
    tolerance=1e-6)
})

test_that("Gaussian residual blocks coexist with scalar, ordinal and nominal likelihoods", {
  set.seed(81)
  original <- expand.grid(time=factor(1:3),group=factor(1:16))
  original$x <- rnorm(nrow(original))
  random <- rnorm(16,sd=.5)[original$group]
  original$y <- .3+.2*original$x+random+
    as.vector(t(matrix(rnorm(nrow(original)),ncol=3) %*%
      chol(.6*.35^abs(outer(1:3,1:3,`-`)))))
  data <- do.call(rbind,lapply(c("gauss","count","ord","nom"),function(trait){
    result <- original
    result$trait <- trait
    result$y <- switch(trait,gauss=result$y,
      count=rpois(nrow(result),exp(.5+.2*result$x+random)),
      ord=as.integer(cut(.2*result$x+random+rlogis(nrow(result)),c(-Inf,-.6,.7,Inf))),
      nom=sample.int(3,nrow(result),replace=TRUE))
    result
  }))
  data$trait <- factor(data$trait,levels=c("gauss","count","ord","nom"))
  family <- familym("trait",gauss=gaussian(),count=poisson(),ord=ordinal(),
    nom=multinomial(reference="3"))
  fit <- sommer:::.mmes_dispatch(y~trait+trait:x-1,
    random=~vsm(ism(group),sigma2=.25,fixedSigma2=TRUE),
    rcov=~vsm(ar1m(time),ism(group)),family=family,data=data,
    control=gaussianBlockControl())
  expect_true(fit$convergence)
  expect_equal(length(fit$laplaceResidual$rows),nrow(original))
  expect_equal(fit$familySpecifications$gauss$free,FALSE)
  expect_equal(fit$laplaceParameterLayout$covarianceCount,2L)
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(unlist(residuals(fit,"pearson")))))
  expect_true(all(is.finite(unlist(residuals(fit,"deviance")))))
  for(mode in c("conditional","full")){
    prediction <- predict(fit,newdata=data[c(1,49,97,145),],type="response",
      se.fit=TRUE,cov.fit=TRUE,uncertainty=mode)
    expect_true(all(is.finite(unlist(prediction$se.fit))))
  }
  expect_equal(postPEV(fit)$uPevList,postPEV(fit,2L)$uPevList,tolerance=1e-7)
  expect_equal(sommer:::.boot_mmes_refit(fit,data$y,1000L)$logLik,fit$llik,
    tolerance=1e-6)
  simulations <- simulate(fit,nsim=2,seed=9)
  expect_equal(dim(simulations),c(nrow(data),2L))
  expect_true(all(is.finite(as.matrix(simulations))))
  expect_identical(simulate(fit,nsim=2,seed=9),simulations)
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(vcov(restored),vcov(fit))
  analyticData <- data[data$trait!="nom",,drop=FALSE]
  set.seed(7381)
  countRows <- analyticData$trait=="count"
  analyticData$y[countRows] <- rnbinom(sum(countRows),size=.8,
    mu=exp(.5+.2*analyticData$x[countRows]+random))
  analyticFamily <- familym("trait",gauss=gaussian(),count=nbinom2(2),ord=ordinal())
  analytic <- sommer:::.mmes_dispatch(y~trait+trait:x-1,
    random=~vsm(ism(group),sigma2=.25,fixedSigma2=TRUE),
    rcov=~vsm(ar1m(time),ism(group)),family=analyticFamily,data=analyticData,
    control=gaussianBlockControl())
  evaluate <- gaussianBlockObjective(analytic)
  proposal <- analytic$optimizer$par+.04
  gradient <- evaluate(proposal)$gradient
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    difference <- (evaluate(upper,FALSE)$logLik-evaluate(lower,FALSE)$logLik)/2e-4
    expect_equal(as.numeric(gradient)[index],difference,tolerance=3e-5)
  }
  incomplete <- data
  gaussianRows <- which(incomplete$trait=="gauss")
  incomplete$time[-gaussianRows] <- NA
  incomplete$time[2] <- NA
  incomplete$x[50] <- NA
  W <- Matrix::Diagonal(nrow(incomplete))
  W[gaussianRows,gaussianRows] <- 2*Matrix::Diagonal(length(gaussianRows))
  prepare <- function(weight)
    sommer:::.mmes_dispatch(y~trait+trait:x-1,random=~group,
      rcov=~vsm(ar1m(time),ism(group)),family=family,data=incomplete,
      W=weight,weights=~group,control=gaussianBlockControl(TRUE))
  setup <- prepare(W)
  expect_equal(length(setup$laplaceResidual$rows),length(gaussianRows)-1L)
  gaussianOnly <- prepare(W[gaussianRows,gaussianRows])
  expect_equal(setup$laplaceResidual,gaussianOnly$laplaceResidual)
  retainedRows <- setdiff(seq_len(nrow(incomplete)),c(2,50))
  retained <- prepare(W[retainedRows,retainedRows])
  expect_equal(setup$laplaceResidual,retained$laplaceResidual)
  retainedGaussian <- setdiff(gaussianRows,2)
  gaussianRetained <- prepare(W[retainedGaussian,retainedGaussian])
  expect_equal(setup$laplaceResidual,gaussianRetained$laplaceResidual)
  W[1,49] <- W[49,1] <- .1
  expect_error(prepare(W),"cross-family")
})

test_that("fixed Gaussian residual descriptors give exact GLS uncertainty and nested tests", {
  data <- gaussianBlockData()
  random <- ~vsm(ism(id),sigma2=.8,fixedSigma2=TRUE)
  residual <- ~vsm(ar1m(time,rho=.35,fixed=TRUE),ism(id),
    sigma2=.6,fixedSigma2=TRUE)
  fit <- sommer:::.mmes_dispatch(y~x,random=random,rcov=residual,family=gaussian(),data=data,
    control=gaussianBlockControl())
  smaller <- sommer:::.mmes_dispatch(y~1,random=random,rcov=residual,family=gaussian(),data=data,
    control=gaussianBlockControl())
  expect_true(fit$convergence)
  expect_equal(fit$npar,2L)
  expect_equal(fit$laplaceParameterLayout$covarianceCount,0L)
  X <- as.matrix(fit$fixedDesign)
  Z <- as.matrix(fit$randomDesign)
  V <- .8*tcrossprod(Z)+sommer:::.laplace_residual_matrix(fit)
  expected <- solve(crossprod(X,solve(V,X)))
  expect_equal(unname(vcov(fit)),unname(expected),tolerance=1e-6)
  expect_true(all(is.na(covparams_mmes_se(fit)$StdError)))
  expect_equal(anova(smaller,fit)$`Chi Df`[2],1L)
  changed <- fit
  changed$laplaceResidual$W <- 2*changed$laplaceResidual$W
  expect_error(anova(smaller,changed),"identical response")
  bootstrap <- boot_mmes(fit,nsim=2,seed=41)
  expect_equal(bootstrap$successful,2L)
  expect_equal(ncol(bootstrap$draws),2L)
  expect_equal(sommer:::.boot_mmes_refit(fit,fit$y,1000L)$logLik,fit$llik,
    tolerance=1e-8)
  section <- sommer:::.mmes_dispatch(y~x,random=random,
    rcov=~dsumm(vsm(ar1m(time,rho=.35,fixed=TRUE),ism(id),
      sigma2=.6,fixedSigma2=TRUE),by=trial),
    data=data,family=gaussian(),control=gaussianBlockControl())
  heterogeneous <- sommer:::.mmes_dispatch(y~x,random=random,
    rcov=~vsm(dsm(trial),ar1m(time,rho=.35,fixed=TRUE),ism(id),
      sigma2=.6,fixedSigma2=TRUE),
    data=data,family=gaussian(),control=gaussianBlockControl())
  expect_true(section$convergence)
  expect_equal(heterogeneous$llik,section$llik,tolerance=1e-8)
  expect_equal(sommer:::.laplace_residual_matrix(heterogeneous),
    sommer:::.laplace_residual_matrix(section),tolerance=1e-10)
  expect_equal(sommer:::.boot_mmes_refit(section,section$y,1000L)$logLik,
    section$llik,tolerance=1e-8)
  estimated <- sommer:::.mmes_dispatch(y~x,random=random,
    rcov=~dsumm(vsm(ar1m(time,rho=.35,fixed=TRUE),ism(id),
      sigma2=.6),by=trial),data=data,family=gaussian(),
    control=gaussianBlockControl())
  expect_true(estimated$convergence)
  expect_equal(sum(estimated$vcParams$free),2L)
  expect_true(all(is.finite(covparams_mmes_se(estimated)$StdError[
    covparams_mmes_se(estimated)$free])))
  evaluate <- gaussianBlockObjective(estimated)
  proposal <- estimated$optimizer$par+.04
  gradient <- evaluate(proposal)$gradient
  for(index in seq_along(proposal)){
    upper <- lower <- proposal
    upper[index] <- upper[index]+1e-4
    lower[index] <- lower[index]-1e-4
    expect_equal(as.numeric(gradient)[index],
      (evaluate(upper,FALSE)$logLik-evaluate(lower,FALSE)$logLik)/2e-4,
      tolerance=2e-5)
  }
})

test_that("Gaussian residual blocks reject incompatible inputs and retain random vcc constraints", {
  data <- gaussianBlockData()
  data$block <- factor(rep(1:6,length.out=nrow(data)))
  fit <- sommer:::.mmes_dispatch(y~x,random=~id+block,rcov=~vsm(ar1m(time),ism(id)),
    vcc=data.frame(parameter=1:2,group=1,scale=c(2,1)),
    family=gaussian(),data=data,control=gaussianBlockControl())
  expect_equal(fit$covPar[[1]],2*fit$covPar[[2]],tolerance=1e-10)
  expect_equal(fit$laplaceParameterLayout$covarianceCount,3L)
  evaluate <- gaussianBlockObjective(fit)
  expect_equal(evaluate(fit$optimizer$par)$logLik,fit$llik)
  setup <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),
    family=gaussian(),data=data,control=gaussianBlockControl(TRUE))
  expect_equal(setup$vcParams$term[2:3],rep("residual",2))
  expect_error(sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),
    family=gaussian(),data=data,
    control=sommer:::.mmes_control_internal(task=list(objective="ml",glmm=list(method="laplace",
      familyParameters=list(dispersion=1,estimate=FALSE))),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE))),"own their dispersion")
  asymmetric <- diag(nrow(data))
  asymmetric[1,2] <- .1
  expect_error(sommer:::.mmes_dispatch(y~x,random=~id,W=asymmetric,data=data,family=gaussian(),
    control=gaussianBlockControl()),"symmetric")
  expect_error(sommer:::.mmes_dispatch(y~x,random=~id,W=-diag(nrow(data)),data=data,family=gaussian(),
    control=gaussianBlockControl()),"positive definite")
  expect_error(sommer:::.mmes_dispatch(y~x,random=~id,weights=~id,data=data,family=gaussian(),
    control=gaussianBlockControl()),"provide W")
  data$time[2] <- NA
  filtered <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),data=data,
    family=gaussian(),control=gaussianBlockControl())
  expect_equal(which(!filtered$obsInfo$included),2L)
  expect_equal(filtered$obsInfo$reason[2],"missing residual coordinate")
})

test_that("ambiguous W subset dimensions cannot silently change observation alignment", {
  data <- gaussianBlockData()
  data$trait <- factor(rep(c("gauss","count"),each=36))
  data$y[data$trait=="count"] <- round(abs(data$y[data$trait=="count"]))+1
  data$x[-c(1:33,37:39)] <- NA
  prepare <- function(W) sommer:::.mmes_dispatch(y~trait+x,random=~id,
    rcov=~vsm(ar1m(time),ism(id)),data=data,W=W,
    family=familym("trait",gauss=gaussian(),count=poisson()),
    control=gaussianBlockControl(TRUE))
  expect_error(prepare(Matrix::Diagonal(36)),"ambiguously match")
  expect_equal(length(prepare(Matrix::Diagonal(72))$laplaceResidual$rows),33L)
})

test_that("Gaussian residual post-fit operations scope and restore saved BLAS controls", {
  configured <- 4L
  updates <- integer()
  controller <- list(provider="mock",get=function() configured,
    set=function(value){ configured <<- value; updates <<- c(updates,value) })
  local_mocked_bindings(.mmes_blas_controller=function() controller,.package="sommer")
  control <- gaussianBlockControl()
  control$execution$blasThreads <- 1L
  fit <- sommer:::.mmes_dispatch(y~x,random=~id,rcov=~vsm(ar1m(time),ism(id)),
    data=gaussianBlockData(),family=gaussian(),control=control)
  expect_identical(configured,4L)
  operations <- list(vcov,postPEV,function(object) simulate(object,nsim=2,seed=1),
    function(object) sommer:::.boot_mmes_refit(object,object$y,1000L))
  for(operation in operations){
    updates <- integer()
    operation(fit)
    expect_identical(configured,4L)
    expect_gt(length(updates),0L)
    expect_identical(tail(updates,1L),4L)
  }
  fit$laplacePostPEVSettings$denseMemoryMB <- 1e-6
  expect_error(vcov(fit),"memory|Memory")
  expect_identical(configured,4L)
})
