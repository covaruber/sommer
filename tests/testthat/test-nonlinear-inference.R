nonlinearInferenceControl <- function(...){
  sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
    start=c(asymptote=1.4,rate=0.65))),
    convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE),...)
}

nonlinearInferenceData <- function(family=gaussian(),randomCurve=FALSE,seed=715){
  set.seed(seed)
  data <- data.frame(x=rep(seq(0,5,length.out=16),12),
    group=factor(rep(1:12,each=16)))
  u <- rep(rnorm(12,sd=0.2),each=16)
  eta <- if(randomCurve) 1.4*exp(-0.65*exp(u)*data$x) else
    1.4*exp(-0.65*data$x)+u
  data$z <- rep(c(-1,1),length.out=nrow(data))
  data$y <- switch(if(grepl("^Negative Binomial",family$family)) "nbinom2" else family$family,
    gaussian=eta+rnorm(nrow(data),sd=0.15),
    binomial=rbinom(nrow(data),1,family$linkinv(eta)),
    poisson=rpois(nrow(data),exp(eta)),
    nbinom2=rnbinom(nrow(data),size=3,mu=exp(eta)),
    Gamma=rgamma(nrow(data),shape=4,scale=exp(eta)/4),
    ziPoisson={ value <- rpois(nrow(data),exp(eta)); value[runif(nrow(data))<0.25] <- 0; value },
    ordinal={
      p <- cbind(family$linkinv(-0.6-eta),
        family$linkinv(0.7-eta)-family$linkinv(-0.6-eta),1-family$linkinv(0.7-eta))
      code <- 1L+rowSums(runif(nrow(data))>t(apply(p,1,cumsum))[,-3,drop=FALSE])
      ordered(c("low","middle","high")[code],levels=c("low","middle","high"))
    },
    multinomial={
      p <- cbind(exp(eta),exp(eta-0.7),1); p <- p/rowSums(p)
      code <- 1L+rowSums(runif(nrow(data))>t(apply(p,1,cumsum))[,-3,drop=FALSE])
      factor(c("red","green","blue")[code],levels=c("red","green","blue"))
    })
  data
}

test_that("post-fit curvature, natural scales and PEV match independent Gaussian algebra", {
  data <- nonlinearInferenceData()
  fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  expect_null(fit$nonlinearInference)
  expect_false(fit$nonlinearDiagnostics$uncertaintyComputed)
  parameters <- fit$optimizer$par
  z <- as.matrix(fit$randomDesign)
  dense <- function(at){
    mean <- exp(at[1])*exp(-exp(at[2])*data$x)
    v <- exp(at[3])*tcrossprod(z)+diag(exp(at[4]),nrow(data))
    r <- chol(v)
    residual <- data$y-mean
    0.5*(nrow(data)*log(2*pi)+2*sum(log(diag(r))))+
      0.5*sum(residual*solve(v,residual))
  }
  expected <- solve(optimHess(parameters,dense,
    control=list(ndeps=rep(1e-4,length(parameters)))))
  post <- post_nonlinear(fit,mode=2)
  expect_equal(unname(vcov(post,type="working")),unname(expected),
    tolerance=2e-5,ignore_attr=TRUE)
  transform <- diag(exp(parameters))
  expect_equal(unname(vcov(post,type="all")),
    unname(transform%*%expected%*%transform),tolerance=2e-5,ignore_attr=TRUE)
  conditional <- solve(crossprod(z)/fit$dispersion+
    diag(1/exp(parameters[3]),ncol(z)))
  expect_equal(vcov(post,type="conditional"),conditional,tolerance=1e-9,ignore_attr=TRUE)
  mode <- function(at){
    h <- crossprod(z)/exp(at[4])+diag(1/exp(at[3]),ncol(z))
    solve(h,crossprod(z,data$y-exp(at[1])*exp(-exp(at[2])*data$x))/exp(at[4]))
  }
  motion <- vapply(seq_along(parameters),function(index){
    upper <- lower <- parameters
    upper[index] <- upper[index]+1e-5
    lower[index] <- lower[index]-1e-5
    as.numeric((mode(upper)-mode(lower))/2e-5)
  },numeric(ncol(z)))
  expect_equal(vcov(post,type="random"),conditional+motion%*%expected%*%t(motion),
    tolerance=2e-5,ignore_attr=TRUE)
  expect_true(all(diag(vcov(post,type="random"))>=diag(conditional)))
  expect_match(post$nonlinearInference$randomMethod,"not an exact frequentist PEV")
  expect_equal(unlist(post$uPevList),diag(vcov(post,type="random")),ignore_attr=TRUE)
  diagonal <- postPEV(fit)
  expect_equal(diagonal$nonlinearInference$randomPEV,diag(vcov(post,type="random")))
  expect_null(postPEV(post,0)$nonlinearInference)
  expect_false(postPEV(post,0)$nonlinearDiagnostics$uncertaintyComputed)
  expect_equal(vcov(fit),vcov(post))
  expect_equal(diag(vcov(post)),diag(vcov(post,type="working"))[1:2]*
    coef(post)^2,ignore_attr=TRUE)
  expect_equal(rownames(confint(post)),names(coef(post)))
  expect_true(all(confint(post)>0))
  expect_equal(nrow(confint(post,type="all")),length(coef(post,type="all")))
  expect_equal(nrow(confint(post,type="working")),length(coef(post,type="working")))
  expect_error(confint(post,parm=1.5),"integer")
  expect_equal(summary(post)$parameterTable$StdError,
    sqrt(diag(vcov(post,type="all"))),ignore_attr=TRUE)
  expect_null(summary(fit)$parameterTable)
  expect_true(all(is.finite(summary(fit,inference=TRUE)$parameterTable$StdError)))
  reported <- covparams_mmes_se(post)
  expect_equal(reported$estimate,covparams_mmes(post)$estimate)
  expect_true(all(reported$StdError>0))
  covarianceMatrix <- covmatrix_mmes(post,1)
  expect_equal(as.numeric(covarianceMatrix$covariance.se),reported$StdError)
  expect_equal(as.numeric(covarianceMatrix$correlation.se),0)
  ratio <- vcDelta(post,h2~V1/(V1+V2))
  variance <- coef(post,type="all")[3:4]
  derivative <- c(variance[2],-variance[1])/sum(variance)^2
  expect_equal(ratio$Estimate,unname(variance[1]/sum(variance)))
  expect_equal(ratio$SE,sqrt(as.numeric(crossprod(derivative,
    vcov(post,type="all")[3:4,3:4]%*%derivative))))
  expect_error(covparams_mmes_se(post,rel_step=1e-4),"inference\\$step")
  expect_equal(AIC(fit),-2*fit$llik+2*fit$npar)
  expect_equal(BIC(fit),-2*fit$llik+log(nobs(fit))*fit$npar)
  expect_equal(AIC(fit,post)$AIC,rep(AIC(fit),2))
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(vcov(restored),vcov(post))
  expect_null(fit$nonlinearInference)
})

test_that("prediction delta covariance includes mode motion and keeps estimands distinct", {
  skip_on_cran()
  data <- nonlinearInferenceData()
  fit <- post_nonlinear(sommer:::.mmes_dispatch(y~0,random=~group,data=data,
    nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl()),mode=2)
  newdata <- data[c(1,5,20),]
  prediction <- predict(fit,newdata,se.fit=TRUE,cov.fit=TRUE,interval="confidence")
  expect_equal(prediction$fit,predict(fit,newdata),tolerance=1e-10)
  expect_equal(prediction$se.fit^2,diag(prediction$cov.fit),ignore_attr=TRUE)
  p <- fit$optimizer$par
  direct <- cbind(fit$nonlinearParameters[1]*exp(-exp(p[2])*newdata$x),
    -fit$nonlinearParameters[1]*exp(-exp(p[2])*newdata$x)*exp(p[2])*newdata$x,
    0,0)
  z <- as.matrix(fit$randomDesign[c(1,5,20),,drop=FALSE])
  full <- direct+z%*%fit$nonlinearInference$modeSensitivity
  expected <- full%*%vcov(fit,type="working")%*%t(full)+
    z%*%vcov(fit,type="conditional")%*%t(z)
  expect_equal(prediction$cov.fit,expected,tolerance=1e-7,ignore_attr=TRUE)
  zero <- predict(fit,newdata,re.form=NA,se.fit=TRUE,cov.fit=TRUE)
  population <- predict(fit,newdata,re.form="population",se.fit=TRUE,cov.fit=TRUE)
  expect_equal(zero,population)
  expect_equal(zero$cov.fit,direct%*%vcov(fit,type="working")%*%t(direct),
    tolerance=1e-7,ignore_attr=TRUE)
  future <- predict(fit,newdata,interval="prediction")
  confidence <- predict(fit,newdata,interval="confidence")
  expect_true(all(future$upper-future$lower>confidence$upper-confidence$lower))
  expect_equal(((future$upper-future$lower)/(2*qnorm(.975)))^2,
    prediction$se.fit^2+fit$dispersion,tolerance=1e-8)
  expect_equal(predict(fit,newdata,re.form="population",interval="prediction")$fit,
    population$fit)
  expect_error(predict(fit,newdata,se.fit=NA),"se.fit")
  expect_error(predict(fit,newdata,interval="confidence",level=1),"level")
  expect_error(predict(fit,newdata,re.form="bad",se.fit=TRUE),"re.form")
})

for(family in list(binomial(),binomial("probit"),binomial("cloglog"),poisson(),
    nbinom2(3),Gamma("log"),gaussian(),zipoisson(.25))){
  test_that(paste(family$family,family$link,"exposes finite post-fit inference"), {
    skip_on_cran()
    data <- nonlinearInferenceData(family)
    control <- nonlinearInferenceControl()
    if(grepl("^Negative Binomial",family$family))
      control$task$glmm$familyParameters <- list(theta=3,estimate=FALSE)
    if(family$family=="Gamma")
      control$task$glmm$familyParameters <- list(dispersion=.25,estimate=FALSE)
    if(family$family=="ziPoisson")
      control$task$glmm$familyParameters <- list(zeroInflation=.25,estimate=FALSE)
    fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
      data=data,family=family,nonlinear=nonlinear_curve("x"),control=control)
    post <- post_nonlinear(fit,mode=2)
    expect_true(all(is.finite(vcov(post,type="all"))),info=family$family)
    expect_equal(coef(post,type="all"),post$nonlinearInference$parameters)
    for(reform in list(NULL,NA,"population")){
      prediction <- predict(post,newdata=data[1:3,],re.form=reform,se.fit=TRUE,cov.fit=TRUE)
      expect_equal(prediction$fit,predict(post,newdata=data[1:3,],re.form=reform),
        tolerance=1e-8)
      expect_true(all(is.finite(prediction$se.fit)))
      expect_equal(prediction$se.fit^2,diag(prediction$cov.fit),ignore_attr=TRUE)
    }
    future <- predict(post,newdata=data[1:2,],re.form="population",interval="prediction")
    expect_true(all(is.finite(future$lower)))
    expect_true(all(future$upper>future$lower))
  })
}

test_that("random curve group uncertainty and prior-population prediction are post-fit", {
  data <- nonlinearInferenceData(poisson(),TRUE)
  fit <- sommer:::.mmes_dispatch(y~0,data=data,family=poisson(),
    nonlinear=nonlinear_curve("x",random=list(rate=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE))),
    control=nonlinearInferenceControl())
  post <- post_nonlinear(fit,mode=2)
  expect_equal(dim(post$nonlinearInference$groupParameterSE),dim(fit$groupParameters))
  expect_true(all(post$nonlinearInference$groupParameterSE>0))
  conditional <- predict(post,data[1:3,],interval="confidence")
  population <- predict(post,data[1:3,],re.form="population",interval="confidence")
  zero <- predict(post,data[1:3,],re.form=NA,interval="confidence")
  expect_gt(max(abs(population$fit-zero$fit)),1e-5)
  expect_true(all(conditional$se.fit>0))
  expect_true(all(population$se.fit>0))
  expect_null(fit$nonlinearInference)
})

for(family in list(ordinal(),ordinal("probit"),multinomial(reference="red"))){
  test_that(paste(family$family,family$link,"retains categorical uncertainty shapes"), {
    skip_on_cran()
    data <- nonlinearInferenceData(family)
    fixed <- if(family$family=="ordinal") y~1 else y~0
    fit <- sommer:::.mmes_dispatch(fixed,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
      data=data,family=family,nonlinear=nonlinear_curve("x"),
      control=nonlinearInferenceControl())
    post <- post_nonlinear(fit)
    for(reform in list(NULL,NA,"population")){
      prediction <- predict(post,data[1:2,],re.form=reform,se.fit=TRUE,cov.fit=TRUE,
        interval="confidence")
      expect_equal(dim(prediction$se.fit),dim(prediction$fit))
      expect_equal(colnames(prediction$se.fit),levels(data$y))
      expect_equal(as.numeric(prediction$se.fit)^2,diag(prediction$cov.fit),
        ignore_attr=TRUE)
      expect_equal(prediction$fit,predict(fit,data[1:2,],re.form=reform),tolerance=1e-8)
      expect_equal(unname(rowSums(prediction$fit)),rep(1,2),tolerance=1e-12)
      expect_lt(max(abs(rowSums(matrix(rowSums(prediction$cov.fit),nrow=2)))),1e-8)
    }
    expect_error(predict(post,data[1:2,],interval="prediction"),"categorical")
    if(family$family=="ordinal"){
      rows <- grep("^family:",names(post$nonlinearInference$parameters))
      expect_true(all(diag(vcov(post,type="all"))[rows]>0))
    }
  })
}

test_that("zero-regression retains uncertainty shapes", {
  skip_on_cran()
  data <- nonlinearInferenceData(zipoisson(.25))
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
    data=data,family=zipoisson(zero=~z),nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  post <- post_nonlinear(fit)
  expect_true(any(grepl("family:",names(post$nonlinearInference$parameters))))
  expect_true(all(is.finite(predict(post,data[1:3,],se.fit=TRUE)$se.fit)))
})

test_that("mixed families retain uncertainty shapes", {
  skip_on_cran()
  left <- nonlinearInferenceData(poisson())
  right <- nonlinearInferenceData(gaussian(),seed=730)
  left$trait <- "count"; right$trait <- "normal"
  data <- rbind(left,right)
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),data=data,
    family=familym("trait",count=poisson(),normal=gaussian()),
    nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl())
  post <- post_nonlinear(fit)
  expect_true(all(is.finite(vcov(post,type="all"))))
  expect_true(all(is.finite(predict(post,data[c(1,193),],se.fit=TRUE)$se.fit)))
})

test_that("inference rejects invalid curvature, controls and uncalibrated tests", {
  data <- nonlinearInferenceData()
  fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  expect_error(post_nonlinear(fit,mode=1.5),"mode")
  expect_error(post_nonlinear(fit,control=list()),"mmes_control")
  bad <- fit; bad$convergence <- FALSE
  expect_error(vcov(bad),"converged")
  bad <- fit; bad$nonlinearDiagnostics$identified <- FALSE
  expect_error(postPEV(bad),"identifiable")
  bad <- fit; bad$nonlinearDiagnostics$boundaryParameters <- "1:sigma2"
  expect_error(vcov(bad),"boundar")
  bad <- fit; bad$llik <- fit$llik+1
  expect_error(vcov(bad),"differs")
  control <- nonlinearInferenceControl(execution=list(denseMemoryMB=.0001))
  expect_error(post_nonlinear(fit,control=control),"workspace")
  expect_error(nonlinearInferenceControl(task=list()),"matched by multiple")
  expect_error(sommer:::.mmes_control_internal(task=list(nonlinear=list(inference=list(step=.2)))),"step<")
  settings <- sommer:::.mme_nonlinear_inference_defaults
  expect_error(sommer:::.nonlinear_positive_inverse(diag(c(1,-1)),settings,"test"),
    "positive definite")
  expect_error(sommer:::.nonlinear_positive_inverse(matrix(1,2,2),settings,"test"),
    "singular")
  expect_error(sommer:::.nonlinear_positive_inverse(diag(1e-320,2),settings,"test"),
    "finite numerical range")
  expect_error(sommer:::.nonlinear_jacobian(function(at) sin(1e4*at),1,settings,"test"),
    "successive step sizes")
  expect_error(anova(fit),"explicit parm and null")
  test <- anova(fit,parm="curve:rate",null=log(.65))
  expect_true(is.finite(test$Chisq))
  expect_error(anova(fit,parm="residual:logVariance",null=0),"fixed/curve")
  expect_error(anova(fit,parm="curve:rate",null=-Inf),"finite")
  expect_error(confint(fit,parm="no"),"parm")
  expect_error(AIC(fit,k=NA),"k")
})

test_that("regular nested fixed comparisons preserve the nonlinear architecture", {
  data <- nonlinearInferenceData()
  data$y <- data$y+.25*data$z
  random <- ~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE)
  small <- sommer:::.mmes_dispatch(y~0,random=random,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  large <- sommer:::.mmes_dispatch(y~0+z,random=random,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  comparison <- anova(small,large)
  expect_equal(comparison$Chisq[2],2*(large$llik-small$llik))
  expect_equal(comparison[["Chi Df"]][2],1)
  expect_error(anova(large,small),"nested")
  altered <- large; altered$input <- altered$input+.01
  expect_error(anova(small,altered),"identical curves")
})

test_that("Gaussian curve intervals have reasonable simulated coverage", {
  covered <- matrix(FALSE,16,2)
  for(index in seq_len(nrow(covered))){
    data <- nonlinearInferenceData(seed=800+index)
    fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
      data=data,nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl())
    interval <- confint(fit)
    truth <- c(1.4,.65)
    covered[index,] <- interval[,1]<=truth & interval[,2]>=truth
  }
  expect_true(all(colMeans(covered)>=.75))
})

test_that("reduced-rank uncertainty restores the full singular covariance basis", {
  set.seed(113)
  data <- expand.grid(x=seq(0,5,length.out=10),environment=factor(1:3),
    subject=factor(1:8))
  loading <- c(1,.7,-.3)
  u <- rnorm(8,sd=.2)
  data$y <- 1.4*exp(-.65*data$x)+loading[data$environment]*u[data$subject]+
    rnorm(nrow(data),sd=.12)
  precision <- solve(.8*diag(8)+.2)
  dimnames(precision) <- list(levels(data$subject),levels(data$subject))
  attr(precision,"inverse") <- TRUE
  random <- ~vsm(rrm(environment,k=1,loadings=matrix(c(1,.7,-.3),3,1),
    fixed=rep(TRUE,3)),ism(subject),Gu=precision,sigma2=.04,fixedSigma2=TRUE)
  fit <- sommer:::.mmes_dispatch(y~0,random=random,data=data,nonlinear=nonlinear_curve("x"),
    control=nonlinearInferenceControl())
  post <- post_nonlinear(fit,mode=2)
  z <- as.matrix(fit$randomDesign)
  sigma <- sommer:::evaluate_covstruct_cpp(fit$covStruct[[1]],fit$covStruct[[1]]$par)
  g <- kronecker(sigma,solve(as.matrix(fit$randomRelationships[[1]])))
  observation <- z%*%g%*%t(z)+diag(fit$dispersion,nrow(data))
  expected <- g-g%*%t(z)%*%solve(observation,z%*%g)
  expect_equal(vcov(post,type="conditional"),expected,tolerance=1e-8,ignore_attr=TRUE)
  expect_equal(qr(vcov(post,type="conditional"))$rank,8L)
  expect_equal(nrow(vcov(post,type="random")),24L)
  expect_true(all(is.finite(predict(post,data[1:3,],se.fit=TRUE)$se.fit)))
  expect_true(all(is.finite(covparams_mmes_se(post)$StdError)))
  expect_true(all(is.finite(covmatrix_mmes(post,1)$covariance.se)))
  map <- sommer:::sommer_nonlinear_coordinate_map(fit$randomDesign,fit$covStruct,
    fit$randomRelationships)
  expect_equal(as.numeric(map%*%fit$latentRandomMode),fit$randomMode,tolerance=1e-12)
  old <- fit
  old$randomBlockWidths <- old$randomDesignBlockWidths
  expect_equal(post_nonlinear(old,mode=2)$nonlinearInference$randomPEV,
    post$nonlinearInference$randomPEV,tolerance=1e-10)
})

test_that("all curve parameterizations and estimated family tails have uncertainty", {
  skip_on_cran()
  data <- nonlinearInferenceData()
  for(curve in c("logistic","saturation")){
    data$y <- if(curve=="logistic") 1.4*plogis(1.2*(data$x-2.5)) else
      1.4*data$x/(.8+data$x)
    set.seed(81)
    data$y <- data$y+rep(rnorm(12,sd=.2),each=16)+rnorm(nrow(data),sd=.12)
    control <- nonlinearInferenceControl()
    control$task$nonlinear$start <- if(curve=="logistic")
      c(asymptote=1.4,rate=1.2,midpoint=2.5) else c(asymptote=1.4,half=.8)
    fit <- sommer:::.mmes_dispatch(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x",curve),
      control=control)
    expect_true(all(is.finite(vcov(fit))))
    expect_true(all(is.finite(confint(fit))))
  }
  for(family in list(nbinom2(3),Gamma("log"))){
    data <- nonlinearInferenceData(family)
    fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),
      data=data,family=family,nonlinear=nonlinear_curve("x"),
      control=nonlinearInferenceControl())
    covariance <- vcov(fit,type="all")
    rows <- grep("^family:",rownames(covariance))
    expect_length(rows,1)
    expect_true(all(diag(covariance)[rows]>0))
  }
})

test_that("grouped-binomial intervals require declared new-data trial counts", {
  skip_on_cran()
  data <- nonlinearInferenceData(binomial())
  set.seed(73)
  data$success <- rbinom(nrow(data),20,plogis(1.4*exp(-.65*data$x)))
  data$failure <- 20-data$success
  fit <- sommer:::.mmes_dispatch(cbind(success,failure)~0,
    random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),data=data,
    family=binomial(),nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl())
  post <- post_nonlinear(fit)
  expect_error(predict(post,data[1:3,],interval="prediction"),"require trials")
  result <- predict(post,data[1:3,],interval="prediction",trials=20)
  expect_true(all(is.finite(result$upper)))
  expect_error(predict(post,data[1:3,],interval="prediction",trials=0),"positive integers")
  expect_error(predict(post,data[1:3,],se.fit=TRUE,trials=20),"only to future")
  expect_true(all(is.finite(predict(post,interval="prediction")$lower)))
})

test_that("covariance constraints and mixed categorical uncertainty remain aligned", {
  skip_on_cran()
  data <- nonlinearInferenceData()
  data$cluster <- factor(rep(1:3,each=64))
  data$y <- data$y+rep(c(-.15,.2,-.05),each=64)
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04)+vsm(ism(cluster),sigma2=.04),
    data=data,vcc=data.frame(parameter=1:2,group=1L,scale=c(1,1)),
    nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl())
  post <- post_nonlinear(fit,mode=2)
  expect_equal(length(coef(post,type="working")),4L)
  covariance <- vcov(post,type="all")
  rows <- grep("^covariance:",rownames(covariance))
  expect_equal(covariance[rows[1],rows],rep(covariance[rows[1],rows[1]],2),
    ignore_attr=TRUE)
  expect_length(post$uPevList,2)
  expect_true(all(is.finite(covparams_mmes_se(post)$StdError)))
  left <- nonlinearInferenceData(ordinal())
  right <- nonlinearInferenceData(multinomial())
  left$y <- as.integer(left$y); right$y <- as.integer(right$y)
  left$trait <- "ordered"; right$trait <- "nominal"
  data <- rbind(left,right)
  control <- nonlinearInferenceControl()
  control$task$glmm$familyParameters <- list(
    ordered=list(thresholds=c(-.6,.7),estimate=FALSE))
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE),data=data,
    family=familym("trait",ordered=ordinal(),nominal=multinomial()),
    nonlinear=nonlinear_curve("x"),control=control)
  post <- post_nonlinear(fit)
  output <- predict(post,data[c(1,2,193,194),],se.fit=TRUE,cov.fit=TRUE,
    interval="confidence")
  expect_named(output$fit,c("ordered","nominal"))
  expect_equal(lapply(output$fit,dim),lapply(output$se.fit,dim))
  expect_equal(length(unlist(output$fit)),nrow(output$cov.fit))
  expect_equal(unname(unlist(output$se.fit))^2,diag(output$cov.fit),
    ignore_attr=TRUE)
})

test_that("Gaussian future-response intervals agree with conditional simulation", {
  data <- nonlinearInferenceData()
  post <- post_nonlinear(sommer:::.mmes_dispatch(y~0,random=~group,data=data,
    nonlinear=nonlinear_curve("x"),control=nonlinearInferenceControl()))
  interval <- predict(post,interval="prediction")
  draws <- as.matrix(simulate(post,nsim=400,seed=127))
  coverage <- mean(draws>=interval$lower & draws<=interval$upper)
  expect_gt(coverage,.9)
  expect_lt(coverage,.999)
})

test_that("Gaussian random curves combine additive effects and exact conditional curvature", {
  skip_on_cran()
  data <- nonlinearInferenceData(gaussian(),TRUE)
  data$cluster <- factor(rep(1:3,each=64))
  data$y <- data$y+rep(c(-.1,.12,-.02),each=64)
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(cluster),sigma2=.0144,fixedSigma2=TRUE),
    data=data,nonlinear=nonlinear_curve("x",
      random=list(rate=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE))),
    control=nonlinearInferenceControl())
  post <- post_nonlinear(fit,mode=2)
  expect_identical(fit$likelihoodMethod,"Laplace Gaussian ML")
  expect_true(all(is.finite(vcov(post,type="all"))))
  expect_length(post$uPevList,2)
  expect_true(all(is.finite(post$nonlinearInference$groupParameterSE)))
  for(reform in list(NULL,NA,"population")){
    result <- predict(post,data[1:3,],re.form=reform,interval="prediction")
    expect_equal(result$fit,predict(fit,data[1:3,],re.form=reform),tolerance=1e-8)
    expect_true(all(result$upper>result$lower))
  }
  old <- fit
  old$randomBlockWidths <- old$randomDesignBlockWidths
  expect_equal(post_nonlinear(old,mode=2)$nonlinearInference$randomPEV,
    post$nonlinearInference$randomPEV,tolerance=1e-10)
})
