specification_control <- function(){
  .mmes_control_build(task=list(objective="ml"),
    convergence=list(maxIters=200L,likelihoodTol=1e-8,parameterTol=1e-8),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

specification_data <- function(){
  set.seed(913)
  data <- expand.grid(x=seq(-1.8,1.8,length.out=7),
    Name=factor(letters[1:3]),Env=factor(c("A","B")))
  data$y <- 1+.3*(data$Env=="B")+1.7*tanh(.9*data$x)+
    rnorm(nrow(data),sd=.18)
  data
}

test_that("parameter declarations and formula markers are strict",{
  expect_s3_class(parameter(fixed=~Env,random=~vsm(ism(Name)),link="log"),
    "sommer_parameter")
  expect_error(parameter(fixed=y~x),"one-sided")
  expect_error(parameter(link=c("log","identity")),"link")
  expect_error(parameter(link="logit"),"link")
  expect_error(graph_term("net"),"formula marker")
  expect_error(graph_term(1),"predictor name")
})

test_that("ordinary fixed-only models match lm and glm exactly",{
  data <- specification_data()
  control <- specification_control()
  fit <- .mmes_specification_fit(y~Env+x,data=data,control=control)
  oracle <- lm(y~Env+x,data)
  expect_equal(coef(fit),coef(oracle))
  expect_equal(fitted(fit),fitted(oracle))
  expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)))
  expect_equal(vcov(fit),vcov(oracle)*(nrow(data)-3)/nrow(data))
  expect_equal(predict(fit,data[1:3,]),predict(oracle,data[1:3,]))
  expect_equal(residuals(fit),residuals(oracle))
  expect_equal(nobs(fit),nrow(data))
  expect_s3_class(summary(fit),"summary.mmes.specification")
  control$task$objective <- "reml"
  reml <- .mmes_specification_fit(y~Env+x,data=data,control=control)
  expect_equal(as.numeric(logLik(reml)),as.numeric(logLik(oracle,REML=TRUE)))
  expect_equal(vcov(reml),vcov(oracle))
  set.seed(711)
  data$z <- rbinom(nrow(data),1,plogis(.2+.8*data$x))
  fit <- .mmes_specification_fit(z~Env+x,data=data,family=binomial(),control=control)
  oracle <- glm(z~Env+x,data=data,family=binomial())
  expect_equal(coef(fit),coef(oracle))
  expect_equal(fitted(fit),fitted(oracle))
  expect_equal(predict(fit,data,type="link"),predict(oracle,data,type="link"))
  expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)))
  expect_equal(vcov(fit),vcov(oracle))
  expect_equal(residuals(fit,type="pearson"),residuals(oracle,type="pearson"))
  expect_equal(dim(simulate(fit,nsim=2,seed=91)),c(nrow(data),2L))
})

test_that("nonlinear parameter regressions agree with an independent nls oracle",{
  data <- specification_data()
  control <- specification_control()
  declarations <- list(a=parameter(),k=parameter(fixed=~Env,link="log"))
  fit <- .mmes_specification_fit(y~Env+a*tanh(k*x),data=data,
    parameters=declarations,start=c(a=1.5,k=1),control=control)
  oracle <- nls(y~b0+b1*(Env=="B")+a*tanh(exp(k0+k1*(Env=="B"))*x),
    data=data,start=list(b0=1,b1=.3,a=1.5,k0=0,k1=0))
  expect_true(fit$convergence)
  expect_equal(unname(coef(fit)),unname(coef(oracle)),tolerance=2e-4)
  expect_equal(fitted(fit),as.numeric(fitted(oracle)),tolerance=2e-5)
  expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)),tolerance=1e-7)
  expect_equal(predict(fit,data[1:4,]),fitted(fit)[1:4],tolerance=1e-10)
  expect_equal(predict(fit,target="zero"),fitted(fit))
  expect_match(summary(fit)$derivativeMethod,"numerical")
  expect_true(all(is.finite(vcov(fit))))
  explicit <- .mmes_specification_fit(y~Env+a*tanh(k*x),data=data,
    parameters=declarations,start=list(coefficients=coef(fit)),control=control)
  expect_equal(as.numeric(logLik(explicit)),as.numeric(logLik(fit)),tolerance=1e-7)
})

test_that("exact additive Gaussian integration preserves relationship DG covariance",{
  data <- specification_data()
  A <- matrix(c(1,.2,.1,.2,1,.3,.1,.3,1),3,3,
    dimnames=list(levels(data$Name),levels(data$Name)))
  Ai <- Matrix::Matrix(solve(A),sparse=TRUE)
  attr(Ai,"inverse") <- TRUE
  control <- specification_control()
  fit <- .mmes_specification_fit(y~Env+a*tanh(k*x),
    random=~vsm(dsm(Env),ism(Name),Gu=Ai,sigma2=.03,fixedSigma2=TRUE),
    data=data,parameters=list(a=parameter(),k=parameter(link="log")),
    start=c(a=1.5,k=1),control=control)
  # The independent oracle uses the retained common covariance descriptors,
  # but independently constructs and factors the marginal observation law.
  z <- fit$randomRecords[[1]]$z
  g <- fit$randomCovariances[[1]]
  zero <- predict(fit,target="zero",type="link")
  v <- z %*% g %*% t(z)+diag(fit$dispersion,nrow(data))
  r <- data$y-zero
  oracle <- -.5*(nrow(data)*log(2*pi)+as.numeric(determinant(v,TRUE)$modulus)+
    sum(r*solve(v,r)))
  expect_equal(as.numeric(logLik(fit)),oracle,tolerance=1e-8)
  expect_equal(g,kronecker(diag(diag(g)[c(1,4)]),A),tolerance=1e-9)
  expect_match(fit$likelihoodMethod,"exact")
  expect_equal(predict(fit,data[1:4,]),fitted(fit)[1:4],tolerance=1e-8)
  expect_equal(predict(fit,target="population"),zero,tolerance=1e-8)
})

test_that("log-rate random effects use genuine conditional-mode Laplace integration",{
  data <- specification_data()
  control <- specification_control()
  fit <- .mmes_specification_fit(y~Env+a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log",
      random=~vsm(ism(Name),sigma2=.025,fixedSigma2=TRUE))),
    start=c(a=1.5,k=1),control=control)
  b <- coef(fit)
  conditional <- function(u){
    eta <- b[1]+b[2]*(data$Env=="B")+
      b[3]*tanh(exp(b[4]+u[data$Name])*data$x)
    -sum(dnorm(data$y,eta,sqrt(fit$dispersion),log=TRUE))+
      sum(u^2)/(.025*2)
  }
  oracle <- optim(rep(0,3),conditional,method="BFGS",
    control=list(reltol=1e-12,maxit=1000))
  h <- optimHess(oracle$par,conditional)
  likelihood <- -oracle$value-3/2*log(.025)-.5*as.numeric(determinant(h,TRUE)$modulus)
  expect_equal(fit$randomMode[[1]],oracle$par,tolerance=2e-5)
  expect_equal(as.numeric(logLik(fit)),unname(likelihood),tolerance=2e-5)
  expect_match(fit$likelihoodMethod,"Laplace")
  expect_equal(predict(fit,data[1:5,]),fitted(fit)[1:5],tolerance=1e-8)
  population <- predict(fit,target="population")
  expect_true(all(is.finite(population)))
  expect_gt(max(abs(population-predict(fit,target="zero"))),1e-5)
  unseen <- data[1:2,]
  unseen$Name <- "unseen"
  expect_error(predict(fit,unseen),"Unseen")
  expect_true(all(is.finite(predict(fit,unseen,target="zero"))))
  expect_error(.mmes_specification_fit(y~Env+a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log",
      random=~vsm(ism(Name),sigma2=.025,fixedSigma2=TRUE))),
    start=c(a=1.5,k=1),control=control,
    specificationMethod=list(integration="exact",approximation="likelihood")),
    "Exact integration")
})

test_that("nonlinear Bernoulli and Poisson likelihoods match independent GLMs",{
  data <- specification_data()
  control <- specification_control()
  set.seed(512)
  data$binary <- rbinom(nrow(data),1,plogis(.1+.6*data$x))
  data$count <- rpois(nrow(data),exp(.2+.3*data$x))
  for(family in list(binomial(),poisson())){
    response <- if(family$family=="binomial") "binary" else "count"
    formula <- as.formula(paste(response,"~a*x"))
    fit <- .mmes_specification_fit(formula,data=data,family=family,
      parameters=list(a=parameter()),start=c(a=.4),control=control)
    oracle <- glm(as.formula(paste(response,"~x")),data=data,family=family)
    expect_equal(unname(coef(fit)),unname(coef(oracle)),tolerance=2e-5)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)),tolerance=1e-8)
  }
})

test_that("nonGaussian additive random effects use Laplace, not a working fit",{
  data <- specification_data()
  set.seed(761)
  data$binary <- rbinom(nrow(data),1,plogis(.1+.7*data$x))
  fit <- .mmes_specification_fit(binary~a*x,data=data,family=binomial(),
    random=~vsm(ism(Name),sigma2=.2,fixedSigma2=TRUE),
    parameters=list(a=parameter()),start=c(a=.5),
    control=specification_control())
  conditional <- function(u){
    eta <- coef(fit)[1]+coef(fit)[2]*data$x+u[data$Name]
    -sum(dbinom(data$binary,1,plogis(eta),log=TRUE))+sum(u^2)/.4
  }
  oracle <- optim(rep(0,3),conditional,method="BFGS",
    control=list(reltol=1e-12,maxit=1000))
  h <- optimHess(oracle$par,conditional)
  likelihood <- -oracle$value-1.5*log(.2)-.5*as.numeric(determinant(h,TRUE)$modulus)
  expect_equal(as.numeric(logLik(fit)),unname(likelihood),tolerance=2e-5)
  expect_match(fit$likelihoodMethod,"Laplace")
})

test_that("Gamma log likelihood and fitted dispersion match independent ML",{
  data <- specification_data()
  set.seed(365)
  data$positive <- rgamma(nrow(data),shape=5,
    scale=exp(.2+.3*data$x)/5)
  control <- specification_control()
  fit <- .mmes_specification_fit(positive~a*x,data=data,family=Gamma("log"),
    parameters=list(a=parameter()),start=c(a=.3),control=control,rcov=~units)
  oracle <- glm(positive~x,data=data,family=Gamma("log"),
    control=glm.control(epsilon=1e-12,maxit=100L))
  phi <- optimize(function(logPhi)
    -sum(dgamma(data$positive,shape=exp(-logPhi),
      scale=fitted(oracle)*exp(logPhi),log=TRUE)),
    interval=c(-10,5),tol=1e-10)
  expect_lt(max(abs(unname(coef(fit))-unname(coef(oracle)))),5e-5)
  expect_equal(fit$dispersion,exp(phi$minimum),tolerance=2e-5)
  expect_equal(as.numeric(logLik(fit)),-phi$objective,tolerance=1e-7)
  ordinary <- .mmes_specification_fit(positive~x,data=data,family=Gamma("log"),
    control=control,rcov=~units)
  expect_equal(ordinary$dispersion,exp(phi$minimum),tolerance=2e-5)
  expect_equal(as.numeric(logLik(ordinary)),-phi$objective,tolerance=1e-7)
  expect_error(.mmes_specification_fit(positive~a*x,data=data,
    family=Gamma("identity"),parameters=list(a=parameter()),
    start=c(a=.3),control=control),"Unsupported native Laplace family/link")
  mixed <- .mmes_specification_fit(positive~a*x,data=data,family=Gamma("log"),
    random=~vsm(ism(Name),sigma2=.025,fixedSigma2=TRUE),
    parameters=list(a=parameter()),start=c(a=.3),control=control,rcov=~units)
  conditional <- function(u){
    mu <- exp(coef(mixed)[1]+coef(mixed)[2]*data$x+u[data$Name])
    -sum(dgamma(data$positive,shape=1/mixed$dispersion,
      scale=mu*mixed$dispersion,log=TRUE))+sum(u^2)/.05
  }
  reference <- optim(rep(0,3),conditional,method="BFGS",
    control=list(reltol=1e-12,maxit=1000L))
  h <- optimHess(reference$par,conditional)
  expect_equal(as.numeric(logLik(mixed)),
    -reference$value-1.5*log(.025)-.5*as.numeric(determinant(h,TRUE)$modulus),
    tolerance=2e-5,ignore_attr=TRUE)
})

test_that("uncertainty is post-fit and its evaluator survives serialization",{
  data <- specification_data()
  fit <- .mmes_specification_fit(y~a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log")),
    start=c(a=1.5,k=1),control=specification_control())
  expect_null(fit$parameterCovariance)
  expect_null(fit$hessian)
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(vcov(restored),vcov(fit),tolerance=1e-9)
  expect_equal(.spec_evaluate(restored$optimizer$par,restored$likelihoodContext),
    -as.numeric(logLik(fit)),tolerance=1e-9)
})

test_that("one shared observation mask respects explicit NA policies",{
  data <- specification_data()
  data$y[1] <- NA
  data$x[2] <- NA
  data$Name[3] <- NA
  control <- specification_control()
  args <- list(fixed=y~a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log",
      random=~vsm(ism(Name),sigma2=.025,fixedSigma2=TRUE))),
    start=c(a=1.5,k=1),control=control)
  fit <- do.call(.mmes_specification_fit,args)
  expect_equal(fit$obsInfo$excluded,1:3)
  expect_equal(nobs(fit),nrow(data)-3L)
  expect_equal(length(fitted(fit)),nrow(data)-3L)
  expect_error(do.call(.mmes_specification_fit,c(args,list(naMethodX="fail"))),"naMethodX")
  expect_error(do.call(.mmes_specification_fit,c(args,list(naMethodY="fail"))),"naMethodY")
  expect_error(do.call(.mmes_specification_fit,c(args,list(naMethodRandom="fail"))),"naMethodRandom")
  expect_error(do.call(.mmes_specification_fit,c(args,list(naMethodR="impute"))),"only")
  expect_error(do.call(.mmes_specification_fit,c(args,list(familyParameters=list(phi=.2)))),
    "familyParameters")
})

test_that("conditional scores certify requested mode tolerance and controls are owned",{
  data <- specification_data()
  control <- specification_control()
  control$task$nonlinear$mode$scoreTol <- 2e-6
  control$task$nonlinear$mode$stepTol <- 3e-10
  fit <- .mmes_specification_fit(y~a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log",
      random=~vsm(ism(Name),sigma2=.025,fixedSigma2=TRUE))),
    start=c(a=1.5,k=1),control=control,rcov=~units)
  expect_lte(fit$modeScore,2e-6)
  expect_equal(fit$likelihoodContext$control$task$nonlinear$mode$stepTol,3e-10)
  expect_equal(fit$control$convergence$likelihoodTol,control$convergence$likelihoodTol)
})

test_that("public wrapper supports homogeneous residuals Gamma and numerical controls",{
  data <- specification_data()
  set.seed(365)
  data$positive <- rgamma(nrow(data),shape=5,
    scale=exp(.2+.3*data$x)/5)
  control <- mmes_control(optimizer=list(gradient="central"),
    convergence=list(maxIters=100L,likelihoodTol=1e-8,parameterTol=1e-8),
    mode=list(tol=2e-6,stepTol=3e-10))
  fit <- mmes(positive~a*x,data=data,family=Gamma("log"),rcov=~units,
    parameters=list(a=parameter()),start=c(a=.3),
    method=mmes_method(objective="ml"),control=control)
  expect_s3_class(fit,"mmes.specification")
  expect_true(is.finite(logLik(fit)))
  expect_equal(predict(fit,data[1:3,]),fitted(fit)[1:3],tolerance=1e-9)
  bad <- mmes_control(optimizer=list(gradient="analytic"))
  expect_error(mmes(positive~a*x,data=data,family=Gamma("log"),rcov=~units,
    parameters=list(a=parameter()),start=c(a=.3),
    method=mmes_method(objective="ml"),control=bad),"analytic")
})

test_that("named graphs share the expression route with namespaced declarations",{
  skip_if_not(is.loaded("_sommer_sommer_nonlinear_graph_values_cpp"))
  data <- specification_data()
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    hidden=list(op="affine",inputs="x",weights=matrix(1L),bias=0L),
    activation=list(op="tanh",input="hidden"),
    scaled=list(op="affine",inputs="activation",weights=matrix(2L),bias=0L),
    output=list(op="identity",input="scaled")),
    output="output",parameters=c("k","a"))
  fit <- .mmes_specification_fit(y~Env+graph_term("net"),data=data,
    predictors=list(net=graph),parameters=list("net::k"=parameter(link="log")),
    start=c("net::k"=1,"net::a"=1.5),control=specification_control())
  oracle <- .mmes_specification_fit(y~Env+a*tanh(k*x),data=data,
    parameters=list(a=parameter(),k=parameter(link="log")),
    start=c(a=1.5,k=1),control=specification_control())
  expect_equal(fitted(fit),fitted(oracle),tolerance=2e-4)
  expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)),tolerance=1e-6)
  expect_true("net::k:(Intercept)" %in% names(coef(fit)))
  expect_error(.mmes_specification_fit(y~Env+graph_term("net"),data=data,
    family=binomial(),predictors=list(net=graph),
    parameters=list("net::k"=parameter(link="log")),
    start=c("net::k"=1,"net::a"=1.5),control=specification_control()),
    "Binomial responses")
})

test_that("invalid options, starts, expressions, and prediction levels fail explicitly",{
  data <- specification_data()
  control <- specification_control()
  fit <- function(...,start=c(a=1.5,k=1),formula=y~a*tanh(k*x))
    .mmes_specification_fit(formula,data=data,control=control,
      parameters=list(a=parameter(),k=parameter(link="log")),start=start,...)
  expect_error(fit(unknown=1),"Unknown")
  expect_error(fit(start=c(a=1,k=-1)),"positive")
  expect_error(fit(start=c(a=1)),"exactly")
  expect_error(fit(formula=y~a*get("x")+k),"Unsupported")
  expect_error(fit(formula=y~a*tanh(k*missing)),"Unsupported")
  expect_error(fit(rcov=~Env),"independent residuals")
  expect_error(fit(predictors=list(net=list())),"exactly")
  control$task$objective <- "reml"
  expect_error(fit(),"require task")
  control$task$objective <- "ml"
  data$x[1] <- NA
  expect_error(fit(naMethodX="fail"),"Missing values")
})
