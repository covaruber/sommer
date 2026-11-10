nonlinearRandomData <- function(curve="exponential", parameter="rate", groups=16L){
  set.seed(917)
  group <- factor(rep(seq_len(groups), each=14L))
  x <- rep(seq(0, 6, length.out=14L), groups)
  theta <- if(curve=="logistic") c(log(4),log(1.2),3) else c(log(4),log(0.6))
  coordinate <- match(parameter, nonlinear_curve("x",curve)$parameters)
  deviations <- rnorm(groups, sd=0.22)
  local <- matrix(rep(theta,each=length(x)),nrow=length(x))
  local[,coordinate] <- local[,coordinate]+deviations[group]
  mean <- switch(curve,
    exponential=exp(local[,1])*exp(-exp(local[,2])*x),
    logistic=exp(local[,1])*plogis(exp(local[,2])*(x-local[,3])),
    saturation=exp(local[,1])*x/(exp(local[,2])+x))
  data.frame(y=mean+rnorm(length(x),sd=0.12),x=x,group=group)
}

nonlinearRandomContext <- function(data, curve="exponential", parameter="rate",
    additive=NULL, random=~vsm(ism(group),sigma2=0.04), fixed=y~0, ...){
  model <- nonlinear_curve("x",curve,random=stats::setNames(list(random),parameter))
  control <- sommer:::.mmes_control_internal(task=list(objective="ml"),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  joint <- sommer:::.nonlinear_joint_formula(additive,model,data)
  prepared <- sommer:::.nonlinear_prepare(fixed,joint,data,control,
    "exclude","exclude","exclude","exclude",NULL,data$x)
  layout <- sommer:::.nonlinear_random_layout(model,prepared)
  context <- sommer:::sommer_nonlinear_laplace_create(prepared$X,
    do.call(cbind,prepared$Z),data$y,numeric(nrow(data)),data$x,curve,
    layout$coordinate,layout$start,prepared$covStruct[seq_along(prepared$rtermss)],
    prepared$Ai,...)
  list(context=context,prepared=prepared,layout=layout,model=model)
}

test_that("curve coordinate derivatives through order three are analytic",{
  for(curve in c("exponential","logistic","saturation")){
    theta <- if(curve=="logistic") c(log(4),log(1.2),2) else c(log(4),log(0.6))
    for(x in c(0,0.5,2,6)){
      value <- sommer:::sommer_nonlinear_curve_derivatives(curve,x,theta)
      expect_equal(as.numeric(value$first),as.numeric(
        sommer:::sommer_nonlinear_curve_eval(curve,x,theta)$jacobian))
      for(step in c(1e-4,1e-5)){
        for(j in seq_along(theta)){
          upper <- lower <- theta
          upper[j] <- upper[j]+step; lower[j] <- lower[j]-step
          up <- sommer:::sommer_nonlinear_curve_derivatives(curve,x,upper)
          down <- sommer:::sommer_nonlinear_curve_derivatives(curve,x,lower)
          expect_equal(value$second[,j],as.numeric((up$first-down$first)/(2*step)),
            tolerance=2e-7)
          expect_equal(value$third[,,j],(up$second-down$second)/(2*step),tolerance=2e-7)
        }
      }
    }
  }
})

test_that("exact conditional Hessians and complete Laplace scores match differences",{
  cases <- list(c("exponential","rate"),c("exponential","asymptote"),
    c("logistic","rate"),c("logistic","midpoint"),c("logistic","asymptote"),
    c("saturation","half"),c("saturation","asymptote"))
  for(case in cases){
    data <- nonlinearRandomData(case[1],case[2],groups=3L)
    fixture <- nonlinearRandomContext(data,case[1],case[2])
    theta <- if(case[1]=="logistic") c(log(4),log(1.2),3) else c(log(4),log(0.6))
    parameters <- c(theta,log(0.04),log(0.12^2))
    evaluate <- function(par) sommer:::sommer_nonlinear_laplace_eval(
      fixture$context,par)
    value <- evaluate(parameters)
    expect_true(value$valid,info=value$reason)
    if(!value$valid) next
    conditional <- function(mode) sommer:::sommer_nonlinear_conditional_eval(
      fixture$context,parameters,mode)
    mode <- as.numeric(value$mode)
    point <- mode+c(0.08,-0.05,0.03)
    reference <- stats::optim(rep(0,3),function(w) conditional(w)$objective,
      function(w) -as.numeric(conditional(w)$score),method="BFGS",
      control=list(reltol=1e-12))
    expect_equal(mode,reference$par,tolerance=2e-6)
    hessian <- as.matrix(conditional(point)$hessian)
    for(step in c(1e-4,1e-5)){
      scoreDifference <- vapply(seq_along(point),function(j){
        up <- down <- point
        up[j] <- up[j]+step; down[j] <- down[j]-step
        -(conditional(up)$objective-conditional(down)$objective)/(2*step)
      },numeric(1))
      expect_equal(as.numeric(conditional(point)$score),scoreDifference,tolerance=2e-6)
      numerical <- vapply(seq_along(point),function(j){
        up <- down <- point
        up[j] <- up[j]+step; down[j] <- down[j]-step
        -as.numeric(conditional(up)$score-conditional(down)$score)/(2*step)
      },numeric(length(point)))
      expect_equal(hessian,numerical,tolerance=2e-6)
      score <- vapply(seq_along(parameters),function(j){
        up <- down <- parameters
        up[j] <- up[j]+step; down[j] <- down[j]-step
        upper <- evaluate(up); lower <- evaluate(down)
        expect_true(upper$valid,info=upper$reason)
        expect_true(lower$valid,info=lower$reason)
        (upper$logLik-lower$logLik)/(2*step)
      },numeric(1))
      expect_equal(as.numeric(value$gradient),score,tolerance=3e-5)
    }
    expect_lt(value$modeScore,1e-7)
    expect_equal(evaluate(parameters)$logLik,value$logLik,tolerance=0)
    conditionalValue <- conditional(mode)
    independent <- -conditionalValue$objective-
      nrow(data)/2*log(2*pi*0.12^2)-3/2*log(0.04)-
      as.numeric(determinant(as.matrix(conditionalValue$hessian),logarithm=TRUE)$modulus)/2
    expect_equal(value$logLik,independent,tolerance=1e-10)
  }
})

test_that("Laplace integration agrees with numerical integration and its Gaussian limit",{
  data <- nonlinearRandomData(groups=2L)
  fixture <- nonlinearRandomContext(data)
  parameters <- c(log(4),log(0.6),log(0.025),log(0.12^2))
  value <- sommer:::sommer_nonlinear_laplace_eval(fixture$context,parameters)
  expect_true(value$valid,info=value$reason)
  exact <- sum(vapply(levels(data$group),function(group){
    rows <- which(data$group==group)
    integrand <- function(b) vapply(b,function(deviation){
      mean <- 4*exp(-0.6*exp(deviation)*data$x[rows])
      mean[data$x[rows]==0] <- 4
      sum(dnorm(data$y[rows],mean,0.12,log=TRUE))+
        dnorm(deviation,0,sqrt(0.025),log=TRUE)
    },numeric(1))
    shift <- max(integrand(seq(-0.8,0.8,length.out=100L)))
    shift+log(integrate(function(b) exp(integrand(b)-shift),-Inf,Inf,
      rel.tol=1e-9)$value)
  },numeric(1)))
  expect_lt(abs(value$logLik-exact),0.015)
  parameters[3] <- log(1e-10)
  limit <- sommer:::sommer_nonlinear_laplace_eval(fixture$context,parameters)
  expect_true(limit$valid,info=limit$reason)
  expect_equal(limit$logLik,sum(dnorm(data$y,4*exp(-0.6*data$x),0.12,log=TRUE)),
    tolerance=2e-6)
})

test_that("random curve fitting separates conditional, zero and population predictions",{
  data <- nonlinearRandomData()
  model <- nonlinear_curve("x",random=list(rate=~vsm(ism(group),sigma2=0.04)))
  control <- sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
    start=c(asymptote=4,rate=0.6))),
    convergence=list(maxIters=300L,likelihoodTol=1e-9,parameterTol=1e-7),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  seed <- .Random.seed
  fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=control)
  expect_true(fit$convergence)
  expect_identical(.Random.seed,seed)
  expect_equal(unname(coef(fit)),c(4,0.6),tolerance=0.1)
  expect_match(fit$executionPlan$objective,"Laplace")
  expect_match(capture.output(print(fit))[1],"Laplace Gaussian ML")
  expect_equal(fitted(fit),predict(fit,newdata=data),tolerance=1e-10)
  expect_equal(predict(fit,re.form=NA),predict(fit,newdata=data,re.form=NA))
  population <- predict(fit,re.form="population")
  expect_equal(population,predict(fit,newdata=data,re.form="population"))
  expect_gt(max(abs(population-predict(fit,re.form=NA))),0.005)
  expect_equal(fit$groupParameters[,"rate"],
    fit$nonlinearParameters["rate"]*exp(fit$curveRandomMode),ignore_attr=TRUE)
  target <- 4*exp(-0.6*data$x)
  expect_lt(mean((fitted(fit)-data$y)^2),mean((target-data$y)^2))
  variance <- fit$covPar[[1]][1]
  for(row in c(2L,5L,14L)){
    reference <- integrate(function(b) fit$nonlinearParameters[1]*
      exp(-fit$nonlinearParameters[2]*exp(b)*data$x[row])*
      dnorm(b,sd=sqrt(variance)),-Inf,Inf,rel.tol=1e-9)$value
    expect_equal(population[row],reference,tolerance=2e-7)
  }
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(predict(restored,newdata=data[1:4,]),fitted(fit)[1:4])
  expect_equal(predict(restored,re.form="population"),population)
  expect_null(fit$Ci)
  expect_null(fit$optimizer$hessian)
  expect_false(fit$nonlinearDiagnostics$uncertaintyComputed)
  expect_null(fit$nonlinearInference)
  expect_lt(fit$nonlinearDiagnostics$modeScore,1e-6)
})

test_that("additive and related curve effects retain analytic covariance scores",{
  data <- nonlinearRandomData(groups=3L)
  data$cluster <- factor(rep(c("a","b","a"),each=14L))
  gu <- Matrix::Matrix(matrix(c(1,0.15,0.05,0.15,1,0.1,0.05,0.1,1),3,3),sparse=TRUE)
  dimnames(gu) <- list(levels(data$group),levels(data$group))
  precision <- solve(gu)
  attr(precision,"inverse") <- TRUE
  fixture <- nonlinearRandomContext(data,additive=~vsm(ism(cluster),sigma2=0.02),
    random=~vsm(ism(group),Gu=precision,sigma2=0.04))
  parameters <- c(log(4),log(0.6),log(0.02),log(0.04),log(0.12^2))
  evaluate <- function(par) sommer:::sommer_nonlinear_laplace_eval(fixture$context,par)
  value <- evaluate(parameters)
  expect_true(value$valid,info=value$reason)
  step <- 1e-5
  numerical <- vapply(seq_along(parameters),function(j){
    up <- down <- parameters
    up[j] <- up[j]+step; down[j] <- down[j]-step
    (evaluate(up)$logLik-evaluate(down)$logLik)/(2*step)
  },numeric(1))
  expect_equal(as.numeric(value$gradient),numerical,tolerance=5e-5)
})

test_that("random curve contracts and failure diagnostics are explicit",{
  expect_identical(nonlinear_curve("x",
    random=list(rate=~group,asymptote=~group))$version,3L)
  expect_error(nonlinear_curve("x",random=list(typo=~group)),"parameter-named")
  expect_error(nonlinear_curve("x",random=list(rate=~group+other)),"one-term")
  expect_error(sommer:::.mmes_control_internal(task=list(nonlinear=list(mode=list(warmStart="on")))),
    "off")
  expect_error(sommer:::.mmes_control_internal(task=list(nonlinear=list(population=list(order=1)))),
    "order>=2")
  data <- nonlinearRandomData(groups=3L)
  fixture <- nonlinearRandomContext(data,maxIters=1L)
  value <- sommer:::sommer_nonlinear_laplace_eval(fixture$context,
    c(log(4),log(0.6),log(0.04),log(0.12^2)))
  expect_false(value$valid)
  expect_match(value$reason,"all deterministic mode starts failed")
  expect_error(nonlinearRandomContext(data,memoryMB=0.001),"workspace")
  expect_error(sommer:::.mmes_dispatch(y~0,random=~group,data=data,nonlinear=nonlinear_curve("x"),
    control=sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
      population=list(order=4L))))),"population controls require random curve")
})

test_that("structured random curves and fixed mean scores use the same coordinate basis",{
  data <- nonlinearRandomData(groups=3L)
  data$occasion <- factor(rep(rep(1:2,each=7L),3L))
  data$z <- rep(c(-1,1),length.out=nrow(data))
  fixture <- nonlinearRandomContext(data,
    random=~vsm(ar1m(occasion),ism(group),sigma2=0.04),fixed=y~0+z)
  structures <- fixture$prepared$covStruct[seq_along(fixture$prepared$rtermss)]
  covariance <- unlist(lapply(structures,function(s) s$par[s$free]),use.names=FALSE)
  parameters <- c(0.03,log(4),log(0.6),covariance,log(0.12^2))
  evaluate <- function(par) sommer:::sommer_nonlinear_laplace_eval(fixture$context,par)
  value <- evaluate(parameters)
  expect_true(value$valid,info=value$reason)
  numerical <- vapply(seq_along(parameters),function(j){
    up <- down <- parameters
    up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
    (evaluate(up)$logLik-evaluate(down)$logLik)/2e-5
  },numeric(1))
  expect_equal(as.numeric(value$gradient),numerical,tolerance=5e-5)
  expect_equal(fixture$layout$count,6L)
})

test_that("additive effects, fixed covariance and missingness survive random curve fitting",{
  data <- nonlinearRandomData(groups=12L)
  data$cluster <- factor(rep(rep(1:3,each=14L),4L))
  data$y <- data$y+c(-0.2,0,0.2)[data$cluster]
  data$x[3] <- NA
  data$y[7] <- NA
  control <- sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
    start=c(asymptote=4,rate=0.6))),
    convergence=list(maxIters=300L),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  fit <- sommer:::.mmes_dispatch(y~0,random=~vsm(ism(cluster),sigma2=0.04,fixedSigma2=TRUE),
    nonlinear=nonlinear_curve("x",random=list(rate=~vsm(ism(group),sigma2=0.04))),
    data=data,control=control)
  expect_true(fit$convergence)
  expect_equal(fit$covPar[[1]][1],0.04)
  expect_length(fit$additiveRandomMode,3L)
  expect_equal(nobs(fit),nrow(data)-2L)
  retained <- data[fit$obsInfo$included,]
  expect_equal(fitted(fit),predict(fit,newdata=retained),tolerance=1e-10)
  expect_equal(predict(fit,re.form="population"),
    predict(fit,newdata=retained,re.form="population"))
  expect_equal(predict(fit,re.form=NA),predict(fit,newdata=retained,re.form=NA))
  expect_null(fit$fixedCovariance)
  bad <- retained[1,]
  bad$group <- "unseen"
  expect_error(predict(fit,newdata=bad),"unseen")
  expect_error(predict(fit,newdata=bad,re.form="population"),"unseen")
  expect_length(predict(fit,newdata=bad,re.form=NA),1L)
  expect_error(sommer:::.mmes_dispatch(y~0,data=retained,nonlinear=nonlinear_curve("x",
    random=list(rate=~vsm(ism(group),sigma2=0,fixedSigma2=TRUE))),
    control=control),"positive")
})

test_that("other traditional curves and relationship predictions are serializable",{
  for(case in list(c("logistic","midpoint"),c("saturation","half"))){
    data <- nonlinearRandomData(case[1],case[2],groups=10L)
    relationship <- 0.85*diag(10)+0.15
    dimnames(relationship) <- list(levels(data$group),levels(data$group))
    precision <- solve(relationship)
    attr(precision,"inverse") <- TRUE
    model <- nonlinear_curve("x",case[1],random=stats::setNames(list(
      ~vsm(ism(group),Gu=precision,sigma2=0.04,fixedSigma2=TRUE)),case[2]))
    start <- if(case[1]=="logistic") c(asymptote=4,rate=1.2,midpoint=3) else
      c(asymptote=4,half=0.6)
    control <- sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(start=start)),
      convergence=list(maxIters=300L),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE))
    fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=control)
    expect_true(fit$convergence)
    expect_equal(unname(coef(fit)),unname(start),tolerance=0.15)
    conditional <- predict(fit,newdata=data[1:4,])
    population <- predict(fit,newdata=data[1:4,],re.form="population")
    expect_equal(conditional,fitted(fit)[1:4],tolerance=1e-10)
    expect_equal(population,predict(fit,re.form="population")[1:4])
    precision <- NULL
    restored <- unserialize(serialize(fit,NULL))
    expect_equal(predict(restored,newdata=data[1:4,]),conditional)
    expect_equal(predict(restored,newdata=data[1:4,],re.form="population"),population)
  }
})

test_that("population quadrature rejects unverified accuracy",{
  data <- nonlinearRandomData(groups=3L)
  fixture <- nonlinearRandomContext(data)
  object <- list(nonlinearRandom=fixture$layout,
    covStruct=fixture$prepared$covStruct[1],
    randomRelationships=fixture$prepared$Ai,nonlinear=fixture$model,
    nonlinearCoordinates=c(log(4),log(0.6)),
    control=sommer:::.mmes_control_internal(task=list(nonlinear=list(
      population=list(order=2L,maxOrder=4L,tol=1e-14)))))
  expect_error(sommer:::.nonlinear_population_prediction(object,data$x,
    fixture$layout$design),"did not meet")
})

test_that("curve covariance equality constraints pull back complete Laplace scores",{
  data <- nonlinearRandomData(groups=3L)
  data$cluster <- factor(rep(c("a","b","a"),each=14L))
  fixture <- nonlinearRandomContext(data,additive=~vsm(ism(cluster),sigma2=0.08))
  prepared <- fixture$prepared
  mapping <- sommer:::.laplace_covariance_map(prepared$covStruct[1:2],
    prepared$rtermss,data.frame(parameter=1:2,group=1L,scale=c(2,1)))
  layout <- list(fixedCount=2L,covarianceCount=length(mapping$map$initial),
    covarianceMap=mapping$map)
  parameters <- c(log(4),log(0.6),mapping$map$initial,log(0.12^2))
  evaluate <- function(par) sommer:::.laplace_constrained_evaluate(par,layout,
    function(expanded) sommer:::sommer_nonlinear_laplace_eval(fixture$context,expanded))
  value <- evaluate(parameters)
  expect_true(value$valid,info=value$reason)
  numerical <- vapply(seq_along(parameters),function(j){
    up <- down <- parameters
    up[j] <- up[j]+1e-5; down[j] <- down[j]-1e-5
    (evaluate(up)$logLik-evaluate(down)$logLik)/2e-5
  },numeric(1))
  expect_equal(as.numeric(value$gradient),numerical,tolerance=5e-5)
})
