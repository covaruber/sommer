phase11SpecificationControl <- function(){
  .mmes_control_build(task=list(objective="ml"),
    convergence=list(maxIters=150L,likelihoodTol=1e-9,parameterTol=1e-9),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
}

test_that("legacy expression links retain their reference likelihood semantics",{
  set.seed(1840)
  data <- data.frame(x=rep(seq(-1,1,length.out=12L),4L))
  data$y <- exp(.3+.2*data$x)+rnorm(nrow(data),sd=.15)
  fit <- .mmes_specification_fit(y~a*x,data=data,family=gaussian("log"),
    parameters=list(a=parameter()),start=c(a=.2),control=phase11SpecificationControl())
  oracle <- glm(y~x,data=data,family=gaussian("log"),
    control=glm.control(epsilon=1e-12))
  expect_equal(as.numeric(predict(fit)),as.numeric(fitted(oracle)),tolerance=2e-6)
  expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)),tolerance=1e-7)
  data$y <- rpois(nrow(data),(1.2+.2*data$x)^2)
  ordinary <- .mmes_specification_fit(y~x,data=data,family=poisson("sqrt"),
    control=phase11SpecificationControl())
  oracle <- glm(y~x,data=data,family=poisson("sqrt"))
  expect_equal(as.numeric(predict(ordinary)),as.numeric(fitted(oracle)),tolerance=1e-12)
})

test_that("specifications reuse grouped-binomial likelihood constants and links",{
  set.seed(1841)
  data <- data.frame(x=rep(seq(-1,1,length.out=12L),4L))
  data$success <- rbinom(nrow(data),5,plogis(.2+.7*data$x))
  data$y <- cbind(data$success,5-data$success)
  for(link in c("logit","probit","cloglog")){
    family <- binomial(link)
    fit <- .mmes_specification_fit(y~a*x,data=data,family=family,
      parameters=list(a=parameter()),start=c(a=.4),control=phase11SpecificationControl())
    oracle <- glm(y~x,data=data,family=family,
      control=glm.control(epsilon=1e-12))
    expect_equal(unname(coef(fit)),unname(coef(oracle)),tolerance=2e-5)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(oracle)),tolerance=1e-7)
    expect_equal(as.numeric(predict(fit)),as.numeric(fitted(oracle)),tolerance=2e-6)
    expect_equal(predict(fit,newdata=data),predict(fit),tolerance=1e-12)
    draw <- simulate(fit,nsim=2,seed=119)
    expect_true(all(round(as.matrix(draw)*5) %in% 0:5))
    expect_lt(max(abs(as.matrix(draw)*5-round(as.matrix(draw)*5))),1e-12)
    expect_error(predict(fit,newdata=data,se.fit=TRUE,interval="prediction"),
      "require trials")
    prediction <- predict(fit,newdata=data,interval="prediction",trials=5)
    expect_true(all(prediction$lower<prediction$upper))
  }
  data$failure <- 5-data$success
  direct <- .mmes_specification_fit(cbind(success,failure)~a*x,data=data,
    family=binomial(),parameters=list(a=parameter()),start=c(a=.4),
    control=phase11SpecificationControl())
  oracle <- glm(cbind(success,failure)~x,data=data,family=binomial())
  expect_equal(as.numeric(logLik(direct)),as.numeric(logLik(oracle)),tolerance=1e-7)
  ordinary <- .mmes_specification_fit(cbind(success,failure)~x,data=data,
    family=binomial(),control=phase11SpecificationControl())
  link <- predict(ordinary,type="link",cov.fit=TRUE)
  x <- model.matrix(oracle)
  expect_equal(unname(link$cov.fit),unname(x %*% vcov(ordinary) %*% t(x)),
    tolerance=1e-12)
  expect_error(predict(ordinary,newdata=data,interval="prediction"),"require trials")
})

test_that("specification NB2 and ZIP nuisances use native family contracts",{
  set.seed(1842)
  data <- data.frame(x=rep(seq(-1,1,length.out=10L),5L))
  data$y <- rnbinom(nrow(data),mu=exp(.3+.4*data$x),size=3)
  nb <- .mmes_specification_fit(y~a*x,data=data,family=nbinom2(3),
    familyParameters=list(theta=3,estimate=FALSE),parameters=list(a=parameter()),
    start=c(a=.2),control=phase11SpecificationControl())
  expect_equal(as.numeric(logLik(nb)),sum(dnbinom(data$y,size=3,
    mu=exp(predict(nb,type="link")),log=TRUE)),tolerance=1e-9)
  expect_true(all(is.finite(residuals(nb,type="pearson"))))
  data$y <- rpois(nrow(data),exp(.2+.5*data$x))
  data$y[seq(1,nrow(data),by=4)] <- 0
  zip <- .mmes_specification_fit(y~a*x,data=data,family=zipoisson(.25),
    familyParameters=list(zeroInflation=.25,estimate=FALSE),
    parameters=list(a=parameter()),start=c(a=.2),control=phase11SpecificationControl())
  mu <- exp(predict(zip,type="link"))
  density <- ifelse(data$y==0,.25+.75*exp(-mu),.75*dpois(data$y,mu))
  expect_equal(as.numeric(logLik(zip)),sum(log(density)),tolerance=1e-8)
  expect_equal(predict(zip),.75*mu,tolerance=1e-12)
  expect_equal(coef(zip,type="all")["family:response:zeroInflation.1"],
    c("family:response:zeroInflation.1"=.25))
  expect_equal(simulate(zip,seed=121),simulate(zip,seed=121))
})

test_that("ordinal specifications identify thresholds and retain probability layouts",{
  set.seed(1843)
  data <- data.frame(x=rep(seq(-1.5,1.5,length.out=15L),5L))
  latent <- .6*data$x+rlogis(nrow(data))
  data$y <- ordered(cut(latent,c(-Inf,-.5,.7,Inf),labels=c("low","mid","high")))
  fit <- .mmes_specification_fit(y~0+a*x,data=data,family=ordinal(),
    parameters=list(a=parameter()),start=c(a=.3),control=phase11SpecificationControl())
  probabilities <- predict(fit)
  expect_equal(colnames(probabilities),levels(data$y))
  expect_equal(unname(rowSums(probabilities)),rep(1,nrow(data)),tolerance=1e-12)
  expect_equal(predict(fit,newdata=data),probabilities,tolerance=1e-12)
  thresholds <- fit$familySpecifications[[1L]]$thresholds
  cumulative <- plogis(outer(predict(fit,type="link"),thresholds,function(eta,t) t-eta))
  oracle <- cbind(cumulative[,1],cumulative[,2]-cumulative[,1],1-cumulative[,2])
  expect_equal(unname(probabilities),unname(oracle),tolerance=1e-12)
  expect_equal(as.numeric(logLik(fit)),sum(log(oracle[cbind(seq_len(nrow(data)),
    as.integer(data$y))])),tolerance=1e-9)
  expect_true(all(diff(thresholds)>0))
  expect_equal(dim(predict(fit,se.fit=TRUE)$se.fit),dim(probabilities))
  expect_error(predict(fit,interval="prediction"),"scalar response")
  expect_true(is.ordered(simulate(fit,seed=120)[[1L]]))
  expect_equal(length(residuals(fit,type="pearson")),nrow(data))
  expect_error(.mmes_specification_fit(y~0+a,data=data,family=ordinal(),
    parameters=list(a=parameter()),start=c(a=.3),control=phase11SpecificationControl()),
    "threshold locations")
})

test_that("nominal parameter regressions preserve reference ordering and block probabilities",{
  set.seed(1844)
  data <- data.frame(x=rep(seq(-1,1,length.out=12L),5L))
  data$y <- factor(sample(c("a","b","c"),nrow(data),replace=TRUE),
    levels=c("a","b","c"))
  declarations <- list(a=parameter(fixed=~0+.familyCoordinate))
  initial <- list(coefficients=c("a:.familyCoordinateresponse:a"=.2,
    "a:.familyCoordinateresponse:c"=-.2))
  fit <- .mmes_specification_fit(y~1+a*x,data=data,family=multinomial(reference="b"),
    parameters=declarations,start=list(coefficients=c(
      "fixed:response:a:Intercept"=0,"fixed:response:c:Intercept"=0,
      initial$coefficients)),control=phase11SpecificationControl())
  probabilities <- predict(fit)
  expect_equal(colnames(probabilities),levels(data$y))
  expect_equal(unname(rowSums(probabilities)),rep(1,nrow(data)),tolerance=1e-12)
  expect_equal(predict(fit,newdata=data),probabilities,tolerance=1e-12)
  expect_equal(as.numeric(logLik(fit)),sum(log(probabilities[
    cbind(seq_len(nrow(data)),as.integer(data$y))])),tolerance=1e-8)
  expect_equal(dim(residuals(fit)),dim(probabilities))
  expect_true(all(simulate(fit,seed=122)[[1L]] %in% levels(data$y)))
  uncertainty <- predict(fit,se.fit=TRUE,cov.fit=TRUE)
  expect_equal(dim(uncertainty$cov.fit),rep(length(probabilities),2L))
})

test_that("specification inference matches a fixed-variance Gaussian linear oracle",{
  set.seed(1845)
  data <- expand.grid(x=seq(-1,1,length.out=8L),group=factor(1:3))
  data$y <- .4+.7*data$x+rep(c(-.3,.15,.25),each=8L)+rnorm(nrow(data),sd=.25)
  fit <- .mmes_specification_fit(y~a*x,data=data,parameters=list(a=parameter()),
    random=~vsm(ism(group),sigma2=.12,fixedSigma2=TRUE),
    familyParameters=list(dispersion=.0625,estimate=FALSE),start=c(a=.5),
    control=phase11SpecificationControl())
  z <- model.matrix(~0+group,data)
  x <- model.matrix(~x,data)
  v <- .12*tcrossprod(z)+diag(.0625,nrow(data))
  expected <- solve(crossprod(x,solve(v,x)))
  expect_equal(unname(vcov(fit)),unname(expected),tolerance=3e-5)
  conditional <- solve(crossprod(z)/.0625+diag(1/.12,3L))
  expect_equal(unname(vcov(fit,type="conditional")),unname(conditional),tolerance=3e-5)
  post <- post_nonlinear(fit,mode=2)
  expect_match(post$nonlinearInference$randomMethod,"not an exact frequentist PEV")
  expect_true(all(diag(vcov(post,type="random"))>=diag(conditional)))
  expect_equal(length(confint(post)),2L*length(coef(fit)))
  expect_equal(predict(fit,se.fit=TRUE)$fit,predict(fit))
  restored <- unserialize(serialize(fit,NULL))
  expect_equal(vcov(restored),vcov(fit),tolerance=1e-9)
  bad <- fit
  bad$convergence <- FALSE
  expect_error(vcov(bad),"converged")
  bad <- fit
  bad$optimizer$par[1L] <- bad$optimizer$par[1L]+.1
  expect_error(vcov(bad),"differs from")
})

test_that("reference prior predictions integrate more than two independent predictor records",{
  data <- data.frame(x=seq(-1,1,length.out=12L),group=factor(rep("all",12L)))
  data$y <- exp(.2+.15*data$x)+.02*sin(2*data$x)
  declaration <- function() parameter(random=~vsm(ism(group),sigma2=.015,
    fixedSigma2=TRUE))
  fit <- .mmes_specification_fit(y~0+exp(a+b*x+c*x^2),data=data,
    parameters=list(a=declaration(),b=declaration(),c=declaration()),
    start=c(a=.2,b=.15,c=.01),
    familyParameters=list(dispersion=.01,estimate=FALSE),
    control=phase11SpecificationControl())
  expect_equal(predict(fit,target="population"),
    exp(coef(fit)[1]+coef(fit)[2]*data$x+coef(fit)[3]*data$x^2+
      .015*(1+data$x^2+data$x^4)/2),tolerance=2e-6)
  expect_equal(predict(fit,target="population",type="link"),
    predict(fit,target="population"),tolerance=1e-12)
  set.seed(42)
  before <- .Random.seed
  draw <- simulate(fit,nsim=2,seed=123,target="population")
  expect_equal(.Random.seed,before)
  expect_equal(dim(draw),c(nrow(data),2L))
})

test_that("mixed ordinal and scalar specifications retain trait layouts",{
  set.seed(1846)
  data <- data.frame(x=rep(seq(-1,1,length.out=30L),2L),
    trait=factor(rep(c("ord","count"),each=30L)))
  data$y <- c(as.integer(cut(.5*data$x[1:30]+rlogis(30L),
    c(-Inf,-.6,.6,Inf))),rpois(30L,exp(.3*data$x[31:60])))
  fit <- .mmes_specification_fit(y~0+a*x,data=data,
    family=familym("trait",ord=ordinal(),count=poisson()),
    parameters=list(a=parameter(fixed=~0+trait)),
    start=list(coefficients=c("a:traitcount"=.2,"a:traitord"=.3)),
    control=phase11SpecificationControl())
  output <- predict(fit)
  expect_named(output,c("ord","count"))
  expect_equal(dim(output$ord),c(30L,3L))
  expect_length(output$count,30L)
  expect_equal(unname(rowSums(output$ord)),rep(1,30L),tolerance=1e-12)
  expect_equal(predict(fit,newdata=data),output,tolerance=1e-12)
  expect_equal(nobs(fit),nrow(data))
  draw <- simulate(fit,seed=128)
  expect_named(draw$sim_1,c("ord","count"))
  expect_true(is.ordered(draw$sim_1$ord))
  expect_true(all(draw$sim_1$count>=0))
  expect_true(all(is.finite(unlist(residuals(fit,type="pearson")))))
})

test_that("nominal population integration retains cross-category prior covariance",{
  set.seed(1847)
  data <- data.frame(x=seq_len(30L),group=factor(rep(1:2,each=15L)))
  data$y <- factor(rep(c("a","b","c","b","a"),6L))
  shape <- matrix(c(1,.4,.4,1.3),2L)
  fit <- .mmes_specification_fit(y~0+a,data=data,family=multinomial(reference="c"),
    parameters=list(a=parameter(fixed=~0+.familyCoordinate,
      random=~vsm(usm(.familyCoordinate,theta=shape,fixed=c(TRUE,TRUE)),
        ism(group),sigma2=.05,fixedSigma2=TRUE))),
    start=list(coefficients=c("a:.familyCoordinateresponse:a"=.2,
      "a:.familyCoordinateresponse:b"=.2)),
    control=phase11SpecificationControl())
  z <- fit$likelihoodContext$randomDesign[[1L]][1:2,,drop=FALSE]
  covariance <- z %*% fit$randomCovariances[[1L]] %*% t(z)
  expect_gt(abs(covariance[1,2]),.001)
  rule <- .laplace_normal_quadrature(40L)
  grid <- expand.grid(a=seq_along(rule$nodes),b=seq_along(rule$nodes))
  nodes <- cbind(rule$nodes[grid$a],rule$nodes[grid$b]) %*% chol(covariance)
  logits <- sweep(nodes,2L,as.numeric(coef(fit)),"+")
  exponential <- exp(cbind(logits,0))
  probabilities <- exponential/rowSums(exponential)
  expected <- as.numeric(crossprod(rule$weights[grid$a]*rule$weights[grid$b],
    probabilities))
  population <- predict(fit,target="population")
  expect_equal(as.numeric(population[1,]),expected,tolerance=1e-9)
  expect_equal(unname(rowSums(population)),rep(1,nrow(data)),tolerance=1e-12)
  draw <- simulate(fit,nsim=500L,seed=129,target="population")
  values <- unlist(lapply(draw,as.character),use.names=FALSE)
  empirical <- as.numeric(table(factor(values,levels=levels(data$y))))/length(values)
  expect_equal(empirical,expected,tolerance=.025)
})

test_that("namespaced graph regressions retain exact nonlinear conditional chain rules",{
  set.seed(1848)
  data <- data.frame(x=seq(-1,1,length.out=8L),group=factor(1L))
  data$y <- tanh(.5*data$x)+rnorm(nrow(data),sd=.1)
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    scaled=list(op="affine",inputs="x",weights=matrix(1L),bias=0L),
    mean=list(op="tanh",input="scaled"),
    output=list(op="identity",input="mean")),
    output="output",parameters="k")
  fit <- .mmes_specification_fit(y~0+graph_term("net"),data=data,
    predictors=list(net=graph),
    parameters=list("net::k"=parameter(link="log",
      random=~vsm(ism(group),sigma2=.025,fixedSigma2=TRUE))),
    familyParameters=list(dispersion=.01,estimate=FALSE),
    start=c("net::k"=.45),control=phase11SpecificationControl())
  conditional <- .spec_conditional(fit$optimizer$par,fit$likelihoodContext,
    fit$randomCovariances)
  mode <- unlist(fit$randomMode,use.names=FALSE)
  expect_equal(conditional$hessian(mode),optimHess(mode,conditional$value,
    gr=conditional$gradient),tolerance=2e-5)
  probe <- mode+.17
  expect_equal(conditional$gradient(probe),.spec_gradient(probe,conditional$value),
    tolerance=2e-6)
  expect_equal(conditional$hessian(probe),optimHess(probe,conditional$value),
    tolerance=2e-5)
  expect_lt(max(abs(conditional$gradient(mode))),1e-8)
  expect_equal(predict(fit,newdata=data),fitted(fit),tolerance=1e-12)
  expect_true(all(is.finite(vcov(fit))))
  k <- exp(coef(fit)[[1L]])
  integrated <- integrate(function(u) vapply(u,function(value)
    exp(sum(dnorm(data$y,tanh(k*exp(value)*data$x),.1,log=TRUE))+
      dnorm(value,sd=sqrt(.025),log=TRUE)),numeric(1)),
    lower=-10,upper=10,rel.tol=1e-10)$value
  expect_lt(abs(as.numeric(logLik(fit))-log(integrated)),.01)
})

test_that("known specification parameters retain conditional uncertainty without outer coordinates",{
  data <- data.frame(y=c(-.2,.1,.3,.2,-.1,.4),
    group=factor(rep(1:2,each=3)))
  fit <- .mmes_specification_fit(y~0,data=data,
    random=~vsm(ism(group),sigma2=.1,fixedSigma2=TRUE),
    familyParameters=list(dispersion=.2,estimate=FALSE),
    control=phase11SpecificationControl())
  inference <- post_nonlinear(fit,mode=2L)
  z <- model.matrix(~0+group,data)
  expected <- unname(solve(crossprod(z)/.2+diag(10,2)))
  expect_equal(dim(vcov(inference)),c(0L,0L))
  expect_equal(unname(vcov(inference,type="conditional")),expected,tolerance=1e-8)
  expect_equal(unname(vcov(inference,type="random")),expected,tolerance=1e-8)
  expect_equal(as.numeric(fit$randomMode[[1L]]),
    as.numeric(expected %*% crossprod(z,data$y)/.2),tolerance=1e-10)
})

for(seed in c(75110L,75123L)){
  test_that(paste("binomial specification information matches a scalar-root oracle",seed),{
    set.seed(seed)
    data <- expand.grid(x=seq(-1.5,1.5,length.out=12L),group=factor(seq_len(4L)))
    eta <- .25+.7*data$x+rnorm(4L,sd=.3)[data$group]
    data$y <- rbinom(nrow(data),1,plogis(eta))
    fit <- .mmes_specification_fit(y~0+a+b*x,data=data,family=binomial(),
      parameters=list(a=parameter(random=~vsm(ism(group),sigma2=.09,
        fixedSigma2=TRUE)),b=parameter()),
      start=c(a=.2,b=.6),control=phase11SpecificationControl())
    rows <- split(seq_len(nrow(data)),data$group)
    reference <- function(par,details=FALSE){
      predictor <- par[1L]+par[2L]*data$x
      mode <- vapply(rows,function(selected)
        uniroot(function(u) sum(plogis(predictor[selected]+u)-data$y[selected])+u/.09,
          c(-10,10),tol=1e-13)$root,numeric(1))
      predictor <- predictor+mode[data$group]
      probability <- plogis(predictor)
      curvature <- vapply(rows,function(selected)
        sum(probability[selected]*(1-probability[selected]))+1/.09,numeric(1))
      value <- -sum(dbinom(data$y,1,probability,log=TRUE))+
        sum(mode^2)/.18+length(mode)*log(.09)/2+sum(log(curvature))/2
      if(details) list(value=value,curvature=curvature) else value
    }
    checked <- post_nonlinear(fit)
    information <- unname(optimHess(fit$optimizer$par,reference))
    accepted <- reference(fit$optimizer$par,TRUE)
    expect_equal(-as.numeric(logLik(fit)),accepted$value,tolerance=1e-9)
    expect_equal(checked$nonlinearInference$information,information,tolerance=5e-5)
    expect_equal(unname(vcov(checked)),solve(information),tolerance=5e-5)
    expect_equal(unname(vcov(checked,type="conditional")),
      diag(1/accepted$curvature),tolerance=1e-8)
  })
}
