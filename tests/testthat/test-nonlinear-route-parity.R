for(case in c("gaussian","logit","probit","cloglog","poisson","Gamma","nbinom2","ziPoisson")){
  local({
    code <- case
    test_that(paste("native and reference graph routes agree for",code),{
      set.seed(1850+match(code,c("gaussian","logit","probit","cloglog","poisson",
        "Gamma","nbinom2","ziPoisson")))
      data <- expand.grid(x=seq(-1.5,1.5,length.out=12L),group=factor(1:2))
      eta <- .2+.4*data$x+c(-.15,.15)[data$group]
      family <- switch(code,gaussian=gaussian(),logit=binomial("logit"),
        probit=binomial("probit"),cloglog=binomial("cloglog"),poisson=poisson(),
        Gamma=Gamma("log"),nbinom2=nbinom2(4),ziPoisson=zipoisson(.2))
      settings <- switch(code,gaussian=list(dispersion=.12,estimate=FALSE),
        Gamma=list(dispersion=.2,estimate=FALSE),
        nbinom2=list(theta=4,estimate=FALSE),
        ziPoisson=list(zeroInflation=.2,estimate=FALSE),list())
      data$y <- switch(code,gaussian=rnorm(nrow(data),eta,sqrt(.12)),
        logit=rbinom(nrow(data),1,family$linkinv(eta)),
        probit=rbinom(nrow(data),1,family$linkinv(eta)),
        cloglog=rbinom(nrow(data),1,family$linkinv(eta)),
        poisson=rpois(nrow(data),exp(eta)),
        Gamma=rgamma(nrow(data),shape=5,scale=exp(eta)/5),
        nbinom2=rnbinom(nrow(data),size=4,mu=exp(eta)),
        ziPoisson={
          value <- rpois(nrow(data),exp(eta))
          value[seq(1L,nrow(data),by=5L)] <- 0L
          value
        })
      graph <- nonlinear_graph(list(
        x=list(op="input",column="x"),
        mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
        output=list(op="identity",input="mean")),
        output="output",parameters=c("b","a"),
        parameterMap=list(b=list(type="fixed"),
          a=list(type="random",group="group")),
        random=list(a=~vsm(ism(group),sigma2=.06,fixedSigma2=TRUE)))
      control <- mmes_control(optimizer=list(maxit=150L),
        convergence=list(likelihoodTol=1e-9,parameterTol=1e-8),
        inference=list(step=5e-4),
        diagnostics=list(verbose=FALSE,dateWarning=FALSE))
      native <- mmes(y~0+graph_term("net"),data=data,family=family,
        predictors=list(net=graph),start=c(a=.15,b=.3),
        method=mmes_method(objective="ml"),familyParameters=settings,control=control)
      reference <- mmes(y~0+a+b*x,data=data,family=family,
        parameters=list(a=parameter(random=~vsm(ism(group),sigma2=.06,
          fixedSigma2=TRUE)),b=parameter()),
        start=c(a=.15,b=.3),method=mmes_method(objective="ml"),
        familyParameters=settings,control=control)
      expect_true(native$convergence)
      expect_true(reference$convergence)
      expect_equal(as.numeric(logLik(native)),as.numeric(logLik(reference)),tolerance=2e-6)
      expect_equal(as.numeric(predict(native)),as.numeric(predict(reference)),tolerance=2e-5)
      for(target in c("zero","population"))
        expect_equal(as.numeric(predict(native,target=target)),
          as.numeric(predict(reference,target=target)),tolerance=2e-5)
      expect_true(all(is.finite(vcov(native))))
      expect_true(all(is.finite(vcov(reference))))
      expect_equal(predict(unserialize(serialize(native,NULL))),
        predict(native),tolerance=1e-12)
      expect_equal(dim(simulate(native,seed=130)),c(nrow(data),1L))
      expect_equal(dim(simulate(reference,seed=130)),c(nrow(data),1L))
    })
  })
}
