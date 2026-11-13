for(case in c("gaussian","logit","probit","cloglog","poisson","Gamma","nbinom2","ziPoisson")){
  local({
    code <- case
    test_that(paste("nonlinear specification supports",code),{
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
      control <- mmes_control(optimizer=list(maxit=150L),
        convergence=list(likelihoodTol=1e-9,parameterTol=1e-8),
        inference=list(step=5e-4),
        diagnostics=list(verbose=FALSE,dateWarning=FALSE))
      reference <- mmes(y~0+a+b*x,data=data,family=family,
        parameters=list(a=parameter(random=~vsm(ism(group),sigma2=.06,
          fixedSigma2=TRUE)),b=parameter()),
        start=c(a=.15,b=.3),method=mmes_method(objective="ml"),
        familyParameters=settings,control=control)
      expect_true(reference$convergence)
      expect_true(is.finite(as.numeric(logLik(reference))))
      for(target in c("zero","population"))
        expect_true(all(is.finite(predict(reference,target=target))))
      expect_true(all(is.finite(vcov(reference))))
      expect_equal(predict(unserialize(serialize(reference,NULL))),
        predict(reference),tolerance=1e-12)
      expect_equal(dim(simulate(reference,seed=130)),c(nrow(data),1L))
    })
  })
}
