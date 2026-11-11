# Run against the current installation, or source after pkgload::load_all().
local({
  if(!"sommer" %in% loadedNamespaces()) library(sommer)
  replications <- getOption("sommer.phase11.replications",100L)
  if(length(replications)!=1L || !is.finite(replications) ||
     replications<1L || replications!=as.integer(replications))
    stop("sommer.phase11.replications must be a positive integer.")
  hadSeed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadSeed) previous <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(hadSeed) assign(".Random.seed",previous,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  control <- mmes_control(optimizer=list(maxit=200L),
    convergence=list(likelihoodTol=1e-9,parameterTol=1e-8),
    execution=list(denseMemoryMB=256),
    diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  graph <- nonlinear_graph(list(
    x=list(op="input",column="x"),
    mean=list(op="affine",inputs="x",weights=matrix(1L),bias=2L),
    output=list(op="identity",input="mean")),
    output="output",parameters=c("slope","intercept"),
    parameterMap=list(slope=list(type="fixed"),
      intercept=list(type="random",group="group")),
    random=list(intercept=~vsm(ism(group),sigma2=.09,fixedSigma2=TRUE)))
  rows <- list()
  profiles <- list()
  for(code in c("gaussian","binomial","poisson")){
    family <- switch(code,gaussian=gaussian(),binomial=binomial(),poisson=poisson())
    settings <- if(code=="gaussian") list(dispersion=.16,estimate=FALSE) else list()
    for(replication in seq_len(replications)){
      seed <- 73100L+replication+1000L*match(code,c("gaussian","binomial","poisson"))
      set.seed(seed)
      data <- expand.grid(x=seq(-1.5,1.5,length.out=12L),group=factor(seq_len(4L)))
      eta <- .25+.7*data$x+rnorm(4L,sd=.3)[data$group]
      data$y <- switch(code,gaussian=rnorm(nrow(data),eta,.4),
        binomial=rbinom(nrow(data),1,plogis(eta)),poisson=rpois(nrow(data),exp(eta)))
      for(route in c("native","reference")){
        start <- proc.time()[["elapsed"]]
        outcome <- tryCatch({
          fit <- if(route=="native")
            mmes(y~0+graph_term("net"),data=data,family=family,
              predictors=list(net=graph),start=c(slope=.6,intercept=.2),
              familyParameters=settings,method=mmes_method(objective="ml"),control=control) else
            mmes(y~0+a+b*x,data=data,family=family,
              parameters=list(a=parameter(random=~vsm(ism(group),sigma2=.09,
                fixedSigma2=TRUE)),b=parameter()),
              start=c(a=.2,b=.6),familyParameters=settings,
              method=mmes_method(objective="ml"),control=control)
          if(!isTRUE(fit$convergence)) stop("Fit failed its convergence certificate.")
          prediction <- predict(fit,newdata=data.frame(x=0,
            group=factor("1",levels=levels(data$group))),
            target="zero",type="link",se.fit=TRUE,interval="confidence")
          if(replication==1L) profiles[[paste(code,route,sep="-")]] <-
            profile_nonlinear(fit,replications=1L,inference=TRUE)
          list(status="ok",reason="",estimate=as.numeric(prediction$fit),
            se=as.numeric(prediction$se.fit),covered=
              as.numeric(prediction$lower)<=.25 && as.numeric(prediction$upper)>=.25,
            logLik=as.numeric(logLik(fit)),bytes=as.numeric(object.size(fit)))
        },error=function(e) list(status="failed",reason=conditionMessage(e),
          estimate=NA_real_,se=NA_real_,covered=NA,logLik=NA_real_,bytes=NA_real_))
        rows[[length(rows)+1L]] <- data.frame(family=code,route=route,
          replication=replication,seed=seed,seconds=proc.time()[["elapsed"]]-start,
          as.data.frame(outcome,stringsAsFactors=FALSE),stringsAsFactors=FALSE)
        message(sprintf("%s/%s replication %d/%d: %s",code,route,
          replication,replications,outcome$status))
      }
    }
  }
  results <- do.call(rbind,rows)
  summary <- do.call(rbind,lapply(split(results,list(results$family,results$route)),function(x)
    data.frame(family=x$family[1L],route=x$route[1L],n=nrow(x),
      failures=sum(x$status!="ok"),bias=mean(x$estimate-.25,na.rm=TRUE),
      rmse=sqrt(mean((x$estimate-.25)^2,na.rm=TRUE)),
      coverage=mean(x$covered,na.rm=TRUE),meanSeconds=mean(x$seconds),
      maxBytes=max(x$bytes,na.rm=TRUE))))
  print(summary,row.names=FALSE)
  if(any(results$status!="ok")) print(results[results$status!="ok",
    c("family","route","replication","seed","reason")],row.names=FALSE)
  metadata <- list(R=R.version.string,platform=R.version$platform,
    BLAS=unname(extSoftVersion()["BLAS"]),
    target="Zero-deviation link predictor at x=0, prior random-intercept variance .09 known",
    control=control,
    inference="Local observed-information/delta; Laplace outside exact Gaussian models",
    caveats=c("Object sizes are not peak resident memory.",
      "Coverage uncertainty is Monte Carlo; inspect failures separately.",
      "One local platform does not qualify other toolchains or BLAS backends."))
  print(metadata)
  output <- list(results=results,summary=summary,profiles=profiles,metadata=metadata)
  destination <- getOption("sommer.phase11.output",NULL)
  if(!is.null(destination)) saveRDS(output,destination)
  if(any(results$status!="ok")) stop("Phase 11 qualification recorded failed fits/inference.")
  if(replications>=30L &&
     any(summary$coverage<.80 | abs(summary$bias)>.25))
    stop("Phase 11 recovery/coverage is outside the declared Monte Carlo envelope.")
  invisible(output)
})
