# Standalone qualification: larger recovery fits are not automated tests.
local({
  library(sommer)
  scope <- sommer:::.mmes_blas_scope(1L)
  on.exit(scope$restore(),add=TRUE)
  hadSeed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadSeed) seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(hadSeed) assign(".Random.seed",seed,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  results <- list()
  covariance <- .01*matrix(c(1,.4,.4,1),2L)
  for(groups in c(8L,16L,32L)){
    set.seed(442L+groups)
    data <- expand.grid(x=seq(0,5,length.out=16L),group=factor(seq_len(groups)))
    deviation <- matrix(rnorm(2L*groups),groups,2L)%*%chol(covariance)
    truth <- 2*exp(deviation[data$group,1L])*
      exp(-.6*exp(deviation[data$group,2L])*data$x)
    data$y <- truth+rnorm(nrow(data),sd=.1)
    model <- nonlinear_curve("x",joint=list(parameters=c("asymptote","rate"),
      random=~vsm(usm(.parameter,theta=matrix(c(1,.4,.4,1),2L)),
        ism(group),sigma2=.01)))
    control <- sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
      start=c(asymptote=2,rate=.6))),
      convergence=list(maxIters=500L,likelihoodTol=1e-10,parameterTol=1e-8),
      execution=list(blasThreads=1L,denseMemoryMB=256),
      diagnostics=list(verbose=FALSE,dateWarning=FALSE,profile="force"))
    elapsed <- system.time(fit <- sommer:::.mmes_dispatch(y~0,data=data,nonlinear=model,control=control))
    if(!isTRUE(fit$convergence) || !isTRUE(fit$nonlinearDiagnostics$identified))
      stop("Joint reference fit did not converge or is unidentified.")
    recovery <- max(abs(fit$nonlinearParameters-c(2,.6)))
    predictionRMSE <- sqrt(mean((as.numeric(fitted(fit))-truth)^2))
    if(recovery>.15 || predictionRMSE>.1)
      stop("Joint recovery/prediction exceeded its declared tolerance.")
    if(!is.null(fit$nonlinearInference) ||
       !inherits(fit$randomDesign,"sparseMatrix"))
      stop("Joint qualification computed fit-time uncertainty or densified its design.")
    results[[as.character(groups)]] <- data.frame(groups=groups,n=nrow(data),
      seconds=unname(elapsed["elapsed"]),
      recovery=recovery,predictionRMSE=predictionRMSE,logLik=fit$llik,
      fitBytes=as.numeric(utils::object.size(fit)),
      asymptote=unname(fit$nonlinearParameters[1L]),
      rate=unname(fit$nonlinearParameters[2L]))
  }
  results <- do.call(rbind,results)
  rownames(results) <- NULL
  print(results,row.names=FALSE)
  invisible(list(results=results,environment=list(R=R.version.string,
    platform=R.version$platform,BLAS=unname(extSoftVersion()["BLAS"])),
    caveats=c("Recovery and predictions are compared with generating truth.",
      "Elapsed times are descriptive; retained R sizes are not native peak RSS.")))
})
