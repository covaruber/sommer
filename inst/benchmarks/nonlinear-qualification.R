# Run against an installed package; timings are descriptive, not CI thresholds.
local({
  library(sommer)
  scope <- sommer:::.mmes_blas_scope(1L)
  on.exit(scope$restore(),add=TRUE)
  hadSeed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadSeed) seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(hadSeed) assign(".Random.seed",seed,envir=.GlobalEnv) else
    if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  cases <- data.frame(groups=c(16L,64L,64L,64L,256L,512L),
    points=c(8L,8L,32L,128L,16L,16L))
  results <- profiles <- list()
  for(path in c("exact-gaussian","curve-gaussian","curve-poisson")){
    for(index in seq_len(nrow(cases))){
      case <- cases[index,]
      set.seed(903L+index)
      data <- expand.grid(x=seq(0,5,length.out=case$points),
        group=factor(seq_len(case$groups)))
      deviation <- rnorm(case$groups,sd=.2)[data$group]
      randomCurve <- path!="exact-gaussian"
      mean <- if(randomCurve) 1.4*exp(-.65*exp(deviation)*data$x) else
        1.4*exp(-.65*data$x)+deviation
      data$y <- if(path=="curve-poisson") rpois(nrow(data),exp(mean)) else
        mean+rnorm(nrow(data),sd=.12)
      settings <- sommer:::.mmes_control_internal(task=list(objective="ml",nonlinear=list(
        start=c(asymptote=1.4,rate=.65))),
        convergence=list(maxIters=500L,likelihoodTol=1e-9,parameterTol=1e-7),
        execution=list(blasThreads=1L,denseMemoryMB=256),
        diagnostics=list(verbose=FALSE,dateWarning=FALSE,profile="force"))
      model <- nonlinear_curve("x",random=if(randomCurve)
        list(rate=~vsm(ism(group),sigma2=.04,fixedSigma2=TRUE)) else list())
      timing <- system.time(fit <- sommer:::.mmes_dispatch(y~0,
        random=if(randomCurve) NULL else ~group,data=data,nonlinear=model,
        family=if(path=="curve-poisson") poisson() else gaussian(),control=settings))
      if(!isTRUE(fit$convergence) || !isTRUE(fit$nonlinearDiagnostics$identified))
        stop("Qualification did not converge: ",path,", case ",index,".")
      if(max(abs(fit$nonlinearParameters-c(1.4,.65)))>.35)
        stop("Qualification curve recovery is outside the declared tolerance.")
      if(!inherits(fit$fixedDesign,"sparseMatrix") ||
         !inherits(fit$randomDesign,"sparseMatrix"))
        stop("Qualification unexpectedly densified a fitted design.")
      if(!is.null(fit$nonlinearInference) ||
         isTRUE(fit$nonlinearDiagnostics$uncertaintyComputed))
        stop("Qualification unexpectedly computed inference during fitting.")
      profile <- profile_nonlinear(fit,replications=3L)
      name <- paste(path,index,sep="-")
      profiles[[name]] <- profile
      measured <- stats::setNames(profile$table$seconds,profile$table$phase)
      results[[name]] <- data.frame(path=path,n=nrow(data),groups=case$groups,
        fixedParameters=length(fit$nonlinearCoordinates),randomFull=length(fit$randomMode),
        randomLatent=length(fit$latentRandomMode),
        randomNonzeros=Matrix::nnzero(fit$randomDesign),
        fitSeconds=unname(timing[["elapsed"]]),
        preparationSeconds=unname(profile$fitProfile$seconds["preparation"]),
        outerSeconds=unname(profile$fitProfile$seconds["outerOptimization"]),
        curveSeconds=unname(measured["curveMean"]),
        curveJacobianSeconds=unname(measured["curveJacobian"]),
        conditionalCurvatureSeconds=unname(measured["conditionalCurvature"]),
        marginalObjectiveSeconds=unname(measured["marginalObjective"]),
        marginalScoreSeconds=unname(measured["marginalScore"]),
        fitBytes=profile$memory$fitBytes,
        designBytes=profile$memory$sparseDesignBytes,
        iterations=fit$optimizer$iterations,
        modeIterations=if(is.null(fit$nonlinearDiagnostics$modeIterations)) 0L else
          fit$nonlinearDiagnostics$modeIterations,
        maxGradient=max(abs(fit$nonlinearDiagnostics$gradient)),
        asymptote=unname(fit$nonlinearParameters[1L]),
        rate=unname(fit$nonlinearParameters[2L]))
    }
  }
  results <- do.call(rbind,results)
  rownames(results) <- NULL
  print(results,row.names=FALSE,digits=5)
  environment <- list(R=R.version.string,platform=R.version$platform,
    BLAS=unname(extSoftVersion()["BLAS"]),compiler=system("R CMD config CXX",intern=TRUE))
  print(environment)
  invisible(list(results=results,profiles=profiles,environment=environment,
    caveats=c("Inclusive overlapping timings; not additive.",
      "R object sizes, not native peak RSS.",
      "One measured Linux/R/BLAS/compiler configuration, not portable timing guarantees.",
      "Traditional curves only; no layers, backpropagation or extra random curve parameters.")))
})
