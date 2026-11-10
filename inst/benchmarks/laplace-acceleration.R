# Run after pkgload::load_all(); verifies parity before reporting CPU timings.
local({
  library(sommer)
  scope <- sommer:::.mmes_blas_scope(1L)
  on.exit(scope$restore(),add=TRUE)
  hadSeed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadSeed) seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(hadSeed) assign(".Random.seed",seed,envir=.GlobalEnv) else
    rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(184)
  control <- function(warmStart="off",gradient="analytic",trace="panels",blocks="off")
    sommer:::.mmes_control_internal(task=list(objective="ml",
      glmm=list(method="laplace",backend="random",warmStart=warmStart,
        gradient=gradient,trace=trace,blocks=blocks,maxit=1000L,modeTol=1e-10)),
      execution=list(blasThreads=1L),diagnostics=list(verbose=FALSE,dateWarning=FALSE))
  cases <- data.frame(case=c("independent-large","chain-large","dense-small"),
    n=c(1536L,1536L,96L),q=c(256L,256L,48L))
  native <- lapply(seq_len(nrow(cases)),function(index){
    case <- cases[index,]
    data <- data.frame(y=1,id=factor(seq_len(case$q)))
    settings <- control()
    settings$output$setup <- TRUE
    prepared <- sommer:::.mmes_dispatch(y~1,random=~vsm(ism(id),sigma2=.4),data=data,
      family=poisson(),control=settings)
    fixed <- Matrix::Matrix(cbind(1,rnorm(case$n)),sparse=TRUE)
    z <- if(case$case=="dense-small")
      Matrix::Matrix(matrix(rnorm(case$n*case$q)/sqrt(case$q),case$n,case$q),sparse=TRUE) else
      Matrix::sparseMatrix(i=seq_len(case$n),j=rep(seq_len(case$q),length.out=case$n),
        x=1,dims=c(case$n,case$q))
    q <- if(case$case=="chain-large")
      Matrix::bandSparse(case$q,k=c(-1,0,1),
        diagonals=list(rep(-.2,case$q-1),rep(1.4,case$q),rep(-.2,case$q-1))) else
      Matrix::Diagonal(case$q)
    response <- rpois(case$n,exp(.3+.2*fixed[,2]))
    create <- function(warmStart="off",trace="panels",blocks="off")
      sommer:::sommer_laplace_objective_create(fixed,z,response,rep(0,case$n),
        "poisson",prepared$covStruct[1],list(q),1000,backend="random",
        warmStart=warmStart,trace=trace,blocks=blocks)
    choices <- list(cold=create(),warm=create(warmStart="auto"),
      selected=create(trace="selected"),
      blocks=create(blocks="on"),blocksSelected=create(blocks="on",trace="selected"))
    proposals <- lapply(seq_len(30),function(step) c(.3+step*.002,.2,log(.4)+step*.003))
    reference <- lapply(proposals,function(par)
      sommer:::sommer_laplace_objective_eval(choices$cold,par,tolerance=1e-10))
    results <- lapply(names(choices),function(name){
      context <- choices[[name]]
      before <- sommer:::sommer_laplace_objective_eval(context,tail(proposals,1)[[1]],
        tolerance=1e-10)$work
      samples <- replicate(3,system.time(for(step in seq_along(proposals)){
        value <- sommer:::sommer_laplace_objective_eval(context,proposals[[step]],tolerance=1e-10)
        stopifnot(value$valid,abs(value$logLik-reference[[step]]$logLik)<1e-8,
          max(abs(value$gradient-reference[[step]]$gradient))<1e-7)
      })[["elapsed"]]/length(proposals))
      work <- sommer:::sommer_laplace_objective_eval(context,tail(proposals,1)[[1]],
        tolerance=1e-10)$work
      data.frame(case=case$case,n=case$n,q=case$q,choice=name,
        seconds=median(samples),
        modeIterations=work$warmStart$modeIterations-before$warmStart$modeIterations,
        selectedBuilds=work$backend$selectedBuilds-before$backend$selectedBuilds,
        panels=work$backend$contractionSolves-before$backend$contractionSolves,
        blocks=work$backend$blocks)
    })
    do.call(rbind,results)
  })
  native <- do.call(rbind,native)
  native$speedup <- ave(native$seconds,native$case,FUN=function(x) x[1]/x)
  print(native,row.names=FALSE,digits=5)

  data <- expand.grid(id=factor(1:40),replicate=1:16)
  data$x <- rnorm(nrow(data))
  eta <- .3+.25*data$x+rnorm(40,sd=.5)[data$id]
  data$y <- rpois(nrow(data),exp(eta))
  full <- function(settings) sommer:::.mmes_dispatch(y~x,random=~vsm(ism(id),sigma2=.3),
    data=data,family=poisson(),control=settings)
  choices <- list(cold=control(),warm=control(warmStart="auto"),
    selected=control(trace="selected"),blocks=control(blocks="on"))
  reference <- full(choices$cold)
  fullTimes <- lapply(names(choices),function(name){
    times <- numeric(3)
    for(batch in seq_len(3)) times[batch] <- system.time({
      value <- full(choices[[name]])
      stopifnot(value$convergence,abs(value$llik-reference$llik)<1e-7,
        max(abs(value$b-reference$b))<1e-5)
    })[["elapsed"]]
    data.frame(case="whole-poisson-fit",choice=name,seconds=median(times),
      modeIterations=value$laplaceDiagnostics$nativeWork$warmStart$modeIterations)
  })
  fullTimes <- do.call(rbind,fullTimes)
  fullTimes$speedup <- fullTimes$seconds[1]/fullTimes$seconds
  print(fullTimes,row.names=FALSE,digits=5)

  base <- data[data$replicate<=4,]
  mixed <- do.call(rbind,lapply(c("count","ord","nom"),function(trait){
    out <- base
    out$trait <- trait
    out$y <- switch(trait,count=base$y,
      ord=as.integer(cut(.3+.25*base$x+rlogis(nrow(base)),c(-Inf,-.5,.7,Inf))),
      nom=sample.int(3,nrow(base),replace=TRUE))
    out
  }))
  mixed$trait <- factor(mixed$trait,levels=c("count","ord","nom"))
  fitMixed <- function(gradient) sommer:::.mmes_dispatch(y~trait+trait:x-1,
    random=~vsm(ism(id),sigma2=.3,fixedSigma2=TRUE),data=mixed,
    family=familym("trait",count=poisson(),ord=ordinal(),nom=multinomial()),
    control=control(gradient=gradient))
  fit <- fitMixed("analytic")
  parameters <- fit$optimizer$par
  context <- sommer:::.laplace_objective_context(fit)
  evaluate <- function(parameters,gradient=TRUE)
    sommer:::.laplace_constrained_evaluate(parameters,fit$laplaceParameterLayout,
      function(expanded) sommer:::sommer_laplace_objective_eval(context,expanded,
        tolerance=1e-10,gradient=gradient))
  proposals <- lapply(seq_len(10),function(step){
    par <- parameters
    par[1] <- par[1]+step*.003
    par
  })
  scores <- lapply(proposals,function(par) evaluate(par,TRUE)$gradient)
  gradientTimes <- lapply(c("analytic","numeric"),function(choice){
    times <- replicate(3,system.time(for(step in seq_len(10)){
      par <- proposals[[step]]
      value <- evaluate(par,choice=="analytic")
      score <- if(choice=="analytic") value$gradient else
        sommer:::.laplace_numeric_gradient(par,function(candidate) evaluate(candidate,FALSE)$logLik)
      stopifnot(max(abs(score-scores[[step]]))<1e-5)
    })[["elapsed"]]/10)
    data.frame(case="mixed-objective-score",choice=choice,seconds=median(times))
  })
  gradientTimes <- do.call(rbind,gradientTimes)
  gradientTimes$speedup <- gradientTimes$seconds[2]/gradientTimes$seconds
  print(gradientTimes,row.names=FALSE,digits=5)
  reference <- fitMixed("analytic")
  mixedTimes <- lapply(c("analytic","numeric"),function(choice){
    times <- replicate(3,system.time({
      value <- fitMixed(choice)
      stopifnot(value$convergence,abs(value$llik-reference$llik)<1e-6,
        max(abs(value$b-reference$b))<1e-4)
    })[["elapsed"]])
    data.frame(case="whole-mixed-fit",choice=choice,seconds=median(times))
  })
  mixedTimes <- do.call(rbind,mixedTimes)
  mixedTimes$speedup <- mixedTimes$seconds[2]/mixedTimes$seconds
  print(mixedTimes,row.names=FALSE,digits=5)
  cat("\nSingle BLAS thread; median of three batches; native batches have thirty proposals.\n")
  cat("Parity is checked for every timed proposal/fit; timings are fixture-specific, not guarantees.\n")
})
