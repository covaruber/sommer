# Run against the working-tree package loaded with pkgload::load_all().
local({
  library(sommer)
  scope <- sommer:::.mmes_blas_scope(1L)
  on.exit(scope$restore(),add=TRUE)
  hadSeed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(hadSeed) seed <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(hadSeed) assign(".Random.seed",seed,envir=.GlobalEnv) else
    rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(238)
  cases <- data.frame(case=c("dense-small-n","dense-large-q","sparse-large-q","small-q"),
    n=c(32L,64L,64L,128L),q=c(256L,512L,512L,32L),dense=c(TRUE,TRUE,FALSE,TRUE))
  results <- lapply(seq_len(nrow(cases)),function(index){
    case <- cases[index,]
    data <- data.frame(y=rep(1,case$q),id=factor(seq_len(case$q)))
    setup <- sommer:::.mmes_dispatch(y~1,random=~vsm(ism(id),sigma2=.3),data=data,
      family=poisson(),control=sommer:::.mmes_control_internal(task=list(objective="ml",
        glmm=list(method="laplace")),output=list(setup=TRUE),
        diagnostics=list(verbose=FALSE,dateWarning=FALSE)))
    design <- if(case$dense)
      Matrix::Matrix(matrix(rnorm(case$n*case$q)/sqrt(case$q),case$n,case$q),sparse=TRUE) else
      Matrix::sparseMatrix(i=seq_len(case$n),j=seq_len(case$n),x=1,
        dims=c(case$n,case$q))
    fixed <- Matrix::Matrix(cbind(1,rnorm(case$n)),sparse=TRUE)
    response <- rpois(case$n,exp(.3+.2*fixed[,2]))
    create <- function(backend)
      sommer:::sommer_laplace_objective_create(fixed,design,response,rep(0,case$n),
        "poisson",setup$covStruct[1],list(Matrix::Diagonal(case$q)),1000,
        backend=backend)
    parameters <- c(.3,.2,log(.3))
    cold <- numeric(3)
    contexts <- lapply(seq_len(3L),function(backend){
      cold[backend] <<- system.time({
        context <- create(c("random","observation","auto")[backend])
        initial <- sommer:::sommer_laplace_objective_eval(context,parameters,tolerance=1e-10)
        stopifnot(initial$valid)
      })[["elapsed"]]
      context
    })
    reference <- lapply(contexts,function(context)
      sommer:::sommer_laplace_objective_eval(context,parameters,tolerance=1e-10))
    stopifnot(all(vapply(reference,`[[`,logical(1),"valid")),
      max(abs(reference[[1]]$logLik-reference[[2]]$logLik))<1e-8,
      max(abs(reference[[1]]$gradient-reference[[2]]$gradient))<1e-7)
    timings <- vapply(seq_along(contexts),function(backend){
      samples <- replicate(3,system.time(for(iteration in seq_len(10L)){
        proposal <- parameters
        proposal[1] <- proposal[1]+iteration*.001
        value <- sommer:::sommer_laplace_objective_eval(contexts[[backend]],proposal,
          tolerance=1e-10)
        stopifnot(value$valid)
      })[["elapsed"]]/10)
      median(samples)
    },numeric(1))
    joint <- vapply(seq_along(contexts),function(backend){
      samples <- replicate(3,system.time(for(iteration in seq_len(10L)){
        proposal <- parameters+c(iteration*.001,0,iteration*.002)
        value <- sommer:::sommer_laplace_objective_eval(contexts[[backend]],proposal,
          tolerance=1e-10)
        stopifnot(value$valid)
      })[["elapsed"]]/10)
      median(samples)
    },numeric(1))
    data.frame(case=case$case,n=case$n,q=case$q,
      auto=reference[[3]]$work$backend$selected,
      randomSeconds=timings[1],observationSeconds=timings[2],
      autoSeconds=timings[3],speedup=timings[1]/timings[2],
      randomJointSeconds=joint[1],observationJointSeconds=joint[2],
      jointSpeedup=joint[1]/joint[2],randomColdSeconds=cold[1],
      observationColdSeconds=cold[2],
      likelihoodDifference=abs(reference[[1]]$logLik-reference[[2]]$logLik),
      gradientDifference=max(abs(reference[[1]]$gradient-reference[[2]]$gradient)))
  })
  print(do.call(rbind,results),row.names=FALSE,digits=5)
  cat("\nTimings: median seconds per warm objective + analytic gradient, three batches of ten.\n")
  cat("One BLAS thread; fixed covariance within each batch; setup and cold covariance build excluded.\n")
  cat("Joint timings change fixed and covariance parameters; cold timings include context construction and first evaluation.\n")
})
