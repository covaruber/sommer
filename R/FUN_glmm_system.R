.laplace_backend <- function(object){
  selected <- object$laplaceDiagnostics$nativeWork$backend$selected
  if(is.null(selected)) "random" else selected
}

.laplace_uses_native_solve <- function(object){
  .laplace_backend(object)=="observation" ||
    identical(object$laplacePostPEVSettings$blocks,"on")
}

.laplace_objective_context <- function(object, response=object$laplaceNativeResponse,
    memoryMB=object$laplacePostPEVSettings$denseMemoryMB){
  settings <- utils::modifyList(.mme_glmm_defaults,object$laplacePostPEVSettings)
  sommer_laplace_objective_create(object$fixedDesign,object$randomDesign,response,
    object$offset,"",object$covStruct,object$randomRelationships,memoryMB,
    familySpecs=object$familySpecifications,familyIndex=object$laplaceFamilyIndex,
    trials=object$priorWeights,residualModel=object$laplaceResidual,
    backend=.laplace_backend(object),warmStart=settings$warmStart,
    trace=settings$trace,blocks=settings$blocks)
}

.laplace_mode_solver <- function(object, memoryMB=object$laplacePostPEVSettings$denseMemoryMB){
  if(!.laplace_uses_native_solve(object)){
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(object$laplaceModeHessian),
      LDL=FALSE,perm=TRUE)
    return(function(rhs) as.matrix(Matrix::solve(factor,rhs,system="A")))
  }
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  settings <- object$laplacePostPEVSettings
  if(.laplace_is_nominal(object)){
    prepared <- object$nominalPreparation
    prepared$control$execution$denseMemoryMB <- memoryMB
    prepared$control$task$glmm$backend <- .laplace_backend(object)
    objective <- .laplace_nominal_objective(prepared)
    value <- objective$evaluate(object$optimizer$par,gradient=FALSE)
    context <- objective$context
    solve <- sommer_multinomial_objective_solve
  }else{
    context <- .laplace_objective_context(object,memoryMB=memoryMB)
    value <- .laplace_constrained_evaluate(object$optimizer$par,
      object$laplaceParameterLayout,function(parameters)
        sommer_laplace_objective_eval(context,parameters,
          maxModeIterations=settings$modeMaxit,tolerance=settings$modeTol,gradient=FALSE))
    solve <- sommer_laplace_objective_solve
  }
  if(!isTRUE(value$valid))
    stop("Post-fit Laplace system could not reconstruct the accepted mode: ",
      value$reason,".",call.=FALSE)
  function(rhs){
    scope <- .mmes_blas_scope(object$control$execution$blasThreads)
    on.exit(scope$restore(),add=TRUE)
    result <- solve(context,as.matrix(rhs))
    colnames(result) <- colnames(rhs)
    result
  }
}

.laplace_joint_solver <- function(object, weightMatrix, memoryMB){
  fixed <- object$fixedDesign
  random <- object$randomDesign
  fixedHessian <- Matrix::crossprod(fixed,weightMatrix %*% fixed)
  cross <- Matrix::crossprod(fixed,weightMatrix %*% random)
  if(!.laplace_uses_native_solve(object)){
    joint <- rbind(cbind(fixedHessian,cross),
      cbind(Matrix::t(cross),object$laplaceModeHessian))
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(joint),LDL=FALSE,perm=TRUE)
    return(function(rhs) as.matrix(Matrix::solve(factor,rhs,system="A")))
  }
  p <- ncol(fixed)
  q <- ncol(random)
  if(8*(6*as.double(p)*q+4*as.double(p)^2)>memoryMB*1048576)
    stop("Laplace fixed-effect Schur workspace exceeds execution$denseMemoryMB.",call.=FALSE)
  solveMode <- .laplace_mode_solver(object,memoryMB)
  applied <- if(p) solveMode(Matrix::t(cross)) else matrix(0,q,0)
  schur <- as.matrix(fixedHessian-cross %*% applied)
  factor <- if(p) chol((schur+t(schur))/2) else NULL
  function(rhs){
    rhs <- as.matrix(rhs)
    if(nrow(rhs)!=p+q)
      stop("Laplace joint solve has incompatible right-hand sides.",call.=FALSE)
    randomSolution <- solveMode(rhs[p+seq_len(q),,drop=FALSE])
    if(!p) return(randomSolution)
    fixedRhs <- rhs[seq_len(p),,drop=FALSE]-cross %*% randomSolution
    fixedSolution <- backsolve(factor,forwardsolve(t(factor),as.matrix(fixedRhs)))
    result <- rbind(fixedSolution,randomSolution-applied %*% fixedSolution)
    colnames(result) <- colnames(rhs)
    result
  }
}
