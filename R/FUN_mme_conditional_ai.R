.conditional_fa_coordinates <- function(model, state){
  blocks <- list()
  add <- function(field, index, transform="identity"){
    blocks[[length(blocks)+1L]] <<- list(field=field, index=index,
      transform=transform)
  }
  if(is.null(model$shared)){
    if(model$free[1L]) add("genetic", seq_along(state$genetic))
    if(model$free[2L] && model$ke) add("residual", seq_along(state$residual))
  }else{
    add("loading", seq_along(state$loading))
    add("weightG", seq_along(state$weightG), "logit")
  }
  if(model$free[1L]) add("specificG",
    if(model$homogeneousG) 1L else seq_len(model$p), "log")
  if(model$free[2L]) add("specificE",
    if(model$homogeneousE) 1L else seq_len(model$p), "log")
  values <- unlist(lapply(blocks, function(block){
    x <- state[[block$field]][block$index]
    switch(block$transform, log=log(x), logit=stats::qlogis(x), x)
  }), use.names=FALSE)
  list(blocks=blocks, values=values)
}

.conditional_fa_decode <- function(model, state, coordinates, values){
  start <- 0L
  for(block in coordinates$blocks){
    count <- length(block$index)
    x <- values[start+seq_len(count)]
    x <- switch(block$transform, log=exp(x), logit=stats::plogis(x), x)
    state[[block$field]][block$index] <- x
    start <- start+count
  }
  if(model$homogeneousG) state$specificG[] <- state$specificG[1L]
  if(model$homogeneousE) state$specificE[] <- state$specificE[1L]
  if(!is.null(model$shared)) state$weightE <- 1-state$weightG
  state
}

.conditional_fa_ai_bounds <- function(model, coordinates, engine){
  tolerance <- max(1e-8,engine$conditionalInverseTol)
  if(tolerance >= .5){
    stop("Conditional FA AI requires execution$inverseTol < 0.5.",call.=FALSE)
  }
  lower <- upper <- numeric()
  for(block in coordinates$blocks){
    count <- length(block$index)
    lower <- c(lower,rep(switch(block$transform,
      log=log(tolerance*model$responseVariance),
      logit=stats::qlogis(tolerance), -Inf),count))
    upper <- c(upper,rep(if(block$transform == "logit")
      stats::qlogis(1-tolerance) else Inf,count))
  }
  list(lower=lower,upper=upper,varianceFloor=tolerance*model$responseVariance,
    allocationFloor=tolerance, feasible=all(coordinates$values >= lower-1e-12 &
      coordinates$values <= upper+1e-12))
}

.conditional_fa_prior_solve <- function(model, state, rhs){
  out <- numeric(model$d)
  specific <- matrix(rhs[model$specific], model$r, model$p)
  genetic <- matrix(rhs[model$genetic], model$r, model$kg)
  solved <- as.matrix(Matrix::solve(model$relationshipFactor, cbind(specific, genetic)))
  out[model$specific] <- sweep(solved[,seq_len(model$p),drop=FALSE],
    2L, state$specificG, `*`)
  dg <- if(is.null(model$shared)) rep(1,model$kg) else state$weightG
  out[model$genetic] <- sweep(solved[,model$p+seq_len(model$kg),drop=FALSE],
    2L, dg, `*`)
  de <- if(is.null(model$shared)) rep(1,model$ke) else state$weightE
  out[model$residual] <- sweep(matrix(rhs[model$residual],model$m,model$ke),
    2L, de, `*`)
  out
}

.conditional_fa_derivative <- function(model, state, system, block, index){
  D <- Matrix::Matrix(0, length(model$y), model$d, sparse=TRUE)
  Q <- Matrix::Matrix(0, model$d, model$d, sparse=TRUE)
  ratio <- numeric(length(model$y))
  logQ <- 0
  field <- block$field
  if(field %in% c("genetic", "residual", "loading")){
    trait <- (index-1L) %% model$p+1L
    factor <- (index-1L) %/% model$p+1L
    S <- switch(field, genetic=model$selectorsG[[factor]],
      residual=model$selectorsE[[factor]],
      loading=model$selectorsG[[factor]]+model$selectorsE[[factor]])
    D <- Matrix::Diagonal(x=as.numeric(model$trait == trait)) %*% S
  }else if(field == "specificE"){
    ratio[] <- if(model$homogeneousE) 1 else as.numeric(model$trait == index)
  }else if(field == "specificG"){
    traits <- if(model$homogeneousG) seq_len(model$p) else index
    for(j in traits){
      indices <- model$specific[,j]
      Q[indices,indices] <- -model$Ai/state$specificG[j]
    }
    logQ <- -model$r*length(traits)
  }else if(field == "weightG"){
    w <- state$weightG[index]
    g <- model$genetic[,index]; e <- model$residual[,index]
    Q[g,g] <- -(1-w)*model$Ai/w
    Q[e,e] <- Matrix::Diagonal(model$m,w/(1-w))
    logQ <- -model$r*(1-w)+model$m*w
  }
  WD <- Matrix::Diagonal(x=system$weight) %*% system$design
  H <- Matrix::crossprod(D,WD)+Matrix::crossprod(WD,D) -
    Matrix::crossprod(system$design,Matrix::Diagonal(x=system$weight*ratio) %*%
      system$design)+Q
  list(design=D, prior=Q, ratio=ratio, logQ=logQ, H=H)
}

.conditional_fa_ai <- function(model, state, posterior, engine){
  coordinates <- .conditional_fa_coordinates(model,state)
  q <- length(coordinates$values); n <- length(model$y)
  # Admit the retained derivative RHS, projection panels and dense curvature.
  bytes <- 8*(4*n*q+4*model$d*min(q,engine$experimental$panelCi2)+
    12*q*q+4*model$r*(model$p+model$kg))
  if(!is.finite(bytes) || bytes > engine$denseMemoryMB*1048576){
    stop("Conditional FA AI workspace exceeds execution$denseMemoryMB; reduce rank/traits or use conditional EM.",call.=FALSE)
  }
  system <- .conditional_fa_system(model,state)
  mean <- as.numeric(posterior$mean)
  residual <- model$y-as.numeric(system$design %*% mean)
  Py <- system$weight*residual
  randomMean <- mean; randomMean[seq_len(model$b)] <- 0
  B <- matrix(0,n,q); score <- numeric(q)
  column <- 0L
  for(block in coordinates$blocks) for(index in block$index){
    column <- column+1L
    derivative <- .conditional_fa_derivative(model,state,system,block,index)
    entries <- summary(as(derivative$H,"generalMatrix"))
    trace <- sum(posterior$covariance[cbind(entries$i,entries$j)]*entries$x)
    quadratic <- -sum(system$weight*derivative$ratio*residual^2) -
      2*sum(Py*as.numeric(derivative$design %*% mean)) +
      sum(mean*as.numeric(derivative$prior %*% mean))
    score[column] <- -.5*(sum(derivative$ratio)-derivative$logQ+trace+quadratic)
    rhs <- as.numeric(Matrix::crossprod(derivative$design,Py)) -
      as.numeric(derivative$prior %*% randomMean)
    B[,column] <- as.numeric(derivative$design %*% randomMean +
      system$design %*% .conditional_fa_prior_solve(model,state,rhs)) +
      state$specificE[model$trait]*derivative$ratio*Py
  }
  information <- matrix(0,q,q)
  for(start in seq.int(1L,q,by=engine$experimental$panelCi2)){
    columns <- start:min(q,start+engine$experimental$panelCi2-1L)
    weighted <- system$weight*B[,columns,drop=FALSE]
    rhs <- as.matrix(Matrix::crossprod(system$design,weighted))
    solved <- mme_conditional_solve_cpp(posterior$workspace,rhs,engine$denseMemoryMB)
    projected <- weighted-system$weight*as.matrix(system$design %*% solved)
    information[,columns] <- .5*crossprod(B,projected)
  }
  if(any(!is.finite(score)) || any(!is.finite(information))){
    stop("Non-finite conditional FA REML score or AI curvature.",call.=FALSE)
  }
  list(score=score, information=(information+t(information))/2,
    coordinates=coordinates, workspaceBytes=bytes)
}

.conditional_fa_ai_step <- function(model, state, posterior, engine, parameterTol){
  ai <- .conditional_fa_ai(model,state,posterior,engine)
  bounds <- .conditional_fa_ai_bounds(model,ai$coordinates,engine)
  values <- ai$coordinates$values
  if(!bounds$feasible){
    stop("Conditional FA AI starting variances/allocations fall below its numerical bounds; revise starting values or execution$inverseTol.",call.=FALSE)
  }
  scaling <- sqrt(pmax(diag(ai$information),.Machine$double.eps))
  active <- (values <= bounds$lower+1e-10 & ai$score<=0) |
    (values >= bounds$upper-1e-10 & ai$score>=0)
  direction <- numeric(length(ai$score))
  rank <- 0L
  working <- active
  for(attempt in seq_len(length(direction)+1L)){
    free <- which(!working)
    if(!length(free)) break
    standardized <- ai$information[free,free,drop=FALSE]/outer(scaling[free],scaling[free])
    decomposition <- eigen(standardized,symmetric=TRUE)
    retained <- decomposition$values > max(1e-10,max(decomposition$values)*1e-9)
    rank <- sum(retained)
    if(any(retained)){
      vectors <- decomposition$vectors[,retained,drop=FALSE]
      direction[free] <- as.numeric(vectors %*%
        (crossprod(vectors,ai$score[free]/scaling[free])/
          decomposition$values[retained]))/scaling[free]
    }
    outward <- (values <= bounds$lower+1e-10 & direction<0) |
      (values >= bounds$upper-1e-10 & direction>0)
    if(!any(outward)) break
    working <- working | outward
    direction[] <- 0
  }
  direction <- direction/max(1,max(abs(direction)))
  projectedScore <- ai$score
  projectedScore[active] <- 0
  scoreNorm <- max(abs(projectedScore)/pmax(1,scaling))
  if(sum(ai$score*direction)<=0){
    direction <- projectedScore/scaling^2
    direction <- direction/max(1,max(abs(direction)))
  }
  accepted <- FALSE; evaluations <- 0L
  if(scoreNorm>parameterTol && sum(ai$score*direction)>0) for(halving in 0:20){
    step <- pmin(bounds$upper,pmax(bounds$lower,values+direction*2^-halving))-values
    updated <- .conditional_fa_decode(model,state,ai$coordinates,
      ai$coordinates$values+step)
    if(any(!is.finite(unlist(updated))) || any(updated$specificG<=0) ||
       any(updated$specificE<=0) ||
       (!is.null(model$shared) && any(updated$weightG<=0 | updated$weightG>=1))) next
    candidate <- .conditional_fa_posterior(model,updated,posterior$workspace,engine,FALSE)
    evaluations <- evaluations+1L
    if(candidate$llik >= posterior$llik &&
       candidate$llik-posterior$llik >= 1e-4*sum(ai$score*step)){
      accepted <- TRUE
      break
    }
  }
  if(!accepted){
    updated <- state; candidate <- posterior
    step <- numeric(length(direction))
  }
  list(state=updated, posterior=candidate, accepted=accepted,
    stationary=scoreNorm<=parameterTol,
    relativeChange=sqrt(sum(step^2))/(1+sqrt(sum(ai$coordinates$values^2))),
    scoreNorm=scoreNorm, rank=rank, evaluations=evaluations,
    activeBounds=which(active),varianceFloor=bounds$varianceFloor,
    allocationFloor=bounds$allocationFloor,workspaceBytes=ai$workspaceBytes)
}
