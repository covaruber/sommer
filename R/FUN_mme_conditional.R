utils::globalVariables(c("par_end", "par_start"))

.conditional_fa_component <- function(cs){
  factor <- cs$factors[[1L]]
  if(identical(factor$model, "fa")){
    decomposition <- .factor_score_decomposition(within(factor, {
      par <- cs$par[par_start:par_end]
    }))
    scale <- exp(cs$par[1L]) / decomposition$normalization
    return(list(loadings=decomposition$loading * sqrt(scale),
      specific=decomposition$specific * scale))
  }
  variance <- exp(cs$par[1L])
  if(factor$model == "diag") variance <- variance * c(1, exp(cs$par[-1L]))
  list(loadings=matrix(0, cs$dim, 0L), specific=variance)
}

.conditional_fa_prepare <- function(X, Z, Ai, y, block, coordinate, structures){
  fail <- function(message) stop("Conditional FA ", message, call.=FALSE)
  if(length(structures) != 2L || length(Ai) != 1L){
    fail("requires exactly one random covariance term.")
  }
  g <- structures[[1L]]; e <- structures[[2L]]
  if(length(g$factors) != 1L || g$factors[[1L]]$model != "fa" ||
     length(e$factors) > 1L ||
     (length(e$factors) && !e$factors[[1L]]$model %in% c("fa", "diag", "identity"))){
    fail("requires one fam() random factor and diagonal or fam() residuals.")
  }
  if(any(vapply(c(g$factors, e$factors), function(f) !is.null(f$section_levels), logical(1)))){
    fail("does not support dsumm() sections.")
  }
  if(ncol(y) != 1L || any(!is.finite(y)) || nrow(y) <= ncol(X) ||
     .sparse_rank(X) != ncol(X)){
    fail("requires one finite response and an estimable full-column-rank fixed design with residual degrees of freedom.")
  }
  responseSd <- stats::sd(as.numeric(y))
  if(!is.finite(responseSd) || responseSd <= 0) fail("requires positive response variance.")
  for(cs in structures){
    if(any(cs$free) && !all(cs$free)){
      fail("currently requires each covariance source to be entirely free or entirely fixed.")
    }
  }
  p <- g$dim; n <- nrow(X); r <- nrow(Ai[[1L]])
  if(length(Z) != p || any(vapply(Z, ncol, integer(1)) != r)){
    fail("requires equal-width trait incidence blocks.")
  }
  incidence <- Reduce(`+`, Z)
  entries <- summary(incidence)
  if(nrow(entries) != n || any(entries$x != 1) ||
     any(Matrix::rowSums(incidence != 0) != 1)){
    fail("requires one individual and one trait per observation (not random regression or overlapping designs).")
  }
  id <- integer(n); id[entries$i] <- entries$j
  trait <- integer(n)
  for(j in seq_len(p)){
    rows <- which(Matrix::rowSums(abs(Z[[j]])) > 0)
    if(any(trait[rows] != 0L)) fail("requires disjoint trait incidence blocks.")
    trait[rows] <- j
  }
  if(any(trait == 0L)) fail("requires a trait coordinate for every retained observation.")
  if(any(tabulate(trait, p) == 0L)) fail("requires observations for every modeled trait.")
  relationshipFactor <- Matrix::Cholesky(Matrix::forceSymmetric(Ai[[1L]]), LDL=FALSE)
  residualFA <- length(e$factors) && e$factors[[1L]]$model == "fa"
  if(e$dim != 1L && (!identical(g$levels, e$levels) || any(coordinate != trait))){
    fail("requires identical trait levels/order in the random and residual factors.")
  }
  block <- match(block, unique(block))
  if(residualFA && anyDuplicated(paste(block, trait, sep=":"))){
    fail("requires at most one observation of each trait per residual pairing key.")
  }
  shared <- g$factors[[1L]]$shared
  otherShared <- if(residualFA) e$factors[[1L]]$shared else NULL
  if(!identical(shared, otherShared)){
    fail("requires matching fam(shared=...) groups on both the random and residual term.")
  }
  if(!is.null(shared) && (!all(g$free) || !all(e$free) ||
      g$factors[[1L]]$order != e$factors[[1L]]$order)){
    fail("requires free genetic and residual FA terms of equal rank for shared loadings.")
  }
  gc <- .conditional_fa_component(g)
  ec <- if(length(e$factors)) .conditional_fa_component(e) else
    list(loadings=matrix(0, p, 0L), specific=exp(e$par[1L]))
  if(!residualFA) ec$loadings <- matrix(0, p, 0L)
  ec$specific <- rep(ec$specific, length.out=p)
  state <- list(genetic=gc$loadings, residual=ec$loadings,
    specificG=gc$specific, specificE=ec$specific)
  if(!is.null(shared)){
    ratios <- colSums(gc$loadings * ec$loadings) / colSums(gc$loadings^2)
    if(any(!is.finite(ratios)) || any(ratios <= 0) ||
       max(abs(ec$loadings - sweep(gc$loadings, 2L, ratios, `*`))) >
         1e-7 * max(1, abs(ec$loadings))){
      fail("requires shared starting loading columns to be positive scalar multiples across sources.")
    }
    state$loading <- sweep(gc$loadings, 2L, sqrt(1 + ratios^2), `*`)
    state$weightG <- 1 / (1 + ratios^2)
    state$weightE <- 1 - state$weightG
  }
  kg <- ncol(state$genetic); ke <- ncol(state$residual)
  m <- max(block); b <- ncol(X)
  specific <- b + seq_len(p * r)
  genetic <- b + p * r + seq_len(kg * r)
  residual <- b + (p + kg) * r + seq_len(ke * m)
  d <- b + (p + kg) * r + ke * m
  selector <- function(columns){
    Matrix::sparseMatrix(i=seq_len(n), j=columns, x=1, dims=c(n, d))
  }
  selectorsG <- lapply(seq_len(kg), function(k) selector(genetic[(k-1L)*r + id]))
  selectorsE <- lapply(seq_len(ke), function(k) selector(residual[(k-1L)*m + block]))
  selectors <- if(is.null(shared)) c(selectorsG, selectorsE) else
    Map(`+`, selectorsG, selectorsE)
  base <- cbind(X, do.call(cbind, Z),
    Matrix::Matrix(0, n, d - b - p*r, sparse=TRUE))
  priorPattern <- Matrix::bdiag(c(list(Matrix::Matrix(0, b, b, sparse=TRUE)),
    rep(list(Ai[[1L]]), p + kg), rep(list(Matrix::Diagonal(m)), ke)))
  structural <- base + Reduce(`+`, c(selectorsG, selectorsE))
  structural@x[] <- 1
  pattern <- as(Matrix::forceSymmetric(Matrix::crossprod(structural) + abs(priorPattern)),
    "generalMatrix")
  pattern@x[] <- 1
  # Only entries occurring in complete-data quadratics are needed from H^-1.
  # These fixed sparse contractions are reused for every EM iteration.
  locations <- pattern
  locations@x <- as.numeric(seq_along(locations@x))
  contract <- function(M){
    entries <- summary(as(M, "generalMatrix"))
    positions <- as.integer(locations[cbind(entries$i, entries$j)])
    if(any(positions == 0L)){
      fail("encountered a conditional quadratic outside the selected moment pattern.")
    }
    list(index=positions, weight=entries$x)
  }
  gramCache <- new.env(parent=emptyenv())
  groups <- lapply(seq_len(p), function(j){
    rows <- which(trait == j)
    B <- base[rows, , drop=FALSE]
    T <- lapply(selectors, function(S) S[rows, , drop=FALSE])
    key <- paste(id[rows], if(ke) block[rows] else 0L, collapse=":")
    gram <- gramCache[[key]]
    if(is.null(gram)){
      gram <- lapply(T, function(S) lapply(T, function(U)
        contract(Matrix::crossprod(S, U))))
      gramCache[[key]] <- gram
    }
    list(rows=rows, baseGram=contract(Matrix::crossprod(B)),
      cross=lapply(T, function(S) contract(Matrix::crossprod(S, B))), gram=gram)
  })
  list(X=X, Z=Z, Ai=Ai[[1L]], y=as.numeric(y), id=id, trait=trait, block=block,
    p=p, r=r, m=m, b=b, d=d, kg=kg, ke=ke, groups=groups, base=base,
    selectorsG=selectorsG, selectorsE=selectorsE, pattern=pattern,
    specific=matrix(specific, r, p), genetic=matrix(genetic, r, kg),
    residual=matrix(residual, m, ke), state=state, shared=shared,
    free=vapply(structures, function(cs) all(cs$free), logical(1)),
    homogeneousG=identical(g$factors[[1L]]$variance, "homogeneous"),
    homogeneousE=e$dim == 1L ||
      (length(e$factors) && e$factors[[1L]]$model == "identity") ||
      (residualFA && identical(e$factors[[1L]]$variance, "homogeneous")),
    gramGroups=length(ls(gramCache)), relationshipFactor=relationshipFactor,
    likelihoodOffset=(n-b)*log(responseSd), responseVariance=responseSd^2,
    logDetAi=as.numeric(determinant(Ai[[1L]], logarithm=TRUE)$modulus))
}

.conditional_fa_system <- function(model, state){
  loadingG <- if(is.null(model$shared)) state$genetic else state$loading
  loadingE <- if(is.null(model$shared)) state$residual else state$loading
  design <- model$base
  for(k in seq_len(model$kg)){
    design <- design + Matrix::Diagonal(x=loadingG[model$trait, k]) %*%
      model$selectorsG[[k]]
  }
  for(k in seq_len(model$ke)){
    design <- design + Matrix::Diagonal(x=loadingE[model$trait, k]) %*%
      model$selectorsE[[k]]
  }
  dg <- if(is.null(model$shared)) rep(1, model$kg) else state$weightG
  de <- if(is.null(model$shared)) rep(1, model$ke) else state$weightE
  prior <- Matrix::bdiag(c(list(Matrix::Matrix(0, model$b, model$b, sparse=TRUE)),
    lapply(c(state$specificG, dg), function(v) model$Ai / v),
    lapply(de, function(v) Matrix::Diagonal(model$m, 1/v))))
  weight <- 1/state$specificE[model$trait]
  weighted <- Matrix::Diagonal(x=weight) %*% design
  H <- as(Matrix::forceSymmetric(Matrix::crossprod(design, weighted) + prior),
    "generalMatrix")
  list(design=design, prior=prior, weight=weight, H=H,
    rhs=as.numeric(Matrix::crossprod(weighted, model$y)),
    logDetQ=(model$p + model$kg) * model$logDetAi -
      model$r * sum(log(c(state$specificG, dg))) - model$m * sum(log(de)))
}

.conditional_fa_posterior <- function(model, state, workspace, engine, moments){
  system <- .conditional_fa_system(model, state)
  posterior <- mme_conditional_moments_cpp(system$H, model$pattern, system$rhs,
    workspace, moments, engine$experimental$panelCi2, engine$denseMemoryMB)
  residual <- model$y - as.numeric(system$design %*% posterior$mean)
  quadratic <- sum(system$weight * residual^2) +
    sum(posterior$mean * as.numeric(system$prior %*% posterior$mean))
  posterior$llik <- -0.5 * (sum(log(state$specificE[model$trait])) -
    system$logDetQ + posterior$logDet + quadratic)
  posterior
}

.conditional_fa_update <- function(model, state, posterior){
  mean <- as.numeric(posterior$mean); covariance <- posterior$covariance
  contraction <- function(M) sum(covariance@x[M$index] * M$weight)
  updated <- state
  loading <- if(is.null(model$shared)) cbind(state$genetic, state$residual) else state$loading
  width <- ncol(loading)
  selectors <- if(is.null(model$shared)) c(model$selectorsG, model$selectorsE) else
    Map(`+`, model$selectorsG, model$selectorsE)
  predictionsAll <- vapply(selectors, function(S) as.numeric(S %*% mean),
    numeric(length(model$y)))
  offsetAll <- model$y - as.numeric(model$base %*% mean)
  sse <- numeric(model$p)
  for(j in seq_len(model$p)){
    group <- model$groups[[j]]
    offset <- offsetAll[group$rows]
    predictions <- predictionsAll[group$rows,,drop=FALSE]
    A <- matrix(0, width, width)
    for(a in seq_len(width)) for(b in seq_len(width)){
      A[a,b] <- sum(predictions[,a] * predictions[,b]) +
        contraction(group$gram[[a]][[b]])
    }
    cross <- as.numeric(crossprod(predictions, offset)) -
      vapply(group$cross, contraction, numeric(1))
    free <- if(is.null(model$shared)) c(rep(model$free[1L], model$kg),
      rep(model$free[2L], model$ke)) else rep(TRUE, width)
    if(any(free)){
      loading[j,free] <- solve(A[free,free,drop=FALSE],
        cross[free] - A[free,!free,drop=FALSE] %*% loading[j,!free])
    }
    sse[j] <- sum(offset^2) + contraction(group$baseGram) -
      2 * sum(loading[j,] * cross) + sum(loading[j,] * (A %*% loading[j,]))
  }
  floor <- .Machine$double.eps * max(1, state$specificG, state$specificE)
  expectedQuadratic <- function(indices, precision){
    values <- mean[indices]
    sum(values * as.numeric(precision %*% values)) +
      sum(covariance[indices,indices,drop=FALSE] * precision)
  }
  if(model$free[1L]){
    values <- vapply(seq_len(model$p), function(j)
      expectedQuadratic(model$specific[,j], model$Ai) / model$r, numeric(1))
    if(model$homogeneousG) values[] <- mean(values)
    updated$specificG <- pmax(values, floor)
  }
  if(model$free[2L]){
    counts <- vapply(model$groups, function(group) length(group$rows), integer(1))
    values <- if(model$homogeneousE) rep(sum(sse)/sum(counts), model$p) else sse/counts
    updated$specificE <- pmax(values, floor)
  }
  if(is.null(model$shared)){
    updated$genetic <- loading[,seq_len(model$kg),drop=FALSE]
    updated$residual <- loading[,model$kg + seq_len(model$ke),drop=FALSE]
  }else{
    dg <- vapply(seq_len(model$kg), function(k)
      expectedQuadratic(model$genetic[,k], model$Ai) / model$r, numeric(1))
    de <- vapply(seq_len(model$ke), function(k)
      expectedQuadratic(model$residual[,k], Matrix::Diagonal(model$m)) / model$m,
      numeric(1))
    dg <- pmax(dg, floor); de <- pmax(de, floor)
    updated$loading <- sweep(loading, 2L, sqrt(dg + de), `*`)
    updated$weightG <- dg/(dg + de)
    updated$weightE <- de/(dg + de)
  }
  updated
}

.conditional_fa_descriptor <- function(cs, loading, specific){
  factor <- if(length(cs$factors)) cs$factors[[1L]] else NULL
  if(!is.null(factor) && factor$model == "fa"){
    k <- ncol(loading)
    rotation <- qr.Q(qr(t(loading[seq_len(k),,drop=FALSE])))
    loading <- loading %*% rotation
    signs <- ifelse(diag(loading[seq_len(k),,drop=FALSE]) < 0, -1, 1)
    loading <- sweep(loading, 2L, signs, `*`)
    if(any(diag(loading[seq_len(k),,drop=FALSE]) <= 0)){
      stop("Conditional FA leading loadings became rank deficient; use a smaller rank or different starting values.", call.=FALSE)
    }
    new <- fam(base::factor(factor$levels, levels=factor$levels), k=k,
      loadings=loading, specific=specific, variance=factor$variance,
      shared=factor$shared)$covFactor
    shape <- new$par
    scale <- sum(loading[1L,]^2) + specific[1L]
    cs$par <- c(log(scale), shape)
    new$par_start <- factor$par_start; new$par_end <- factor$par_end
    new$free <- factor$free
    cs$factors[[1L]] <- new
  }else{
    cs$par <- c(log(specific[1L]), if(cs$dim > 1L &&
      !identical(factor$model, "identity")) log(specific[-1L]/specific[1L]))
    if(length(cs$factors)) cs$factors[[1L]]$par <- cs$par[-1L]
  }
  cs
}

.mmes_conditional_fa <- function(model, structures, engine, iterations,
                                 likelihoodTol, parameterTol, verbose){
  state <- model$state
  optimizer <- if(identical(engine$updateMethod,"ai")) "ai" else "em"
  if(optimizer == "ai" && !.conditional_fa_ai_bounds(model,
      .conditional_fa_coordinates(model,state),engine)$feasible){
    stop("Conditional FA AI starting variances/allocations fall below its numerical bounds; revise starting values or execution$inverseTol.",call.=FALSE)
  }
  posterior <- .conditional_fa_posterior(model, state, NULL, engine, TRUE)
  history <- posterior$llik
  converged <- !any(model$free)
  updates <- character(); aiDiagnostics <- list()
  for(iteration in seq_len(if(converged) 0L else iterations)){
    update <- "em"
    if(optimizer == "ai"){
      proposal <- .conditional_fa_ai_step(model,state,posterior,engine,
        parameterTol)
      aiDiagnostics <- proposal[c("scoreNorm","rank","evaluations","workspaceBytes",
        "activeBounds","varianceFloor","allocationFloor")]
      if(proposal$stationary){
        converged <- TRUE
        break
      }
      if(proposal$accepted){
        updated <- proposal$state
        update <- "ai"
      }else{
        updated <- .conditional_fa_update(model,state,posterior)
        coordinates <- .conditional_fa_coordinates(model,updated)
        bounds <- .conditional_fa_ai_bounds(model,coordinates,engine)
        updated <- .conditional_fa_decode(model,updated,coordinates,
          pmin(bounds$upper,pmax(bounds$lower,coordinates$values)))
        update <- "em-recovery"
      }
    }else{
      updated <- .conditional_fa_update(model, state, posterior)
    }
    candidate <- .conditional_fa_posterior(model, updated, posterior$workspace, engine, TRUE)
    improvement <- candidate$llik - posterior$llik
    if(!is.finite(improvement) || improvement < -1e-7 * (1 + abs(posterior$llik))){
      stop("Conditional FA update decreased the observed restricted likelihood; no candidate was accepted.", call.=FALSE)
    }
    change <- sqrt(sum((unlist(updated) - unlist(state))^2)) /
      (1 + sqrt(sum(unlist(state)^2)))
    state <- updated; posterior <- candidate
    updates <- c(updates,update)
    history <- c(history, posterior$llik)
    if(verbose) cat("Conditional FA", update, iteration, "logLik", posterior$llik,
      "relative change", change, "\n")
    if(optimizer == "em" && abs(improvement) <= likelihoodTol && change <= parameterTol){
      converged <- TRUE
      break
    }
  }
  if(!converged) warning("Conditional FA ", toupper(optimizer), " did not converge within convergence$maxIters; inspect conditionalFA$likelihood and increase the limit or revise the rank.", call.=FALSE)
  mean <- as.numeric(posterior$mean)
  loadingG <- if(is.null(model$shared)) state$genetic else
    sweep(state$loading, 2L, sqrt(state$weightG), `*`)
  loadingE <- if(is.null(model$shared)) state$residual else
    sweep(state$loading, 2L, sqrt(state$weightE), `*`)
  scoresG <- matrix(mean[model$genetic], model$r, model$kg)
  genetic <- matrix(mean[model$specific], model$r, model$p) +
    tcrossprod(scoresG, if(is.null(model$shared)) loadingG else state$loading)
  for(j in seq_along(structures)){
    structures[[j]] <- .conditional_fa_descriptor(structures[[j]],
      if(j == 1L) loadingG else loadingE,
      if(j == 1L) state$specificG else state$specificE)
  }
  natural <- lapply(structures, function(cs) vapply(seq_along(cs$par), function(j){
    kind <- .vc_kind(cs, j)
    .vc_to_natural(cs$par[j], kind$kind, kind$lower, kind$upper)
  }, numeric(1)))
  bu <- matrix(c(mean[seq_len(model$b)], as.vector(genetic)), ncol=1L)
  count <- if(is.null(model$shared)) sum(vapply(structures, function(cs) sum(cs$free), integer(1))) else
    model$p*model$kg + model$kg +
      (if(model$homogeneousG) 1L else model$p) + (if(model$homogeneousE) 1L else model$p)
  reportedHistory <- history + model$likelihoodOffset
  list(b=bu[seq_len(model$b),,drop=FALSE], u=matrix(as.vector(genetic), ncol=1L),
    bu=bu, uList=list(genetic), W=do.call(cbind, c(list(model$X), model$Z)),
    partitions=list(cbind(model$b + (seq_len(model$p)-1L)*model$r + 1L,
      model$b + seq_len(model$p)*model$r)), theta=vector("list", 2L),
    covStruct=structures, covPar=lapply(natural, function(x) matrix(x, ncol=1L)),
    covParWorking=lapply(structures, `[[`, "par"), monitor=matrix(unlist(natural), ncol=1L),
    llik=matrix(reportedHistory, nrow=1L), convergence=converged, CiMode=0L, C=NULL, Ci=NULL,
    theta_se=NULL, InfMat=NULL, covarianceInformation="none",
    retainCoefficientMatrix=FALSE, AIC=-2*tail(reportedHistory,1L)+2*model$b,
    BIC=-2*tail(reportedHistory,1L)+log(length(model$y))*model$b,
    solver="ldlt", conditionalFA=list(shared=model$shared, state=state,
      likelihood=reportedHistory, rawLikelihood=history,
      likelihoodOffset=model$likelihoodOffset,
      factorScoresG=scoresG, parameterCount=count, optimizer=optimizer,
      updates=updates, ai=aiDiagnostics),
    engineDiagnostics=list(conditionalFA=list(active=TRUE, fallbackReason="",
      symbolicAnalyses=posterior$symbolicAnalyses,
      numericFactorizations=posterior$numericFactorizations,
      selectedEntries=posterior$selectedEntries, panelColumns=posterior$panelColumns,
      panelBytes=posterior$panelBytes, equations=model$d,
      gramGroups=model$gramGroups,
      optimizer=optimizer, updates=updates, ai=aiDiagnostics,
      method=paste("exact sparse conditional-moment REML",toupper(optimizer)))))
}
