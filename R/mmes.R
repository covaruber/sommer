# Refactored mmes front end: unified formula environments, centralized
# observation filtering, and language-object parsing for random/residual terms.
mmes <- function(fixed, random, rcov, data, W,
                 nIters=30, tolParConvLL=1e-04,
                 tolParConvNorm=1e-04, tolParInv=1e-06,
                 naMethodX="exclude", naMethodY="exclude",
                 naMethodRandom="exclude", naMethodR="exclude",
                 returnParam=FALSE, dateWarning=TRUE,
                 verbose=TRUE, stepWeight=NULL, emWeight=NULL,
                 contrasts=NULL, getPEV=TRUE, henderson=TRUE,
                 computeCi=0, solver="ldlt", pcgTol=1.0e-8,
                 pcgMaxIters=0, pcgTraceProbes=8,
                 pcgLanczosSteps=20){
  
  if(!isTRUE(henderson)){
    stop("This mmes() interface is Henderson-only. Use the separate MNR/direct-inversion mmer interface for henderson=FALSE.",
         call.=FALSE)
  }
  
  desc <- utils::packageDescription("sommer")
  my.date <- as.Date(desc$Date) + 90
  if(dateWarning && Sys.Date() > my.date){
    cat("Version out of date. Please update sommer to the newest version using:\n",
        "install.packages('sommer') in a new session\n",
        "Use the 'dateWarning' argument to disable the warning message.\n", sep="")
  }
  
  # ---- Helpers ---------------------------------------------------------
  formula_env <- function(f, fallback){
    e <- if(inherits(f, "formula")) environment(f) else NULL
    if(is.null(e)) fallback else e
  }
  
  # Split only top-level additions. '+' inside I(), vsm(), etc. is untouched.
  split_plus <- function(expr){
    if(is.call(expr) && identical(expr[[1L]], as.name("+"))){
      c(split_plus(expr[[2L]]), split_plus(expr[[3L]]))
    } else list(expr)
  }
  
  expr_label <- function(x) paste(deparse(x, width.cutoff=500L), collapse="")

  call_name <- function(expr){
    if(!is.call(expr)) return(NULL)
    head <- expr[[1L]]
    if(is.symbol(head)) return(as.character(head))
    if(is.call(head) && length(head) == 3L &&
       as.character(head[[1L]]) %in% c("::", ":::")){
      return(as.character(head[[3L]]))
    }
    NULL
  }
  
  has_call <- function(expr, names){
    if(!is.call(expr)) return(FALSE)
    head <- call_name(expr)
    if(length(head) == 1L && head %in% names) return(TRUE)
    any(vapply(as.list(expr)[-1L], has_call, logical(1), names=names))
  }
  
  eval_model_expr <- function(expr, data_full, enclos){
    if(is.null(data_full)) eval(expr, envir=enclos)
    else eval(expr, envir=data_full, enclos=enclos)
  }
  
  # Observation-level variables are symbols resolving to vectors of length n,
  # or matrices/data.frames with n rows. Objects such as Gu/Ai are therefore
  # not mistaken for observation variables unless they actually have n rows.
  observation_ok <- function(expr, data_full, enclos, n){
    vars <- unique(all.vars(expr))
    if(!length(vars)) return(rep(TRUE, n))
    ok <- rep(TRUE, n)
    for(v in vars){
      val <- tryCatch({
        if(!is.null(data_full) && v %in% names(data_full)) data_full[[v]]
        else get(v, envir=enclos, inherits=TRUE)
      }, error=function(e) NULL)
      if(is.null(val)) next
      if(is.data.frame(val) || is.matrix(val)){
        if(nrow(val) == n) ok <- ok & stats::complete.cases(val)
      } else if(length(val) == n){
        ok <- ok & !is.na(val)
      }
    }
    ok
  }
  
  method_keep <- function(ok, method, what){
    method <- tolower(method)
    if(method %in% c("exclude", "omit")) return(ok)
    if(method %in% c("include", "pass")) return(rep(TRUE, length(ok)))
    if(method == "fail" && any(!ok))
      stop("Missing values found in ", what, ".", call.=FALSE)
    if(method == "fail") return(rep(TRUE, length(ok)))
    stop("Unknown missing-data method '", method, "' for ", what, ".", call.=FALSE)
  }
  
  # ---- Unified evaluation context -------------------------------------
  callEnv <- parent.frame()
  fixedEnv <- formula_env(fixed, callEnv)
  if(missing(rcov)){
    rcov <- stats::as.formula("~units", env=fixedEnv)
  } else if(is.null(environment(rcov))){
    environment(rcov) <- fixedEnv
  }
  if(!missing(random) && is.null(environment(random))) environment(random) <- fixedEnv
  
  dataSupplied <- !missing(data)
  data_full <- if(dataSupplied) as.data.frame(data) else NULL
  
  # Build the fixed model frame before filtering. model.frame follows normal
  # R lookup rules: data first, then the formula environment.
  mf_full <- try(stats::model.frame(fixed, data=data_full,
                                    na.action=stats::na.pass,
                                    drop.unused.levels=FALSE), silent=TRUE)
  if(inherits(mf_full, "try-error")){
    stop("Unable to evaluate the fixed formula. Variables may be supplied in 'data' or in the formula/calling environment.\n",
         as.character(mf_full), call.=FALSE)
  }
  nObs <- nrow(mf_full)
  if(nObs < 1L) stop("No observations are available for model fitting.", call.=FALSE)
  if(!is.null(data_full) && nrow(data_full) != nObs){
    stop("The fixed formula and 'data' do not describe the same number of observations.", call.=FALSE)
  }
  
  # If data was omitted, create a row scaffold. Formula variables remain
  # available through the formula environment and are not copied unnecessarily.
  if(is.null(data_full)) data_full <- data.frame(.sommer_row=seq_len(nObs))
  data_full$.sommer_row <- seq_len(nObs)
  data_full$units <- factor(paste0("u", seq_len(nObs)),
                            levels=paste0("u", seq_len(nObs)))
  
  # ---- Parse/evaluate random and residual expressions on full rows -----
  randomExprs <- if(missing(random)) list() else split_plus(random[[2L]])
  randomLabels <- vapply(randomExprs, expr_label, character(1))
  randomFits <- vector("list", length(randomExprs))
  
  if(length(randomExprs)){
    randomEnv <- formula_env(random, fixedEnv)
    for(u in seq_along(randomExprs)){
      ex <- randomExprs[[u]]
      if(!has_call(ex, c("vsm", "covm", "spl2Dc"))){
        ex <- as.call(list(as.name("vsm"), as.call(list(as.name("ism"), ex))))
      }
      randomExprs[[u]] <- ex
      randomLabels[u] <- expr_label(ex)
      randomFits[[u]] <- eval_model_expr(ex, data_full, randomEnv)
      ff <- randomFits[[u]]
      if(is.null(ff$covStruct) || !identical(ff$covStruct$type, "kron") ||
         is.null(ff$covStruct$descriptor_version) || ff$covStruct$descriptor_version < 2L){
        stop("All random covariance terms must use the CovarianceFactor v2 vsm() descriptor interface.",
             call.=FALSE)
      }
    }
  }
  
  rcovExprs <- split_plus(rcov[[2L]])
  if(length(rcovExprs) != 1L){
    stop("The Henderson interface accepts one residual vsm() term. Use arbitrary Kronecker products inside that vsm() term instead of summing residual terms.",
         call.=FALSE)
  }
  residualExpr <- rcovExprs[[1L]]
  if(!has_call(residualExpr, c("vsm", "gvs", "spl2Da", "spl2Db"))){
    residualExpr <- as.call(list(as.name("vsm"), as.call(list(as.name("ism"), residualExpr))))
  }
  residualLabel <- expr_label(residualExpr)
  residualEnv <- formula_env(rcov, fixedEnv)
  rf_full <- eval_model_expr(residualExpr, data_full, residualEnv)
  if(is.null(rf_full$covStruct) || !identical(rf_full$covStruct$type, "kron") ||
     is.null(rf_full$covStruct$descriptor_version) || rf_full$covStruct$descriptor_version < 2L){
    stop("The residual covariance term must use the CovarianceFactor v2 vsm() descriptor interface.",
         call.=FALSE)
  }
  if(is.null(rf_full$residualLocalIndex)) rf_full$residualLocalIndex <- rep(1L, nObs)
  if(length(rf_full$residualLocalIndex) != nObs){
    stop("Residual local-index vector has incompatible length.", call.=FALSE)
  }
  
  # ---- Centralized observation map ------------------------------------
  # Response and fixed RHS are separated so naMethodY and naMethodX retain
  # their historical meaning.
  responseMF <- mf_full[, 1L, drop=FALSE]
  responseOK <- stats::complete.cases(responseMF)
  fixedOK <- if(ncol(mf_full) > 1L) stats::complete.cases(mf_full[, -1L, drop=FALSE]) else rep(TRUE, nObs)
  
  randomOK <- rep(TRUE, nObs)
  if(length(randomExprs)){
    for(ex in randomExprs) randomOK <- randomOK & observation_ok(ex, data_full, formula_env(random, fixedEnv), nObs)
  }
  residualOK <- observation_ok(residualExpr, data_full, residualEnv, nObs) & !is.na(rf_full$residualLocalIndex)
  
  keepY <- method_keep(responseOK, naMethodY, "the response")
  keepX <- method_keep(fixedOK, naMethodX, "fixed-effect variables")
  keepRandom <- method_keep(randomOK, naMethodRandom, "random-effect variables")
  keepResidual <- method_keep(residualOK, naMethodR, "residual covariance variables")
  keep <- keepY & keepX & keepRandom & keepResidual
  
  reason <- rep("included", nObs)
  reason[!keepY] <- "missing response"
  reason[keepY & !keepX] <- "missing fixed covariate"
  reason[keepY & keepX & !keepRandom] <- "missing random-effect variable"
  reason[keepY & keepX & keepRandom & !keepResidual] <- "missing residual coordinate"
  
  obsInfo <- data.frame(originalRow=seq_len(nObs), responseOK=responseOK,
                        fixedOK=fixedOK, randomOK=randomOK,
                        residualOK=residualOK, included=keep,
                        reason=reason, stringsAsFactors=FALSE)
  if(!any(keep)) stop("No observations remain after applying the missing-data rules.", call.=FALSE)
  
  # The model frame is the authoritative fixed/response representation.
  mf <- mf_full[keep, , drop=FALSE]
  data <- data_full[keep, , drop=FALSE]
  dataor <- data_full
  
  # Response matrix: preserve multivariate/model.frame response behavior.
  yobj <- stats::model.response(mf)
  if(is.null(dim(yobj))) yobj <- matrix(yobj, ncol=1L)
  yvar <- Matrix::Matrix(yobj, sparse=TRUE)
  responseNames <- colnames(yobj)
  if(is.null(responseNames)) responseNames <- as.character(fixed[[2L]])
  if(ncol(yvar) == 1L) colnames(yvar) <- responseNames[1L]
  
  # ---- Random structures, now subset exactly once ---------------------
  Z <- list(); Ai <- list(); covStruct <- list(); Zind <- numeric()
  rTermsNames <- list(); rtermss <- randomLabels
  
  if(length(randomFits)){
    for(u in seq_along(randomFits)){
      ff <- randomFits[[u]]
      Zi <- lapply(ff$Z, function(x){
        if(nrow(x) == nObs) x[keep, , drop=FALSE]
        else if(nrow(x) == sum(keep)) x
        else stop("Random-effect design has incompatible number of rows in term: ", rtermss[u], call.=FALSE)
      })
      Z <- c(Z, Zi)
      pGu <- to_sparse(ff$Gu)
      attr(pGu, "inverse") <- TRUE
      Ai[[u]] <- pGu
      covStruct[[u]] <- ff$covStruct
      Zind <- c(Zind, rep(u, length(Zi)))
      s2 <- paste(all.vars(randomExprs[[u]]), collapse=":")
      rTermsNames[[u]] <- paste(s2, ff$covStruct$par_names, sep=":")
    }
  }
  nRandomStruct <- length(covStruct)
  
  # ---- Residual structure ---------------------------------------------
  rf <- rf_full
  if(isTRUE(rf$covStruct$free[1])) rf$covStruct$par[1] <- rf$covStruct$par[1] + log(5)
  residualStructIndex <- nRandomStruct + 1L
  covStruct[[residualStructIndex]] <- rf$covStruct
  localIndex <- as.integer(rf$residualLocalIndex[keep])
  
  # Preserve the established residual-block pairing semantics for this
  # refactor. The important change here is that it is applied after the
  # single centralized observation mask, so R, X, Z and y see identical rows.
  # .sommer_row is internal and must never participate in pairing.
  residualVars <- unique(all.vars(residualExpr))
  residualVars <- setdiff(residualVars, "units")
  pairingVariables <- setdiff(names(data),
                              c(responseNames, residualVars, "units", ".sommer_row"))
  if(length(pairingVariables)){
    baseKey <- do.call(paste, c(data[pairingVariables], sep="\r"))
  } else {
    baseKey <- rep("all", nrow(data))
  }
  
  pairLocal <- paste(baseKey, localIndex, sep="\r")
  occurrence <- ave(seq_along(pairLocal), pairLocal, FUN=seq_along)
  blockKey <- paste(baseKey, occurrence, sep="\r")
  residualBlock <- match(blockKey, unique(blockKey))
  if(anyDuplicated(paste(residualBlock, localIndex, sep=":"))){
    stop("Internal residual-layout error: a block contains duplicate local covariance coordinates.", call.=FALSE)
  }
  s2 <- paste(all.vars(residualExpr), collapse=":")
  rTermsNames[[residualStructIndex]] <- paste(s2, rf$covStruct$par_names, sep=":")
  
  # ---- Fixed-effect design using terms()/assign -----------------------
  X <- Matrix::sparse.model.matrix(fixed, data=mf, contrasts.arg=contrasts)
  tt <- attr(mf, "terms")
  fixedTerms <- attr(tt, "term.labels")
  assignX <- attr(X, "assign")
  if(is.null(assignX)) assignX <- attr(stats::model.matrix(tt, mf, contrasts.arg=contrasts), "assign")
  partitionsX <- list()
  if(attr(tt, "intercept") == 1L){
    ii <- which(assignX == 0L)
    if(length(ii)) partitionsX[["1"]] <- matrix(ii, nrow=1L)
  }
  for(ix in seq_along(fixedTerms)){
    ii <- which(assignX == ix)
    if(length(ii)) partitionsX[[fixedTerms[ix]]] <- matrix(ii, nrow=1L)
  }
  if("(Intercept)" %in% colnames(X)) colnames(X)[colnames(X) == "(Intercept)"] <- "Intercept"
  
  # ---- Weights ---------------------------------------------------------
  if(missing(W)){
    W <- Matrix::Diagonal(n=nrow(yvar), x=1)
    useH <- FALSE
  }else{
    # W may be supplied for all original rows or already-filtered rows.
    if(nrow(W) == nObs && ncol(W) == nObs) W <- W[keep, keep, drop=FALSE]
    else if(nrow(W) != sum(keep) || ncol(W) != sum(keep))
      stop("W must have dimensions equal to either the original or retained number of observations.", call.=FALSE)
    W <- as(as(as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    useH <- TRUE
  }
  
  if(is.null(emWeight)){
    # EM-heavy warm start: early REML iterations are intentionally more
    # conservative and rely on the EM information block to stabilize the
    # variance-component update; later iterations transition smoothly toward
    # AI-dominated updates as the estimate enters the asymptotic regime.
    if(nIters <= 1L){
      emWeight <- 1
    } else {
      emWeight <- exp(seq(log(1), log(0.05), length.out = nIters))
      emWeight[1L] <- 1
      emWeight[length(emWeight)] <- 0.05
    }
  }
  if(length(emWeight) == 1L) emWeight <- rep(emWeight, nIters)
  if(length(emWeight) != nIters) emWeight <- rep(emWeight, length.out = nIters)
  if(any(!is.finite(emWeight)) || any(emWeight < 0 | emWeight > 1))
    stop("emWeight must contain finite values between 0 and 1.", call.=FALSE)
  
  if(is.null(stepWeight)){
    w <- which(emWeight <= .5)
    stepWeight <- rep(.9, nIters)
    if(nIters > 1){
      if(length(w) > 1) stepWeight[w[1:2]] <- c(.5,.7)
      else stepWeight[seq_len(min(2L,nIters))] <- c(.5,.7)[seq_len(min(2L,nIters))]
    }
  }
  if(length(stepWeight) == 1L) stepWeight <- rep(stepWeight, nIters)
  if(length(stepWeight) != nIters) stepWeight <- rep(stepWeight, length.out = nIters)
  if(any(!is.finite(stepWeight)) || any(stepWeight <= 0))
    stop("stepWeight must contain finite positive values.", call.=FALSE)
  
  if(length(Ai)){
    nInverses <- sum(vapply(Ai, function(x) isTRUE(attr(x,"inverse")), logical(1)))
    if(nInverses != length(Ai)){
      stop("The Henderson algorithm requires every Gu relationship matrix to be supplied as an inverse matrix with attr(Gu,'inverse')=TRUE.", call.=FALSE)
    }
  }
  
  if(returnParam){
    return(list(yvar=yvar, X=X, Z=Z, Zind=Zind, Ai=Ai,
                W=W, useH=useH, residualBlock=residualBlock,
                residualIndex=localIndex, nIters=nIters,
                tolParConvLL=tolParConvLL, tolParConvNorm=tolParConvNorm,
                tolParInv=tolParInv, verbose=verbose, covStruct=covStruct,
                stepWeight=stepWeight, emWeight=emWeight,
                rtermss=rtermss, partitionsX=partitionsX,
                getPEV=getPEV, rTermsNames=rTermsNames,
                obsInfo=obsInfo))
  }
  
  res <- .Call("_sommer_ai_mme_sp2", PACKAGE="sommer",
               X, Z, Zind, Ai, yvar, W, useH,
               residualBlock, localIndex,
               nIters, tolParConvLL, tolParConvNorm,
               tolParInv, covStruct, emWeight, stepWeight,
               verbose, computeCi, solver, pcgTol, pcgMaxIters,
               pcgTraceProbes, pcgLanczosSteps)
  
  rownames(res$b) <- colnames(X)
  if(length(randomFits) && length(res$u)) rownames(res$u) <- unlist(lapply(Z, colnames))
  rownames(res$bu) <- c(rownames(res$b), rownames(res$u))
  rownames(res$monitor) <- unlist(rTermsNames)
  if(!is.null(res$monitorOriginalScale)) rownames(res$monitorOriginalScale) <- unlist(rTermsNames)
  
  res$data <- data
  res$dataOriginal <- dataor
  res$obsInfo <- obsInfo
  res$y <- yvar
  res$partitionsX <- partitionsX
  res$covStruct <- covStruct
  
  if(length(randomFits) && length(rtermss)){
    names(res$theta) <- c(rtermss, residualLabel)
    names(res$covPar) <- c(rtermss, residualLabel)
    names(res$partitions) <- rtermss
    names(res$uList) <- rtermss
    if(getPEV && computeCi > 0) names(res$uPevList) <- rtermss
    for(i in seq_along(res$partitions)){
      colnames(res$uList[[i]]) <- covStruct[[i]]$levels
      rr <- res$partitions[[i]][1,1]:res$partitions[[i]][1,2]
      rownames(res$uList[[i]]) <- rownames(res$bu)[rr]
      if(getPEV && computeCi > 0){
        colnames(res$uPevList[[i]]) <- covStruct[[i]]$levels
        rownames(res$uPevList[[i]]) <- rownames(res$bu)[rr]
      }
    }
    res$args <- list(fixed=fixed, random=random, rcov=rcov)
    res$Dtable <- data.frame(type=c(rep("fixed",length(res$partitionsX)),rep("random",length(res$partitions))),
                             term=c(names(res$partitionsX),names(res$partitions)),
                             include=FALSE, average=FALSE)
  }else{
    names(res$theta) <- residualLabel
    res$args <- list(fixed=fixed, rcov=rcov)
    res$Dtable <- data.frame(type=rep("fixed",length(res$partitionsX)),
                             term=names(res$partitionsX), include=FALSE, average=FALSE)
  }
  
  class(res) <- "mmes"
  res
}
