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
                 computeCi=0, solver="auto", pcgTol=1.0e-8,
                 pcgMaxIters=0, pcgTraceProbes=8,
                 pcgLanczosSteps=20, REML=TRUE,
                 family=stats::gaussian(), pqlControl=list(),
                 .pqlInner=FALSE, .pqlFixedDispersion=FALSE,
                 .pqlWorkingPrecision=NULL, .pqlBaseW=NULL,
                 .pqlBaseFactor=NULL){

  WWasMissing <- missing(W)

  if(length(henderson) != 1L || !is.logical(henderson) || is.na(henderson)){
    stop("henderson must be a single TRUE/FALSE value.", call.=FALSE)
  }

  if(!inherits(family, "family")){
    stop("family must be a family object, such as stats::binomial() or stats::poisson().",
         call.=FALSE)
  }
  isGaussianIdentity <- identical(family$family, "gaussian") &&
    identical(family$link, "identity")
  if(!.pqlInner && !isGaussianIdentity){
    return(get(".mmes_pql", mode="function")(
      fixed=fixed,
      random=if(missing(random)) NULL else random,
      rcov=if(missing(rcov)) NULL else rcov,
      data=if(missing(data)) NULL else data,
      W=if(missing(W)) NULL else W,
      family=family,
      pqlControl=pqlControl,
      mmesArgs=list(
        nIters=nIters, tolParConvLL=tolParConvLL,
        tolParConvNorm=tolParConvNorm, tolParInv=tolParInv,
        naMethodX=naMethodX, naMethodY=naMethodY,
        naMethodRandom=naMethodRandom, naMethodR=naMethodR,
        dateWarning=dateWarning, verbose=verbose, stepWeight=stepWeight,
        emWeight=emWeight, contrasts=contrasts, getPEV=getPEV,
        henderson=henderson, computeCi=computeCi, solver=solver,
        pcgTol=pcgTol, pcgMaxIters=pcgMaxIters,
        pcgTraceProbes=pcgTraceProbes, pcgLanczosSteps=pcgLanczosSteps,
        REML=REML
      )
    ))
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
  if(isTRUE(.pqlFixedDispersion)){
    rf$covStruct$par[1L] <- 0
    rf$covStruct$free[1L] <- FALSE
  }
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

  # ---- Optional Lee-van der Werf observation rotation -----------------
  rotationTerms <- which(vapply(
    randomFits,
    function(z) !is.null(z$rotation),
    logical(1)
  ))
  rotationInfo <- NULL
  responsePrepared <- FALSE
  preparedMean <- 0
  preparedSd <- 1
  preparedIntercept <- FALSE

  if(length(rotationTerms)){
    if(length(rotationTerms) != 1L){
      stop("Only one random-effect term may request rotation.", call.=FALSE)
    }
    if(!isGaussianIdentity || isTRUE(.pqlInner)){
      stop("rotation=TRUE currently requires a Gaussian identity-link model.",
           call.=FALSE)
    }
    if(!isTRUE(henderson) && ncol(yvar) != 1L){
      stop("The direct rotation path currently requires one response column.",
           call.=FALSE)
    }
    if(!WWasMissing || !is.null(.pqlWorkingPrecision) ||
       !is.null(.pqlBaseW) || !is.null(.pqlBaseFactor)){
      stop("rotation=TRUE currently requires the default identity W matrix.",
           call.=FALSE)
    }
    if(rf$covStruct$dim != 1L || any(localIndex != 1L)){
      stop("rotation=TRUE currently requires an identity residual structure (rcov=~units).",
           call.=FALSE)
    }
    if(computeCi != 0L){
      stop("rotation=TRUE currently requires computeCi=0; rotated PEV support is not yet available.",
           call.=FALSE)
    }

    rotationTerm <- rotationTerms[[1L]]
    focal <- randomFits[[rotationTerm]]$rotation
    U <- focal$vectors
    focalBlocks <- which(Zind == rotationTerm)
    focalZ <- Z[focalBlocks]
    nLevels <- nrow(U)

    if(!length(focalZ) || any(vapply(focalZ, ncol, integer(1)) != nLevels)){
      stop("The rotated random term has incompatible incidence dimensions.",
           call.=FALSE)
    }

    rowsByLevelList <- list()
    observationBlockTerm <- integer()
    for(j in seq_along(focalZ)){
      ss <- Matrix::summary(focalZ[[j]])
      counts <- tabulate(ss$j, nbins=nLevels)
      if(!nrow(ss) || any(abs(ss$x - 1) > 1e-12) ||
         any(counts == 0L) || length(unique(counts)) != 1L){
        stop(
          paste0(
            "rotation=TRUE requires each covariance coordinate to contain ",
            "every Gu level equally often with unit incidence."
          ),
          call.=FALSE
        )
      }
      rows <- split(ss$i, factor(ss$j, levels=seq_len(nLevels)))
      rowsByLevelList[[j]] <- do.call(rbind, rows)
      observationBlockTerm <- c(
        observationBlockTerm,
        rep(j, counts[[1L]])
      )
    }
    rowsByLevel <- do.call(cbind, rowsByLevelList)
    if(length(rowsByLevel) != nrow(yvar) ||
       !identical(sort(as.integer(rowsByLevel)), seq_len(nrow(yvar)))){
      stop(
        paste0(
          "rotation=TRUE requires a complete balanced layout: every retained ",
          "observation must belong to exactly one relationship-level block."
        ),
        call.=FALSE
      )
    }

    rotateRows <- function(M){
      out <- as.matrix(M)
      for(j in seq_len(ncol(rowsByLevel))){
        rr <- rowsByLevel[,j]
        out[rr,] <- crossprod(U, out[rr,,drop=FALSE])
      }
      to_sparse(Matrix::Matrix(out, sparse=TRUE))
    }

    yOriginal <- yvar
    XOriginal <- X
    ZOriginal <- Z

    yvar <- rotateRows(yvar)
    X <- rotateRows(X)
    for(j in seq_along(Z)){
      if(!j %in% focalBlocks) Z[[j]] <- rotateRows(Z[[j]])
    }
    for(j in seq_along(focalBlocks)){
      blockColumns <- which(observationBlockTerm == j)
      selectedRows <- as.integer(rowsByLevel[,blockColumns,drop=FALSE])
      selectedModes <- rep(seq_len(nLevels), length(blockColumns))
      Z[[focalBlocks[j]]] <- Matrix::sparseMatrix(
        i=selectedRows,
        j=selectedModes,
        x=1,
        dims=c(nrow(yvar), nLevels),
        dimnames=list(NULL, focal$modes)
      )
    }
    Ai[[rotationTerm]] <- randomFits[[rotationTerm]]$GuRot
    attr(Ai[[rotationTerm]], "inverse") <- TRUE

    rotationInfo <- list(
      term=rotationTerm,
      termName=rtermss[[rotationTerm]],
      vectors=U,
      precision=focal$precision,
      covariance=focal$covariance,
      levels=focal$levels,
      modes=focal$modes,
      rowsByLevel=rowsByLevel,
      formulation=if(isTRUE(henderson)) "henderson-eigen-coefficients" else "direct-observation-covariance",
      yOriginal=yOriginal,
      XOriginal=XOriginal,
      ZOriginal=ZOriginal
    )
    responsePrepared <- TRUE
    preparedMean <- mean(as.numeric(yOriginal))
    preparedSd <- stats::sd(as.numeric(yOriginal))
    preparedIntercept <- "Intercept" %in% colnames(XOriginal)
    if(!is.finite(preparedSd) || preparedSd <= 0){
      stop("rotation=TRUE requires a response with positive finite variance.",
           call.=FALSE)
    }
  }
  
  # ---- Data-driven starting values for variance-component scales ------
  # Replaces vsm()'s flat sigma2 default (0.15 random / 0.75 residual, both
  # scale-blind) with values informed by the data. Only touches par[1]
  # (log_sigma2) entries still flagged sigma2_is_default with free[1]=TRUE;
  # user-supplied or fixedSigma2=TRUE values are always left untouched. Any
  # failure here is silently ignored and the old flat defaults are kept,
  # since this only affects the optimization starting point, never the
  # converged answer.
  tryCatch({
    if(ncol(X) >= 1L && ncol(X) <= 2000L){
      Xd <- as.matrix(X)
      Yd <- as.matrix(yvar)
      qrX <- qr(Xd)
      dfResid <- nrow(Xd) - qrX$rank
      if(dfResid >= 1L){
        beta <- qr.coef(qrX, Yd)
        if(!anyNA(beta)){
          resid <- Yd - Xd %*% beta
          varResid0 <- mean(colSums(resid^2)) / dfResid
          if(is.finite(varResid0) && varResid0 > 0){
            floorVar <- 1e-6 * varResid0
            residShare <- if(nRandomStruct > 0L) 0.5 * varResid0 else varResid0
            randomPoolShare <- 0.5 * varResid0
            r <- rowMeans(resid)
            
            # Phase 2: for plain ism()-only random terms with an identity
            # relationship matrix, use a one-way ANOVA method-of-moments
            # variance-component estimate (classical unequal-n formula)
            # from the fixed-effects-only residuals, instead of an
            # arbitrary equal split.
            anovaEst <- rep(NA_real_, nRandomStruct)
            if(nRandomStruct > 0L){
              for(u in seq_len(nRandomStruct)){
                ff <- randomFits[[u]]
                isSimpleGrouping <-
                  length(ff$covStruct$factors) == 0L &&
                  Matrix::isDiagonal(ff$Gu) &&
                  isTRUE(all(Matrix::diag(ff$Gu) == 1)) &&
                  length(all.vars(randomExprs[[u]])) == 1L
                if(isSimpleGrouping){
                  gvar <- all.vars(randomExprs[[u]])[1]
                  if(gvar %in% names(data)){
                    grp <- as.factor(data[[gvar]])
                    k <- nlevels(grp)
                    if(k > 1L && k < length(r)){
                      grpMeans <- tapply(r, grp, mean)
                      grpN <- as.numeric(table(grp))
                      grandMean <- mean(r)
                      msBetween <- sum(grpN * (grpMeans - grandMean)^2) / (k - 1L)
                      dfWithin <- length(r) - k
                      if(dfWithin > 0L){
                        msWithin <- sum((r - grpMeans[as.character(grp)])^2) / dfWithin
                        n0 <- (length(r) - sum(grpN^2) / length(r)) / (k - 1L)
                        if(is.finite(n0) && n0 > 0){
                          vcEst <- (msBetween - msWithin) / n0
                          if(is.finite(vcEst) && vcEst > 0) anovaEst[u] <- vcEst
                        }
                      }
                    }
                  }
                }
              }
            }
            
            # Phase 1 fallback: split the remaining pool evenly across
            # every random term that Phase 2 could not estimate directly.
            nFallback <- sum(is.na(anovaEst))
            fallbackShare <- if(nFallback > 0L) randomPoolShare / nFallback else NA_real_
            
            if(isTRUE(rf$covStruct$sigma2_is_default) && isTRUE(rf$covStruct$free[1])){
              rf$covStruct$par[1] <- log(max(residShare, floorVar))
              covStruct[[residualStructIndex]] <- rf$covStruct
            }
            if(nRandomStruct > 0L){
              for(u in seq_len(nRandomStruct)){
                cs <- covStruct[[u]]
                if(isTRUE(cs$sigma2_is_default) && isTRUE(cs$free[1])){
                  share <- if(!is.na(anovaEst[u])) anovaEst[u] else fallbackShare
                  if(is.finite(share)){
                    cs$par[1] <- log(max(share, floorVar))
                    covStruct[[u]] <- cs
                  }
                }
              }
            }
          }
        }
      }
    }
  }, error=function(e) NULL)

  if(responsePrepared){
    standardizedResponse <- (rotationInfo$yOriginal - preparedMean) / preparedSd
    yvar <- rotateRows(standardizedResponse)
  }
  
  # ---- Weights ---------------------------------------------------------
  if(!is.null(.pqlWorkingPrecision)){
    if(length(.pqlWorkingPrecision) == nObs){
      workingPrecision <- .pqlWorkingPrecision[keep]
    }else if(length(.pqlWorkingPrecision) == sum(keep)){
      workingPrecision <- .pqlWorkingPrecision
    }else{
      stop("PQL working precision must have one value per original or retained observation.",
           call.=FALSE)
    }
    if(any(!is.finite(workingPrecision)) || any(workingPrecision <= 0)){
      stop("PQL working precision must be finite and positive.", call.=FALSE)
    }

    if(!is.null(.pqlBaseFactor)){
      if(nrow(.pqlBaseFactor) != length(workingPrecision) ||
         ncol(.pqlBaseFactor) != length(workingPrecision)){
        stop("Cached PQL base factor has incompatible dimensions.", call.=FALSE)
      }
      baseFactor <- .pqlBaseFactor
    }else if(is.null(.pqlBaseW)){
      baseW <- Matrix::Diagonal(n=length(workingPrecision), x=1)
    }else if(nrow(.pqlBaseW) == nObs && ncol(.pqlBaseW) == nObs){
      baseW <- .pqlBaseW[keep, keep, drop=FALSE]
    }else if(nrow(.pqlBaseW) == sum(keep) && ncol(.pqlBaseW) == sum(keep)){
      baseW <- .pqlBaseW
    }else{
      stop("PQL base W must have dimensions equal to either the original or retained number of observations.",
           call.=FALSE)
    }
    if(is.null(.pqlBaseFactor)){
      baseW <- as(as(as(baseW, "dMatrix"), "generalMatrix"), "CsparseMatrix")
      if(!isSymmetric(baseW)){
        stop("PQL base W must be symmetric positive definite.", call.=FALSE)
      }
      baseFactor <- tryCatch(Matrix::chol(baseW), error=function(e) NULL)
      if(is.null(baseFactor)){
        stop("PQL base W must be positive definite.", call.=FALSE)
      }
    }
    weightedFactor <- Matrix::Diagonal(n=length(workingPrecision),
                                       x=sqrt(workingPrecision)) %*% baseFactor
    W <- Matrix::crossprod(weightedFactor)
    W <- as(as(as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    useH <- TRUE
  }else if(missing(W)){
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
  
  if (is.null(emWeight)) {
    taperIters <- min(nIters, 18L)

    if (taperIters <= 1L) {
      emWeight <- 1
    } else {
      emWeight <- rep(0.03, nIters)
      emWeight[seq_len(taperIters)] <- exp(seq(log(1), log(0.03), length.out = taperIters))
    }
  }
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

  if(length(REML) != 1L || !is.logical(REML) || is.na(REML)){
    stop("REML must be a single TRUE/FALSE value.", call.=FALSE)
  }

  if(isTRUE(henderson)){
    # ---- Solver selection ------------------------------------------------
    # "auto" (the default) picks a solver based on the density of the random-
    # effect relationship matrices actually supplied: pedigree-style Ai
    # matrices are typically sparse (a handful of nonzeros per row), while
    # genomic/marker-based relationship matrices are essentially fully dense.
    # The supernodal CHOLMOD factorization amortizes dense fill-in with
    # threaded BLAS-3 kernels and tends to outperform the simplicial LDLT path
    # once any random effect has a dense Gu; otherwise LDLT stays the default.
    solverChoices <- c("auto", "ldlt", "pcg", "cholmod")
    if(length(solver) != 1L || !is.character(solver) || is.na(solver) ||
       !(tolower(solver) %in% solverChoices)){
      stop("solver must be one of 'auto', 'ldlt', 'pcg', or 'cholmod'.", call.=FALSE)
    }
    solver <- tolower(solver)
    if(solver == "auto"){
      hasDenseGu <- length(Ai) > 0L && any(vapply(Ai, function(a){
        n <- nrow(a)
        if(n <= 1L) return(FALSE)
        (Matrix::nnzero(a) / (as.double(n) * as.double(n))) > 0.2
      }, logical(1)))
      solver <- if(hasDenseGu) "cholmod" else "ldlt"
    }
    if(!REML && !(solver %in% c("ldlt", "cholmod"))){
      stop("REML=FALSE (maximum likelihood) currently requires solver='ldlt' or ",
           "solver='cholmod' (solver='auto' resolves to one of these already).",
           call.=FALSE)
    }
    message(crayon::blue(paste("Solver selected:", solver)))
  }else{
    # The direct-inversion engine has no Henderson-style sparse solver
    # choice; it always inverts the n x n phenotypic covariance directly.
    solver <- "direct"
    if(!(computeCi %in% c(0L, 2L))){
      stop("With henderson=FALSE (the direct-inversion engine), computeCi must be 0 (no PEV) or 2 (full PEV, small models only); computeCi=1 (Takahashi selected inverse) is Henderson-only.",
           call.=FALSE)
    }
    message(crayon::blue("Engine selected: direct inversion (henderson=FALSE)"))
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
                obsInfo=obsInfo, solver=solver, REML=REML,
                henderson=henderson, rotation=rotationInfo))
  }

  if(isTRUE(henderson)){
    res <- .Call("_sommer_ai_mme_sp2", PACKAGE="sommer",
                 X, Z, Zind, Ai, yvar, W, useH,
                 residualBlock, localIndex,
                 nIters, tolParConvLL, tolParConvNorm,
                 tolParInv, covStruct, emWeight, stepWeight,
                 verbose, computeCi, solver, pcgTol, pcgMaxIters,
                 pcgTraceProbes, pcgLanczosSteps, REML,
                 responsePrepared, preparedMean, preparedSd,
                 preparedIntercept)

  }else{
    res <- .Call("_sommer_ai_reml_direct_sp2", PACKAGE="sommer",
                 X, Z, Zind, Ai, yvar, W, useH,
                 residualBlock, localIndex,
                 nIters, tolParConvLL, tolParConvNorm,
                 tolParInv, covStruct, emWeight, stepWeight,
                 verbose, computeCi, REML,
                 responsePrepared, preparedMean, preparedSd,
                 preparedIntercept)
  }
  res$engine <- if(isTRUE(henderson)) "henderson" else "direct"

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
  res$REML <- REML
  
  if(length(randomFits) && length(rtermss)){
    names(res$theta) <- c(rtermss, residualLabel)
    names(res$covPar) <- c(rtermss, residualLabel)
    names(res$covStruct) <- c(rtermss, residualLabel)
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
    names(res$covPar) <- residualLabel
    names(res$covStruct) <- residualLabel
    res$args <- list(fixed=fixed, rcov=rcov)
    res$Dtable <- data.frame(type=rep("fixed",length(res$partitionsX)),
                             term=names(res$partitionsX), include=FALSE, average=FALSE)
  }

  if(!is.null(rotationInfo)){
    res$buEngine <- res$bu
    res$uEngine <- res$u
    res$uListEngine <- res$uList

    rotationTerm <- rotationInfo$term
    rotated <- rotationInfo$vectors %*% res$uList[[rotationTerm]]
    rownames(rotated) <- rotationInfo$levels
    colnames(rotated) <- colnames(res$uList[[rotationTerm]])
    res$uList[[rotationTerm]] <- rotated

    publicUNames <- unlist(lapply(rotationInfo$ZOriginal, colnames))
    publicU <- unlist(lapply(res$uList, as.vector), use.names=FALSE)
    res$u <- matrix(publicU, ncol=1L, dimnames=list(publicUNames, NULL))
    res$bu <- rbind(res$b, res$u)
    res$W <- do.call(cbind, c(list(rotationInfo$XOriginal), rotationInfo$ZOriginal))
    res$y <- rotationInfo$yOriginal
    res$rotation <- rotationInfo
  }
  
  class(res) <- "mmes"
  res$covParNative <- get(".covparams_mmes_se", mode="function")(res)
  res
}
