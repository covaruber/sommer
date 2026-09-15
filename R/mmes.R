mmes <- function(fixed, random, rcov, data, W,
                 nIters=30, tolParConvLL=1e-04,
                 tolParConvNorm=1e-04, tolParInv=1e-06,
                 naMethodX="exclude", naMethodY="exclude",
                 returnParam=FALSE, dateWarning=TRUE,
                 verbose=TRUE, stepWeight=NULL, emWeight=NULL,
                 contrasts=NULL, getPEV=TRUE, henderson=TRUE, computeCi=0){

  # This interface is intentionally Henderson-only. Direct-inversion/MNR
  # methods should use a separate front end.
  if(!isTRUE(henderson)){
    stop("This mmes() interface is Henderson-only. Use the separate MNR/direct-inversion mmer interface for henderson=FALSE.",
         call. = FALSE)
  }

  desc <- utils::packageDescription("sommer")
  my.date <- as.Date(desc$Date) + 90
  if(dateWarning && Sys.Date() > my.date){
    cat("Version out of date. Please update sommer to the newest version using:\n",
        "install.packages('sommer') in a new session\n",
        "Use the 'dateWarning' argument to disable the warning message.\n", sep="")
  }

  if(missing(data)){
    data <- environment(fixed)
    cat("data argument not provided\n")
  }else{
    data <- as.data.frame(data)
  }

  if(missing(rcov)) rcov <- as.formula("~units")

  # Handle missingness with sommer's existing helper.
  dataor <- data
  provdat <- subdata(data, fixed=fixed,
                     na.method.Y=naMethodY,
                     na.method.X=naMethodX)
  data <- provdat$datar
  nonMissing <- provdat$good

  # Dedicated observation identity used by residual structures.
  data$units <- levels(as.factor(paste0("u", seq_len(nrow(data)))))

  response <- strsplit(as.character(fixed[2]), split="[+]")[[1]]
  responsef <- as.formula(paste(response, "~1"))
  mfna <- try(model.frame(responsef, data=data, na.action=na.pass), silent=TRUE)
  if(inherits(mfna, "try-error")){
    stop("Please provide the 'data' argument for all variables used in the model.",
         call. = FALSE)
  }
  mfna <- eval(mfna, data, parent.frame())
  yvar <- sparse.model.matrix(as.formula(paste("~", response, "-1")), data)
  if(ncol(yvar) == 1) colnames(yvar) <- response

  # ------------------------------------------------------------------
  # Random structures
  # ------------------------------------------------------------------
  Z <- list()
  Ai <- list()
  covStruct <- list()
  Zind <- numeric()
  rTermsNames <- list()
  rtermss <- character()

  if(!missing(random)){
    yuyu <- strsplit(as.character(random[2]), split="[+]")[[1]]
    rtermss <- vapply(yuyu, function(x){
      strsplit(as.character((as.formula(paste("~", x)))[2]), split="[+]")[[1]]
    }, character(1))

    for(u in seq_along(rtermss)){
      term <- rtermss[u]
      checkvs <- intersect(all.names(as.formula(paste0("~", term))), c("vsm","spl2Dc"))
      if(length(checkvs) == 0L){
        term <- paste0("sommer::vsm(sommer::ism(", term, "))")
        rtermss[u] <- term
      }

      ff <- eval(parse(text=term), data, parent.frame())
      if(is.null(ff$covStruct) || !identical(ff$covStruct$type, "kron") ||
         is.null(ff$covStruct$descriptor_version) || ff$covStruct$descriptor_version < 2L){
        stop("All random covariance terms must use the CovarianceFactor v2 vsm() descriptor interface.",
             call. = FALSE)
      }

      Zi <- lapply(ff$Z, function(x){
        if(nrow(x) != length(nonMissing)) x[nonMissing,,drop=FALSE] else x
      })

      Z <- c(Z, Zi)
      
      pGu <- to_sparse(ff$Gu) 
      attr(pGu, "inverse") <- TRUE
      Ai[[u]] <- pGu# ff$Gu 
      covStruct[[u]] <- ff$covStruct
      Zind <- c(Zind, rep(u, length(Zi)))

      s2 <- paste(all.vars(as.formula(paste("~", term))), collapse=":")
      rTermsNames[[u]] <- paste(s2, ff$covStruct$par_names, sep=":")
    }
  }

  nRandomStruct <- length(covStruct)

  # ------------------------------------------------------------------
  # Residual structure
  # ------------------------------------------------------------------
  rcov_terms <- strsplit(as.character(rcov[2]), split="[+]")[[1]]
  if(length(rcov_terms) != 1L){
    stop("The new Henderson interface currently accepts one residual vsm() term. ",
         "Use arbitrary Kronecker products inside that vsm() term instead of summing residual terms.",
         call. = FALSE)
  }

  rterm <- strsplit(as.character((as.formula(paste("~", rcov_terms[1])))[2]),
                    split="[+]")[[1]]
  checkvs <- intersect(all.names(as.formula(paste0("~", rterm))),
                       c("vsm","gvs","spl2Da","spl2Db"))
  if(length(checkvs) == 0L){
    rterm <- paste0("sommer::vsm(sommer::ism(", rterm, "))")
  }

  rf <- eval(parse(text=rterm), data, parent.frame())
  if(is.null(rf$covStruct) || !identical(rf$covStruct$type, "kron") ||
     is.null(rf$covStruct$descriptor_version) || rf$covStruct$descriptor_version < 2L){
    stop("The residual covariance term must use the CovarianceFactor v2 vsm() descriptor interface.",
         call. = FALSE)
  }
  if(is.null(rf$residualLocalIndex)){
    # A bare ~units model has no covariance-shaping factor. Its local index is 1.
    rf$residualLocalIndex <- rep(1L, nrow(data))
  }

  # Preserve the historical more generous residual starting scale, but only
  # on the single product-level variance. The parameter is log(sigma2).
  if(isTRUE(rf$covStruct$free[1])){
    rf$covStruct$par[1] <- rf$covStruct$par[1] + log(5)
  }

  residualStructIndex <- nRandomStruct + 1L
  covStruct[[residualStructIndex]] <- rf$covStruct

  residualVariables <- setdiff(all.vars(as.formula(paste("~", rterm))), "units")
  pairingVariables <- setdiff(names(data), c(response, residualVariables, "units"))

  localIndex <- as.integer(rf$residualLocalIndex)
  if(length(localIndex) != nrow(data)){
    stop("Residual local-index vector has incompatible length.", call. = FALSE)
  }

  if(length(pairingVariables)){
    baseKey <- do.call(paste, c(data[pairingVariables], sep="\r"))
  }else{
    baseKey <- rep("all", nrow(data))
  }

  # If replicated rows share the same block/local coordinate, split replicate
  # occurrences into independent blocks rather than silently overwriting them.
  pairLocal <- paste(baseKey, localIndex, sep="\r")
  occurrence <- ave(seq_along(pairLocal), pairLocal, FUN=seq_along)
  blockKey <- paste(baseKey, occurrence, sep="\r")
  residualBlock <- match(blockKey, unique(blockKey))

  if(anyDuplicated(paste(residualBlock, localIndex, sep=":"))){
    stop("Internal residual-layout error: a block contains duplicate local covariance coordinates.",
         call. = FALSE)
  }

  s2 <- paste(all.vars(as.formula(paste("~", rterm))), collapse=":")
  rTermsNames[[residualStructIndex]] <-
    paste(s2, rf$covStruct$par_names, sep=":")

  # ------------------------------------------------------------------
  # Fixed-effect design
  # ------------------------------------------------------------------
  data$`1` <- 1
  newfixed <- fixed
  fixedTerms <- gsub(" ", "", strsplit(as.character(fixed[3]), split="[+-]")[[1]])
  mf <- try(model.frame(newfixed, data=data, na.action=na.pass), silent=TRUE)
  mf <- eval(mf, parent.frame())
  X <- Matrix::sparse.model.matrix(newfixed, mf, contrasts.arg=contrasts)

  partitionsX <- list()
  for(ix in seq_along(fixedTerms)){
    effs <- colnames(Matrix::sparse.model.matrix(
      as.formula(paste("~", fixedTerms[ix], "-1")), mf))
    effs2 <- colnames(Matrix::sparse.model.matrix(
      as.formula(paste("~", fixedTerms[ix])), mf))
    partitionsX[[ix]] <- matrix(which(colnames(X) %in% c(effs,effs2)), nrow=1)
  }
  names(partitionsX) <- fixedTerms

  classColumns <- lapply(data, class)
  for(ix in seq_along(fixedTerms)){
    colnamesBase <- colnames(X)[partitionsX[[ix]]]
    colnamesBaseList <- strsplit(colnamesBase, ":")
    toRemoveList <- strsplit(fixedTerms[ix], ":")[[1]]
    if(!("1" %in% toRemoveList)){
      toRemoveList <- all.vars(as.formula(paste("~", paste(toRemoveList, collapse="+"))))
    }
    for(j in seq_along(toRemoveList)){
      if(toRemoveList[[j]] %in% names(classColumns) &&
         classColumns[[toRemoveList[[j]]]] != "numeric"){
        nc <- nchar(gsub(" ", "", toRemoveList[[j]], fixed=TRUE))
        colnamesBaseList <- lapply(colnamesBaseList, function(h){
          if(is.na(h[j])) return(h)
          if(length(grep(toRemoveList[[j]], h[j])) == 1 && nchar(h[j]) > nc){
            h[j] <- substr(h[j], 1+nc, nchar(h[j]))
          }
          h
        })
      }
    }
    colnames(X)[partitionsX[[ix]]] <-
      unlist(lapply(lapply(colnamesBaseList, na.omit), paste, collapse=":"))
  }

  step1 <- gsub(" ", "", strsplit(as.character(fixed[3]), split="[-]")[[1]])
  step2 <- unlist(lapply(step1, function(x) strsplit(x, split="[+]")[[1]]))
  if(length(intersect(c("1","-1"), step2)) == 0L) colnames(X)[1] <- "Intercept"

  # ------------------------------------------------------------------
  # Weights / information weights
  # ------------------------------------------------------------------
  if(missing(W)){
    x <- data.frame(d=as.factor(seq_len(length(yvar))))
    W <- sparse.model.matrix(~d-1, x)
    useH <- FALSE
  }else{
    W <- as(as(as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    useH <- TRUE
  }

  if(is.null(emWeight)){
    emWeight <- stan(logspace(seq(1,-1,-2/nIters), p=3))
  }
  if(is.null(stepWeight)){
    w <- which(emWeight <= .5)
    stepWeight <- rep(.9, nIters)
    if(nIters > 1){
      if(length(w) > 1) stepWeight[w[1:2]] <- c(.5,.7)
      else stepWeight[1:2] <- c(.5,.7)
    }
  }

  # Henderson requires inverse relationship matrices.
  if(length(Ai)){
    nInverses <- sum(vapply(Ai, function(x) isTRUE(attr(x, "inverse")), logical(1)))
    if(nInverses != length(Ai)){
      stop("The Henderson algorithm requires every Gu relationship matrix to be supplied as an inverse matrix with attr(Gu,'inverse')=TRUE.",
           call. = FALSE)
    }
  }

  if(returnParam){
    return(list(
      yvar=yvar, X=X, Z=Z, Zind=Zind, Ai=Ai,
      W=W, useH=useH,
      residualBlock=residualBlock,
      residualIndex=localIndex,
      nIters=nIters,
      tolParConvLL=tolParConvLL,
      tolParConvNorm=tolParConvNorm,
      tolParInv=tolParInv,
      verbose=verbose,
      covStruct=covStruct,
      stepWeight=stepWeight,
      emWeight=emWeight,
      rtermss=rtermss,
      partitionsX=partitionsX,
      getPEV=getPEV,
      rTermsNames=rTermsNames
    ))
  }

  # if(verbose) message("henderson's formulation used")

  # At this boundary all covariance models have been compiled by vsm() to
  # universal CovarianceFactor descriptors. ai_mme_sp2() receives only generic
  # evaluator/derivative/report/trust specifications and arbitrary Kronecker
  # layouts; it does not dispatch on covariance-model names.
  res <- .Call("_sommer_ai_mme_sp2", PACKAGE="sommer",
               X, Z, Zind,
               Ai, yvar,
               W, useH,
               residualBlock, localIndex,
               nIters, tolParConvLL, tolParConvNorm,
               tolParInv, covStruct,
               emWeight, stepWeight,
               verbose, computeCi)

  rownames(res$b) <- colnames(X)
  if(!missing(random) && length(res$u)){
    rownames(res$u) <- unlist(lapply(Z, colnames))
  }
  rownames(res$bu) <- c(rownames(res$b), rownames(res$u))
  rownames(res$monitor) <- unlist(rTermsNames)

  res$data <- data
  res$y <- yvar
  res$partitionsX <- partitionsX
  res$covStruct <- covStruct

  if(!missing(random) && length(rtermss)){
    names(res$theta) <- c(rtermss, rterm)
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
    res$Dtable <- data.frame(
      type=c(rep("fixed", length(res$partitionsX)),
             rep("random", length(res$partitions))),
      term=c(names(res$partitionsX), names(res$partitions)),
      include=FALSE, average=FALSE
    )
  }else{
    names(res$theta) <- rterm
    res$args <- list(fixed=fixed, rcov=rcov)
    res$Dtable <- data.frame(
      type=rep("fixed", length(res$partitionsX)),
      term=names(res$partitionsX),
      include=FALSE, average=FALSE
    )
  }

  class(res) <- "mmes"
  res
}
