.mmes_prepare_weights <- function(W, WWasMissing, nObs, keep, nRetained,
                                  weightBlocks, pqlPrecision=NULL,
                                  pqlBaseW=NULL, pqlBaseFactor=NULL){
  if(!is.null(pqlPrecision)){
    if(length(pqlPrecision) == nObs){
      workingPrecision <- pqlPrecision[keep]
    }else if(length(pqlPrecision) == sum(keep)){
      workingPrecision <- pqlPrecision
    }else{
      stop("PQL working precision must have one value per original or retained observation.",
           call.=FALSE)
    }
    if(any(!is.finite(workingPrecision)) || any(workingPrecision <= 0)){
      stop("PQL working precision must be finite and positive.", call.=FALSE)
    }

    if(!is.null(pqlBaseFactor)){
      if(nrow(pqlBaseFactor) != length(workingPrecision) ||
         ncol(pqlBaseFactor) != length(workingPrecision)){
        stop("Cached PQL base factor has incompatible dimensions.", call.=FALSE)
      }
      baseFactor <- pqlBaseFactor
    }else if(is.null(pqlBaseW)){
      baseW <- Matrix::Diagonal(n=length(workingPrecision), x=1)
    }else if(nrow(pqlBaseW) == nObs && ncol(pqlBaseW) == nObs){
      baseW <- pqlBaseW[keep, keep, drop=FALSE]
    }else if(nrow(pqlBaseW) == sum(keep) && ncol(pqlBaseW) == sum(keep)){
      baseW <- pqlBaseW
    }else{
      stop("PQL base W must have dimensions equal to either the original or retained number of observations.",
           call.=FALSE)
    }
    if(is.null(pqlBaseFactor)){
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
  }else if(WWasMissing){
    W <- Matrix::Diagonal(n=nRetained, x=1)
    useH <- FALSE
  }else{
    if(nrow(W) == nObs && ncol(W) == nObs) W <- W[keep, keep, drop=FALSE]
    else if(nrow(W) != sum(keep) || ncol(W) != sum(keep))
      stop("W must have dimensions equal to either the original or retained number of observations.", call.=FALSE)
    W <- as(as(as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    useH <- TRUE
  }
  if(!is.null(weightBlocks)){
    if(WWasMissing && is.null(pqlPrecision)){
      stop("weights block formulas describe a supplied W matrix; provide W as well.", call.=FALSE)
    }
    if(length(weightBlocks) != nrow(W)){
      stop("The weights block formula does not match the retained W dimensions.", call.=FALSE)
    }
    nonzeroWeights <- Matrix::summary(W)
    if(nrow(nonzeroWeights) && any(weightBlocks[nonzeroWeights$i] !=
                                   weightBlocks[nonzeroWeights$j])){
      stop("W contains nonzero entries across groups declared independent by weights.", call.=FALSE)
    }
  }
  list(W=W, useH=useH)
}
