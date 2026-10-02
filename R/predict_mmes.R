#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction

.mmes_engine_contrast <- function(object, D){
  if(!is.null(object$factorScoreInfo)){
    mappings <- object$factorScoreInfo$mappings
    nFixed <- length(object$b)
    nAugmented <- nrow(object$C)
    out <- Matrix::Matrix(0, nrow(D), nAugmented, sparse=TRUE)
    if(nFixed) out[, seq_len(nFixed)] <- D[, seq_len(nFixed), drop=FALSE]
    for(mapping in mappings){
      publicRanges <- mapping$publicRanges
      if(mapping$augmented){
        factorCount <- ncol(mapping$loading)
        q <- nrow(mapping$loading)
        for(level in seq_len(q)){
          publicColumns <- publicRanges[level,1]:publicRanges[level,2]
          for(component in seq_len(factorCount)){
            augmentedRange <- mapping$augmentedRanges[[component]][1L,]
            columns <- augmentedRange[1L]:augmentedRange[2L]
            out[, columns] <- out[, columns, drop=FALSE] +
              D[, publicColumns, drop=FALSE] * mapping$loading[level, component]
          }
          augmentedRange <- mapping$augmentedRanges[[factorCount + level]][1L,]
          columns <- augmentedRange[1L]:augmentedRange[2L]
          out[, columns] <- out[, columns, drop=FALSE] +
            D[, publicColumns, drop=FALSE] * sqrt(mapping$specific[level])
        }
      }else{
        ranges <- mapping$augmentedRanges[[1L]]
        for(level in seq_len(nrow(publicRanges))){
          publicColumns <- publicRanges[level,1]:publicRanges[level,2]
          augmentedRange <- ranges[level,]
          columns <- augmentedRange[1L]:augmentedRange[2L]
          out[, columns] <- D[, publicColumns, drop=FALSE]
        }
      }
    }
    return(out)
  }
  if(is.null(object$rotation)) return(D)

  term <- object$rotation$term
  U <- object$rotation$vectors
  ranges <- object$partitions[[term]]
  for(j in seq_len(nrow(ranges))){
    rr <- ranges[j,1]:ranges[j,2]
    D[,rr] <- D[,rr,drop=FALSE] %*% U
  }
  D
}


.dt_strip_labels <- function(labels, vars){
  vapply(strsplit(labels, ":", fixed=TRUE), function(parts){
    if(length(parts) == length(vars)){
      for(k in seq_along(parts)){
        if(startsWith(parts[k], vars[k])) parts[k] <- substring(parts[k], nchar(vars[k]) + 1L)
      }
    }
    paste(parts, collapse=":")
  }, character(1))
}

.predict_df <- function(object, df){
  if(is.numeric(df) && length(df) == 1L && !is.na(df) && df > 0) return(df)
  if(identical(df, "residual")){
    return(nrow(object$W) - .fixed_rank(object))
  }
  stop("df must be a positive number (Inf for normal-based tests), \"residual\", \"satterthwaite\" or \"kr\".",
       call.=FALSE)
}

.predict_comparisons <- function(pvals, vcov, ids, sed, pairwise, adjust, df, object, D){
  out <- list()
  v <- as.matrix(vcov)
  dv <- diag(v)
  S <- sqrt(pmax(outer(dv, dv, "+") - 2 * v, 0))
  diag(S) <- 0
  dimnames(S) <- list(ids, ids)
  k <- length(ids)
  if(sed){
    out$sed <- S
    ut <- S[upper.tri(S)]
    out$avsed <- if(length(ut)) c(mean=mean(ut), min=min(ut), max=max(ut)) else
      c(mean=NA_real_, min=NA_real_, max=NA_real_)
  }
  if(!isFALSE(pairwise)){
    if(isTRUE(pairwise)){
      if(k > 2000L){
        stop("pairwise=TRUE would produce ", k * (k - 1) / 2, " comparisons; ",
             "compare against a reference level (pairwise=\"<level>\") or use sed=TRUE.",
             call.=FALSE)
      }
      i <- rep(seq_len(k), times=rev(seq_len(k)) - 1L)
      j <- unlist(lapply(seq_len(k), function(a) if(a < k) (a + 1L):k else integer(0)))
    }else if(is.character(pairwise) && length(pairwise) == 1L){
      ref <- match(pairwise, ids)
      if(is.na(ref)) stop("Reference level '", pairwise, "' is not among the predictions.", call.=FALSE)
      i <- setdiff(seq_len(k), ref)
      j <- rep(ref, length(i))
    }else{
      stop("pairwise must be TRUE, FALSE or a single reference level.", call.=FALSE)
    }
    difference <- pvals$predicted.value[i] - pvals$predicted.value[j]
    SED <- S[cbind(i, j)]
    if(is.character(df) && df %in% c("satterthwaite", "kr")){
      # Small-sample df only for contrasts among fixed effects; BLUP contrasts stay normal-based.
      mach <- .mmes_df_machinery(object, second=(df == "kr"))
      nb <- nrow(object$b)
      Dm <- as.matrix(D)
      ddf <- rep(Inf, length(i))
      for(r in seq_along(i)){
        Lfull <- Dm[i[r], ] - Dm[j[r], ]
        if(any(abs(Lfull[-seq_len(nb)]) > 0)) next
        L <- matrix(Lfull[seq_len(nb)], nrow=1)
        if(all(L == 0)) next
        if(df == "kr"){
          ddf[r] <- .kenward_roger(L, as.numeric(object$b), mach)$ddf
          SED[r] <- sqrt(drop(L %*% mach$PhiA %*% t(L)))
        }else{
          ddf[r] <- .satterthwaite_ddf(L, mach)
        }
      }
    }else{
      ddf <- rep(.predict_df(object, df), length(i))
    }
    statistic <- difference / SED
    p <- ifelse(is.finite(ddf), 2 * stats::pt(-abs(statistic), ddf), 2 * stats::pnorm(-abs(statistic)))
    out$pairwise <- data.frame(level1=ids[i], level2=ids[j], difference=difference,
                               SED=SED, statistic=statistic, df=ddf, p.value=p,
                               p.adjusted=stats::p.adjust(p, method=adjust),
                               stringsAsFactors=FALSE)
  }
  out
}

"predict.mmes" <- function(object, Dtable=NULL, D, levels=NULL, sed=FALSE,
                           pairwise=FALSE, adjust="none", df=Inf, ...){
  if(is.character(D)){classify <- D}else{classify="id"} # save a copy before D is overwriten
  # Prediction variances are obtained by solving C %*% X = t(D) for the
  # handful of rows in D (see predict_mmes_vcov_cpp), so the complete or
  # selected-inverse Ci is not required; only the stored C/Cscale are used.
  if(is.null(object$C) || nrow(object$C) == 0){
    stop("The predict function requires the mixed-model coefficient matrix C to be available in the object.", call.=FALSE)
  }
  if(length(levels) && !is.character(D)){
    stop("levels can only be used when D is a character classify term.", call.=FALSE)
  }
  # complete the Dtable withnumber of effects in each term
  xEffectN <- lapply(object$partitionsX, as.vector)
  nz <- unlist(lapply(object$uList,function(x){nrow(x)*ncol(x)}))
  # add a value but if there's no intercept consider it
  lenEff <- length(xEffectN)
  toAdd <- xEffectN[[lenEff]];
  if(length(toAdd) > 0){ # there's a value in the last element of xEffectN
    add <- max(toAdd) + 1 # this is our last column of fixed effects
  }else{ # there's not a value in the last element of xEffectN
    if(lenEff > 1){
      toAdd2 <- xEffectN[[lenEff-1]]
      add <- max(toAdd2) + 1
    }
  }
  if(!is.null(nz)){
    zEffectsN <- list()
    for(i in 1:length(nz)){
      end= add + nz[i] - 1
      zEffectsN[[i]] <- add:end
      add = end + 1
    }
    names(zEffectsN) <- names(nz)
    effectsN = c(xEffectN,zEffectsN)
  }else{
    effectsN = xEffectN
  }
  
  # fill the Dt table for rules
  if(is.null(Dtable) & is.character(D) ){ # if user didn't provide the Dtable but D is character
    Dtable <- object$Dtable # we extract it from the model
    termsInDtable <- apply(data.frame(Dtable$term),1,function(xx){all.names(as.formula(paste0("~",xx)))})
    termsInDtable <- lapply(termsInDtable, function(x){intersect(x,colnames(object$data))})
    termsInDtable <- lapply(termsInDtable, function(x){return(unique(c(x, paste(x, collapse=":"))))})
    # term identified
    termsInDtableN <- unlist(lapply(termsInDtable,length))
    pickTerm <- which( unlist(lapply(termsInDtable, function(xxx){ifelse(length(which(xxx == D)) > 0, 1, 0)})) > 0)
    if(length(pickTerm) == 0){
      isInDf <- which(colnames(object$data) %in% D)
      if(length(isInDf) > 0){ # is in data frame but not in model
        stop(paste("Predict:",classify,"not in the model but present in the original dataset. You may need to provide the
                   Dtable argument to know how to predict", classify), call. = FALSE)
      }else{
        stop(paste("Predict:",classify,"not in the model and not present in the original dataset. Please correct D."), call. = FALSE)
      }
    }
    # check if the term to predict is fixed or random
    pickTermIsFixed = ifelse("fixed" %in% Dtable[pickTerm,"type"], TRUE,FALSE)
    ## 1) when we predict a random effect, fixed effects are purely "average"
    if(!pickTermIsFixed){Dtable[which(Dtable$type %in% "fixed"),"average"]=TRUE; Dtable[pickTerm,"include"]=TRUE}
    ## 2) when we predict a fixed effect, random effects are ignored and the fixed effect is purely "include"
    if(pickTermIsFixed){Dtable[pickTerm,"include"]=TRUE}
    ## 3) for predicting a random effect, the interactions are ignored and only main effect is "included", then we follow 1)
    ## 4) for a model with pure interaction trying to predict a main effect of the interaction we "include" and "average" the interaction and follow 1)
    if(length(pickTerm) == 1){ # only one effect identified
      if(termsInDtableN[pickTerm] > 1){ # we are in #4 (is a pure interaction model)
        Dtable[pickTerm,"average"]=TRUE
      }
    }else{# more than 1, there's main effect and interactions, situation #3
      main <- which(termsInDtableN[pickTerm] == min(termsInDtableN[pickTerm])[1])
      Dtable[pickTerm,"include"]=FALSE;  Dtable[pickTerm,"average"]=FALSE # reset
      Dtable[pickTerm[main],"include"]=TRUE
    }
  }
  ## if user has provided D as a classify then we create the D matrix
  if(is.character(D)){
    if(is.null(Dtable$levels)) Dtable$levels <- vector("list", nrow(Dtable))
    if(length(levels)){
      if(!is.list(levels) || is.null(names(levels)) || any(names(levels) == "")){
        stop("levels must be a named list whose names are terms of the Dtable.", call.=FALSE)
      }
      unknown <- setdiff(names(levels), Dtable$term)
      if(length(unknown)){
        stop("levels names not found in Dtable$term: ", paste(unknown, collapse=", "),
             ". Valid terms are: ", paste(Dtable$term, collapse=", "), call.=FALSE)
      }
      for(nm in names(levels)){
        if(!is.null(levels[[nm]])) Dtable$levels[[match(nm, Dtable$term)]] <- levels[[nm]]
      }
    }
    hasLevels <- !vapply(Dtable$levels, is.null, logical(1))
    termVars <- lapply(Dtable$term, function(xx){
      v <- intersect(all.names(as.formula(paste0("~", xx))), colnames(object$data))
      unique(c(v, paste(v, collapse=":")))
    })
    classifyVars <- strsplit(D, ":", fixed=TRUE)[[1]]
    # create model matrices to form D
    P <- sparse.model.matrix(as.formula(paste0("~",D,"-1")), data=object$data)
    colnames(P) <- .dt_strip_labels(colnames(P), classifyVars)
    tP <- t(P)
    W <- object$W
    D = tP %*% W
    colnames(D) <- c(rownames(object$b),rownames(object$u))
    cd <- colnames(D)
    # levels of included classify terms select (and may add) the prediction rows
    classifyRows <- which(Dtable$include & hasLevels & vapply(seq_len(nrow(Dtable)), function(i){
      if(Dtable$type[i] == "fixed") setequal(strsplit(Dtable$term[i], ":", fixed=TRUE)[[1]], classifyVars)
      else classify %in% termVars[[i]]
    }, logical(1)))
    if(length(classifyRows)){
      rowLevels <- unique(as.character(unlist(Dtable$levels[classifyRows])))
      newRows <- setdiff(rowLevels, rownames(D))
      unknown <- setdiff(newRows, cd)
      if(length(unknown)){
        stop("Levels not present in the data nor among the fitted coefficients of '", classify,
             "': ", paste(unknown, collapse=", "), call.=FALSE)
      }
      if(length(newRows)){
        D <- rbind(D, Matrix::Matrix(0, nrow=length(newRows), ncol=ncol(D), sparse=TRUE,
                                     dimnames=list(newRows, cd)))
      }
      D <- D[rowLevels, , drop=FALSE]
    }
    rd <- rownames(D)
    coocIndex <- match(rd, rownames(tP))
    noData <- is.na(coocIndex)
    cooc <- tP[ifelse(noData, 1L, coocIndex), , drop=FALSE]
    interceptColumn <- unique(c(grep("Intercept",rownames(object$b) ))) # ,which(rownames(object$b)=="1")
    for(jRow in 1:nrow(D)){ # for each effect add 1's where missing
      myMatch <- which(cd == rd[jRow])
      if(length(myMatch) > 0){D[jRow,myMatch]=1}
    }
    # apply rules in Dtable
    for(iRow in 1:nrow(Dtable)){
      w <- effectsN[[iRow]]
      if(!length(w)) next
      term <- Dtable[iRow,"term"]
      isFixed <- Dtable[iRow,"type"] == "fixed" && term != "1"
      vars <- strsplit(term, ":", fixed=TRUE)[[1]]
      lv <- Dtable$levels[[iRow]]
      full <- NULL
      if(isFixed){
        full <- sparse.model.matrix(reformulate(term, intercept=FALSE), data=object$data)
        fullLabels <- .dt_strip_labels(colnames(full), vars)
        colLabels <- .dt_strip_labels(cd[w], vars)
      }else{
        colLabels <- cd[w]
      }
      isCovariate <- isFixed && all(vars %in% colnames(object$data)) &&
        all(vapply(vars, function(x) is.numeric(object$data[[x]]), logical(1)))
      if(!is.null(lv) && isCovariate){
        if(!is.numeric(lv) || !(length(lv) %in% c(1L, length(w)))){
          stop("levels for covariate term '", term, "' must be numeric values at which to predict.",
               call.=FALSE)
        }
        D[,w] <- matrix(rep(lv, each=nrow(D)), nrow=nrow(D))
        next
      }
      if(!is.null(lv)){
        available <- if(isFixed) fullLabels else colLabels
        unknown <- setdiff(as.character(lv), available)
        if(length(unknown)){
          stop("Unknown levels for term '", term, "': ", paste(unknown, collapse=", "), call.=FALSE)
        }
        if(!Dtable[iRow,"include"] && !Dtable[iRow,"average"]){
          warning("levels given for term '", term, "' are ignored because it is neither included nor averaged.",
                  call.=FALSE)
        }
      }
      selected <- if(is.null(lv)) rep(TRUE, length(w)) else colLabels %in% as.character(lv)
      fullSelected <- if(isFixed){
        if(is.null(lv)) rep(TRUE, ncol(full)) else fullLabels %in% as.character(lv)
      }else NULL
      if(Dtable[iRow,"include"]){
        subD <- D[,w,drop=FALSE]
        subD[which(subD > 0, arr.ind = TRUE)] = 1
        if(any(!selected)) subD[, which(!selected)] <- 0
        if(Dtable[iRow,"average"]){
          if(isFixed){
            # cells of the full term (reference levels included) present for each prediction row
            present <- as.matrix(cooc %*% full[, fullSelected, drop=FALSE]) > 0
            averageN <- rowSums(present)
            averageN[noData] <- sum(fullSelected)
          }else{
            averageN <- Matrix::rowSums(subD > 0)
          }
          subD <- Matrix::Diagonal(x=1/pmax(averageN, 1)) %*% subD
        }
        D[,w] <- subD
      }else{
        if(Dtable[iRow,"average"]){
          averageN <- if(isFixed) sum(fullSelected) else sum(selected)
          D[,w] <- matrix(rep(ifelse(selected, 1/averageN, 0), each=nrow(D)), nrow=nrow(D))
        }else{
          D[,w] <- D[,w] * 0
        }
      }
    }
    if(length(interceptColumn) > 0){D[,interceptColumn] = 1}
  }else{ }# user has provided D as a matrix to do direct multiplication
  ## calculate predictions and standard errors
  bu <- object$bu
  predicted.value <- D %*% bu
  vcov <- predict_mmes_vcov_cpp(object, .mmes_engine_contrast(object, D))
  std.error <- sqrt(diag(vcov))
  pvals <- data.frame(id=rownames(D),predicted.value=predicted.value[,1], std.error=std.error)
  if(is.character(classify)){colnames(pvals)[1] <- classify}
  out <- list(pvals=pvals,D=D,vcov=vcov, Dtable=Dtable)
  if(sed || !isFALSE(pairwise)){
    out <- c(out, .predict_comparisons(pvals, vcov, rownames(D), sed, pairwise,
                                       adjust, df, object, D))
  }
  return(out)
}

"print.predict.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
                 The predictions are obtained by averaging/aggregating across
                 the hypertable calculated from model terms constructed solely
                 from factors in the include sets. You can customize the model
                 terms used with the 'Dtable' argument.\n")
  ))
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}
