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

.mmes_prediction_memory <- function(object, predictionRows, fixedOnly, VarU,
                                    sed, pairwise, maxMemoryMB){
  engineColumns <- nrow(object$C)
  covarianceCopies <- if(fixedOnly || !VarU) 3 else 10
  workspaceMB <- 8e-6 * (3 * engineColumns * predictionRows +
    (covarianceCopies + if(sed || !isFALSE(pairwise)) 4 else 0) * predictionRows^2)
  if(!is.finite(workspaceMB) || workspaceMB > maxMemoryMB){
    stop("Prediction for ", predictionRows, " rows needs approximately ",
         format(round(workspaceMB, 1), scientific=FALSE),
         " MB of dense solve/covariance workspace, exceeding maxMemoryMB=",
         maxMemoryMB, ". Each full covariance alone needs ",
         format(round(8e-6 * predictionRows^2, 1), scientific=FALSE),
         " MB. Request fewer levels or smaller explicit contrast batches",
         if(VarU && !fixedOnly) ", or use VarU=FALSE" else "",
         "; increase maxMemoryMB only if sufficient memory is available. ",
         "PEV=FALSE does not remove the full vcov output.",
         call.=FALSE)
  }
  invisible(NULL)
}

.mmes_observed_predictions <- function(object, classify, Dtable, levels,
                                       VarU, sed, pairwise, maxMemoryMB){
  if(length(classify) != 1L || is.na(classify) || !nzchar(classify)){
    stop("observed=TRUE requires one categorical fixed-effect classify term in D.",
         call.=FALSE)
  }
  predictionRules <- if(is.null(Dtable)) object$Dtable else Dtable
  requiredColumns <- c("type", "term", "include", "average")
  if(!is.data.frame(predictionRules) ||
     !all(requiredColumns %in% names(predictionRules)) ||
     nrow(predictionRules) != nrow(object$Dtable) ||
     anyDuplicated(predictionRules$term) ||
     !setequal(predictionRules$term, object$Dtable$term)){
    stop("With observed=TRUE, Dtable must contain each fitted term exactly once, ",
         "with type, term, include and average columns; start from object$Dtable.",
         call.=FALSE)
  }
  predictionRules <- predictionRules[match(object$Dtable$term, predictionRules$term), , drop=FALSE]
  if(!identical(as.character(predictionRules$type), as.character(object$Dtable$type)) ||
     !is.logical(predictionRules$include) || !is.logical(predictionRules$average) ||
     anyNA(predictionRules$include) || anyNA(predictionRules$average)){
    stop("Dtable types must match the fitted terms, and include/average must be nonmissing logical flags.",
         call.=FALSE)
  }
  predictionRules$levels <- NULL
  fixedRules <- predictionRules$type == "fixed"
  if(any(predictionRules$include[fixedRules] != object$Dtable$include[fixedRules]) ||
     any(predictionRules$average[fixedRules] != object$Dtable$average[fixedRules])){
    stop("With observed=TRUE, fixed Dtable flags must remain as in object$Dtable: ",
         "all fixed terms use retained-record averages. Only random-term flags may be changed. ",
         "Use explicit estimable contrasts for other fixed-effect rules.", call.=FALSE)
  }
  classifyVariables <- strsplit(classify, ":", fixed=TRUE)[[1L]]
  fixedVariables <- all.vars(stats::delete.response(stats::terms(object$args$fixed)))
  if(anyDuplicated(classifyVariables) ||
     !all(classifyVariables %in% intersect(fixedVariables, names(object$data)))){
    stop("observed=TRUE requires D to name categorical variables in the fitted fixed formula.",
         call.=FALSE)
  }
  classifyData <- object$data[, classifyVariables, drop=FALSE]
  if(!all(vapply(classifyData, function(values) is.factor(values) || is.character(values),
                 logical(1))) || anyNA(classifyData)){
    stop("observed=TRUE requires nonmissing categorical classify variables in the retained fitting data.",
         call.=FALSE)
  }
  if(nrow(classifyData) != nrow(object$W)){
    stop("The retained prediction data do not match the fitted design rows.", call.=FALSE)
  }
  groupCodes <- lapply(classifyData, function(values) match(values, unique(values)))
  groupKeys <- do.call(paste, c(groupCodes, sep="\r"))
  firstRows <- which(!duplicated(groupKeys))
  groupIndices <- match(groupKeys, groupKeys[firstRows])
  predictionLabels <- do.call(paste, c(lapply(classifyData[firstRows, , drop=FALSE],
                                              as.character), sep=":"))
  if(anyDuplicated(predictionLabels)){
    stop("Observed classify labels are ambiguous when joined with ':'. ",
         "Rename the conflicting levels or supply explicit contrasts.", call.=FALSE)
  }
  selectedGroups <- seq_along(firstRows)
  if(length(levels)){
    if(!is.list(levels) || is.null(names(levels)) || anyNA(names(levels)) ||
       anyDuplicated(names(levels)) ||
       !all(names(levels) %in% c(classify, predictionRules$term[!fixedRules]))){
      stop("With observed=TRUE, levels must be a named list for D and/or fitted random terms.",
           call.=FALSE)
    }
  }
  classifyLevels <- levels[[classify]]
  if(!is.null(classifyLevels)){
    if(!is.character(classifyLevels) || !length(classifyLevels) ||
       anyNA(classifyLevels) || anyDuplicated(classifyLevels)){
      stop("With observed=TRUE, classify levels must contain unique observed labels.",
           call.=FALSE)
    }
    selectedGroups <- match(classifyLevels, predictionLabels)
    if(anyNA(selectedGroups)){
      stop("Requested levels were not observed in the retained fitting data: ",
           paste(head(classifyLevels[is.na(selectedGroups)], 5L), collapse=", "),
           ". observed=TRUE does not add unobserved combinations.", call.=FALSE)
    }
  }
  activeRandom <- !fixedRules & (predictionRules$include | predictionRules$average)
  .mmes_prediction_memory(object, length(selectedGroups), !any(activeRandom), VarU, sed,
                          pairwise, maxMemoryMB)
  selectedRows <- match(groupIndices, selectedGroups)
  usedRows <- which(!is.na(selectedRows))
  groupSizes <- tabulate(groupIndices, nbins=length(firstRows))
  averages <- Matrix::sparseMatrix(i=selectedRows[usedRows], j=usedRows,
    x=1 / groupSizes[groupIndices[usedRows]],
    dims=c(length(selectedGroups), nrow(object$W)))
  fixedColumns <- nrow(object$b)
  fixedContrasts <- averages %*% object$W[, seq_len(fixedColumns), drop=FALSE]
  contrasts <- cbind(fixedContrasts, Matrix::Matrix(0, nrow(fixedContrasts),
    nrow(object$bu) - fixedColumns, sparse=TRUE))
  for(ruleIndex in which(!fixedRules)){
    randomIndex <- match(predictionRules$term[ruleIndex], names(object$partitions))
    ranges <- object$partitions[[randomIndex]]
    randomColumns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate)
      seq.int(ranges[coordinate, 1L], ranges[coordinate, 2L])), use.names=FALSE)
    selectedColumns <- rep(TRUE, length(randomColumns))
    requestedLevels <- levels[[predictionRules$term[ruleIndex]]]
    if(!is.null(requestedLevels)){
      effects <- object$uList[[randomIndex]]
      effectRows <- rep(rownames(effects), times=ncol(effects))
      effectCoordinates <- rep(colnames(effects), each=nrow(effects))
      effectLabels <- list(rownames(object$bu)[randomColumns], effectRows,
        paste(effectCoordinates, effectRows, sep=":"),
        paste(effectRows, effectCoordinates, sep=":"), effectCoordinates)
      if(!is.character(requestedLevels) || !length(requestedLevels) ||
         anyNA(requestedLevels) || anyDuplicated(requestedLevels)){
        stop("Random-term levels must be unique nonmissing character labels for '",
             predictionRules$term[ruleIndex], "'.", call.=FALSE)
      }
      unknownLevels <- setdiff(requestedLevels, unlist(effectLabels, use.names=FALSE))
      if(length(unknownLevels)){
        stop("Unknown levels for term '", predictionRules$term[ruleIndex], "': ",
             paste(head(unknownLevels, 5L), collapse=", "), call.=FALSE)
      }
      selectedColumns <- Reduce(`|`, lapply(effectLabels, function(labels)
        labels %in% requestedLevels))
    }
    if(!activeRandom[ruleIndex]){
      if(!is.null(requestedLevels)){
        warning("levels given for term '", predictionRules$term[ruleIndex],
                "' are ignored because it is neither included nor averaged.", call.=FALSE)
      }
      next
    }
    termColumns <- randomColumns[selectedColumns]
    if(!length(termColumns)){
      stop("No fitted coefficients selected for random term '",
           predictionRules$term[ruleIndex], "'.", call.=FALSE)
    }
    if(predictionRules$include[ruleIndex]){
      randomDesign <- object$W[, termColumns, drop=FALSE]
      if(predictionRules$average[ruleIndex]){
        if(any(randomDesign@x != 1)){
          stop("include=TRUE with average=TRUE requires an indicator random design for term '",
               predictionRules$term[ruleIndex],
               "'. Use include=TRUE alone for empirical numeric-design averages.",
               call.=FALSE)
        }
        weights <- averages %*% randomDesign
        weights <- as(weights != 0, "dMatrix")
        presentCount <- Matrix::rowSums(weights)
        if(any(presentCount == 0)){
          stop("No selected random levels co-occur with some observed prediction rows for term '",
               predictionRules$term[ruleIndex], "'.", call.=FALSE)
        }
        weights <- Matrix::Diagonal(x=1 / presentCount) %*% weights
      }else{
        weights <- averages %*% randomDesign
      }
      contrasts[, termColumns] <- weights
    }else{
      contrasts[, termColumns] <- 1 / length(termColumns)
    }
  }
  dimnames(contrasts) <- list(predictionLabels[selectedGroups], rownames(object$bu))
  list(D=contrasts, observations=groupSizes[selectedGroups],
       classify=classifyData[firstRows[selectedGroups], , drop=FALSE],
       Dtable=predictionRules, levels=levels,
       randomTerms=predictionRules$term[activeRandom])
}

"predict.mmes" <- function(object, Dtable=NULL, D, levels=NULL, sed=FALSE,
                           pairwise=FALSE, adjust="none", df=Inf,
                           PEV=TRUE, VarU=TRUE, observed=FALSE,
                           maxMemoryMB=NULL, ...){
  .mmes_require_coefficient_matrix(object)
  if(length(PEV) != 1L || is.na(PEV) || !is.logical(PEV)){
    stop("PEV must be TRUE or FALSE.", call.=FALSE)
  }
  if(length(VarU) != 1L || is.na(VarU) || !is.logical(VarU)){
    stop("VarU must be TRUE or FALSE.", call.=FALSE)
  }
  if(length(observed) != 1L || is.na(observed) || !is.logical(observed)){
    stop("observed must be TRUE or FALSE.", call.=FALSE)
  }
  if(!is.null(Dtable) && "levels" %in% names(Dtable)){
    stop("Dtable no longer accepts a levels column; supply selections and numeric settings through the levels argument.",
         call.=FALSE)
  }
  if(!is.null(levels) &&
     (!is.list(levels) || (length(levels) &&
      (is.null(names(levels)) || anyNA(names(levels)) ||
       any(names(levels) == "") || anyDuplicated(names(levels)))))){
    stop("levels must be a named list with unique nonmissing term names.", call.=FALSE)
  }
  resolvedLevels <- levels
  if(is.null(maxMemoryMB)){
    maxMemoryMB <- object$control$execution$denseMemoryMB
    if(is.null(maxMemoryMB)) maxMemoryMB <- 1000
  }
  maxMemoryMB <- .mme_positive(maxMemoryMB, "maxMemoryMB")
  if(inherits(object, "mmes.glmm") && VarU){
    stop("VarU is currently available only for Gaussian mmes fits.", call.=FALSE)
  }
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
  observedInfo <- NULL
  if(observed){
    if(is.null(object$fixedRank) || !is.character(D)){
      stop("observed=TRUE requires a rank-reduced fit and a categorical fixed classify term. ",
           "For explicit D matrices use observed=FALSE.", call.=FALSE)
    }
    observedInfo <- .mmes_observed_predictions(object, D, Dtable, levels,
      VarU, sed, pairwise, maxMemoryMB)
    D <- observedInfo$D
    Dtable <- observedInfo$Dtable
    resolvedLevels <- observedInfo$levels
  }else if(!is.null(object$fixedRank) && is.character(D)){
    stop("Classify prediction for a rank-reduced fixed model is disabled by default ",
         "before constructing the Cartesian level grid. Zero-constrained coefficients ",
         "do not make unobserved combinations estimable. Use observed=TRUE for ",
         "observation-weighted fixed predictions over combinations present in the retained ",
         "fitting data, optionally selecting random contributions with Dtable, ",
         "or supply an explicit estimable D matrix. ",
         "Dense covariance allocation is checked against maxMemoryMB.",
         call.=FALSE)
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
    Dtable$levels <- NULL
    termsInDtable <- apply(data.frame(Dtable$term),1,function(xx){all.names(as.formula(paste0("~",xx)))})
    termsInDtable <- lapply(termsInDtable, function(x){intersect(x,colnames(object$data))})
    termsInDtable <- lapply(termsInDtable, function(x){return(unique(c(x, paste(x, collapse=":"))))})
    # term identified
    termsInDtableN <- unlist(lapply(termsInDtable,length))
    d_vars <- strsplit(D, ":", fixed=TRUE)[[1]]
    pickTerm <- which(vapply(termsInDtable, function(xxx){
      any(xxx == D) || (length(d_vars) > 1 && any(vapply(xxx, function(cand){
        setequal(strsplit(cand, ":", fixed=TRUE)[[1]], d_vars)
      }, logical(1))))
    }, logical(1)))
    if(length(pickTerm) == 0){
      isInDf <- intersect(d_vars, colnames(object$data))
      if(length(isInDf) == length(d_vars)){ # is in data frame but not in model
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
    termLevels <- stats::setNames(vector("list", nrow(Dtable)), Dtable$term)
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
        if(!is.null(levels[[nm]])) termLevels[[nm]] <- levels[[nm]]
      }
    }
    hasLevels <- !vapply(termLevels, is.null, logical(1))
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
      rowLevels <- unique(as.character(unlist(termLevels[classifyRows])))
      newRows <- setdiff(rowLevels, rownames(D))
      # Collect all available levels for classify across data, coefficients, and random effect dimensions
      availableLevels <- c(cd, rownames(D))
      if(length(classifyVars) > 1){
        for(k in seq_along(object$uList)){
          uMat <- object$uList[[k]]
          if(!is.null(colnames(uMat)) && !is.null(rownames(uMat))){
            availableLevels <- c(availableLevels,
              as.vector(outer(colnames(uMat), rownames(uMat), paste, sep=":")),
              as.vector(outer(rownames(uMat), colnames(uMat), paste, sep=":")))
          }
        }
      }
      unknown <- setdiff(newRows, availableLevels)
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
      lv <- termLevels[[iRow]]

      if(isFixed){
        numVars <- vars[vapply(vars, function(x) is.numeric(object$data[[x]]), logical(1))]
        facVars <- vars[vapply(vars, function(x) !is.numeric(object$data[[x]]), logical(1))]
        hasCovariate <- length(numVars) > 0

        if(!hasCovariate){
          full <- sparse.model.matrix(reformulate(term, intercept=FALSE), data=object$data)
          fullLabels <- .dt_strip_labels(colnames(full), vars)
          colLabels <- .dt_strip_labels(cd[w], vars)
          if(!is.null(lv)){
            unknown <- setdiff(as.character(lv), fullLabels)
            if(length(unknown)){
              stop("Unknown levels for term '", term, "': ", paste(unknown, collapse=", "), call.=FALSE)
            }
            if(!Dtable[iRow,"include"] && !Dtable[iRow,"average"]){
              warning("levels given for term '", term, "' are ignored because it is neither included nor averaged.",
                      call.=FALSE)
            }
          }
          selected <- if(is.null(lv)) rep(TRUE, length(w)) else colLabels %in% as.character(lv)
          fullSelected <- if(is.null(lv)) rep(TRUE, ncol(full)) else fullLabels %in% as.character(lv)
          if(Dtable[iRow,"include"]){
            subD <- D[,w,drop=FALSE]
            subD[which(subD > 0, arr.ind = TRUE)] = 1
            if(any(!selected)) subD[, which(!selected)] <- 0
            if(Dtable[iRow,"average"]){
              present <- as.matrix(cooc %*% full[, fullSelected, drop=FALSE]) > 0
              averageN <- rowSums(present)
              averageN[noData] <- sum(fullSelected)
              subD <- Matrix::Diagonal(x=1/pmax(averageN, 1)) %*% subD
            }
            D[,w] <- subD
          }else{
            if(Dtable[iRow,"average"]){
              averageN <- sum(fullSelected)
              D[,w] <- matrix(rep(ifelse(selected, 1/averageN, 0), each=nrow(D)), nrow=nrow(D))
            }else{
              D[,w] <- D[,w] * 0
            }
          }
        }else{
          # Fixed term with covariate(s)
          if(is.null(lv)){
            refValues <- stats::setNames(lapply(numVars, function(v) mean(object$data[[v]], na.rm=TRUE)), numVars)
            termLevels[[iRow]] <- refValues
            lv <- refValues
          }else if(is.numeric(lv)){
            if(any(!is.finite(lv))){
              stop("levels for covariate term '", term, "' must be finite numeric values at which to predict.", call.=FALSE)
            }
            if(length(lv) == 1){
              refValues <- stats::setNames(as.list(rep(lv, length(numVars))), numVars)
              termLevels[[iRow]] <- refValues
              lv <- refValues
            }else if(length(lv) == length(w) && length(facVars) == 0){
              termLevels[[iRow]] <- lv
            }else{
              refValues <- stats::setNames(as.list(rep(lv[1], length(numVars))), numVars)
              termLevels[[iRow]] <- refValues
              lv <- refValues
            }
          }else if(is.list(lv)){
            unknownVars <- setdiff(names(lv), vars)
            if(length(unknownVars) > 0){
              stop("Unknown variables in levels for term '", term, "': ", paste(unknownVars, collapse=", "), call.=FALSE)
            }
            for(nv in intersect(names(lv), numVars)){
              val <- lv[[nv]]
              if(!is.numeric(val) || any(!is.finite(val))){
                stop("levels for covariate term '", term, "' must be finite numeric values at which to predict.", call.=FALSE)
              }
            }
            for(nv in numVars){
              if(is.null(lv[[nv]])){
                lv[[nv]] <- mean(object$data[[nv]], na.rm=TRUE)
              }
            }
            termLevels[[iRow]] <- lv
          }else{
            stop("levels for covariate term '", term, "' must be numeric values or a named list at which to predict.", call.=FALSE)
          }

          if(!Dtable[iRow,"include"] && !Dtable[iRow,"average"]){
            if(!is.null(termLevels[[iRow]])){
              warning("levels given for term '", term, "' are ignored because it is neither included nor averaged.", call.=FALSE)
            }
            D[,w] <- 0
            next
          }

          if(length(facVars) == 0){
            if(is.numeric(lv) && length(lv) == length(w)){
              D[,w] <- matrix(rep(lv, each=nrow(D)), nrow=nrow(D))
            }else{
              prod_val <- prod(vapply(numVars, function(v) as.numeric(lv[[v]]), numeric(1)))
              D[,w] <- prod_val
            }
          }else{
            contrastsArg <- if(!is.null(object$inputArgs$contrasts)) object$inputArgs$contrasts else if(!is.null(object$args$contrasts)) object$args$contrasts else NULL
            fixedFormula <- if(!is.null(object$args$fixed)) object$args$fixed else if(!is.null(object$inputArgs$fixed)) object$inputArgs$fixed else NULL

            allFacVars <- facVars
            levsList <- lapply(allFacVars, function(fv){
              l <- levels(object$data[[fv]])
              if(is.null(l)) unique(as.character(object$data[[fv]])) else l
            })
            names(levsList) <- allFacVars
            gridDf <- expand.grid(levsList, stringsAsFactors=TRUE)
            for(nv in numVars){
              gridDf[[nv]] <- as.numeric(lv[[nv]])
            }
            otherFixedVars <- setdiff(all.vars(stats::delete.response(terms(fixedFormula))), c(names(gridDf), names(object$data)[1]))
            for(ov in otherFixedVars){
              if(ov %in% colnames(object$data)){
                if(is.numeric(object$data[[ov]])){
                  gridDf[[ov]] <- mean(object$data[[ov]], na.rm=TRUE)
                }else{
                  l <- levels(object$data[[ov]])
                  gridDf[[ov]] <- factor(if(is.null(l)) object$data[[ov]][1] else l[1], levels=if(is.null(l)) unique(object$data[[ov]]) else l)
                }
              }
            }
            mmGrid <- model.matrix(stats::delete.response(terms(fixedFormula)), gridDf, contrasts.arg=contrastsArg)

            specifiedFac <- intersect(names(lv), allFacVars)
            if(length(specifiedFac) > 0){
              rowMatches <- rep(TRUE, nrow(gridDf))
              for(sf in specifiedFac){
                rowMatches <- rowMatches & (gridDf[[sf]] == as.character(lv[[sf]]))
              }
              selRow <- which(rowMatches)[1]
              if(!is.na(selRow)){
                colVals <- mmGrid[selRow, w]
                D[,w] <- matrix(rep(colVals, each=nrow(D)), nrow=nrow(D))
              }
            }else{
              if(Dtable[iRow, "include"]){
                gridRowNames <- if(length(allFacVars) == 1) as.character(gridDf[[allFacVars[1]]]) else apply(gridDf[, allFacVars, drop=FALSE], 1, paste, collapse=":")
                for(jRow in seq_len(nrow(D))){
                  matchIdx <- match(rd[jRow], gridRowNames)
                  if(is.na(matchIdx) && length(allFacVars) == 1){
                    rdParts <- strsplit(rd[jRow], ":", fixed=TRUE)[[1]]
                    matchIdx <- match(rdParts[1], gridRowNames)
                    if(is.na(matchIdx) && length(rdParts) >= 2) matchIdx <- match(rdParts[2], gridRowNames)
                  }
                  if(!is.na(matchIdx)){
                    D[jRow, w] <- mmGrid[matchIdx, w]
                  }
                }
              }else if(Dtable[iRow, "average"]){
                expCols <- colMeans(mmGrid[, w, drop=FALSE])
                D[,w] <- matrix(rep(expCols, each=nrow(D)), nrow=nrow(D))
              }else{
                D[,w] <- 0
              }
            }
          }
        }
      }else if(Dtable[iRow, "type"] == "random"){
        # Random effect term
        r <- match(iRow, which(Dtable$type == "random"))
        effects <- object$uList[[r]]
        ranges <- object$partitions[[r]]
        if(ncol(effects) > 1){
          available <- c(cd[w],
            as.vector(outer(colnames(effects), rownames(effects), paste, sep=":")),
            as.vector(outer(rownames(effects), colnames(effects), paste, sep=":")),
            colnames(effects), rownames(effects))
        }else{
          available <- c(cd[w], rownames(effects))
        }
        if(!is.null(lv)){
          unknown <- setdiff(as.character(lv), available)
          if(length(unknown)){
            stop("Unknown levels for term '", term, "': ", paste(unknown, collapse=", "), call.=FALSE)
          }
          if(!Dtable[iRow,"include"] && !Dtable[iRow,"average"]){
            warning("levels given for term '", term, "' are ignored because it is neither included nor averaged.", call.=FALSE)
          }
        }
        if(!Dtable[iRow,"include"] && !Dtable[iRow,"average"]){
          D[,w] <- 0
          next
        }
        if(Dtable[iRow,"include"]){
          if(ncol(effects) > 1){
            # Structured random term
            for(jRow in seq_len(nrow(D))){
              D[jRow, w] <- 0
              rd_parts <- strsplit(rd[jRow], ":", fixed=TRUE)[[1]]
              if(length(rd_parts) >= 2){
                col_env <- match(rd_parts[1], colnames(effects))
                row_name <- match(rd_parts[2], rownames(effects))
                if(is.na(col_env) || is.na(row_name)){
                  col_env <- match(rd_parts[2], colnames(effects))
                  row_name <- match(rd_parts[1], rownames(effects))
                }
                if(!is.na(col_env) && !is.na(row_name)){
                  target_col <- ranges[col_env, 1L] + row_name - 1L
                  D[jRow, target_col] <- 1
                }
              }else if(length(rd_parts) == 1){
                row_name <- match(rd[jRow], rownames(effects))
                col_env <- match(rd[jRow], colnames(effects))
                if(!is.na(row_name)){
                  target_cols <- ranges[, 1L] + row_name - 1L
                  D[jRow, target_cols] <- 1 / length(target_cols)
                }else if(!is.na(col_env)){
                  target_cols <- seq.int(ranges[col_env, 1L], ranges[col_env, 2L])
                  D[jRow, target_cols] <- 1 / length(target_cols)
                }
              }
            }
          }else{
            # Unstructured random term
            for(jRow in seq_len(nrow(D))){
              idx <- match(rd[jRow], rownames(effects))
              if(is.na(idx)){
                rd_parts <- strsplit(rd[jRow], ":", fixed=TRUE)[[1]]
                if(length(rd_parts) >= 2){
                  idx <- match(rd_parts[1], rownames(effects))
                  if(is.na(idx)) idx <- match(rd_parts[2], rownames(effects))
                }
              }
              if(!is.na(idx)){
                D[jRow, w] <- 0
                D[jRow, w[idx]] <- 1
              }
            }
            if(Dtable[iRow,"average"]){
              D[,w] <- 1 / length(w)
            }
          }
        }else{
          # !include, but average
          if(Dtable[iRow,"average"]){
            D[,w] <- 1 / length(w)
          }else{
            D[,w] <- 0
          }
        }
      }
    }
    if(length(interceptColumn) > 0){D[,interceptColumn] = 1}
    resolvedLevels <- Filter(Negate(is.null), termLevels)
  }else{ }# user has provided D as a matrix to do direct multiplication
  ## calculate predictions and standard errors
  .mmes_check_estimable_fixed(object, D)
  nFixed <- nrow(object$b)
  fixedOnly <- nFixed == ncol(D) ||
    !any(D[, seq.int(nFixed + 1L, ncol(D)), drop=FALSE] != 0)
  if(!is.null(object$fixedRank)){
    .mmes_prediction_memory(object, nrow(D), fixedOnly, VarU, sed, pairwise,
                            maxMemoryMB)
  }
  bu <- object$bu
  predicted.value <- D %*% bu
  solveVcov <- .mmes_vcov_solver(object)
  vcov <- solveVcov(D)
  pev <- if(PEV) vcov else NULL
  varU <- NULL
  samplingVcov <- NULL
  if(VarU){
    .mmes_varu_check(object)
    if(fixedOnly){
      varU <- matrix(0, nrow(D), nrow(D))
      samplingVcov <- vcov
    }else{
      randomD <- D
      if(nFixed) randomD[, seq_len(nFixed)] <- 0
      prior <- .mmes_prior_vcov(object, randomD)
      randomError <- solveVcov(randomD)
      varU <- prior - randomError
      varU <- (varU + t(varU))/2
      fixedD <- D
      if(nFixed < ncol(D)) fixedD[, seq.int(nFixed + 1L, ncol(D))] <- 0
      fixedSampling <- solveVcov(fixedD)
      samplingVcov <- fixedSampling + varU
    }
  }
  std.error <- sqrt(diag(vcov))
  predictionIds <- rownames(D)
  if(is.null(predictionIds)) predictionIds <- as.character(seq_len(nrow(D)))
  pvals <- data.frame(id=predictionIds, predicted.value=predicted.value[,1],
                      std.error=std.error, pev=diag(vcov),
                      var.u=if(is.null(varU)) rep(NA_real_, nrow(D)) else diag(varU))
  if(is.character(classify)){colnames(pvals)[1] <- classify}
  out <- list(pvals=pvals,D=D,vcov=vcov,PEV=pev,VarU=varU,
              sampling.vcov=samplingVcov,Dtable=Dtable,levels=resolvedLevels)
  if(!is.null(observedInfo)){
    out$observed <- list(classify=observedInfo$classify,
      observations=observedInfo$observations,
      randomTerms=observedInfo$randomTerms,
      weighting=paste0("observation-weighted fixed predictor within retained observed combinations; ",
        if(length(observedInfo$randomTerms)) "random contributions follow the returned Dtable" else
          "random contributions excluded"))
  }
  if(sed || !isFALSE(pairwise)){
    out <- c(out, .predict_comparisons(pvals, vcov, predictionIds, sed, pairwise,
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
