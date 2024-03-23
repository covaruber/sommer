#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction

"predict.mmer" <- function(object, Dtable=NULL, D, ...){
  if(is.character(D)){classify <- D}else{classify="id"} # save a copy before D is overwriten
  # complete the Dtable withnumber of effects in each term
  xEffectN <- object$xEffectsN
  zEffectsN <- object$zEffectsN
  effectsN <- c(xEffectN,zEffectsN)
  start = end = numeric(); add <- 1
  for(i in 1:length(effectsN)){
    start[i] = add
    end[i] = start[i] + effectsN[i] - 1
    add = end[i] + 1
  }
  # fill the Dt table
  if(is.null(Dtable) & is.character(D)){ # if user didn't provide the Dtable
    Dtable <- object$Dtable # we extract it from the model
    termsInDtable <- apply(data.frame(Dtable$term),1,function(xx){all.vars(as.formula(paste0("~",xx)))})
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
    interceptColumn <- grep("Intercept",object$Beta$Effect )
    Dtable$start <- start
    Dtable$end <- end
    # create model matrices to form D
    P <- sparse.model.matrix(as.formula(paste0("~",D,"-1")), data=object$data)
    colnames(P) <- gsub(D,"",colnames(P))
    tP <- t(P)
    W <- object$W
    D = tP %*% W
    toRemove <- names(object$U)
    for(ii in 1:length(toRemove)){
      names(object$U[[ii]][[1]]) <- gsub(toRemove[ii],"",names(object$U[[ii]][[1]]))
    }
    colnames(D) <- c(as.character(object$Beta$Effect),unlist(lapply(object$U,function(y){names(y[[1]])})))
    rd <- rownames(D)
    cd <- colnames(D)
    for(jRow in 1:nrow(D)){ # for each effect add 1's where missing
      myMatch <- which(cd == rd[jRow])
      if(length(myMatch) > 0){D[jRow,myMatch]=1}
    }
    # apply rules in Dtable
    for(iRow in 1:nrow(Dtable)){
      s <- Dtable[iRow,"start"]; e <- Dtable[iRow,"end"]
      # include/exclude rule
      if(Dtable[iRow,"include"]){ # set to 1
        subD <- D[,s:e,drop=FALSE]
        subD[which(subD > 0, arr.ind = TRUE)] = 1
        D[,s:e] <- subD
        # average rule
        if(Dtable[iRow,"average"]){ # set to 1
          # average the include set
          for(o in 1:nrow(subD)){
            v <- which(subD[o,] > 0);  subD[o,v] <- subD[o,v]/(length(v) + length(interceptColumn))
          }
          D[,s:e] <- subD
        }
      }else{ # set to zero
        if(Dtable[iRow,"average"]){ # set to 1
          subD <- D[,s:e,drop=FALSE] + 1
          subD <- subD/subD
          subD[which(subD > 0, arr.ind = TRUE)] = subD[which(subD > 0, arr.ind = TRUE)]/(ncol(subD) + length(interceptColumn))
          D[,s:e] <- subD
        }else{
          D[,s:e] <- D[,s:e] * 0
        }
      }
    }
    if(length(interceptColumn) > 0){D[,interceptColumn] = 1}
    if(length(which(Dtable$term == "1")) > 0){Dtable[which(Dtable$term == "1"),"include"]=TRUE}
  }else{ }# user has provided D as a matrix to do direct multiplication
  ## calculate predictions and standard errors
  bu <- c(object$Beta$Estimate,unlist(object$U)) #object$bu
  predicted.value <- D %*% bu
  ## standard error
  Vi <- object$Vi
  ViW = Vi %*% W; # ViW (nxp)
  WtViW = t(W) %*% ViW; # W'ViW (pxp)
  # inverse of W'ViW
  WtViWi <- try(
    solve(WtViW),
    silent = TRUE
  )
  if(inherits(WtViWi,"try-error") ){
    WtViW = WtViW + diag(mean(diag(Vi)), nrow(WtViW),nrow(WtViW))
    WtViWi <- try(
      solve(WtViW),
      silent = TRUE
    )
  }
  # get G (not sure if this is the right G) var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
  # H <- do.call(adiag1, lapply(object$VarU, function(x){do.call(adiag1,x)}))
  # H <- adiag1(object$VarBeta*0,H)
  # G <- H

  # Gs <- lapply(as.list(1:length(object$K)), function(x){object$K[[x]]*as.vector(object$sigma[[x]])})
  # Gs2 <- adiag1(object$VarBeta*0,do.call(adiag1,Gs))
  # tWG <- W %*% Gs2
  # G <- t(tWG) %*% object$P %*% tWG

  PEV <- do.call(adiag1, lapply(object$PevU, function(x){do.call(adiag1,x)}))
  PEV <- adiag1(object$VarBeta,PEV)
  G <- PEV # it should ideally be G but we don't extract that

  # build the projection matrix
  P <- object$P # Vi - (Vi%*%X%*%(XtViX)%*%t(X)%*%Vi)
  # project W using P
  WtPW <- t(W) %*% (P) %*% W
  # add G to complete the coefficient matrix
  Ci <- WtPW + G # G because is the inverse of Gi which is what goes into the M matrix
  vcov <- D %*% Ci %*% t(D)
  std.error <- sqrt(diag(vcov))
  pvals <- data.frame(id=rownames(D),predicted.value=predicted.value[,1], std.error=std.error)
  if(is.character(classify)){colnames(pvals)[1] <- classify}
  return(list(pvals=pvals,D=D,vcov=vcov, Dtable=Dtable))
}


"print.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
    The predictions are obtained by averaging/aggregating across
    the hypertable calculated from model terms constructed solely
    from factors in the include sets. You can customize the model
    terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  # print(x$hypertable)
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}
