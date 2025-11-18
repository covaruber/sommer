

corImputation <- function(wide, Gu=NULL, nearest=10, roundR=FALSE){
  if(is.null(rownames(wide))){stop("Rownames of the input matrix cannot be NULL. Please add them", call. = FALSE)}
  if(is.null(Gu)){
    X <- apply(wide, 2, imputev)
    Gu <- cor(t(X))
  }
  wide2 <- wide
  rowNamesWide <-  rownames(wide)
  # for each feature
  for(iEnv in 1:ncol(wide)){ # iEnv=10
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    toPredict <- 1:nrow(wide)
    # if(length(toPredict) > 1){
    likelihood=Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData], drop=FALSE]
    # }else{
    #   likelihood=matrix(Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData]], nrow=1)
    #   rownames(likelihood) <- as.character(rowNamesWide)[toPredict]
    #   colnames(likelihood) <- as.character(rowNamesWide)[withData]
    # }
    replacement <- numeric()
    for(iInd in 1:nrow(likelihood)){ # iInd=1
      # averaging only the positively correlated to avoid decrease   # wide[iInd, iEnv]
      indLik <- sort(abs(likelihood[iInd,]), decreasing = TRUE)
      toAverage <- indLik[1:min(c(nearest,length(withData)))]
      indLikToAverage <- likelihood[iInd,names(toAverage)]
      replacement[iInd] <- mean(wide[names(which(indLikToAverage > 0)),iEnv]) 
    }
    names(replacement) <- rownames(likelihood)
    # time to replace the missing data
    dd=data.frame(replacement=replacement, index=1:length(replacement),
                  imputed=1, id=rownames(likelihood),orVal=wide[,iEnv])
    dd$imputed[which(dd$id %in% names(withoutData))]=0
    head(dd)
    if(roundR){
      wide2[toPredict,iEnv] <- round(replacement)
    }else{
      wide2[toPredict,iEnv] <- replacement
    }
  }
  # for each individual
  for(jRow in 1:nrow(wide)){
    miss <- which(is.na(wide[jRow,]))
    if(length(miss) > 0){
      dd=data.frame(full=as.vector(unlist(wide2[jRow,])), partial=as.vector(unlist(wide[jRow,])))
      # model <- fastLm(partial~full,data=dd[which(!is.na(dd$partial)),])
      # model <- lm(partial~full,data=dd[which(!is.na(dd$partial)),])
      model <- mmes(partial~full,data=dd[which(!is.na(dd$partial)),], verbose = FALSE)
      pp=as.vector(model$b[1,1])+(dd[which(is.na(dd$partial)),"full"]*as.vector(model$b[2,1]))
      if(roundR){
        wide[jRow,miss] <- round(pp)#round(predict(model,newdata = dd[which(is.na(dd$partial)),]))
      }else{
        wide[jRow,miss] <- pp#predict(model,newdata = dd[which(is.na(dd$partial)),])
      }
    }
  }
  stillEmpty <- which(is.na(wide), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide[stillEmpty] <- mean(wide, na.rm=TRUE)}
  stillEmpty <- which(is.na(wide2), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide2[stillEmpty] <- mean(wide2, na.rm=TRUE)}
  
  return(list(imputed=wide, corImputed=wide2))
}


r2 <- function(object, object2=NULL){
  if(!inherits(object, "mmes")){
    stop("This function is only available for models fitted with the mmes() function.", call. = FALSE)
  }
  result <- list()
  for(iPart in 1:length(object$uPevList)){
    pev <- object$uPevList[[iPart]]
    variance0 <- as.numeric(diag(object$theta[[iPart]]))
    variance <- apply(data.frame(variance0),1,function(x){rep(x,nrow(pev))})
    if(!is.null(object2)){
      if(!is.null(object2$Ai)){
        for(iVar in 1:length(variance0)){ # if user provided A matrices we subsitute with more accurate values
          variance[,iVar] <- variance0[iVar] * diag(solve(object2$Ai[[iPart]]))
        }
      }else{
        stop("object2 needs to be a model that sets the argument 'returnParam=TRUE' so we can extract the relationship matrices. Please correct. ", call. = FALSE)
      }
      
    }
    result[[iPart]] <- (variance - pev)/variance
  }
  names(result) <- names(object$uList)
  return(result)
}




matrix.trace <- function(x){
  if (!is.square.matrix(x))
    stop("argument x is not a square matrix")
  return(sum(diag(x)))
}

is.diagonal.matrix <- function (x, tol = 1e-08){
  y <- x
  diag(y) <- rep(0, nrow(y))
  return(all(abs(y) < tol))
}

is.square.matrix <-function(x){
  return(nrow(x) == ncol(x))
}

hadamard.prod <-function (x, y){
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (!is.numeric(y)) {
    stop("argument y is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (is.matrix(y)) {
    Ymat <- y
  }
  else {
    if (is.vector(y)) {
      Ymat <- matrix(y, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (nrow(Xmat) != nrow(Ymat))
    stop("argumentx x and y do not have the same row order")
  if (ncol(Xmat) != ncol(Ymat))
    stop("arguments x and y do not have the same column order")
  return(Xmat * Ymat)
}




