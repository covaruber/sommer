propMissing <- function(x){
  length(which(is.na(x)))/length(x)
}

corImputation <- function(wide, Gu=NULL, nearest=10, roundR=FALSE){
  if(is.null(rownames(wide))){stop("Rownames of the input matrix cannot be NULL. Please add them", call. = FALSE)}
  if(is.null(Gu)){
    X <- apply(wide, 2, sommer::imputev)
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

logspace <- function (n, start, end) {
  exp(seq(log(start), log(end), length.out = n))
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


wald.test <- function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL,
                       verbose = FALSE) {
  if (is.null(Terms) & is.null(L))
    stop("One of the arguments Terms or L must be used.")
  if (!is.null(Terms) & !is.null(L))
    stop("Only one of the arguments Terms or L must be used.")
  if (is.null(Terms)) {
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
  }
  else w <- length(Terms)
  if (is.null(H0))
    H0 <- rep(0, w)
  if (w != length(H0))
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
  if (is.null(L)) {
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for (i in 1:w) L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),
                            sep = ""), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if (is.null(df))
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else {
    fstat <- stat/nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p), Ftest = c(Fstat = fstat,
                                                                df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
  }
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
                 L = L, result = res, verbose = verbose, df = df), class = "wald.test")
}

print.wald.test <- function(x, digits = 2, ...){
  Terms <- x[["Terms"]]; b <- x[["b"]]; H0 <- x[["H0"]]; v <- x[["result"]][["chi2"]]; df <- x[["df"]]
  verbose <- x[["verbose"]]
  namb <- names(b)[Terms]
  cat("Wald test:\n", "----------\n", sep = "")
  if(verbose){
    cat("\nCoefficients:\n")
    print(format(b, digits = digits), quote = FALSE)
    cat("\nVar-cov matrix of the coefficients:\n")
    print(format(x[["Sigma"]], digits = digits), quote = FALSE)
    cat("\nTest-design matrix:\n")
    print(x[["L"]])
    cat("\nPositions of tested coefficients in the vector of coefficients:", paste(Terms, collapse = ", "), "\n")
    if(is.null(namb))
      cat("\nH0: ", paste(paste(format(b[Terms], digits), format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
    else{
      cat("\nH0: ", paste(paste(namb, format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
    }
    #    cat("\nTest results:\n")
  }
  cat("\nChi-squared test:\n")
  cat("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), ", df = ", v["df"],
      ", P(> X2) = ", format(v["P"], digits = digits, nsmall = 1), "\n", sep = "")
  if(!is.null(df)){
    v <- x[["result"]][["Ftest"]]
    cat("\nF test:\n")
    cat("W = ", format(v["Fstat"], digits = digits, nsmall = 1),
        ", df1 = ", v["df1"],
        ", df2 = ", v["df2"],
        ", P(> W) = ", format(v["P"], digits = digits), "\n", sep = "")
  }
}

leg <- function(x,n=1,u=-1,v=1, intercept=TRUE, intercept1=FALSE){
  
  init0 <- as.character(substitute(list(x)))[-1L]
  
  if(system.file(package = "orthopolynom") == ""){
    stop("Please install the orthopolynom package to use this function.",call. = FALSE)
  }
  requireNamespace("orthopolynom",quietly=TRUE)
  (leg4coef <- orthopolynom::legendre.polynomials(n=n, normalized=TRUE))
  leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                  x=orthopolynom::scaleX(x, u=u, v=v))))
  colnames(leg4) <- paste("leg",0:(ncol(leg4)-1),sep="")
  if(!intercept){
    leg4 <- leg4[, 2:ncol(leg4), drop = FALSE]
  }
  if(intercept1){
    leg4 <- leg4*sqrt(2)
    # leg4[,1] <- leg4[,1]*sqrt(2)
  }
  attr(leg4,"variables") <- c(init0)
  return(leg4)
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

adiag1 <- function (..., pad = as.integer(0), do.dimnames = TRUE){
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(list(args[[1]]), list(jj),
                               list(pad = pad))))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) & is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) & length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) & length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) +
                  dim.a[[i]])
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

imputev <- function(x, method="median"){
  if(is.numeric(x)){
    if(method=="mean"){
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }else if(method=="median"){
      x[which(is.na(x))] <- median(x,na.rm=TRUE)
    }else{
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }
  }else{
    if(method=="mean"){
      stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
    }else if(method=="median"){
      tt <- table(x)
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }else{
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }
  }
  return(x)
}

fdr <- function(p, fdr.level=0.05){
  ##### transform to real p-values
  # if maximum value is grater than 1 means that is in -log10 or LOD scale
  # if maximum value is less than one means that the user is using already raw  p.values
  if(max(p, na.rm = TRUE) > 1){ # is in -lod 10 scale
    pval <- 10^-p
  }else{
    pval <- p
  }
  ########## make sure there is a value
  #ro1 <- c(pval)
  #ro2 <- p.adjust(ro1, method="fdr")
  #ro3 <- -log(c(ro2,0.05), 10)
  #ro4 <- 10^-ro3
  #ro5 <-p.adjust(ro4, method="fdr")
  ##### adjust for FDR ---- ADJUSTED IN P.VAL SCALE -----
  pvalA <- p.adjust(pval, method="fdr")
  #plot(pvalA)
  ##### ---- VALS IN LOG.10 SCALE -----
  pvalA.l10 <- -log(pvalA, 10)
  #plot(pvalA.l10)
  ##### ---- FDR IN P.VAL SCALE FOR ADJUSTED ----
  fdr.p <- fdr.level
  ## FDR in original scale
  #1) transform the values to p-values
  # pvalA
  #2) find which value adjusted is equal to 0.05 and go back to the original value
  sortedd <- sort(pvalA, decreasing = TRUE)
  closer <- sortedd[which(sortedd < fdr.level)[1]] # closest value found to the fdr.level indicated by the user
  vv <- which(pvalA == closer)[1]
  
  #vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0){
    fdr.10 <- p[vv]#fdr.or <- min(p[vv])
    #fdr <- 0.05
  }else{
    fdr.10 <- NULL
  }
  ######
  result <- list(p.ad=pvalA, fdr.p=fdr.p, p.log10=pvalA.l10, fdr.10=fdr.10 )
  return(result)
}

fdr2 <- function(p, fdr.level=0.05){
  
  
  ##### transform to real p-values
  # if maximum value is grater than 1 means that is in -log10 or LOD scale
  # if maximum value is less than one means that the user is using already raw  p.values
  if(max(p, na.rm = TRUE) > 1){ # is in -log 10 scale
    pval <- 10^-p
  }else{
    pval <- p
  }
  
  ##for endelmans function
  pen <- -log(pval,10)
  ########## make sure there is a value
  #ro1 <- c(pval)
  #ro2 <- p.adjust(ro1, method="fdr")
  #ro3 <- -log(c(ro2,0.05), 10)
  #ro4 <- 10^-ro3
  #ro5 <-p.adjust(ro4, method="fdr")
  ##### adjust for FDR ---- ADJUSTED IN P.VAL SCALE -----
  pvalA <- p.adjust(pval, method="fdr")
  #plot(pvalA)
  ##### ---- VALS IN LOG.10 SCALE -----
  pvalA.l10 <- -log(pvalA, 10)
  #plot(pvalA.l10)
  ##### ---- FDR IN P.VAL SCALE FOR ADJUSTED ----
  fdr.p <- fdr.level
  ## FDR in original scale
  #1) transform the values to p-values
  # pvalA
  #2) find which value adjusted is equal to 0.05 and go back to the original value
  sortedd <- sort(pvalA, decreasing = TRUE)
  closer <- sortedd[which(sortedd < fdr.level)[1]] # closest value found to the fdr.level indicated by the user
  vv <- which(pvalA == closer)[1]
  
  #vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0 & !is.na(vv)){
    fdr.10 <- p[vv]#fdr.or <- min(p[vv])
    #fdr <- 0.05
  }else{
    
    fdrendel <- function(dd, fdr.level=0.05){
      qvalue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
          print("ERROR: p-values not in valid range.")
          return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
          pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        #print(pi0)
        if (pi0 <= 0) {
          #pi0 <- abs(pi0)
          print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
          return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
          idx <- sort.list(x)
          fc <- factor(x)
          nl <- length(levels(fc))
          bin <- as.integer(fc)
          tbl <- tabulate(bin)
          cs <- cumsum(tbl)
          tbl <- rep(cs, tbl)
          tbl[idx] <- tbl
          return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
          qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]],
                              1)
        }
        return(qvalue)
      }
      ## go back to p values
      q.ans <- qvalue(10^-dd)
      temp <- cbind(q.ans, dd)
      temp <- temp[order(temp[, 1]), ]
      
      temp2 <- tapply(temp[, 2], temp[, 1], mean)
      qvals <- as.numeric(rownames(temp2))
      x <- which.min(abs(qvals - fdr.level))
      first <- max(1, x - 2)
      last <- min(x + 2, length(qvals))
      if ((last - first) < 4) {
        last <- first + 3
      }
      #print(qvals[first:last])
      #print(temp2[first:last])
      splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last],
                             df = 3)
      popo <- predict(splin, x = fdr.level)$y
      return(popo)
    }
    
    fdr.10 <- fdrendel( pen,fdr.level = fdr.level)
  }
  ######
  result <- list(p.ad=pvalA, fdr.p=fdr.p, p.log10=pvalA.l10, fdr.10=fdr.10 )
  return(result)
}


