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

vpredict <- function (object, transform){
  
  # if(object$method %in% c("EMMA","EM")){
  #   stop("The pin function only works for 'NR' and 'AI' methods.",call. = FALSE)
  # }
  # 
  # if(object$method %in% c("MNR","MAI","MEMMA")){
  #   pframe <- as.list(object$sigma)#as.list(summary(object)[[3]][,1])
  # }else{
  #   pframe <- as.list(object$var.comp[,1])
  # }
  pframe <- as.list(summary(object)$varcomp[,1])
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  ## deriv creates a derivative for a simple expression
  # i.e. dx2x <- deriv(~ x^2, "x") ; dx2x
  ## eval evaluates a derivative expression for certain values specified in the deriv() expresion
  
  ## 1) get derivatives of the expression provided, 
  ## i.e. V1/(V1+V4) with respect to each of the variables existing, i.e. v1, v2, v3, v4
  ## 2) Evaluate the expression (derivatives) for values provided for each variable
  ## i.e. with the actual values of the variance components
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
                 pframe)
  X <- as.vector(attr(tvalue, "gradient")) # just make it a sinpe vector of derivatives
  tname <- if (length(transform) == 3) 
    transform[[2]]
  else ""
  n <- length(pframe) ## number of parameters available, i.e. V1,V2,V3,V4
  i <- rep(1:n, 1:n) ## repeat each parameter by its own
  j <- sequence(1:n) ## makes a sequence from 1 to the number provided, i.e. if sequence(1:2) = 1 1 2, because it makes the sequence for 1:1 and then 1:2
  k <- 1 + (i > j) # all where i <= j get a 1, all i > j get a 2
  Vmat <- object$sigmaSE
  toext <- upper.tri(Vmat)
  diag(toext) <- TRUE
  Vmat <- Vmat[which(toext,arr.ind = TRUE)] ## extract the upper triangular
  se <- sqrt(abs(sum(Vmat * X[i] * X[j] * k))) ## only X[i] and X[j] that match with the Vi's indicated are != than zero
  
  ### Vmat are the second derivatives
  ### X are the derivatives of the expression with respect to each parameter of interest
  ## X[i] * X[j] * k  multiplies once the var(var.i) and var(var.j) and twice the covar(covar.ij)
  ## then takes the sqrt(abs( sum[c(var(i),covar(i,j),var(j)] ))) , that's the SE
  ## Vmat * X[i] * X[j] * k --> 2nd derivatives * derivatives of i.e. h2 with respect to each term
  ## d''(x) * d'(x) * d'
  ## those var(vc.i) and covar(covar.ij) from the variance comp. come from the inverse if the second derivatives (Fisher's)
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
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


