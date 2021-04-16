GWAS <- function(fixed, random, rcov, data, weights, 
                 iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                 init=NULL, constraints=NULL, method="NR", 
                 getPEV=TRUE,
                 na.method.X="exclude",
                 na.method.Y="exclude",
                 return.param=FALSE, 
                 date.warning=TRUE,
                 verbose=FALSE,
                 stepweight=NULL, emupdate=NULL,
                 M=NULL, gTerm=NULL, n.PC = 0, min.MAF = 0.05, 
                 P3D = TRUE){
  
  if(is.null(gTerm)){
    stop("Please provide the name of the genetic term in the model in the gTerm argument", call. = FALSE)
  }
  
  if(length(which(is.na(M))) > 0){
    stop("Please provide an imputed marker matrix (M).", call. = FALSE)
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  
  
  ## return all parameters for a mixed model
  res <- mmer(fixed=fixed, random=random, rcov=rcov, data=data, weights=weights, 
              iters=iters, tolpar=tolpar, tolparinv=tolparinv, 
              init=init, constraints=constraints, method=method, 
              getPEV=getPEV,
              na.method.X=na.method.X,
              na.method.Y=na.method.Y,
              return.param=TRUE, 
              date.warning=date.warning,
              verbose=verbose)
  # print(str(res))
  
  if(return.param){
    lastmodel <- res
  }else{
    ## fit initial models if P3D
    if(P3D){
      lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, res$X,res$Gx,
                         res$Z,res$K,res$R,res$GES,res$GESI,res$ws,
                         res$iters,res$tolpar,res$tolparinv, res$selected,res$getPEV,res$verbose,
                         TRUE, res$stepweight, res$emupdate)
      Vinv <- lastmodel$Vi
    }else{
      lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, res$X,res$Gx,
                         res$Z,res$K,res$R,res$GES,res$GESI,res$ws,
                         res$iters,res$tolpar,res$tolparinv, res$selected,res$getPEV,FALSE,
                         TRUE, res$stepweight, res$emupdate)
      Vinv <- lastmodel$Vi
    }
    ## get names of random effects
    re_names <- unlist(res$re_names)
    re_names <- gsub("\\(Intercept):","",re_names)
    gTermi <- which(re_names %in% gTerm)
    if(length(gTermi) < 1){
      stop("No match between the Gterm and the random effects present.")
    }
    ## get input matrices
    Y <- scale(res$yvar) 
    Z <- do.call(cbind,res$Z[gTermi])
    
    if (n.PC > 0) {
      Kb <- do.call(adiag1,res$K[gTermi])
      eig.vec <- eigen(Kb)$vectors
      X <- do.call(cbind,res$X)
      zbeig <- Z %*% eig.vec[,1:n.PC]
      X <- make.full(cbind(X,zbeig))
    }else {
      X <- do.call(cbind,res$X)
      X <- make.full(X)
    }
    
    if(nrow(M) != ncol(Z)){
      stop(paste("Marker matrix M needs to have same numbers of rows(",nrow(M),") than columns of the gTerm incidence matrix(",ncol(Z),")."),call. = FALSE)
    }
    
    if(length(which(rownames(M) != colnames(Z)))){
      M <- M[colnames(Z),]
    }
    m <- ncol(M); colnamesM <- colnames(M)
    ######################
    ## scorecalc function
    ######################
    cat(red("Performing GWAS evaluation\n"))
    preScores <- .Call("_sommer_gwasForLoop",PACKAGE = "sommer",
                       M,Y,as.matrix(Z),X,Vinv,min.MAF,TRUE
                       ) # we need a different function for P3D that uses MNR
    # preScores <- gwasForLoop(
    #                    M,Y,as.matrix(Z),X,Vinv,min.MAF
    # ) # we need a different function for P3D that uses MNR
    v2 <- length(Y) - ((ncol(X)+1)*ncol(Y)) # ncol(XZMi)
    scores <- -log10(pbeta(preScores, v2/2, 1/2))
    ########################
    rownames(scores) <- colnames(M)
    colnames(scores) <- colnames(Y)
    
    lastmodel$scores <- scores
    lastmodel$method <- method
    lastmodel$constraints <- res[[8]]
    class(lastmodel)<-c("mmergwas")
  }
  
  return(lastmodel)
}
