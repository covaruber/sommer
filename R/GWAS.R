GWAS <- function(fixed, random, rcov, data, weights, W,
                 nIters=20, tolParConvLL = 1e-03, tolParInv = 1e-06,
                 init=NULL, constraints=NULL, method="NR",
                 getPEV=TRUE,
                 naMethodX="exclude",
                 naMethodY="exclude",
                 returnParam=FALSE,
                 dateWarning=TRUE,date.warning=TRUE,
                 verbose=FALSE,
                 stepWeight=NULL, emWeight=NULL,
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

  dwToUse <- c(dateWarning,date.warning)
  v<-which(!dwToUse)
  if(length(v) == 0){v<-1}
  # print("hey")
  ## return all parameters for a mixed model (NOT ACTUAL FITTING)
  res <- mmer(fixed=fixed, random=random, rcov=rcov, data=data, W=W,  # silence weights and W when testing
              nIters=nIters, tolParConvLL=tolParConvLL, tolParInv=tolParInv,
              init=init, constraints=constraints, method=method,
              getPEV=getPEV,
              naMethodX=naMethodX,
              naMethodY=naMethodY,
              returnParam=TRUE,
              dateWarning=dwToUse[v],
              verbose=verbose,
              stepWeight=stepWeight, emWeight=emWeight
              )
  # print("hey")
  # print(str(res))

  if(returnParam){ # if user only wants the initial parameters and objects
    lastmodel <- res
  }else{  # if actual fitting is required
    ## fit initial models if P3D
    if(P3D){ # same VC for all markers
      lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, res$X,res$Gx,
                         res$Z,res$K,res$R,res$GES,res$GESI,res$W, res$isInvW,
                         res$nIters,res$tolParConvLL,res$tolParInv, res$selected,res$getPEV,res$verbose,
                         TRUE, res$stepWeight, res$emWeight)
      Vinv <- lastmodel$Vi
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
      )
      v2 <- length(Y) - ((ncol(X)+1)*ncol(Y)) # ncol(XZMi)
      pvals <- pbeta(preScores, v2/2, 1/2)
      scores <- -log10(pvals)
      ########################
      rownames(scores) <- colnames(M)
      colnames(scores) <- colnames(Y)
      
      lastmodel$pvals <- pvals
      lastmodel$scores <- scores
      lastmodel$shape1 <- v2/2
      lastmodel$shape2 <- 1/2
      lastmodel$method <- method
      lastmodel$constraints <- res[[8]]
      class(lastmodel)<-c("mmergwas")
    }else{ # if different variance for each marker

      ## get names of random effects
      re_names <- unlist(res$re_names)
      re_names <- gsub("\\(Intercept):","",re_names)
      gTermi <- which(re_names %in% gTerm)
      if(length(gTermi) < 1){
        stop("No match between the Gterm and the random effects present.")
      }
      Y <- scale(res$yvar)
      # print(head(Y))
      Z <- do.call(cbind,res$Z[gTermi])
      if(nrow(M) != ncol(Z)){
        stop(paste("Marker matrix M needs to have same numbers of rows(",nrow(M),") than columns of the gTerm incidence matrix(",ncol(Z),")."),call. = FALSE)
      }

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
      #
      preScores <- bs <- list()
      for(iMarker in 1:ncol(M)){
        # print(iMarker)
        mi <- M[,iMarker]
        if(var(mi) > 0){ # if var > 0
          Xmi <- cbind(X,mi)
          lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, list(Xmi),res$Gx,
                             res$Z,res$K,res$R,res$GES,res$GESI,res$W, res$isInvW,
                             res$nIters,res$tolParConvLL,res$tolParInv, res$selected,res$getPEV,FALSE,
                             TRUE, res$stepWeight, res$emWeight)
          b = lastmodel$Beta
          var.b = diag(as.matrix(lastmodel$VarBeta))
          se.b = sqrt(abs(var.b))
          t.val  = b/se.b
          # print(t.val)
          preScores[[iMarker]] = 1 - pnorm(t.val[(nrow(t.val)-ncol(Y)+1):nrow(t.val),])
          # print(b)
          bs[[iMarker]] = b[(nrow(b)-ncol(Y)+1):nrow(b),] # last fixed effect where we put the marker
        }else{ # if var == 0
          preScores[[iMarker]] <- rep(1,ncol(Y))
          bs[[iMarker]] = rep(0,ncol(Y)) # effect
        }
      } # for loop for each marker
      preScores <- do.call(rbind,preScores)
      bs <- do.call(rbind, bs)
      scores <- apply(preScores, 2, function(x){-log10(x)})
      
      # scores <- as.matrix(-log10(preScores))
      
      rownames(scores) <- colnames(M) # marker names
      # print(head(scores));print(colnames(Y))
      colnames(scores) <- colnames(Y) # trait names
     
      lastmodel$effects <- bs
      lastmodel$scores <- scores
      lastmodel$pvals <- preScores
      lastmodel$method <- method
      lastmodel$constraints <- res[[8]]
      class(lastmodel)<-c("mmergwas")

    }

  } # if actual fitting is required

  return(lastmodel)
}
