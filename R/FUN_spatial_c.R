spl2Dc <-  function(x.coord,y.coord,at.var=NULL,at.levels=NULL, type="PSANOVA", 
                    nsegments = c(10,10), penaltyord = c(2,2), degree = c(3,3), 
                    nestorder = c(1,1), thetaC=NULL, theta=NULL ) {
  
  ##
  if(length(degree) == 1){degree <- rep(degree,2) }
  if(length(nsegments) == 1){nsegments <- rep(nsegments,2) }
  if(length(penaltyord) == 1){penaltyord <- rep(penaltyord,2) }
  if(length(nestorder) == 1){nestorder <- rep(nestorder,2) }
  
  x.coord.name <- as.character(substitute(list(x.coord)))[-1L]
  y.coord.name <- as.character(substitute(list(y.coord)))[-1L]
  
  # x.coord.name <- "col"
  # y.coord.name <- "row"
  
  if(!is.numeric(x.coord)){
    stop("x.coord argument in spl2D() needs to be numeric.", call. = FALSE)
  }
  if(!is.numeric(y.coord)){
    stop("y.coord argument in spl2D() needs to be numeric.", call. = FALSE)
  }
  
  interpret.covarrubias.formula <-
    function(formula) {
      env <- environment(formula) 
      if(inherits(formula, "character"))          
        formula <- as.formula(formula)
      tf <- terms.formula(formula, specials = c("SAP", "PSANOVA"))
      terms <- attr(tf, "term.labels")
      nt <- length(terms)
      if(nt != 1)
        stop("Error in the specification of the spatial effect: only a sigle bidimensional function is allowed")
      
      res <- eval(parse(text = terms[1]), envir = env)
      res
    }
  
  bbase <- function(X., XL., XR., NDX., BDEG.) {
    # Function for B-spline basis
    dx <- (XR. - XL.)/NDX.
    knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
    P <- outer(X., knots, tpower, BDEG.)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = BDEG. + 1) / (gamma(BDEG. + 1) * dx ^ BDEG.)
    B <- (-1) ^ (BDEG. + 1) * P %*% t(D)
    res <- list(B = B, knots = knots)
    res 
  }
  
  tpower <-
    function(x, t, p) {
      # Function for truncated p-th power function
      return((x - t) ^ p * (x > t))
    }
  
  Rten2 <-
    function(X1,X2) {
      one.1 <- matrix(1,1,ncol(X1))
      one.2 <- matrix(1,1,ncol(X2))
      kronecker(X1,one.2)*kronecker(one.1,X2)
    }
  
  MM.basis <-
    function (x, xl, xr, ndx, bdeg, penaltyord, decom = 1) {
      Bb = bbase(x,xl,xr,ndx,bdeg)
      knots <- Bb$knots
      B = Bb$B
      m = ncol(B)
      n = nrow(B)
      D = diff(diag(m), differences=penaltyord)
      P.svd = svd(crossprod(D))
      U.Z = (P.svd$u)[,1:(m-penaltyord)] # eigenvectors
      d = (P.svd$d)[1:(m-penaltyord)]  # eigenvalues
      Z = B%*%U.Z
      U.X = NULL
      if(decom == 1) {
        U.X = ((P.svd$u)[,-(1:(m-penaltyord))])
        X = B%*%U.X
      } else if (decom == 2){
        X = NULL
        for(i in 0:(penaltyord-1)){
          X = cbind(X,x^i)
        }
      } else if(decom == 3) {
        U.X = NULL
        for(i in 0:(penaltyord-1)){
          U.X = cbind(U.X,knots[-c((1:penaltyord),(length(knots)- penaltyord + 1):length(knots))]^i)
        }
        X = B%*%U.X
      } else if(decom == 4) { # Wood's 2013
        X = B%*%((P.svd$u)[,-(1:(m-penaltyord))])
        id.v <- rep(1, nrow(X))
        D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
        Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
        X <- X%*%Xf
        U.X = ((P.svd$u)[,-(1:(m-penaltyord)), drop = FALSE])%*%Xf
      }
      list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
    }
  
  ####################
  ### if we want to use at.var and at.levels
  if(is.null(at.var)){
    at.var <- rep("A",length(x.coord))
    at.name <- "FIELDINST"
    at.levels <- "A"
  }else{
    at.name <- as.character(substitute(list(at)))[-1L]
    if(length(at.var) != length(x.coord)){stop("at.var has different length than x.coord and y.coord, please fix.", call. = FALSE)}
    if(is.null(at.levels)){at.levels <- levels(as.factor(at.var))}
  }
  #######################
  index <- 1:length(x.coord)
  if(!is.numeric(x.coord)){stop("x.coord argument in spl2D() needs to be numeric.", call. = FALSE)}
  if(!is.numeric(y.coord)){stop("y.coord argument in spl2D() needs to be numeric.", call. = FALSE)}
  #######################
  ## split data by the "at.name" argument
  dat <- data.frame(x.coord, y.coord, at.var, index); colnames(dat) <- c(x.coord.name,y.coord.name,at.name,"index")
  missby <- which(is.na(dat[,at.name]))
  if(length(missby)>0){stop("We will split using the at.name argument and you have missing values in this column.\nPlease correct.", call. = FALSE)}
  dat[,at.name] <- as.factor(dat[,at.name])
  data0 <- split(dat, dat[,at.name])
  names(data0) <- levels(dat[,at.name])
  #######################################
  # make sure there's no missing data in coordinate variables
  nasx <- which(is.na(dat[,x.coord.name]))
  nasy <- which(is.na(dat[,y.coord.name]))
  if(length(nasx) > 0 | length(nasy) >0){
    stop("x.coord and y.coord columns cannot have NA's", call. = FALSE)
  }
  #######################################
  ## now apply the same to all environments
  multires <- lapply(data0, function(dxy){
    
    
    x1 <- dxy[ ,x.coord.name]
    x2 <- dxy[ ,y.coord.name]
    # print(x2)
    #type = type
    
    MM1 = MM.basis(x1, min(x1), max(x1), nsegments[1], degree[1], penaltyord[1], 4)
    MM2 = MM.basis(x2, min(x2), max(x2), nsegments[2], degree[2], penaltyord[2], 4)
    
    X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
    X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B
    
    c1 = ncol(B1); c2 = ncol(B2)
    
    # Nested bases
    if(nestorder[1] == 1) {
      MM1n <- MM1
      Z1n <- Z1
      c1n <- c1
      d1n <- d1	
    } else {
      MM1n = MM.basis(x1, min(x1), max(x1), nsegments[1]/nestorder[1], degree[1], penaltyord[1], 4)
      Z1n <- MM1n$Z
      d1n <- MM1n$d
      c1n <-  ncol(MM1n$B)  					
    }
    if(nestorder[2] == 1) {
      MM2n <- MM2
      Z2n <- Z2
      c2n <- c2
      d2n <- d2	
    } else {
      MM2n = MM.basis(x2, min(x2), max(x2), nsegments[2]/nestorder[2], degree[2], penaltyord[2], 4)
      Z2n <- MM2n$Z
      d2n <- MM2n$d
      c2n <-  ncol(MM2n$B)  					
    }
    
    x.fixed <- y.fixed <- ""
    for(i in 0:(penaltyord[1]-1)){
      if(i == 1) 
        x.fixed <- c(x.fixed, x.coord.name)
      else if( i > 1)
        x.fixed <- c(x.fixed, paste(x.coord.name, "^", i, sep = ""))
    }
    for(i in 0:(penaltyord[2]-1)){
      if(i == 1) 
        y.fixed <- c(y.fixed, y.coord.name)
      else if( i > 1)
        y.fixed <- c(y.fixed, paste(y.coord.name, "^", i, sep = ""))
    }
    xy.fixed <- NULL
    for(i in 1:length(y.fixed)) {
      xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
    }
    xy.fixed <- xy.fixed[xy.fixed != ""]
    names.fixed <- xy.fixed
    
    smooth.comp <- paste("f.", x.coord.name,".", y.coord.name,"", sep = "")
    
    if(type == "SAP") {
      names.random <- paste(smooth.comp, c(x.coord.name, y.coord.name), sep = "|")				
      X = Rten2(X2, X1)		
      # Delete the intercept
      X <- X[,-1,drop = FALSE]
      Z = cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2n, Z1n))
      
      dim.random <- c((c1 -penaltyord[1])*penaltyord[2] , (c2 - penaltyord[2])*penaltyord[1], (c1n - penaltyord[1])*(c2n - penaltyord[2]))		
      dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
      names(dim$fixed) <- names.fixed
      names(dim$random) <- paste(smooth.comp, "Global")
      
      # Variance/Covariance components
      g1u <- rep(1, penaltyord[2])%x%d1
      g2u <- d2%x%rep(1, penaltyord[1])
      g1b <- rep(1, c2n - penaltyord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - penaltyord[1])
      
      g <- list()	
      g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
      g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)
      
      names(g) <- names.random
      
    } else {		
      one1. <- X1[,1, drop = FALSE]
      one2. <- X2[,1, drop = FALSE]
      
      x1. <- X1[,-1, drop = FALSE]
      x2. <- X2[,-1, drop = FALSE]
      
      # Fixed and random matrices
      X <- Rten2(X2, X1)
      # Delete the intercept
      X <- X[,-1,drop = FALSE]
      Z <- cbind(Rten2(one2., Z1), Rten2(Z2, one1.), Rten2(x2., Z1), Rten2(Z2, x1.), Rten2(Z2n, Z1n))
      
      dim.random <- c((c1-penaltyord[1]), (c2-penaltyord[2]), (c1-penaltyord[1])*(penaltyord[2]-1), (c2-penaltyord[2])*(penaltyord[1]-1), (c1n-penaltyord[2])*(c2n-penaltyord[2]))
      
      # Variance/Covariance components		
      g1u <- d1
      g2u <- d2
      
      g1v <- rep(1, penaltyord[2] - 1)%x%d1
      g2v <- d2%x%rep(1,penaltyord[1] - 1)
      
      g1b <- rep(1, c2n - penaltyord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - penaltyord[1])
      
      g <- list()
      
      if(type == "SAP.ANOVA") {
        g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
        g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
        g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, dim.random[4]), g1b)
        g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, g2b)
        
        names.random <- c(paste("f.", x.coord.name,"", sep = ""), paste("f.", y.coord.name,"", sep = ""), paste(smooth.comp, c(x.coord.name, y.coord.name), sep = "."))			
        dim <- list(fixed = rep(1, ncol(X)), random = c(dim.random[1:2], sum(dim.random[-(1:2)])))		
        names(dim$fixed) <- names.fixed
        names(dim$random) <- c(names.random[1:2], paste(smooth.comp, "Global"))
        names(g) <- names.random
      } else {
        g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
        g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
        g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, sum(dim.random[4:5])))
        g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, rep(0, dim.random[5]))
        g[[5]] <- c(rep(0, sum(dim.random[1:4])), g1b + g2b)
        
        names.random <- c(paste("f.", x.coord.name,"", sep = ""), paste("f.", y.coord.name,"", sep = ""),
                          paste("f.", x.coord.name,".", y.coord.name, sep = ""),
                          paste(x.coord.name,".f.", y.coord.name,"", sep = ""),
                          paste("f.", x.coord.name,".f.", y.coord.name,"", sep = ""))
        # print(names.random)
        
        dim <- list(fixed = rep(1, ncol(X)), random = dim.random)		
        names(dim$fixed) <- names.fixed
        names(dim$random) <- names.random
        names(g) <- names.random
      }		
    }
    colnames(X) <- names.fixed
    colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
    
    attr(dim$fixed, "random") <- attr(dim$fixed, "sparse") <- rep(FALSE, length(dim$fixed))
    attr(dim$fixed, "spatial") <- rep(TRUE, length(dim$fixed))
    
    attr(dim$random, "random") <- attr(dim$random, "spatial") <- rep(TRUE, length(dim$random)) 
    attr(dim$random, "sparse") <- rep(FALSE, length(dim$random))
    
    terms <- list()
    terms$MM <- list(MM1 = MM1, MM2 = MM2)
    terms$MMn <- list(MM1 = MM1n, MM2 = MM2n)
    #terms$terms.formula <- res
    
    # attr(terms, "term") <- smooth.comp
    
    # Initialize variance components
    init.var <- rep(1, length(g))
    
    res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var)	
    M <- list(cbind(res$X,res$Z))
    
    return(M)
  })
  
  names(multires) <- at.levels
  # print(str(multires))
  # print(lapply(multires,colnames))
  toFill <- lapply(data0,function(x){x$index})
  #############################################
  ## CAPTURE COLNAMES AND BUILD A MATRIX WITH THOSE NAMES
  myNames <- list()
  for(k in 1:1){
    # print(lapply(multires,function(x){colnames(x[[k]])}))
    myNames[[k]] <- lapply(multires,function(x){colnames(x[[k]])})
  };  names(myNames) <- c("all")
  # print(myNames)
  #################################
  ## move matrices to the right size
  Zup <- list() # store incidence matrices
  Zup2 <- list() # store incidence matrices
  Kup <- list() # store relationship matrices between levels in Z
  typevc <- numeric() # store wheter is a variance (1) or covariance (2;allowed to be negative) component
  re_name <- character() # store the name of the random effect
  counter <- 1
  counter2 <- 1
  for(k in 1:1){ # for each tensor product (single matrix in this case)
    for(j in 1:length(multires)){ # for each environment or by.level
      if(names(multires)[j] %in% at.levels){
        # print("yes")
        Z <- matrix(0,nrow=nrow(dat),ncol=length(myNames[[k]][[j]]))
        colnames(Z) <- myNames[[k]][[j]]
        # print(dim(Z))
        
        prov <- multires[[j]][[k]]
        # print(dim(prov))
        Z[toFill[[j]],colnames(prov)] <- as.matrix(prov) 
        attr(Z,"variables") <- c(x.coord.name, y.coord.name)
        Zup[[counter]] <- Z
        names(Zup)[counter] <- paste0(names(multires)[j],":",names(myNames)[k])
        Gu1 <- diag(ncol(Z)); colnames(Gu1) <- rownames(Gu1) <- colnames(Z)
        Kup[[counter]] <- Gu1
        typevc[counter] <- 1
        re_name[counter] <- names(Zup)[counter]
        counter <- counter + 1
      } # else don't fill that portion of the matrix
    }
  }
  Gti=NULL
  Gtc=NULL
  vcs <- diag(length(Zup)); rownames(vcs) <- colnames(vcs) <- names(Zup)
  #################################
  namess2 <- c(x.coord.name, y.coord.name)
  if(is.null(theta)){
    theta <- diag(length(Zup))*.15
    colnames(theta) <- rownames(theta) <- names(Zup)
  }
  if(is.null(thetaC)){
    thetaC <- diag(length(Zup))
    colnames(thetaC) <- rownames(thetaC) <- names(Zup)
  }
  thetaF <- diag(length(Zup))
  partitionsR <- list() ## only meaningful for residuals
  Zup <- lapply(Zup,function(x){as(x,Class="sparseMatrix")})
  Kup <- as(Kup[[1]],Class="sparseMatrix")
  S3 <- list(Z=Zup,Gu=Kup,theta=theta,thetaC=thetaC,thetaF=thetaF,partitionsR=partitionsR)
  return(S3)
}