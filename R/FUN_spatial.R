spl2D <-  function(x.coord,y.coord,at,at.levels, type="PSANOVA", nseg = c(10,10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1) ) {
  
  ##
  init0 <- as.character(substitute(list(x.coord)))[-1L]
  init1 <- as.character(substitute(list(y.coord)))[-1L]
  
  x.coord.name <- x.coord
  y.coord.name <- y.coord
  
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
  
  bbase <-
    function(X., XL., XR., NDX., BDEG.) {
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
    function (x, xl, xr, ndx, bdeg, pord, decom = 1) {
      Bb = bbase(x,xl,xr,ndx,bdeg)
      knots <- Bb$knots
      B = Bb$B
      m = ncol(B)
      n = nrow(B)
      D = diff(diag(m), differences=pord)
      P.svd = svd(crossprod(D))
      U.Z = (P.svd$u)[,1:(m-pord)] # eigenvectors
      d = (P.svd$d)[1:(m-pord)]  # eigenvalues
      Z = B%*%U.Z
      U.X = NULL
      if(decom == 1) {
        U.X = ((P.svd$u)[,-(1:(m-pord))])
        X = B%*%U.X
      } else if (decom == 2){
        X = NULL
        for(i in 0:(pord-1)){
          X = cbind(X,x^i)
        }
      } else if(decom == 3) {
        U.X = NULL
        for(i in 0:(pord-1)){
          U.X = cbind(U.X,knots[-c((1:pord),(length(knots)- pord + 1):length(knots))]^i)
        }
        X = B%*%U.X
      } else if(decom == 4) { # Wood's 2013
        X = B%*%((P.svd$u)[,-(1:(m-pord))])
        id.v <- rep(1, nrow(X))
        D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
        Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
        X <- X%*%Xf
        U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
      }
      list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
    }
  
  ####################
  ### if we want to use at and at.levels
  if(!missing(at)){
    col1 <- deparse(substitute(at))
    # print(col1)
    dat <- data.frame(x.coord, y.coord, at); colnames(dat) <- c("x.coord","y.coord",col1)
    by <- col1
    if(!missing(at.levels)){
      by.levels=at.levels
    }else{by.levels=NULL}
  }else{
    by=NULL
    by.levels=NULL
    dat <- data.frame(x.coord,y.coord); colnames(dat) <- c("x.coord","y.coord")
  }
  #######################
  
  x.coord <- "x.coord"
  y.coord <- "y.coord"
  
  if(is.null(by)){
    dat$FIELDINST <- "FIELD1"
    by="FIELDINST"
    dat[,by] <- as.factor(dat[,by])
    data0 <- split(dat, dat[,by])
    if(!is.null(by.levels)){
      keep <- names(data0)[which(names(data0) %in% by.levels)]
      
      if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
      
      data0 <- data0[[keep]]
      if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
    }
  }else{
    check <- which(colnames(dat)==by)
    if(length(check)==0){stop("by argument not found in the dat provided", call. = FALSE)}else{
      
      missby <- which(is.na(dat[,by]))
      if(length(missby)>0){stop("We will split using the by argument and you have missing values in this column.\nPlease correct.", call. = FALSE)}
      
      dat[,by] <- as.factor(dat[,by])
      data0 <- split(dat, dat[,by])
      
      if(!is.null(by.levels)){
        keep <- names(data0)[which(names(data0) %in% by.levels)]
        
        if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
        
        data0 <- data0[[keep]]
        if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
        #print(str(data0))
      }
      
    }
  }
  
  nasx <- which(is.na(dat[,x.coord]))
  nasy <- which(is.na(dat[,y.coord]))
  if(length(nasx) > 0 | length(nasy) >0){
    stop("x.coord and y.coord columns cannot have NA's", call. = FALSE)
  }
  #res <- interpret.covarrubias.formula(formula)
  
  ####
  #### now apply the same to all environments
  multires <- lapply(data0, function(data){
    
    
    x1 <- data[ ,x.coord]
    x2 <- data[ ,y.coord]
    
    #type = type
    
    MM1 = MM.basis(x1, min(x1), max(x1), nseg[1], degree[1], pord[1], 4)
    MM2 = MM.basis(x2, min(x2), max(x2), nseg[2], degree[2], pord[2], 4)
    
    X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
    X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B
    
    c1 = ncol(B1); c2 = ncol(B2)
    
    # Nested bases
    if(nest.div[1] == 1) {
      MM1n <- MM1
      Z1n <- Z1
      c1n <- c1
      d1n <- d1	
    } else {
      MM1n = MM.basis(x1, min(x1), max(x1), nseg[1]/nest.div[1], degree[1], pord[1], 4)
      Z1n <- MM1n$Z
      d1n <- MM1n$d
      c1n <-  ncol(MM1n$B)  					
    }
    if(nest.div[2] == 1) {
      MM2n <- MM2
      Z2n <- Z2
      c2n <- c2
      d2n <- d2	
    } else {
      MM2n = MM.basis(x2, min(x2), max(x2), nseg[2]/nest.div[2], degree[2], pord[2], 4)
      Z2n <- MM2n$Z
      d2n <- MM2n$d
      c2n <-  ncol(MM2n$B)  					
    }
    
    x.fixed <- y.fixed <- ""
    for(i in 0:(pord[1]-1)){
      if(i == 1) 
        x.fixed <- c(x.fixed, x.coord)
      else if( i > 1)
        x.fixed <- c(x.fixed, paste(x.coord, "^", i, sep = ""))
    }
    for(i in 0:(pord[2]-1)){
      if(i == 1) 
        y.fixed <- c(y.fixed, y.coord)
      else if( i > 1)
        y.fixed <- c(y.fixed, paste(y.coord, "^", i, sep = ""))
    }
    xy.fixed <- NULL
    for(i in 1:length(y.fixed)) {
      xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
    }
    xy.fixed <- xy.fixed[xy.fixed != ""]
    names.fixed <- xy.fixed
    
    smooth.comp <- paste("f(", x.coord,",", y.coord,")", sep = "")
    
    if(type == "SAP") {
      names.random <- paste(smooth.comp, c(x.coord, y.coord), sep = "|")				
      X = Rten2(X2, X1)		
      # Delete the intercept
      X <- X[,-1,drop = FALSE]
      Z = cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2n, Z1n))
      
      dim.random <- c((c1 -pord[1])*pord[2] , (c2 - pord[2])*pord[1], (c1n - pord[1])*(c2n - pord[2]))		
      dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
      names(dim$fixed) <- names.fixed
      names(dim$random) <- paste(smooth.comp, "Global")
      
      # Variance/Covariance components
      g1u <- rep(1, pord[2])%x%d1
      g2u <- d2%x%rep(1, pord[1])
      g1b <- rep(1, c2n - pord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - pord[1])
      
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
      
      dim.random <- c((c1-pord[1]), (c2-pord[2]), (c1-pord[1])*(pord[2]-1), (c2-pord[2])*(pord[1]-1), (c1n-pord[2])*(c2n-pord[2]))
      
      # Variance/Covariance components		
      g1u <- d1
      g2u <- d2
      
      g1v <- rep(1, pord[2] - 1)%x%d1
      g2v <- d2%x%rep(1,pord[1] - 1)
      
      g1b <- rep(1, c2n - pord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - pord[1])
      
      g <- list()
      
      if(type == "SAP.ANOVA") {
        g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
        g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
        g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, dim.random[4]), g1b)
        g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, g2b)
        
        names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""), paste(smooth.comp, c(x.coord, y.coord), sep = "|"))			
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
        
        names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""),
                          paste("f(", x.coord,"):", y.coord, sep = ""),
                          paste(x.coord,":f(", y.coord,")", sep = ""),
                          paste("f(", x.coord,"):f(", y.coord,")", sep = ""))
        
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
    M <- cbind(res$X,res$Z)
    
    return(M)
  })
  
  nrowss <- (unlist(lapply(multires,nrow)))
  nranges <- (unlist(lapply(multires,ncol)))
  
  names(multires) <- gsub(" ",".",names(multires))
  names(multires) <- gsub("#",".",names(multires))
  names(multires) <- gsub("-",".",names(multires))
  names(multires) <- gsub("/",".",names(multires))
  names(multires) <- gsub("%",".",names(multires))
  names(multires) <- gsub("\\(",".",names(multires))
  names(multires) <- gsub(")",".",names(multires))
  
  
  st <- 1 # for number of rows
  st2 <- 1 # for number of column
  end2 <- numeric() # for end of number of column
  dataflist <- list()
  #glist <- list()
  for(u in 1:length(multires)){
    prov <- multires[[u]]
    mu <- as.data.frame(matrix(0,nrow = sum(nrowss), ncol = ncol(prov)))
    colnames(mu) <- paste(names(multires)[u],colnames(prov), sep="_")
    
    nam <- paste("at",names(multires)[u],"2Dspl", sep="_")
    end <- as.numeric(unlist(st+(nrowss[u]-1)))
    
    mu[st:end,] <- prov
    dataflist[[nam]] <- as(as.matrix(mu), Class="sparseMatrix")
    st <- end+1
    ## for keeping track of the inits
    # end2 <- as.numeric(unlist(st2+(ncol(prov)-1)))
    # glist[[nam]] <- st2:end2
    # st2 <- end2+1
  }
  
  ## now build the last dataframe and adjust the glist
  # newdatspl <- as.data.frame(do.call(cbind,dataflist))
  # nn <- ncol(dat) # to add to the glist
  # glist <- lapply(glist, function(x){x+nn})
  # newdat <- data.frame(dat,newdatspl)
  
  ## now make the formula
  
  #funny <- paste(paste("grp(",names(dataflist),")",sep=""), collapse=" + ")
  
  ## important
  # newdat: is the neew data frame with original data and splines per location matrices
  # glist: is the argument to provide in group in asreml to indicate where each grouping starts and ends
  # funny: formula to add to your random formula
  dataflist <- lapply(dataflist,as.matrix)
  fin <- Reduce("+",dataflist)
  attr(fin,"variables") <- c(init0, init1)
  # fin <-dataflist#list(newdat=dataflist, funny=funny) # 
  return(fin)
}


nna<- function (pheno, trait = "y", rown = "row", coln = "col", nrows = 1, 
                ncols = 2) 
{
  bases <- c(rown, coln, trait)
  into <- intersect(colnames(pheno), bases)
  if (length(into) < 3) {
    stop("Please provide the arguments 'rown', 'coln' and 'trait'.", 
         call. = FALSE)
  }
  newd <- pheno[, c(rown, coln, trait)]
  rn <- rownames(newd)
  newd <- apply(newd, 2, as.numeric)
  rownames(newd) <- rn
  newd <- data.frame(newd)
  pheno[, c(rown, coln, trait)] <- newd
  pheno$nnx <- NA
  nnx <- apply(newd, 1, function(x) {
    s1 <- as.numeric(x[1]-nrows)
    s2 <- as.numeric(x[1] - 1)
    s3 <- as.numeric(x[1] + 1)
    s4 <- as.numeric(x[1] + nrows)
    ns <- c(s1:s2, s3:s4)
    
    s5 <- as.numeric(x[2] - ncols)
    s6 <- as.numeric(x[2] - 1)
    s7 <- as.numeric(x[2] + 1)
    s8 <- as.numeric(x[2] + ncols)
    we <- c(s5:s6, s7:s8)
    
    
    ns <- ns[which(ns > 0)]
    we <- we[which(we > 0)]
    v1 <- which((newd[, rown] %in% ns) & (newd[, coln] %in% x[2]))
    v2 <- which(newd[, coln] %in% we & newd[, rown] %in% x[1])
    if (length(v1) > 0 & length(v2) > 0) {
      yy <- x[3] - apply(newd[c(v1, v2), ], 2, mean, na.rm = TRUE)[3]
    }  else if (length(v1) > 0 & length(v2) == 0) {
      yy <- x[3] - apply(newd[c(v1), ], 2, mean, na.rm = TRUE)[3]
    }  else if (length(v1) > 0 & length(v2) == 0) {
      yy <- x[3] - apply(newd[c(v2), ], 2, mean, na.rm = TRUE)[3]
    } else{
      yy<-x[3]
    }
    return(yy)
  })
  pheno$nnx <- as.numeric(unlist(nnx))
  return(pheno)
}


fill.design <- function(x,rows="ROW",ranges="RANGE", by, extra){
  
  checo <- which(colnames(x) %in% c(rows,ranges))
  if(length(checo)!=2){
    stop("Please double check the rows and ranges argument that you provided. We did not find such columns.\n")
  }
  
  roro <- x[,which(colnames(x)==rows)]
  raro <- x[,which(colnames(x)==ranges)]
  
  if(!is.numeric(roro) | !is.numeric(raro)){
    stop("Please make sure that the columns; ",rows," and ",ranges," are both numeric.\n",call. = FALSE)
  }
  
  if(!missing(extra)){
    extras=TRUE
  }else{extras=FALSE}
  
  fac <- which(unlist(lapply(x,class))=="factor")
  if(length(fac)>0){
    for(u in 1:length(fac)){
      fac2 <- fac[u]
      x[,fac2] <- as.character(x[,fac2])
    }
  }
  
  
  # if user has multiple field
  if(!missing(by)){ # if multiple fields
    y <- split(x,x[,by])
    xnew <- lapply(y,function(x, extra, extras){ # the data, the additional element, and the T/F
      
      x <- x[order(x[,rows], x[,ranges]), ] # order by rows and columns before starint the whole process
      
      roro <- x[,which(colnames(x)==rows)]
      raro <- x[,which(colnames(x)==ranges)]
      
      bybo <- na.omit(unique(x[,by])) # level of by
      ## rows needed
      ro1 <- seq(min(roro,na.rm=TRUE),max(roro,na.rm=TRUE))
      ## ranges needed
      ra1 <- seq(min(raro,na.rm=TRUE),max(raro,na.rm=TRUE))
      ## order the needed one
      needed <- expand.grid(ro1,ra1); colnames(needed) <- c(rows,ranges)
      needed <- needed[ order(needed[,rows], needed[,ranges]), ]
      head(needed)
      ## order the provided
      x <- x[ order(x[,rows], x[,ranges]), ]
      head(x)
      
      ## create empty frame
      dis <- dim(x)
      dis2 <- dim(needed)
      newf <- data.frame(matrix(NA,dis2[1],dis[2]))
      colnames(newf) <- colnames(x)
      head(newf)
      ## create rownames in both so is easier to merge
      
      #1) store rownames of original dataset
      sto <- rownames(x)
      #2) assign new rownames
      rownames(x) <- paste(x[,rows],x[,ranges],sep=".")
      
      #3) do the same for the new frame
      rownames(needed) <- rownames(newf) <- paste(needed[,rows],needed[,ranges],sep=".")
      
      #4) complete the merge
      newf[rownames(x),] <- x
      newf[,c(rows,ranges)] <- needed
      head(newf)
      rownames(newf) <- NULL
      newf[,by] <- bybo
      
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      if(extras){
        for(u in 1:length(extra)){
          pextra <- extra[u]
          pox <- table(newf[,c(rows,ranges,pextra)]) # all zero should be filled in
          leve <- dim(pox)[3] # levels for the extra factor
          levelnames <- na.omit(unique(newf[,pextra]))
          ################
          for(o in 1:leve){ # for each level of the xtra column figure out which rows and ranges are there
            ## let's see until which row extends
            init1 <- which(apply(as.matrix(pox[,,o]),1,sum) > 0) #end in row direction
            stend1 <- c(init1[1],init1[length(init1)]) # start and end in row direction
            
            init2 <- which(apply(as.matrix(pox[,,o]),2,sum) > 0) #end in row direction
            stend2 <- c(init2[1],init2[length(init2)]) # start and end in range direction
            
            kk <- which(newf[,rows] >= stend1[1] & newf[,rows] <= stend1[2] & newf[,ranges] >= stend2[1] & newf[,ranges] <= stend2[2] )
            newf[kk,pextra] <- levelnames[o]
          }
          ################
        }
      }
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      ### if user wants to fill additional factor columns
      return(newf)
    }, extra=extra, extras=extras)
    
    xnew <- do.call(rbind,xnew)
  }else{
    
    roro <- x[,which(colnames(x)==rows)]
    raro <- x[,which(colnames(x)==ranges)]
    
    cat("Argument 'by' not provided. Single field assumed.\n")
    ro1 <- seq(min(roro,na.rm=TRUE),max(roro,na.rm=TRUE))
    ## ranges needed
    ra1 <- seq(min(raro,na.rm=TRUE),max(raro,na.rm=TRUE))
    ## order the needed one
    needed <- expand.grid(ro1,ra1); colnames(needed) <- c(rows,ranges)
    needed <- needed[ order(needed[,rows], needed[,ranges]), ]
    head(needed)
    ## order the provided
    x <- x[ order(x[,rows], x[,ranges]), ]
    head(x)
    
    ## create empty frame
    dis <- dim(x)
    dis2 <- dim(needed)
    newf <- data.frame(matrix(NA,dis2[1],dis[2]))
    colnames(newf) <- colnames(x)
    head(newf)
    ## create rownames in both so is easier to merge
    
    #1) store rownames of original dataset
    sto <- rownames(x)
    #2) assign new rownames
    rownames(x) <- paste(x[,rows],x[,ranges],sep=".")
    
    #3) do the same for the new frame
    rownames(needed) <- rownames(newf) <- paste(needed[,rows],needed[,ranges],sep=".")
    
    #4) complete the merge
    newf[rownames(x),] <- x
    newf[,c(rows,ranges)] <- needed
    head(newf)
    rownames(newf) <- NULL
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    if(extras){
      for(u in 1:length(extra)){
        pextra <- extra[u]
        pox <- table(newf[,c(rows,ranges,pextra)]) # all zero should be filled in
        leve <- dim(pox)[3] # levels for the extra factor
        levelnames <- na.omit(unique(newf[,pextra]))
        ################
        for(o in 1:leve){ # for each level of the xtra column figure out which rows and ranges are there
          ## let's see until which row extends
          init1 <- which(apply(as.matrix(pox[,,o]),1,sum) > 0) #end in row direction
          stend1 <- c(init1[1],init1[length(init1)]) # start and end in row direction
          
          init2 <- which(apply(as.matrix(pox[,,o]),2,sum) > 0) #end in row direction
          stend2 <- c(init2[1],init2[length(init2)]) # start and end in range direction
          
          kk <- which(newf[,rows] >= stend1[1] & newf[,rows] <= stend1[2] & newf[,ranges] >= stend2[1] & newf[,ranges] <= stend2[2] )
          newf[kk,pextra] <- levelnames[o]
        }
        ################
      }
    }
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    ### if user wants to fill additional factor columns
    
    xnew <- newf
  }
  
  return(xnew)
}


