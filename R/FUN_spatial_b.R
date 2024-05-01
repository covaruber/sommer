spl2Db <-  function(x.coord,y.coord, at.var=NULL,at.levels=NULL, nsegments = c(10,10),
                    degree = c(3,3), penaltyord = c(2,2),nestorder = c(1,1), 
                    minbound=NULL, maxbound=NULL, method="Lee", what="bits" ) {
  
  if(length(degree) == 1){degree <- rep(degree,2) }
  if(length(nsegments) == 1){nsegments <- rep(nsegments,2) }
  if(length(penaltyord) == 1){penaltyord <- rep(penaltyord,2) }
  if(length(nestorder) == 1){nestorder <- rep(nestorder,2) }
  
  x.coord.name <- as.character(substitute(list(x.coord)))[-1L]
  y.coord.name <- as.character(substitute(list(y.coord)))[-1L]
  
  if(is.null(at.var)){
    at.var <- rep("A",length(x.coord))
    at.name <- "FIELDINST"
    at.levels <- "A"
  }else{
    at.name <- as.character(substitute(list(at)))[-1L]
    if(length(at.var) != length(x.coord)){stop("at.var has different length than x.coord and y.coord, please fix.", call. = FALSE)}
    if(is.null(at.levels)){at.levels <- levels(as.factor(at.var))}
  }
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
  ##########################################
  #### now calculate TP design matrices for each by.level
  multires <- lapply(data0, function(dxy){ # for each environment
    # use function to extract incidence matrices
    TPXZg <- tpsmmbwrapper(columncoordinates=x.coord.name, rowcoordinates=y.coord.name,
                           maxbound=maxbound, minbound=minbound, penaltyord=penaltyord,
                           data=dxy, nsegments=nsegments, nestorder=nestorder, asreml="grp", method=method)
    
    # extract the incidence matrices
    fC <- TPXZg$data[,TPXZg$grp$TP.R.1_fcol]
    fR <- TPXZg$data[,TPXZg$grp$TP.C.1_frow]
    fC.R <- TPXZg$data[,TPXZg$grp$TP.R.2_fcol]
    C.fR <- TPXZg$data[,TPXZg$grp$TP.C.2_frow]
    fC.fR <- TPXZg$data[,TPXZg$grp$TP_fcol_frow]
    rest <- TPXZg$data[,min(c(which(colnames(TPXZg$data)=="TP.col"),which(colnames(TPXZg$data)=="TP.row"))):(min(c(TPXZg$grp$TP.R.1_fcol,TPXZg$grp$TP.C.1_frow))-1)]
    # all <- TPXZg$data[,TPXZg$grp$All]
    # return output
    return(list(fC,fR,fC.R,C.fR,fC.fR,rest))
  })
  names(multires) <- at.levels
  # print(str(multires))
  toFill <- lapply(data0,function(x){x$index})
  #############################################
  ## CAPTURE COLNAMES AND BUILD A MATRIX WITH THOSE NAMES
  myNames <- list()
  for(k in 1:6){
    myNames[[k]] <- lapply(multires,function(x){colnames(x[[k]])})
  };  names(myNames) <- c("fC","fR","fC.R", "C.fR","fC.fR","rest")
  #################################
  ## move matrices to the right size
  Zup <- list() # store incidence matrices
  Zup2 <- list() # store incidence matrices
  Kup <- list() # store relationship matrices between levels in Z
  typevc <- numeric() # store wheter is a variance (1) or covariance (2;allowed to be negative) component
  re_name <- character() # store the name of the random effect
  counter <- 1
  counter2 <- 1
  for(k in 1:6){ # for each tensor product
    for(j in 1:length(multires)){ # for each environment or by.level
      if(names(multires)[j] %in% at.levels){
        Z <- matrix(0,nrow=nrow(dat),ncol=length(myNames[[k]][[j]]))
        colnames(Z) <- myNames[[k]][[j]]
        prov <- multires[[j]][[k]]
        Z[toFill[[j]],colnames(prov)] <- as.matrix(prov) 
        attr(Z,"variables") <- c(x.coord.name, y.coord.name)
        if(k < 6){
          Zup[[counter]] <- Z
          names(Zup)[counter] <- paste0(names(multires)[j],":",names(myNames)[k])
          Gu1 <- diag(ncol(Z)); colnames(Gu1) <- rownames(Gu1) <- colnames(Z)
          Kup[[counter]] <- Gu1
          typevc[counter] <- 1
          re_name[counter] <- names(Zup)[counter]
          counter <- counter + 1
        }else{
          Zup2[[counter2]] <- Z
          names(Zup2)[counter2] <- paste0(names(multires)[j],":",names(myNames)[k])
          counter2 <- counter2 + 1
        }
      } # else don't fill that portion of the matrix
    }
  }
  Gti=NULL
  Gtc=NULL
  vcs <- diag(length(Zup)); rownames(vcs) <- colnames(vcs) <- names(Zup)
  #################################
  if(what=="bits"){
    namess2 <- c(x.coord.name, y.coord.name)
    S3 <- list(Z=Zup,K=Kup,Gti=Gti,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=vcs, terms=namess2)
  }else if(what == "base"){
    S3 <- do.call(cbind,Zup2)
  }else{stop("method not recognized.",call. = FALSE)}
  
  return(S3)
}

spl2Dmats <-  function(x.coord.name,
                       y.coord.name,
                       data,                       
                       at.name,
                       at.levels,
                       nsegments = NULL,
                       minbound=NULL, 
                       maxbound=NULL,
                       degree = c(3,3), 
                       penaltyord = c(2,2), 
                       nestorder = c(1,1), 
                       method="Lee" ) {
  
  # x.coord.name <- as.character(substitute(list(x.coord)))[-1L]
  # y.coord.name <- as.character(substitute(list(y.coord)))[-1L]
  
  x.coord <- data[,x.coord.name]
  y.coord <- data[,y.coord.name]
  
  data[,paste0(x.coord.name,"f")] <- as.factor(data[,x.coord.name])
  data[,paste0(y.coord.name,"f")] <- as.factor(data[,y.coord.name])
  
  if(!is.numeric(x.coord)){stop("x.coord argument in spl2D() needs to be numeric.", call. = FALSE)}
  if(!is.numeric(y.coord)){stop("y.coord argument in spl2D() needs to be numeric.", call. = FALSE)}
  #######################
  ## split data by the "at.name" argument
  if(missing(at.name)){ # if user doesn't provide the at.name argument
    dat <- data.frame(x.coord, y.coord); colnames(dat) <- c(x.coord.name,y.coord.name)
    dat$FIELDINST <- "FIELD1"
    at.name="FIELDINST"
    dat[,at.name] <- as.factor(dat[,at.name])
    data0 <- split(dat, dat[,at.name])
    at.levels="FIELD1"
    # for the actual dataset
    data[,at.name] <- "FIELD1"
  }else{
    check <- which(colnames(data)==at.name)
    if(length(check)==0){stop("at.name column not found in the data provided", call. = FALSE)}else{
      at <- data[,at.name]
      dat <- data.frame(x.coord, y.coord, at); colnames(dat) <- c(x.coord.name,y.coord.name,at.name)
      missby <- which(is.na(dat[,at.name]))
      if(length(missby)>0){stop("We will split using the at.name argument and you have missing values in this column.\nPlease correct.", call. = FALSE)}
      dat[,at.name] <- as.factor(dat[,at.name])
      data0 <- split(dat, dat[,at.name])
      names(data0) <- levels(dat[,at.name])
    }
  }
  if(missing(at.levels)){
    at.levels <- levels(dat[,at.name])
  }
  #######################################
  # make sure there's no missing data in coordinate variables
  nasx <- which(is.na(dat[,x.coord.name]))
  nasy <- which(is.na(dat[,y.coord.name]))
  if(length(nasx) > 0 | length(nasy) >0){
    stop("x.coord and y.coord columns cannot have NA's", call. = FALSE)
  }
  ##########################################
  #### now calculate TP design matrices for each by.level
  multires <- lapply(data0, function(dxy){
    # use function to extract incidence matrices
    
    TPXZg <- tpsmmbwrapper(columncoordinates=x.coord.name, rowcoordinates=y.coord.name,
                             maxbound=maxbound, minbound=minbound, penaltyord=penaltyord,
                             data=dxy, nsegments=nsegments, nestorder=nestorder, asreml="grp", method=method)

    # extract the incidence matrices
    fC <- TPXZg$data[,TPXZg$grp$TP.R.1_fcol]
    fR <- TPXZg$data[,TPXZg$grp$TP.C.1_frow]
    fC.R <- TPXZg$data[,TPXZg$grp$TP.R.2_fcol]
    C.fR <- TPXZg$data[,TPXZg$grp$TP.C.2_frow]
    fC.fR <- TPXZg$data[,TPXZg$grp$TP_fcol_frow]
    rest <- TPXZg$data[,min(c(which(colnames(TPXZg$data)=="TP.col"),which(colnames(TPXZg$data)=="TP.row"))):(min(c(TPXZg$grp$TP.R.1_fcol,TPXZg$grp$TP.C.1_frow))-1)]
    all <- TPXZg$data[,TPXZg$grp$All]
    # return output
    return(list(fC,fR,fC.R,C.fR,fC.fR,all,rest))
  })
  # print(str(multires))
  nrows <- unlist(lapply(data0,nrow))
  end <- numeric(); for(l in 1:length(nrows)){end[l] <- sum(nrows[1:l]) }
  start <- numeric(); for(l in 1:length(nrows)){start[l] <- end[l]-nrows[l]+1 }
  #############################################
  ## CAPTURE COLNAMES AND BUILD A MATRIX WITH THOSE NAMES
  nColList <- list()
  for(k in 1:6){ # for each environment or by.level
    nColList[[k]] <- unlist(lapply(multires,function(x){colnames(x[[k]])}))
  }
  uniqueNames <- lapply(nColList,unique) # 5 element in a list with names for matrices
  # build the matrices
  Zl <- list()
  for(k in 1:6){ # for each tensor product
    Z <- matrix(0,nrow=nrow(dat),ncol=length(uniqueNames[[k]]))
    colnames(Z) <- uniqueNames[[k]]
    for(j in 1:length(multires)){ # for each environment or by.level
      if(names(multires)[j] %in% at.levels){
        prov <- multires[[j]][[k]]
        Z[start[j]:end[j],colnames(prov)] <- as.matrix(prov) 
      } # else don't fill that portion of the matrix
    }
    attr(Z,"variables") <- c(x.coord.name, y.coord.name)
    Zl[[k]] <- Z
  }
  names(Zl) <- c("fC","fR","fC.R", "C.fR","fC.fR","all")
  
  data[,at.name] <- as.factor(data[,at.name])
  dataToreturn <- split(data, data[at.name])
  dataToreturn <- do.call(rbind,dataToreturn)
  rest <- lapply(multires,function(x){x[[7]]})
  rest <- do.call(rbind,rest)
  dataToreturn <- cbind(dataToreturn,rest)
  Zl$data <- dataToreturn
  return(Zl)
}

tpsmmbwrapper <- function (columncoordinates, rowcoordinates, data, nsegments=NULL, 
                           minbound=NULL, maxbound=NULL, degree = c(3, 3), penaltyord = c(2, 2), 
                           nestorder = c(1, 1), asreml = "mbf", eigenvalues = "include", 
                           method = "Lee", stub = NULL) 
{
  if (missing(columncoordinates)) 
    stop("columncoordinates argument must be set")
  if (missing(rowcoordinates)) 
    stop("rowcoordinates argument must be set")
  if (missing(data)) 
    stop("data argument must be set")
  col <- sort(unique(data[[columncoordinates]]))
  nuc <- length(col)
  col.match <- match(data[[columncoordinates]], col)
  row <- sort(unique(data[[rowcoordinates]]))
  nur <- length(row)
  row.match <- match(data[[rowcoordinates]], row)
  nv <- length(data[[columncoordinates]])
  if (is.null(minbound)) {
    cminval <- min(col)
    rminval <- min(row)
  } else {
    cminval <- min(c(minbound[1], min(col)))
    if (length(minbound) < 2) {
      rminval <- min(c(minbound[1], min(row)))
    }
    else {
      rminval <- min(c(minbound[2], min(row)))
    }
  }
  if (is.null(maxbound)) {
    cmaxval <- max(col)
    rmaxval <- max(row)
  }
  else {
    cmaxval <- max(c(maxbound[1], max(col)))
    if (length(maxbound) < 2) {
      rmaxval <- max(c(maxbound[1], max(row)))
    }
    else {
      rmaxval <- max(c(maxbound[2], max(row)))
    }
  }
  if (is.null(nsegments)) {
    nsegcol <- nuc - 1
    nsegrow <- nur - 1
  }
  else {
    nsegcol <- max(c(nsegments[1], 2))
  }
  if (length(nsegments) < 2) {
    nsegrow <- max(c(nsegments[1], 2))
  }
  else {
    nsegrow <- max(c(nsegments[2], 2))
  }
  nestcol <- floor(nestorder[1])
  if (length(nestorder) < 2) 
    nestrow <- floor(nestorder[1])
  else nestrow <- floor(nestorder[2])
  nsncol <- 0
  if (nestcol > 1) {
    if (nsegcol%%nestcol != 0) 
      warning("Column nesting ignored: number of column segments must be a multiple of nesting order")
    else nsncol <- nsegcol/nestcol
  }
  nsnrow <- 0
  if (nestrow > 1) {
    if (nsegrow%%nestrow != 0) 
      warning("Row nesting ignored: number of row segments must be a multiple of nesting order")
    else nsnrow <- nsegrow/nestrow
  }
  Bc <- bbasis(col, cminval, cmaxval, nsegcol, degree[1])
  nc <- ncol(Bc)
  if (length(degree) < 2) 
    degr <- degree[1]
  else degr <- degree[2]
  Br <- bbasis(row, rminval, rmaxval, nsegrow, degr)
  nr <- ncol(Br)
  if (nsncol > 0) {
    Bcn <- bbasis(col, cminval, cmaxval, nsncol, degree[1])
    ncn <- ncol(Bcn)
  }
  else ncn <- nc
  if (nsnrow > 1) {
    Brn <- bbasis(row, rminval, rmaxval, nsnrow, degr)
    nrn <- ncol(Brn)
  }
  else nrn <- nr
  diff.c <- penaltyord[[1]]
  Dc <- diff(diag(nc), diff = diff.c)
  svd.c <- svd(crossprod(Dc))
  nbc <- nc - diff.c
  U.Zc <- svd.c$u[, c(1:nbc)]
  U.Xc <- svd.c$u[, -c(1:nbc)]
  L.c <- sqrt(svd.c$d[c(1:nbc)])
  diagc <- L.c^2
  BcU <- Bc %*% U.Zc
  BcX <- Bc %*% U.Xc
  BcULi <- BcU %*% diag(1/L.c)
  if ("include" %in% eigenvalues) {
    BcZmat.df <- as.data.frame(BcULi)
    BcZmat <- BcULi
  }
  else {
    BcZmat.df <- as.data.frame(BcU)
    BcZmat <- BcU
  }
  BcZmat.df$TP.col <- col
  mat1c <- matrix(rep(1, nuc), nrow = nuc)
  BcXadj <- BcX - mat1c %*% t(mat1c) %*% BcX/nuc
  Xfc <- (svd(crossprod(BcXadj)))$u[, c(ncol(BcXadj):1)]
  BcX <- BcX %*% Xfc
  if (BcX[1, 1] < 0) 
    BcX[, 1] <- -1 * BcX[, 1]
  if (BcX[1, 2] > 0) 
    BcX[, 2] <- -1 * BcX[, 2]
  if (nsncol > 0) {
    Dcn <- diff(diag(ncn), diff = diff.c)
    svd.cn <- svd(crossprod(Dcn))
    nbcn <- ncn - diff.c
    U.Zcn <- svd.cn$u[, c(1:nbcn)]
    U.Xcn <- svd.cn$u[, -c(1:nbcn)]
    L.cn <- sqrt(svd.cn$d[c(1:nbcn)])
    BcnU <- Bcn %*% U.Zcn
    BcnX <- Bcn %*% U.Xcn
  }
  else {
    nbcn <- nbc
    BcnU <- BcU
    L.cn <- L.c
  }
  if (length(penaltyord) < 2) {
    diff.r <- penaltyord[1]
  }
  else {
    diff.r <- penaltyord[2]
  }
  Dr <- diff(diag(nr), diff = diff.r)
  svd.r <- svd(crossprod(Dr))
  nbr <- nr - diff.r
  U.Zr <- svd.r$u[, c(1:nbr)]
  U.Xr <- svd.r$u[, -c(1:nbr)]
  L.r <- sqrt(svd.r$d[c(1:nbr)])
  diagr <- L.r^2
  BrU <- Br %*% U.Zr
  BrX <- Br %*% U.Xr
  BrULi <- BrU %*% diag(1/L.r)
  if ("include" %in% eigenvalues) {
    BrZmat.df <- as.data.frame(BrULi)
    BrZmat <- BrULi
  }
  else {
    BrZmat.df <- as.data.frame(BrU)
    BrZmat <- BrU
  }
  BrZmat.df$TP.row <- row
  mat1r <- matrix(rep(1, nur), nrow = nur)
  BrXadj <- BrX - mat1r %*% t(mat1r) %*% BrX/nur
  Xfr <- (svd(crossprod(BrXadj)))$u[, c(ncol(BrXadj):1)]
  BrX <- BrX %*% Xfr
  if (BrX[1, 1] < 0) 
    BrX[, 1] <- -1 * BrX[, 1]
  if (BrX[1, 2] > 0) 
    BrX[, 2] <- -1 * BrX[, 2]
  if (nsnrow > 0) {
    Drn <- diff(diag(nrn), diff = diff.r)
    svd.rn <- svd(crossprod(Drn))
    nbrn <- nrn - diff.r
    U.Zrn <- svd.rn$u[, c(1:nbrn)]
    U.Xrn <- svd.rn$u[, -c(1:nbrn)]
    L.rn <- sqrt(svd.rn$d[c(1:nbrn)])
    BrnU <- Brn %*% U.Zrn
    BrnX <- Brn %*% U.Xrn
  }
  else {
    nbrn <- nbr
    BrnU <- BrU
    L.rn <- L.r
  }
  A <- 10^(floor(log10(max(row))) + 1)
  row.index <- rep(row, times = nuc)
  col.index <- rep(col, each = nur)
  index <- A * col.index + row.index
  C.R <- A * data[[columncoordinates]] + data[[rowcoordinates]]
  BcrZ1 <- BcnU[col.match, ] %x% matrix(rep(1, nbrn), nrow = 1, 
                                        ncol = nbrn)
  BcrZ2 <- matrix(rep(1, nbcn), nrow = 1, ncol = nbcn) %x% 
    BrnU[row.match, ]
  BcrZ <- BcrZ1 * BcrZ2
  diagrx <- rep(L.cn^2, each = nbrn)
  diagcx <- rep(L.rn^2, times = nbcn)
  if ("Lee" %in% method) {
    diagcr <- diagrx + diagcx
  }
  if ("Wood" %in% method) {
    diagcr <- diagrx * diagcx
  }
  if (!("Lee" %in% method) & !("Wood" %in% method)) {
    stop("Invalid setting of method argument")
  }
  BcrZLi <- BcrZ %*% diag(1/sqrt(diagcr))
  if ("include" %in% eigenvalues) {
    BcrZmat.df <- as.data.frame(BcrZLi)
    BcrZmat <- BcrZLi
  }
  else {
    BcrZmat.df <- as.data.frame(BcrZ)
    BcrZmat <- BcrZ
  }
  BcrZmat.df$TP.CxR <- C.R
  tracelist <- list()
  for (i in 1:diff.c) {
    nm <- paste0("Xc", i, ":Zr")
    tempmat <- (BcX[col.match, i] %x% matrix(rep(1, nbr), 
                                             nrow = 1)) * BrZmat[row.match, ]
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagr), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  for (i in 1:diff.r) {
    nm <- paste0("Zc:Xr", i)
    tempmat <- BcZmat[col.match, ] * (matrix(rep(1, nbc), 
                                             nrow = 1) %x% BrX[row.match, i])
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagc), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  if ("include" %in% eigenvalues) 
    tracelist["Zc:Zr"] <- sum(BcrZmat * BcrZmat)
  else {
    tempmatsc <- BcrZmat * (rep(1, nv) %*% matrix((1/diagcr), 
                                                  nrow = 1))
    tracelist["Zc:Zr"] <- sum(tempmatsc * BcrZmat)
  }
  outdata <- as.data.frame(data)
  outdata$TP.col <- data[[columncoordinates]]
  outdata$TP.row <- data[[rowcoordinates]]
  outdata$TP.CxR <- C.R
  BcrX1 <- BcX[col.match, ] %x% matrix(rep(1, diff.r), nrow = 1)
  BcrX2 <- matrix(rep(1, diff.c), nrow = 1) %x% BrX[row.match, 
                                                    ]
  BcrX <- BcrX1 * BcrX2
  fixed <- list()
  fixed$col <- data.frame(row.names = C.R)
  for (i in 1:diff.c) {
    c.fixed <- paste("TP.C", ".", i, sep = "")
    outdata[c.fixed] <- BcX[col.match, i]
    fixed$col[c.fixed] <- BcX[col.match, i]
  }
  fixed$row <- data.frame(row.names = C.R)
  for (i in 1:diff.r) {
    r.fixed <- paste("TP.R", ".", i, sep = "")
    outdata[r.fixed] <- BrX[row.match, i]
    fixed$row[r.fixed] <- BrX[row.match, i]
  }
  ncolX <- diff.c * diff.r
  fixed$int <- data.frame(row.names = C.R)
  for (i in 1:ncolX) {
    cr.fixed <- paste("TP.CR", ".", i, sep = "")
    outdata[cr.fixed] <- BcrX[, i]
    fixed$int[cr.fixed] <- BcrX[, i]
  }
  if (!missing(stub)) {
    cname <- paste0("BcZ", stub, ".df")
    rname <- paste0("BrZ", stub, ".df")
    crname <- paste0("BcrZ", stub, ".df")
  }
  else {
    cname <- "BcZ.df"
    rname <- "BrZ.df"
    crname <- "BcrZ.df"
  }
  mbftext <- paste0("list(TP.col=list(key=c(\"TP.col\",\"TP.col\"),cov=\"", 
                    cname, "\"),")
  mbftext <- paste0(mbftext, "TP.row=list(key=c(\"TP.row\",\"TP.row\"),cov=\"", 
                    rname, "\"),")
  mbftext <- paste0(mbftext, "TP.CxR=list(key=c(\"TP.CxR\",\"TP.CxR\"),cov=\"", 
                    crname, "\"))")
  mbflist <- eval(parse(text = mbftext))
  if ("grp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    start0 <- start
    scale <- 1
    j <- 1
    for (i in 1:diff.c) {
      nm0 <- paste0(names(fixed$col[i]), "_frow")
      listnames[j] <- nm0
      for (k in 1:nbr) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$col[[i]] * BrZmat[row.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbr, 
                      by = 1)
      start <- start + nbr
      j <- j + 1
    }
    for (i in 1:diff.r) {
      nm0 <- paste0(names(fixed$row[i]), "_fcol")
      listnames[j] <- nm0
      for (k in 1:nbc) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$row[[i]] * BcZmat[col.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbc, 
                      by = 1)
      start <- start + nbc
      j <- j + 1
    }
    m <- 0
    nm0 <- "TP_fcol_frow"
    listnames[j] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- scale * BcrZmat[, k]
    }
    grp[[j]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    end <- start + (nbcn * nbrn)
    j <- j + 1
    listnames[j] <- "All"
    grp[[j]] <- seq(from = start0 + 1, to = end, by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("sepgrp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    nm0 <- "TP_C"
    listnames[1] <- nm0
    for (i in 1:diff.c) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$col[[i]]
    }
    grp[[1]] <- seq(from = start + 1, to = start + diff.c, 
                    by = 1)
    start <- start + diff.c
    nm0 <- "TP_R"
    listnames[2] <- nm0
    for (i in 1:diff.r) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$row[[i]]
    }
    grp[[2]] <- seq(from = start + 1, to = start + diff.r, 
                    by = 1)
    start <- start + diff.r
    nm0 <- "TP_fcol"
    listnames[3] <- nm0
    for (k in 1:nbc) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcZmat[col.match, k]
    }
    grp[[3]] <- seq(from = start + 1, to = start + nbc, by = 1)
    start <- start + nbc
    nm0 <- "TP_frow"
    listnames[4] <- nm0
    for (k in 1:nbr) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BrZmat[row.match, k]
    }
    grp[[4]] <- seq(from = start + 1, to = start + nbr, by = 1)
    start <- start + nbr
    grp <- structure(grp, names = listnames)
    nm0 <- "TP_fcol_frow"
    listnames[5] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcrZmat[, k]
    }
    grp[[5]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("own" %in% asreml) {
    grp <- list()
    listnames <- list()
    listnames[1] <- "All"
    start <- length(outdata)
    nm0 <- "Xc_Zr"
    Xc_Zr <- (BcX[col.match, ] %x% matrix(rep(1, nbr), nrow = 1)) * 
      (matrix(rep(1, diff.c), nrow = 1) %x% BrZmat[row.match, 
                                                   ])
    nXc_Zr <- ncol(Xc_Zr)
    for (i in 1:nXc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Xc_Zr[, i]
    }
    nm0 <- "Zc_Xr"
    Zc_Xr <- (BcZmat[col.match, ] %x% matrix(rep(1, diff.r), 
                                             nrow = 1)) * (matrix(rep(1, nbc), nrow = 1) %x% BrX[row.match, 
                                                                                                 ])
    nZc_Xr <- ncol(Zc_Xr)
    for (i in 1:nZc_Xr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Xr[, i]
    }
    nm0 <- "Zc_Zr"
    Zc_Zr <- BcrZmat
    nZc_Zr <- ncol(Zc_Zr)
    for (i in 1:nZc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Zr[, i]
    }
    grp[[1]] <- seq(from = start + 1, to = start + nXc_Zr + 
                      nZc_Xr + nZc_Zr, by = 1)
    grp <- structure(grp, names = listnames)
  }
  res <- list()
  res$data <- outdata
  res$mbflist <- mbflist
  res[["BcZ.df"]] <- BcZmat.df
  res[["BrZ.df"]] <- BrZmat.df
  res[["BcrZ.df"]] <- BcrZmat.df
  res$dim <- c(diff.c = diff.c, nbc = nbc, nbcn = nbcn, diff.r = diff.r, 
               nbr = nbr, nbrn = nbrn)
  res$trace <- tracelist
  if ("grp" %in% asreml) 
    res$grp <- grp
  if ("sepgrp" %in% asreml) 
    res$grp <- grp
  if ("own" %in% asreml) 
    res$grp <- grp
  if ("mbf" %in% asreml) 
    res$grp <- NULL
  if (!("include" %in% eigenvalues)) 
    res$eigen <- list(diagc = diagc, diagr = diagr, diagcr = diagcr)
  res
}

bbasis <- function (x, xl, xr, ndx, deg) 
{
  tpower <- function(x, t, p) {
    (x - t)^p * (x > t)
  }
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}

tps <- function (columncoordinates, rowcoordinates, nsegments=NULL, 
                      minbound=NULL, maxbound=NULL, degree = c(3, 3), penaltyord = c(2, 2), 
                      nestorder = c(1, 1), asreml = "grp", eigenvalues = "include", 
                      method = "Lee", stub = NULL) 
{
  if (missing(columncoordinates)) 
    stop("columncoordinates argument must be set")
  if (missing(rowcoordinates)) 
    stop("rowcoordinates argument must be set")
  col <- columncoordinates
  nuc <- length(col)
  col.match <- match(columncoordinates, col)
  row <- sort(unique(rowcoordinates))
  nur <- length(row)
  row.match <- match(rowcoordinates, row)
  nv <- length(columncoordinates)
  if (is.null(minbound)) {
    cminval <- min(col)
    rminval <- min(row)
  } else {
    cminval <- min(c(minbound[1], min(col)))
    if (length(minbound) < 2) {
      rminval <- min(c(minbound[1], min(row)))
    }
    else {
      rminval <- min(c(minbound[2], min(row)))
    }
  }
  if (is.null(maxbound)) {
    cmaxval <- max(col)
    rmaxval <- max(row)
  }
  else {
    cmaxval <- max(c(maxbound[1], max(col)))
    if (length(maxbound) < 2) {
      rmaxval <- max(c(maxbound[1], max(row)))
    }
    else {
      rmaxval <- max(c(maxbound[2], max(row)))
    }
  }
  if (is.null(nsegments)) {
    nsegcol <- nuc - 1
    nsegrow <- nur - 1
  }
  else {
    nsegcol <- max(c(nsegments[1], 2))
  }
  if (length(nsegments) < 2) {
    nsegrow <- max(c(nsegments[1], 2))
  }
  else {
    nsegrow <- max(c(nsegments[2], 2))
  }
  nestcol <- floor(nestorder[1])
  if (length(nestorder) < 2) 
    nestrow <- floor(nestorder[1])
  else nestrow <- floor(nestorder[2])
  nsncol <- 0
  if (nestcol > 1) {
    if (nsegcol%%nestcol != 0) 
      warning("Column nesting ignored: number of column segments must be a multiple of nesting order")
    else nsncol <- nsegcol/nestcol
  }
  nsnrow <- 0
  if (nestrow > 1) {
    if (nsegrow%%nestrow != 0) 
      warning("Row nesting ignored: number of row segments must be a multiple of nesting order")
    else nsnrow <- nsegrow/nestrow
  }
  Bc <- bbasis(col, cminval, cmaxval, nsegcol, degree[1])
  nc <- ncol(Bc)
  if (length(degree) < 2) 
    degr <- degree[1]
  else degr <- degree[2]
  Br <- bbasis(row, rminval, rmaxval, nsegrow, degr)
  nr <- ncol(Br)
  if (nsncol > 0) {
    Bcn <- bbasis(col, cminval, cmaxval, nsncol, degree[1])
    ncn <- ncol(Bcn)
  }
  else ncn <- nc
  if (nsnrow > 1) {
    Brn <- bbasis(row, rminval, rmaxval, nsnrow, degr)
    nrn <- ncol(Brn)
  }
  else nrn <- nr
  diff.c <- penaltyord[[1]]
  Dc <- diff(diag(nc), diff = diff.c)
  svd.c <- svd(crossprod(Dc))
  nbc <- nc - diff.c
  U.Zc <- svd.c$u[, c(1:nbc)]
  U.Xc <- svd.c$u[, -c(1:nbc)]
  L.c <- sqrt(svd.c$d[c(1:nbc)])
  diagc <- L.c^2
  BcU <- Bc %*% U.Zc
  BcX <- Bc %*% U.Xc
  BcULi <- BcU %*% diag(1/L.c)
  if ("include" %in% eigenvalues) {
    BcZmat.df <- as.data.frame(BcULi)
    BcZmat <- BcULi
  }
  else {
    BcZmat.df <- as.data.frame(BcU)
    BcZmat <- BcU
  }
  BcZmat.df$TP.col <- col
  mat1c <- matrix(rep(1, nuc), nrow = nuc)
  BcXadj <- BcX - mat1c %*% t(mat1c) %*% BcX/nuc
  Xfc <- (svd(crossprod(BcXadj)))$u[, c(ncol(BcXadj):1)]
  BcX <- BcX %*% Xfc
  if (BcX[1, 1] < 0) 
    BcX[, 1] <- -1 * BcX[, 1]
  if (BcX[1, 2] > 0) 
    BcX[, 2] <- -1 * BcX[, 2]
  if (nsncol > 0) {
    Dcn <- diff(diag(ncn), diff = diff.c)
    svd.cn <- svd(crossprod(Dcn))
    nbcn <- ncn - diff.c
    U.Zcn <- svd.cn$u[, c(1:nbcn)]
    U.Xcn <- svd.cn$u[, -c(1:nbcn)]
    L.cn <- sqrt(svd.cn$d[c(1:nbcn)])
    BcnU <- Bcn %*% U.Zcn
    BcnX <- Bcn %*% U.Xcn
  }
  else {
    nbcn <- nbc
    BcnU <- BcU
    L.cn <- L.c
  }
  if (length(penaltyord) < 2) {
    diff.r <- penaltyord[1]
  }
  else {
    diff.r <- penaltyord[2]
  }
  Dr <- diff(diag(nr), diff = diff.r)
  svd.r <- svd(crossprod(Dr))
  nbr <- nr - diff.r
  U.Zr <- svd.r$u[, c(1:nbr)]
  U.Xr <- svd.r$u[, -c(1:nbr)]
  L.r <- sqrt(svd.r$d[c(1:nbr)])
  diagr <- L.r^2
  BrU <- Br %*% U.Zr
  BrX <- Br %*% U.Xr
  BrULi <- BrU %*% diag(1/L.r)
  if ("include" %in% eigenvalues) {
    BrZmat.df <- as.data.frame(BrULi)
    BrZmat <- BrULi
  }
  else {
    BrZmat.df <- as.data.frame(BrU)
    BrZmat <- BrU
  }
  BrZmat.df$TP.row <- row
  mat1r <- matrix(rep(1, nur), nrow = nur)
  BrXadj <- BrX - mat1r %*% t(mat1r) %*% BrX/nur
  Xfr <- (svd(crossprod(BrXadj)))$u[, c(ncol(BrXadj):1)]
  BrX <- BrX %*% Xfr
  if (BrX[1, 1] < 0) 
    BrX[, 1] <- -1 * BrX[, 1]
  if (BrX[1, 2] > 0) 
    BrX[, 2] <- -1 * BrX[, 2]
  if (nsnrow > 0) {
    Drn <- diff(diag(nrn), diff = diff.r)
    svd.rn <- svd(crossprod(Drn))
    nbrn <- nrn - diff.r
    U.Zrn <- svd.rn$u[, c(1:nbrn)]
    U.Xrn <- svd.rn$u[, -c(1:nbrn)]
    L.rn <- sqrt(svd.rn$d[c(1:nbrn)])
    BrnU <- Brn %*% U.Zrn
    BrnX <- Brn %*% U.Xrn
  }
  else {
    nbrn <- nbr
    BrnU <- BrU
    L.rn <- L.r
  }
  A <- 10^(floor(log10(max(row))) + 1)
  row.index <- rep(row, times = nuc)
  col.index <- rep(col, each = nur)
  index <- A * col.index + row.index
  C.R <- A * columncoordinates + rowcoordinates
  BcrZ1 <- BcnU[col.match, ] %x% matrix(rep(1, nbrn), nrow = 1, 
                                        ncol = nbrn)
  BcrZ2 <- matrix(rep(1, nbcn), nrow = 1, ncol = nbcn) %x% 
    BrnU[row.match, ]
  BcrZ <- BcrZ1 * BcrZ2
  diagrx <- rep(L.cn^2, each = nbrn)
  diagcx <- rep(L.rn^2, times = nbcn)
  if ("Lee" %in% method) {
    diagcr <- diagrx + diagcx
  }
  if ("Wood" %in% method) {
    diagcr <- diagrx * diagcx
  }
  if (!("Lee" %in% method) & !("Wood" %in% method)) {
    stop("Invalid setting of method argument")
  }
  BcrZLi <- BcrZ %*% diag(1/sqrt(diagcr))
  if ("include" %in% eigenvalues) {
    BcrZmat.df <- as.data.frame(BcrZLi)
    BcrZmat <- BcrZLi
  }
  else {
    BcrZmat.df <- as.data.frame(BcrZ)
    BcrZmat <- BcrZ
  }
  BcrZmat.df$TP.CxR <- C.R
  tracelist <- list()
  for (i in 1:diff.c) {
    nm <- paste0("Xc", i, ":Zr")
    tempmat <- (BcX[col.match, i] %x% matrix(rep(1, nbr), 
                                             nrow = 1)) * BrZmat[row.match, ]
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagr), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  for (i in 1:diff.r) {
    nm <- paste0("Zc:Xr", i)
    tempmat <- BcZmat[col.match, ] * (matrix(rep(1, nbc), 
                                             nrow = 1) %x% BrX[row.match, i])
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagc), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  if ("include" %in% eigenvalues) 
    tracelist["Zc:Zr"] <- sum(BcrZmat * BcrZmat)
  else {
    tempmatsc <- BcrZmat * (rep(1, nv) %*% matrix((1/diagcr), 
                                                  nrow = 1))
    tracelist["Zc:Zr"] <- sum(tempmatsc * BcrZmat)
  }
  # outdata <- as.data.frame(data)
  outdata <- data.frame(TP.col=columncoordinates)
  outdata$TP.row <- rowcoordinates
  outdata$TP.CxR <- C.R
  BcrX1 <- BcX[col.match, ] %x% matrix(rep(1, diff.r), nrow = 1)
  BcrX2 <- matrix(rep(1, diff.c), nrow = 1) %x% BrX[row.match, 
  ]
  BcrX <- BcrX1 * BcrX2
  fixed <- list()
  fixed$col <- data.frame(row.names = C.R)
  for (i in 1:diff.c) {
    c.fixed <- paste("TP.C", ".", i, sep = "")
    outdata[c.fixed] <- BcX[col.match, i]
    fixed$col[c.fixed] <- BcX[col.match, i]
  }
  fixed$row <- data.frame(row.names = C.R)
  for (i in 1:diff.r) {
    r.fixed <- paste("TP.R", ".", i, sep = "")
    outdata[r.fixed] <- BrX[row.match, i]
    fixed$row[r.fixed] <- BrX[row.match, i]
  }
  ncolX <- diff.c * diff.r
  fixed$int <- data.frame(row.names = C.R)
  for (i in 1:ncolX) {
    cr.fixed <- paste("TP.CR", ".", i, sep = "")
    outdata[cr.fixed] <- BcrX[, i]
    fixed$int[cr.fixed] <- BcrX[, i]
  }
  if (!missing(stub)) {
    cname <- paste0("BcZ", stub, ".df")
    rname <- paste0("BrZ", stub, ".df")
    crname <- paste0("BcrZ", stub, ".df")
  }
  else {
    cname <- "BcZ.df"
    rname <- "BrZ.df"
    crname <- "BcrZ.df"
  }
  mbftext <- paste0("list(TP.col=list(key=c(\"TP.col\",\"TP.col\"),cov=\"", 
                    cname, "\"),")
  mbftext <- paste0(mbftext, "TP.row=list(key=c(\"TP.row\",\"TP.row\"),cov=\"", 
                    rname, "\"),")
  mbftext <- paste0(mbftext, "TP.CxR=list(key=c(\"TP.CxR\",\"TP.CxR\"),cov=\"", 
                    crname, "\"))")
  mbflist <- eval(parse(text = mbftext))
  if ("grp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    start0 <- start
    scale <- 1
    j <- 1
    for (i in 1:diff.c) {
      nm0 <- paste0(names(fixed$col[i]), "_frow")
      listnames[j] <- nm0
      for (k in 1:nbr) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$col[[i]] * BrZmat[row.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbr, 
                      by = 1)
      start <- start + nbr
      j <- j + 1
    }
    for (i in 1:diff.r) {
      nm0 <- paste0(names(fixed$row[i]), "_fcol")
      listnames[j] <- nm0
      for (k in 1:nbc) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$row[[i]] * BcZmat[col.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbc, 
                      by = 1)
      start <- start + nbc
      j <- j + 1
    }
    m <- 0
    nm0 <- "TP_fcol_frow"
    listnames[j] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- scale * BcrZmat[, k]
    }
    grp[[j]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    end <- start + (nbcn * nbrn)
    j <- j + 1
    listnames[j] <- "All"
    grp[[j]] <- seq(from = start0 + 1, to = end, by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("sepgrp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    nm0 <- "TP_C"
    listnames[1] <- nm0
    for (i in 1:diff.c) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$col[[i]]
    }
    grp[[1]] <- seq(from = start + 1, to = start + diff.c, 
                    by = 1)
    start <- start + diff.c
    nm0 <- "TP_R"
    listnames[2] <- nm0
    for (i in 1:diff.r) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$row[[i]]
    }
    grp[[2]] <- seq(from = start + 1, to = start + diff.r, 
                    by = 1)
    start <- start + diff.r
    nm0 <- "TP_fcol"
    listnames[3] <- nm0
    for (k in 1:nbc) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcZmat[col.match, k]
    }
    grp[[3]] <- seq(from = start + 1, to = start + nbc, by = 1)
    start <- start + nbc
    nm0 <- "TP_frow"
    listnames[4] <- nm0
    for (k in 1:nbr) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BrZmat[row.match, k]
    }
    grp[[4]] <- seq(from = start + 1, to = start + nbr, by = 1)
    start <- start + nbr
    grp <- structure(grp, names = listnames)
    nm0 <- "TP_fcol_frow"
    listnames[5] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcrZmat[, k]
    }
    grp[[5]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("own" %in% asreml) {
    grp <- list()
    listnames <- list()
    listnames[1] <- "All"
    start <- length(outdata)
    nm0 <- "Xc_Zr"
    Xc_Zr <- (BcX[col.match, ] %x% matrix(rep(1, nbr), nrow = 1)) * 
      (matrix(rep(1, diff.c), nrow = 1) %x% BrZmat[row.match, 
      ])
    nXc_Zr <- ncol(Xc_Zr)
    for (i in 1:nXc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Xc_Zr[, i]
    }
    nm0 <- "Zc_Xr"
    Zc_Xr <- (BcZmat[col.match, ] %x% matrix(rep(1, diff.r), 
                                             nrow = 1)) * (matrix(rep(1, nbc), nrow = 1) %x% BrX[row.match, 
                                             ])
    nZc_Xr <- ncol(Zc_Xr)
    for (i in 1:nZc_Xr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Xr[, i]
    }
    nm0 <- "Zc_Zr"
    Zc_Zr <- BcrZmat
    nZc_Zr <- ncol(Zc_Zr)
    for (i in 1:nZc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Zr[, i]
    }
    grp[[1]] <- seq(from = start + 1, to = start + nXc_Zr + 
                      nZc_Xr + nZc_Zr, by = 1)
    grp <- structure(grp, names = listnames)
  }
  res <- list()
  res$data <- outdata
  res$mbflist <- mbflist
  res[["BcZ.df"]] <- BcZmat.df
  res[["BrZ.df"]] <- BrZmat.df
  res[["BcrZ.df"]] <- BcrZmat.df
  res[["All"]] <- as.matrix(outdata[,grp$All])
  res$dim <- c(diff.c = diff.c, nbc = nbc, nbcn = nbcn, diff.r = diff.r, 
               nbr = nbr, nbrn = nbrn)
  res$trace <- tracelist
  if ("grp" %in% asreml) 
    res$grp <- grp
  if ("sepgrp" %in% asreml) 
    res$grp <- grp
  if ("own" %in% asreml) 
    res$grp <- grp
  if ("mbf" %in% asreml) 
    res$grp <- NULL
  if (!("include" %in% eigenvalues)) 
    res$eigen <- list(diagc = diagc, diagr = diagr, diagcr = diagcr)
  res
}