atcg1234 <- function(data, ploidy=2, format="ATCG", maf=0, multi=TRUE, silent=FALSE, by.allele=FALSE, imp=TRUE, ref.alleles=NULL){

  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  ##### START GBS.TO.BISNP DATA ######
  gbs.to.bisnp <- function(x) {
    y <- rep(NA,length(x))
    y[which(x=="A")] <- "AA"
    y[which(x=="T")] <- "TT"
    y[which(x=="C")] <- "CC"
    y[which(x=="G")] <- "GG"
    y[which(x=="R")] <- "AG"
    y[which(x=="Y")] <- "CT"
    y[which(x=="S")] <- "CG"
    y[which(x=="W")] <- "AT"
    y[which(x=="K")] <- "GT"
    y[which(x=="M")] <- "AC"
    y[which(x=="+")] <- "++"
    y[which(x=="0")] <- "NN"
    y[which(x=="-")] <- "--"
    y[which(x=="N")] <- NA
    return(y)
  }
  ##### END GBS.TO.BISNP DATA ######
  imputeSNP <- function(data){
    #######
    data2 <- apply(data,2,function(x){
      areNA <- which(is.na(x))
      if(length(areNA)>0){
        pos.all <- table(data[,1])
        totake <- names(pos.all)[which(pos.all == max(pos.all))]
        x[areNA] <- totake
      }
      return(x)
    })
    #######
    return(data2)
  }
  #### apply with progress bar ######
  apply_pb <- function(X, MARGIN, FUN, ...){
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  ###### zero.one function
  zero.one <- function(da){
    # this function takes a matrix of markers in biallelic format and returns a matrix of
    # presense/absense of alleles
    mar.nam <- colnames(da)#unique(gsub("\\.\\d","", names(da))) # find a dot and a number after the dot
    mat.list <- list(NA) # list of matrices for each marker
    wi=0 # counter
    if(!silent){
      count <- 0
      tot <- length(mar.nam)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    for(i in 1:length(mar.nam)){ # for each marker
      wi=wi+1
      if(!silent){
        count <- count + 1
      }
      
      v <- which(colnames(da)==mar.nam[i])#grep(mar.nam[i], colnames(da))
      
      if(length(v)==0){
        qqqqq <- grep(mar.nam[i-1],names(da))
        qqqqq2 <- names(da)[qqqqq[length(qqqqq)] + 1]
        
        stop(paste("Marker",qqqqq2,"has a problem"), call.=FALSE)
      }else if(length(v) == 1){ # for markers with a single column
        prov <- matrix(da[,v])
      }else{prov <- da[,v]}
      ##################################
      alls <- unique(unlist(strsplit(prov,"")))
      alls <- alls[which(!is.na(alls))]
      ninds <- dim(prov)[1]
      fff <- apply(data.frame(alls),1,function(h){
        temp <- numeric(length = ninds)
        temp[grep(h,prov)]<-1
        #make sure is full rank
        
        return(temp)
      })#1 # assigning 1's
      #if(FULL){ # if user want to make sure only get the columns that will ensure full rank
      #  fff <- t(unique(t(fff)))
      #}
      colnames(fff) <- paste(mar.nam[i],alls, sep="/")
      
      mat.list[[i]] <- fff
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      
    }
    
    fin.mat <- do.call(cbind,mat.list)
    rownames(fin.mat) <- rownames(da)
    
    # make full rank
    #q <- qr(fin.mat)
    #chas <- q$pivot[seq(q$rank)]
    #fin.mat <- as.matrix(fin.mat[,chas])
    #############
    return(fin.mat)
  }
  
  ## remove all markers or columns that are all missing data
  all.na <- apply(data,2,function(x){length(which(is.na(x)))/length(x)})
  bad.na <- which(all.na==1)
  if(length(bad.na) > 0){
    data <- data[,-bad.na]
  }
  
  if(is.null(ref.alleles)){
  #############################
  if(by.allele){ ####&&&&&&&&&&&&&&&&&&&&&& use zero.one function
    user.code <- apply(data[,c(1:(round(dim(data)[2]/20)))], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
    AA <- sum(user.code, na.rm = TRUE)/length(user.code)
    if(AA > .9){ # means user is using single letter
      rnd <- rownames(data)
      data <- apply(data,2,gbs.to.bisnp);#W2[1:5,1:5]
      rownames(data) <- rnd
    }
    M <- zero.one(data)
    
  }else{ ###&&&&&&&&&&&&&&&&&&&&&&&&
    n.g <- apply(data,2,function(x){length(table(x))})
    bad <- which(n.g > 3)
    if(length(bad) == dim(data)[2]){
      cat("Error. All your markers are multiallelic. This function requires at least one bi-allelic marker\n")
    }
    
    # tells you which markers have double letter code, i.e. TT instead of T
    # 1: has only one letter
    # 0: has two letters
    user.code <- apply(data[,c(1:(round(dim(data)[2]/20)))], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
    AA <- sum(user.code, na.rm = TRUE)/length(user.code)
    if(AA > .9){
      rrn <- rownames(data)
      
      cat("Converting GBS or single-letter code to biallelic code\n")
      if(silent){
        data <- apply(data, 2,gbs.to.bisnp)
      }else{
        data <- apply_pb(data, 2,gbs.to.bisnp) 
      }
      rownames(data) <- rrn
      data <- as.data.frame(data)
    }
    #### apply with progress bar ######
    s1 <- rownames(data)
    s2 <- colnames(data)
    data <- as.data.frame(t(data))
    rownames(data) <- s2
    colnames(data) <- s1
    bases <- c("A", "C", "G", "T","l","m","n","p","h","k","-","+","e","f","g","a","b","c","d")
    ## get reference allele function
    get.ref <- function(x, format) {
      if (format == "numeric") {
        ref.alt <- c(0, 1)
      }
      if (format == "AB") {
        ref.alt <- c("A", "B")
      }
      if (format == "ATCG") {
        y <- paste(na.omit(x), collapse = "")
        ans <- apply(array(bases), 1, function(z, y) {
          length(grep(z, y, fixed = T))
        }, y)
        if (sum(ans) > 2) {
          ref.alt <- (bases[which(ans == 1)])[1:2]
          #stop("Error in genotype matrix: More than 2 alleles")
        }
        if (sum(ans) == 2) {
          ref.alt <- bases[which(ans == 1)]
        }
        if (sum(ans) == 1) {
          ref.alt <- c(bases[which(ans == 1)], NA)
        }
      }
      return(ref.alt)
    }
    
    get.multi <- function(x, format) {
      if (format == "numeric") {
        ref.alt <- c(0, 1)
      }
      if (format == "AB") {
        ref.alt <- c("A", "B")
      }
      if (format == "ATCG") {
        y <- paste(na.omit(x), collapse = "")
        ans <- apply(array(bases), 1, function(z, y) {
          length(grep(z, y, fixed = T))
        }, y)
        if (sum(ans) > 2) {
          ref.alt <- TRUE
        }
        if (sum(ans) == 2) {
          ref.alt <- FALSE
        }
        if (sum(ans) == 1) {
          ref.alt <- FALSE
        }
      }
      return(ref.alt)
    }
    
    ####################################
    ## convert to matrix format
    ####################################
    markers <- as.matrix(data)
    ####################################
    # get reference alleles
    ####################################
    cat("Obtaining reference alleles\n")
    if(silent){
      tmp <- apply(markers, 1, get.ref, format=format)
    }else{
      tmp <- apply_pb(markers, 1, get.ref, format=format) 
    }
    
    if(multi){ # if markers with multiple alleles should be removed
      cat("Checking for markers with more than 2 alleles. If found will be removed.\n")
      if(silent){
        tmpo <- apply(markers, 1, get.multi, format = format)
      }else{
        tmpo <- apply_pb(markers, 1, get.multi, format = format) 
      }
      ###&&&&&&&&&&&& HERE WE MUST INSERT THE NEW FUNCTIONALITY, WHERE WE DETECTED MULTIPLE ALLELES
      multi.allelic <- which(!tmpo) # good markers
      markers <- markers[multi.allelic,]
      tmp <- tmp[, multi.allelic]
    }
    
    Ref <- tmp[1, ]
    Alt <- tmp[2, ]
    ####################################
    ## bind reference allele and markers and convert to numeric format based on the 
    # reference/alternate allele found
    ####################################
    cat("Converting to numeric format\n")
    if(silent){
      M <- apply(cbind(Ref, markers), 1, function(x) {
        y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
        ans <- as.integer(lapply(y, function(z) {
          ifelse(z[1] < 0, ploidy, ploidy - length(z))
        }))
        return(ans)
      })
    }else{
      M <- apply_pb(cbind(Ref, markers), 1, function(x) {
        y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
        ans <- as.integer(lapply(y, function(z) {
          ifelse(z[1] < 0, ploidy, ploidy - length(z))
        }))
        return(ans)
      })
    }
    
    gid.geno <- s1 #colnames(geno)
    rownames(M) <- gid.geno
    ####################################
    # identify bad markers
    ####################################
    bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
    if (bad > 0) {
      stop("Invalid marker calls.")
    }
    
  }
  #rownames(M) <- rownames(data)
  ####################################
  rownames(tmp) <- c("Alt","Ref")
  }else{# user provides reference alleles and just want a conversion
    
    common.mark <- intersect(colnames(data), colnames(ref.alleles))
    data <- data[,common.mark]
    tmp <- ref.alleles[,common.mark]; #rownames(refa) <- c("Alt","Ref")
    cat("Converting to numeric format\n")
    M <- apply_pb(data.frame(1:ncol(data)),1,function(k){
      x <- as.character(data[,k])
      x2 <- strsplit(x,"")
      x3 <- unlist(lapply(x2,function(y){length(which(y == tmp[2,k]))}))
      return(x3)
    })
    #M <- M-1
    colnames(M) <- colnames(data)
    
  }
  
  ####################################
  # by column or markers calculate MAF
  ####################################
  cat("Calculating minor allele frequency (MAF)\n")
  if(silent){
    MAF <- apply(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }else{
    MAF <- apply_pb(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }
  ####################################
  # which markers have MAF > 0, JUST GET THOSE
  ####################################
  polymorphic <- which(MAF > maf)
  M <- M[, polymorphic]
  ####################################
  # function to impute markers with the mode
  ####################################
  
  # time to impute
  if(imp){
    missing <- which(is.na(M))
    if (length(missing) > 0) {
      cat("Imputing missing data with mode \n")
      if(silent){
        M <- apply(M, 2, impute.mode)
      }else{
        M <- apply_pb(M, 2, impute.mode)
      }
    }
  }else{
    cat("Imputation not required. Be careful using non-imputed matrices in mixed model solvers\n")
  }
  ## ploidy 2 needs to be adjusted to -1,0,1
  if(ploidy == 2){
    M <- M - 1
  }
  
  return(list(M=M,ref.alleles=tmp))
}

build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE){
  # build hybrid marker matrix
  
  if(!is.null(custom.hyb)){
    pheno <- custom.hyb
    found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
    if(found != 3){
      stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
    }
    return.combos.only=FALSE
  }else{
    a <- rownames(M1)
    b <- rownames(M2)
    pheno <- expand.grid(a,b)
    pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
    pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=":")
  }
  
  if(!return.combos.only){
    # check that marker matrices are in -1,0,1 format
    checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
    checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))
    
    checkM1[which(checkM1 > 0)] <- 1
    checkM2[which(checkM2 > 0)] <- 1
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    
    
    ## add markers coming from parents M1
    Z1 <- model.matrix(~Var1-1,pheno);dim(Z1); 
    colnames(Z1) <- gsub("Var1","",colnames(Z1))
    M1 <- M1[colnames(Z1),]
    #M1[1:4,1:4]; Z1[1:4,1:4]; 
    ## add markers coming from parents M2
    Z2 <- model.matrix(~Var2-1,pheno);dim(Z2); 
    colnames(Z2) <- gsub("Var2","",colnames(Z2))
    M2 <- M2[colnames(Z2),]
    #M2[1:4,1:4]; Z2[1:4,1:4];  
    
    ## create the 
    # Z3 <- model.matrix(~hybrid-1,pheno);dim(Z3);
    # colnames(Z3) <- gsub("hybrid","",colnames(Z3))
    # hyb.names <- colnames(Z3)[as.vector(apply(Z3,1,function(x){which(x==1)}))] # names of hybrids
    hyb.names <- pheno$hybrid
    ## marker matrix for hybrids one for each parent
    cat(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))
    
    # M1 <- as(M1, Class="sparseMatrix")
    # M2 <- as(M2, Class="sparseMatrix")
    # Z1 <- as(Z1, Class="sparseMatrix")
    # Z2 <- as(Z2, Class="sparseMatrix")
    
    cat("Extracting M1 contribution\n")
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Md <- Z1 %*% M1;  # was already converted to -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      Md <- 2*Z1 %*% M1 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      Md <- Z1 %*% M1 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    cat("Extracting M2 contribution\n")
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      Mf <- Z2 %*% M2;  # was already converted to -1,1
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      Mf <- 2*Z2 %*% M2 - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      Mf <- Z2 %*% M2 - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
    }
    
    ## marker matrix coded as additive -1,0,1
    Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
    rownames(Mdf) <- hyb.names
    #hist(Mdf)
    
    ## dominance matrix for hybrids (0,1 coded)
    Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
    rownames(Delta) <- hyb.names
    #hist(Delta)
    cat("Done!!\n")
    return(list(HMM.add=Mdf, HMM.dom=Delta, data.used=pheno))
    
  }else{
    return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
  }
}

h2.fun <- function(object, data, gTerm=NULL, eTerm=NULL, md=NULL) {
  
  if(missing(object)){
    stop("Please provide a model object of type mmer.\n", call. = FALSE)
  }
  if(is.null(gTerm)){
    stop("Please specify the gTerm in your model.\n", call. = FALSE)
  }
  if(missing(data)){
    stop("Please provide the dataset used to fit the model.\n", call. = FALSE)
  }
  
  prov <- mmer(fixed=object$call$fixed,
               #random=object$call$random,
               rcov=object$call$rcov,
               data=object$data, return.param = TRUE,#reshape.results=TRUE,
               na.method.Y = object$call$na.method.Y, 
               na.method.X = object$call$na.method.X)
  
  if(!is.null(eTerm)){
    elevels <- as.character(na.omit(unique(data[,eTerm])))
  }else{
    data[,"FIELD"] <- "FIELD1"
    eTerm <- "FIELD"
    elevels <- "FIELD1"
  }
  
  h2s <- list()
  for(e in elevels){ # e <- elevels[1]
    if(length(elevels)==1){
      geTerm <- paste(gTerm,sep=":")
    }else{
      geTerm <- paste(e,gTerm,sep=":")
    }
    geTerm
    ## now get the heritability for each location
    v <- which(names(object$PevU) %in% geTerm)
    if(length(v)==0){
      stop("The gTerm provided was not found in the model. You may not be providing the entire name.\n")
    }
    
    ## average sample size in such environment
    expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
    expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
    subdata <- droplevels(data[which(data[,eTerm]==e),])
    
    check <- grep("\\(",gTerm)
    if(length(check) > 0){
      gTerm2 <- expi(gTerm)
      nn <- median(table(subdata[,gTerm2]))
    }else{
      nn <- median(table(subdata[,gTerm]))
    }
    #print(nn)
    
    if(length(elevels)==1){ # single field
      f2 <- grep("units",names(object$sigma))
      ve <- object$sigma[[f2]]
    }else{
      f1 <- grep(e,names(object$sigma))
      f2 <- grep("units",names(object$sigma))
      ve <- object$sigma[[intersect(f1,f2)]]
    }
    ve
    ##
    G <- 0
    pev.mat <- 0
    vcs <- numeric()
    counter <- 0
    meandiag <- numeric() # store the diagonal of the covariance matrix
    co <- 0
    for(gs in geTerm){# gs <- geTerm[1]
      co = co+1
      counter <- counter + 1
      if(counter==1){
        pev.mat <- pev.mat + object$PevU[[gs]][[1]]
      }; pev.mat[1:4,1:4]
      if(length(prov$K)==0){A <- diag(nrow(pev.mat))}else{A <- prov$K[[v]]}
      if(is.null(md)){
        meandiag[co] <- mean(diag(A))
      }else{meandiag[co] <- md}
      vc <- as.numeric(object$sigma[[gs]])
      vcs[[gs]] <- vc
      G <- G + A * vc
    }
    meandiag <- mean(meandiag) # make sure we take the mean if more than one genetic term was fitted although that is not allowed
    
    G[1:3,1:3]
    try(ginv <- solve(G), silent = TRUE)
    if(class(ginv)=="try-error"){
      cat("Adding a small amount to the diagonal of A to make it positive-definite.\n")
      ginv <- solve(G+diag(1e-3,nrow(G)))
    }
    ginv[1:4,1:4]
    nv <- nrow(G)
    vc <- sum(as.numeric(vcs))
    
    id.mat<-diag(nv) #create identity matrix
    esh2<-1-(sum(pev.mat*ginv)/nv)
    #library(Matrix)
    M<-id.mat-(((1/meandiag)*ginv)%*%pev.mat)
    #M <- make.full(M)
    eM<-eigen(M)# eigenvalues of M
    sm<-rep(NA,nv)
    sm<-ifelse((Im(eM$values) ==0 |Re(eM$values) ==0 )==T, 1, 0)#number of non-zero eigenvectors
    neM<-sm*Re(eM$values) #  eigen values get a zero, and non-zero are multiplied by 1
    seM<-sum(neM) # add all eigen values, full heritability
    h2<- seM/sum(sm)
    
    pevs <- mean(diag(pev.mat))
    ve <- as.vector(ve)
    
    h2.cullis <- 1 - (mean(diag(pev.mat))/((meandiag)*vc))
    h2.stdrd <- as.numeric(vc/(vc+(ve/nn)))
    
    h2s[[e]] <- data.frame(PEV=pevs, Vg=vc, Ve=ve,N.rep=nn, H2.stdrd=h2.stdrd,H2.cullis=h2.cullis,H2.oakey.eigen=h2)
    # print(ve)
  }
  h2d <- as.data.frame(do.call(rbind,h2s))
  h2d[,eTerm] <- names(h2s)
  h2d[,"Trait"] <- colnames(object$Y)
  #h2d <- h2d[,c("Trait","Env",setdiff(colnames(h2d),c("Trait","Env")))]
  return(h2d)
  # pevs <- diag(mix2$PEV.u.hat$id$Yield)
  # vg <- mix2$var.comp$id[1,1]
  # ve <- mix2$var.comp$units[1,1]
  # pevm <- mix2$PEV.u.hat$id$Yield
  # m <- nrow(A)
  # Gi <- solve(diag(1,m)*vg)
  # 1 - mean(pevs/(2*vg)); mean(1 - pevs/(2*vg))
  # vg/(vg+ve)
  # uu <- (0.5*Gi %*% pevm)/m
  # 1 - matrix.trace(as.matrix(uu)); 1 - sum(diag(uu))
  # 1 - mean(diag(0.5*Gi %*% pevm))
}