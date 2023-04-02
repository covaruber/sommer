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
      ncolsData <- dim(data)[2]
      ncolsData <- max(2,round(ncolsData/20))
      # print(ncolsData)
      user.code <- apply(data[,c(1:ncolsData)], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
      print(user.code)
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
      ncolsData <- dim(data)[2]
      ncolsData <- max(2,round(ncolsData/20))
      # print(ncolsData)
      user.code <- apply(data[,c(1:ncolsData)], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
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
    
    # M1 <- as(M1, Class="dgCMatrix")
    # M2 <- as(M2, Class="dgCMatrix")
    # Z1 <- as(Z1, Class="dgCMatrix")
    # Z2 <- as(Z2, Class="dgCMatrix")
    
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

LD.decay <- function(markers,map,silent=FALSE,unlinked=FALSE,gamma=.95){
  
  #if(is.null(rownames(map))){
  good <- which(!duplicated(map$Locus))
  map <- map[good,]
  rownames(map) <- map$Locus 
  #}
  ## clean markers and map from sd=0 and markers with no position
  markers <- markers[,which(apply(markers,2,sd)>0)]
  map <- map[which(!is.na(map$LG)),] # which(apply(map,1,function(x){length(which(is.na(x)))}) ==0)
  #############
  fullmap <- map
  lgs <- unique(map$LG)
  ############################
  dat.list <- list(NA) # will contain the LD (r) and genetic distance (d) for each LG
  ############################
  if (!silent) {
    count <- 0
    tot <- length(lgs)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!unlinked){
    
    LDM <- list()
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      #ords2 <- ords
      head(ords2);tail(ords2)
      #################
      #intersect markers and map
      #################  
      inn <- intersect(ords2$Locus,colnames(markers))
      
      if(length(inn)>0){
        ords2 <- ords2[inn,]
        
        cor2.mat <- cor(as.matrix(markers[,inn]),use="pairwise.complete.obs")^2
        LDM[[k]] <- cor2.mat
        
        N <- dim(markers[,inn])[1]
        cor2.matp <- 1-pchisq(cor2.mat*N,df=1)
        
        #cor(fofo[,which(apply(fofo,2,sd) > 0)])^2
        # trnsforming to double haploid data
        # this is a mtrix of genetic distances
        mat.dist <- matrix(0,dim(ords2)[1], dim(ords2)[1])
        for(i in 1:dim(ords2)[1]){
          for(j in 1:dim(ords2)[1]){
            mat.dist[i,j] <- ords2$Position[i] - ords2$Position[j]
          }
        }
        # get absolute values
        mat.dist <- abs(mat.dist) 
        # make zero the valuies of upper and lower triangular
        cor2.mat[lower.tri(cor2.mat)] <- 0
        cor2.matp[lower.tri(cor2.mat)] <- 0
        mat.dist[lower.tri(mat.dist)] <- 0
        
        y.dist <- as.vector(mat.dist) #los deshace columna por column
        y.cor <- as.vector(cor2.mat) 
        p.cor <- as.vector(cor2.matp)
        
        dat <- data.frame(d=y.dist,r2=y.cor,p=p.cor)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
        if (!silent) {
          setTxtProgressBar(pb, (count/tot))
        }
      }else{
        dat <- data.frame(d=0,r2=0,p=0)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
      }
      
    }
    # make a big data frame
    if (!silent) {
      setTxtProgressBar(pb, (count/tot))
    }
    
    big <- do.call("rbind",dat.list)
    
    # big contains all the LD and and genetic distances for all groups
    resp <- list(by.LG=dat.list, all.LG=big, LDM=LDM)
    
  }else{ # if WANTS TO KNOW THE THRESHOLD OF UNLINKED
    doo <- intersect(colnames(markers), fullmap$Locus)
    cor2.mat <- cor(as.matrix(markers[,doo]),use="pairwise.complete.obs")^2
    
    dat.list <- numeric()#list()#store the r2 values of unlinked markers with chromosome k
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      # markers in kth chromosome
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      # markers not in kth chromosome
      ords3 <- fullmap[which(fullmap[,"LG"] != k),]
      # markers in kth and present
      sik <- intersect(ords2$Locus,colnames(cor2.mat))
      # markers not in kth and present
      nok <- intersect(ords3$Locus,colnames(cor2.mat))
      #ords2 <- ords
      #head(ords2);tail(ords2)
      #into <- setdiff(nok,ords2$Locus)
      
      step1 <- cor2.mat[sik,nok]
      #boxplot(unlist(step1))
      dat.list[k] <- quantile(sqrt(step1[upper.tri(step1)]),gamma)^2
      
      if (!silent) {
        setTxtProgressBar(pb, (count/tot))
      }
    }
    # make a big data frame
    resp <- list(by.LG=dat.list, all.LG=mean(dat.list))
  }
  return(resp)
}