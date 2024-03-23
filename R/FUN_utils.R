
vcsExtract <- function(object){
  nre <- length(object$sigma)
  namesre <- names(object$sigma)
  vcs <- list()
  for(i in 1:nre){
    toextract <- which(object$constraints[[i]] > 0,arr.ind = TRUE)
    vcs[[i]] <- object$sigma[[i]][toextract]
    names1 <- apply(toextract,1,function(x){paste(colnames(object$sigma[[i]])[x[1]],colnames(object$sigma[[i]])[x[2]], sep="-")})
    names2 <- paste(namesre[i],names1, sep=".")
    names(vcs[[i]]) <- names2
  }
  vcs <- unlist(vcs)
  return(vcs)
}

#### =========== ####
## SUMMARY FUNCTION mmer #
#### =========== ####
"summary.mmer" <- function(object, ...) {

  replace.values <- function(Values,Search,Replace){
    dd0 <- data.frame(Values)
    vv <- which(Values%in%Search)
    dd <- data.frame(Search,Replace)
    rownames(dd) <- Search
    dd0[vv,"Values"] <- as.character(dd[Values[vv],"Replace"])
    return(dd0[,1])
  }

  if(!object$reshapeOutput){stop("summary function only works for reshaped output.", call. = FALSE)}
  #dim(object$u.hat)
  digits = max(3, getOption("digits") - 3)
  #forget <- length(object)

  groupss.nn <- lapply(object$U,function(x){
    unlist(lapply(x,length))
  })
  groupss.nn <- do.call(rbind,groupss.nn)

  lll <- object$monitor[1,]
  lll <- lll[which(lll > 1e-300 | lll < 0)]
  lll2 <- lll[length(lll)]

  LLAIC <- data.frame(as.numeric(lll2), as.numeric(object$AIC),
                      as.numeric(object$BIC), object$method, object$convergence)
  colnames(LLAIC) = c("logLik","AIC","BIC","Method","Converge")
  rownames(LLAIC) <- "Value"

  method=object$method
  #extract fixed effects
  coef <- as.data.frame((object$Beta))#, Std.Error=(matrix(sqrt(diag(object$Var.beta.hat)),ncol=1)), t.value=(matrix((object$beta.hat-0)/sqrt(diag(object$Var.beta.hat)), ncol=1)))
  # if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}

  ## se and t values for fixed effects
  ts <- ncol(object$sigma[[1]])
  s2.beta <- diag(as.matrix(object$VarBeta))
  coef$Std.Error <- sqrt(abs(s2.beta))
  coef$t.value <- coef$Estimate/coef$Std.Error
  # print(coef)
  # nse.beta <- length(s2.beta)/ts
  # inits <- seq(1,length(s2.beta),nse.beta)
  # ends <- inits+nse.beta-1
  # seti <- list() # stardard errors partitioned by trait
  # for(u in 1:ts){
  #   prox <- data.frame(coef[,u],sqrt(abs(s2.beta[inits[u]:ends[u]])))
  #   prox$`t value` <- prox[,1]/prox[,2]
  #   colnames(prox) <- c("Estimate","Std. Error","t value")
  #   rownames(prox) <- rownames(coef)
  #   seti[[u]] <- prox
  # }
  # names(seti) <- colnames(object$sigma[[1]])

  vcsl <- list()
  consl <- list()
  for(k in 1:length(object$sigma)){
    x <- object$sigma[[k]]
    y <- object$constraints[[k]]
    xn <- names(object$sigma)[k]
    vcs <- numeric()
    cons <- numeric()
    counter <-1
    for(i in 1:ncol(x)){
      for(j in i:ncol(x)){
        # print(y[i,j])
        if( y[i,j] != 0 ){
          vcs[counter] <- x[i,j]
          cons[counter] <- y[i,j]
          names(vcs)[counter] <- paste(colnames(x)[i],colnames(x)[j],sep="-" )
          counter <- counter+1
        }
      }
    }
    vcsl[[xn]] <- as.data.frame(vcs)
    consl[[xn]] <- as.data.frame(cons)
  }
  mys2 <- do.call(rbind,vcsl)
  mycons <- do.call(rbind,consl)

  rrr <- lapply(vcsl,rownames)
  rrr <- rrr[which(unlist(lapply(rrr, length)) > 0)]
  for(o in 1:length(rrr)){rrr[[o]] <- paste(names(rrr)[o],rrr[[o]],sep=".")}
  rownames(mys2) <- as.vector(unlist(rrr))

  varcomp <- as.data.frame(cbind(mys2,sqrt(abs(diag(object$sigmaSE)))))
  varcomp[,3] <- varcomp[,1]/varcomp[,2]
  colnames(varcomp) <- c("VarComp","VarCompSE","Zratio")
  varcomp$Constraint <- replace.values(mycons$cons, 1:3, c("Positive","Unconstr","Fixed"))

  output <- list(groups=groupss.nn, varcomp=varcomp, betas=coef, method=method,logo=LLAIC)
  attr(output, "class")<-c("summary.mmer", "list")
  return(output)
}

"print.summary.mmer"<-function (x, digits = max(3, getOption("digits") - 3),  ...){

  nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$varcomp)),1,nchar))),na.rm = TRUE)

  if(nmaxchar0 < 26){
    nmaxchar0 <- 26
  } # + 26 spaces we have nmaxchar0+26  spaces to put the title

  nmaxchar <- nmaxchar0+34 ## add spaces from the 3 columns
  nmaxchar2 <- nmaxchar0+18
  nmaxchar3 <- nmaxchar0+34-46 #round(nmaxchar0/2)
  rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
  rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
  digits = max(3, getOption("digits") - 3)
  ################################################
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by REML",rlt,"\n", collapse = ""))
  cat(paste(rlh," sommer 4.3 ",rlh, "\n", collapse = ""))
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\n")
  cat("")
  print(x$logo)#, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nVariance-Covariance components:\n")
  print(x$varcomp, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nFixed effects:\n")
  if(nrow(x$betas) > 8){
    print(x$betas[1:8,], digits = digits)
    cat("   ... please access the object to see more\n")
  }else{
    print(x$betas, digits = digits)
  }
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nGroups and observations:\n")
  print(x$groups, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nUse the '$' sign to access results and parameters")#\nArguments set to FALSE for multiresponse models:\n'draw', and 'gwas.plots'\n")
  ################################################
}

"summary.mmec" <- function(object, ...) {

  replace.values <- function(Values,Search,Replace){
    dd0 <- data.frame(Values)
    vv <- which(Values%in%Search)
    dd <- data.frame(Search,Replace)
    rownames(dd) <- Search
    dd0[vv,"Values"] <- as.character(dd[Values[vv],"Replace"])
    return(dd0[,1])
  }

  digits = max(3, getOption("digits") - 3)

  lll <- object$llik
  lll2 <- lll[length(lll)]

  LLAIC <- data.frame(as.numeric(lll2), as.numeric(object$AIC),
                      as.numeric(object$BIC), "AI", object$convergence)
  colnames(LLAIC) = c("logLik","AIC","BIC","Method","Converge")
  rownames(LLAIC) <- "Value"
  method="AI"
  coef <- data.frame(Estimate=object$b)

  ## se and t values for fixed effects
  nX <- length(object$b)
  VarBeta <- object$Ci[1:nX,1:nX]
  s2.beta <- diag(as.matrix(VarBeta))
  coef$Std.Error <- sqrt(abs(s2.beta))
  coef$t.value <- coef$Estimate/coef$Std.Error

  mys2 <- object$monitor[,which(object$llik[1,] == max(object$llik[1,]))]
  varcomp <- as.data.frame(cbind(mys2,sqrt(diag(ginv(object$avInf/2)))))
  varcomp[,3] <- varcomp[,1]/varcomp[,2]
  colnames(varcomp) <- c("VarComp","VarCompSE","Zratio")

  varcomp$Constraint <- replace.values(object$constraints, 1:3, c("Positive","Unconstr","Fixed"))

  output <- list(varcomp=varcomp, betas=coef, method=method,logo=LLAIC)
  attr(output, "class")<-c("summary.mmec", "list")
  return(output)
}

"print.summary.mmec"<-function (x, digits = max(3, getOption("digits") - 3),  ...){

  nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$varcomp)),1,nchar))),na.rm = TRUE)

  if(nmaxchar0 < 26){
    nmaxchar0 <- 26
  } # + 26 spaces we have nmaxchar0+26  spaces to put the title

  nmaxchar <- nmaxchar0+34 ## add spaces from the 3 columns
  nmaxchar2 <- nmaxchar0+18
  nmaxchar3 <- nmaxchar0+34-46 #round(nmaxchar0/2)
  rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
  rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
  digits = max(3, getOption("digits") - 3)
  ################################################
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by REML",rlt,"\n", collapse = ""))
  cat(paste(rlh," sommer 4.3 ",rlh, "\n", collapse = ""))
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\n")
  cat("")
  print(x$logo)#, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nVariance-Covariance components:\n")
  print(x$varcomp, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nFixed effects:\n")
  if(nrow(x$betas) > 8){
    print(x$betas[1:8,], digits = digits)
    cat("   ... please access the object to see more\n")
  }else{
    print(x$betas, digits = digits)
  }
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nUse the '$' sign to access results and parameters")#\nArguments set to FALSE for multiresponse models:\n'draw', and 'gwas.plots'\n")
  ################################################
}

#### =========== ####
## FITTED FUNCTION ##
#### =========== ####

"fitted.mmer" <- function(object,...){

  if(is.null(object$call$random)){
    originalModelForMatrices <- mmer(fixed=object$call$fixed,
                                     # random=object$call$random,
                                     rcov=object$call$rcov,
                                     data=object$dataOriginal, returnParam = TRUE,#reshapeOutput =FALSE,
                                     init = object$sigma_scaled, constraints = object$constraints,
                                     naMethodY = object$call$naMethodY,
                                     naMethodX = object$call$naMethodX,...)

    originalModelForParameters <- mmer(fixed=object$call$fixed,
                                       # random=object$call$random,
                                       rcov=object$call$rcov, nIters=1, verbose=FALSE,
                                       data=object$dataOriginal, reshapeOutput =FALSE,
                                       init = object$sigma_scaled, constraints = object$constraints,
                                       naMethodY = object$call$naMethodY,
                                       naMethodX = object$call$naMethodX,...)
  }else{
    originalModelForMatrices <- mmer(fixed=object$call$fixed,
                                     random=object$call$random,
                                     rcov=object$call$rcov,
                                     data=object$dataOriginal, returnParam = TRUE,#reshapeOutput =FALSE,
                                     init = object$sigma_scaled, constraints = object$constraints,
                                     naMethodY = object$call$naMethodY,
                                     naMethodX = object$call$naMethodX,...)

    originalModelForParameters <- mmer(fixed=object$call$fixed,
                                       random=object$call$random,
                                       rcov=object$call$rcov, nIters=1, verbose=FALSE,
                                       data=object$dataOriginal, reshapeOutput =FALSE,
                                       init = object$sigma_scaled, constraints = object$constraints,
                                       naMethodY = object$call$naMethodY,
                                       naMethodX = object$call$naMethodX,...)
  }
  ys <- object$terms$response[[1]]
  nt <- length(ys) # number of traits
  TT <- diag(nt) # diagonal matrix

  Xo <- do.call(cbind,originalModelForMatrices$X)
  X.mv.original <- kronecker(Xo,TT)

  Xb=X.mv.original%*%originalModelForParameters$Beta

  Zu=NULL
  if(!is.null(object$call$random)){
    nz <- length(originalModelForMatrices$Z)
    Zu <- vector(mode = "list", length = nz) # list for Zu
    for(ir in 1:nz){ # for each random effect
      Z <- originalModelForMatrices$Z[[ir]] # provisional Z
      Zu[[ir]] <- kronecker(Z,TT) %*% originalModelForParameters$U[[ir]] # calculate Zu
    }
  }
  names(Zu) <- names(object$U)

  if(!is.null(object$call$random)){
    y.hat <- Xb + Reduce("+",Zu) # y.hat = Xb + Zu.1 + ... + Zu.n
  }else{
    y.hat <- Xb # y.hat = Xb
  }

  Zudf <- as.matrix(do.call(cbind,Zu)); colnames(Zudf) <- paste0(names(Zu),".fitted")
  Xbdf <- as.matrix(Xb); colnames(Xbdf) <- "Xb.fitted"
  y.hat.df <- matrix(y.hat[,1],byrow = TRUE, ncol=nt)
  colnames(y.hat.df) <- paste0(object$terms$response[[1]],".fitted")
  dataWithFitted <- cbind(object$data,y.hat.df,Xbdf,Zudf)
  # build summary table
  nLevels <- c(unlist(lapply(originalModelForMatrices$X,ncol)), unlist(lapply(originalModelForMatrices$Z,ncol)))
  namesLevels <- c(unlist(object$terms$fixed),names(object$U))
  formLevels <- c(rep("fixed",length(unlist(object$terms$fixed))),rep("random",length(names(object$U))))
  id <- 1:length(formLevels)
  used <- rep(TRUE,length(id))
  fittedSummary <- data.frame(namesLevels,formLevels,nLevels,id,used)
  colnames(fittedSummary) <- c("term","type","nLevels","id","used")

  return(list(dataWithFitted=dataWithFitted,Xb=Xb, Zu=Zu,fittedSummary=fittedSummary))
}

"print.fitted.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("\n  The fitted values are obtained by adding Xb + Zu.1 + ... + Zu.n
                 containing: \n")
  ))
  print(x$fittedSummary)
  cat(blue(paste("\n  head of fitted values: \n")
  ))
  head(x$dataWithFitted,...)
}

"fitted.mmec" <- function(object,...){

  ff <- object$W %*% object$bu

  return(ff)
}

"print.fitted.mmec"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("\n  The fitted values are obtained by adding Xb + Zu.1 + ... + Zu.n
                 containing: \n")
  ))
  cat(blue(paste("\n  head of fitted values: \n")
  ))
  head(x,...)
}

#### =========== ######
## RESIDUALS FUNCTION #
#### =========== ######
"residuals.mmer" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)

  pp <- fitted.mmer(object)
  responses <- object$terms$response[[1]]
  dat <- pp$dataWithFitted
  for(ir in 1:length(responses)){
    dat[,paste0(responses[ir],".residuals")] = dat[,responses[ir]] - dat[,paste0(responses[ir],".fitted")]
  }

  return(dat)
}

"print.residuals.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

"residuals.mmec" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  ff <- fitted.mmec(object)
  e <- object$y - ff
  return(e)
}

"print.residuals.mmec"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}
#### =========== ######
## RANEF FUNCTION #
#### =========== ######

"randef" <- function(object) {
  output<- object$U
  return(output)
}

#### =========== ####
## COEF FUNCTION ####
#### =========== ####

"coef.mmer" <- function(object, ...){
  object$Beta
}

"print.coef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

"coef.mmec" <- function(object, ...){
  object$b
}

"print.coef.mmec"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

#### =========== ####
## ANOVA FUNCTION ###
#### =========== ####
anova.mmer <- function(object, object2=NULL, type=1, ...) {
  sequential.fit <- function(object, type=1){

    fixed <- object$call$fixed
    yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
    fixedtermss <- apply(data.frame(yuyuf),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fixedtermss <- fixedtermss[which(fixedtermss != "-1")]
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
    fixedtermss <- c("1",fixedtermss)
    response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    # responsef <- as.formula(paste(response,"~1"))
    # fixedtermss <- as.list(fixedtermss)
    # fixedtermss[[length(fixedtermss)+1]] <- unlist(fixedtermss)
    SSM <- numeric(length = length(fixedtermss)+1)
    df <- numeric(length = length(fixedtermss)+1)
    models <- character(length = length(fixedtermss)+1)
    # if(length(fixedtermss) < 0){
    #   until <- 1
    # }else{
    until <- length(fixedtermss)
    # }

    namesl <- character(length = length(fixedtermss)+1)
    for(i in (until+1):1){

      if(i == (length(fixedtermss)+1)){ # first the full model
        myf <- paste(response,"~",paste(fixedtermss,collapse = " + "))
        fixedi <- as.formula(myf)
        usef <- fixedtermss
      }else if(i == 1){
        myf <- paste(response,"~1")
        fixedi <- as.formula(myf)
      }else{
        if(type == 1){
          usef <- setdiff(fixedtermss,fixedtermss[i])
          myf <- paste(response,"~",paste(usef,collapse = " + "))
          fixedi <- as.formula(myf)
        }else if(type == 2){
          usef <- fixedtermss[1:i]
          myf <- paste(response,"~",paste(usef,collapse = " + "))
          fixedi <- as.formula(myf)
        }else{
          stop("anova type not enabled.", call. = FALSE)
        }
      }
      namesl[i] <- paste(usef,collapse = ",")
      models[i] <- myf
      cat("\n")
      cat(paste("Fitting:", myf))
      # cat("\n")
      # print(fixedi)
      if(is.null(object$call$random)){
        prov0 <- mmer(fixed=fixedi,
                      # random = object$call$random,
                      rcov=object$call$rcov, nIters=1,
                      data=object$data, returnParam = FALSE,reshapeOutput=FALSE,
                      verbose = FALSE,
                      init = object$sigma_scaled, constraints = object$constraints,
                      naMethodY = object$call$naMethodY,
                      naMethodX = object$call$naMethodX)
        prov <- mmer(fixed=fixedi,
                     # random = object$call$random,
                     rcov=object$call$rcov,
                     data=object$data, returnParam = TRUE,
                     naMethodY = object$call$naMethodY,
                     naMethodX = object$call$naMethodX)
      }else{
        prov0 <- mmer(fixed=fixedi,
                      random = object$call$random,
                      rcov=object$call$rcov, nIters=1,
                      data=object$data, returnParam = FALSE,reshapeOutput=FALSE,
                      verbose = FALSE,
                      init = object$sigma_scaled, constraints = object$constraints,
                      naMethodY = object$call$naMethodY,
                      naMethodX = object$call$naMethodX)
        prov <- mmer(fixed=fixedi,
                     random = object$call$random,
                     rcov=object$call$rcov,
                     data=object$data, returnParam = TRUE,
                     naMethodY = object$call$naMethodY,
                     naMethodX = object$call$naMethodX)
      }
      Y <- prov[[1]]
      n= nrow(Y)
      J = matrix(1, nrow=n, ncol=n)
      # Total sum of Square
      SSTO = t(Y) %*% Y - (1/n)*t(Y)%*%J%*%Y
      Xlist <- list()
      # ncolsx <- numeric()
      for(o in 1:length(prov[[3]])){
        Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
        # ncolsx[o] <- ncol(Xlist[[o]])
      }
      X <- do.call(cbind,Xlist)
      head(X)
      ncolsx <- unlist(lapply(Xlist,ncol))
      ncolsx <- ncolsx[which(ncolsx>0)]
      xdf <- ncolsx[length(ncolsx)]
      #X <- as.matrix(X[,1:2])
      b = as.matrix(prov0$Beta)
      # regression sum of square
      ss <- (SSR = t(b)%*%t(X)%*%Y - (1/n)%*%t(Y)%*%J%*%Y)
      if(i == (length(fixedtermss)+1)){
        SSM[i] = as.numeric(ss)
      }else{
        SSM[i] = SSM[length(fixedtermss)+1] - as.numeric(ss)
      }
      df[i] <- xdf#ncol(X)
    }

    myanova <- data.frame(Df=rev(df), Sum.Sq=rev(SSM))

    # print(myanova)
    # print(c("Full",rev(namesl[-c(length(namesl))])))
    # rownames(myanova) <- c("Full",rev(namesl[-c(length(namesl))]))

    #else{
    #   rownames(myanova) <- c("Full", "None", rev(fixedtermss[-c(1)]))
    # }
    myanova$Mean.Sq <- myanova$Sum.Sq/myanova$Df

    vc <- object$sigma
    vare <- as.numeric(vc[[length(vc)]])
    dfe <- nrow(Y)-sum(myanova$Df)
    myanova2 <- data.frame(Df=dfe, Sum.Sq=NA,Mean.Sq=vare)
    rownames(myanova2) <- "Residuals"
    myanovaf <- rbind(myanova,myanova2)

    myanovaf$F.value <- myanovaf$Mean.Sq/vare
    # print(c(myanovaf,dfe))
    myanovaf$`Pr(>F)` <- round(apply(myanovaf,1,function(x){1-pf(x[4],x[1],dfe)}),4)
    myanovaf[nrow(myanovaf),4:5] <- NA
    myanovaf$Models <- c(rev(models),"")
    if(type==1){
      rownames(myanovaf) <- c("Full", rev(fixedtermss), "Residuals")
      myanovaf <- myanovaf[-c(nrow(myanovaf)-1),]
    }else{
      rownames(myanovaf) <- c("Full","all", rev(fixedtermss[-c(1)]), "Residuals")
      # print(rev(fixedtermss))
    }
    cat("\n\n")
    return(myanovaf)
  }
  signifo <- function(x){
    if(x >= 0 & x < 0.001){y="***"}
    if(x >= 0.001 & x < 0.01){y="**"}
    if(x >= 0.01 & x < 0.05){y="*"}
    if(x >= 0.05 & x < 0.1){y="."}
    if(x > 0.1){y=""}
    return(y)
  }
  ########################################
  digits = max(3, getOption("digits") - 3)
  if(is.null(object2)){
    # stop("The 'anova' function for the sommer package only works to compare mixed models by likelihood ratio tests (LRT), was not intended to provide regular sum of squares output.")
    result <- sequential.fit(object,type=type)
  }else{
    #if(object$maxim){ # user used REML=TRUE, not possible to do LRT
    #  stop("Please fit the models using ML instead of REML by setting the argument REML=FALSE and try again")
    #}else{ #if user used REML=FALSE, then proceed
    if(object$method != object2$method){
      stop("Error! When comparing models please use the same method for the fitted models.")
    }else{
      yu <- summary(object)
      yu2 <- summary(object2)
      dis=c(dim(yu$varcomp)[1]+dim(object$Beta)[1]*dim(object$Beta)[2],
            dim(yu2$varcomp)[1]+dim(object2$Beta)[1]*dim(object2$Beta)[2]) # dimensions
      mods=c("mod1","mod2")
      lls=c( yu$logo[1,1],  yu2$logo[1,1]) # likelihoods
      aics=c(object$AIC, object2$AIC) # AIC's
      bics=c(object$BIC, object2$BIC) # AIC's
      vv=which(dis == max(dis))[1] # which has more variance components BIGGER
      vv2=c(1:2)[which(c(1:2)!= vv)] # SMALLER
      LR = (lls[vv] - lls[vv2])
      r.stat= abs(-2*((LR))) # -2(LL1 - LL2)
      df=dis[vv]-dis[vv2]
      chichi=pchisq((r.stat), df, lower.tail=FALSE)
      if(chichi > 1e-5){
        chichi <- round(chichi,5)
      }
      chichi2=paste(as.character(chichi),signifo(chichi), sep=" ")
      ### construct the table
      cat("Likelihood ratio test for mixed models\n")
      cat("==============================================================\n")
      result=data.frame(Df=c(dis[vv],dis[vv2]), AIC=c(aics[vv],aics[vv2]),
                        BIC=c(bics[vv],bics[vv2]), loLik=c(lls[vv],lls[vv2]),
                        Chisq=c("",as.character(round(r.stat,5))),
                        ChiDf=c("",as.character(df)), PrChisq=c("",chichi2 ))
      rownames(result) <- c(mods[vv],mods[vv2])
      print(result)
      # print(result)
      cat("==============================================================\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
    #}
  }
  return(result)
}

anova.mmec <- function(object, object2=NULL, ...) {
  signifo <- function(x){
    if(x >= 0 & x < 0.001){y="***"}
    if(x >= 0.001 & x < 0.01){y="**"}
    if(x >= 0.01 & x < 0.05){y="*"}
    if(x >= 0.05 & x < 0.1){y="."}
    if(x > 0.1){y=""}
    return(y)
  }
  ########################################
  digits = max(3, getOption("digits") - 3)
  if(is.null(object2)){
    stop("The 'anova' function for the sommer package only works to compare mixed models by likelihood ratio tests (LRT), was not intended to provide regular sum of squares output.")
    # result <- sequential.fit(object,type=type)
  }else{
    dis=c(
          nrow(object$monitor)+nrow(object$b),
          nrow(object2$monitor)+nrow(object2$b)
          ) # dimensions
    mods=c("mod1","mod2")
    lls=c( object$llik[ncol(object$llik)],  object2$llik[ncol(object2$llik)] ) # likelihoods
    aics=c(object$AIC, object2$AIC) # AIC's
    bics=c(object$BIC, object2$BIC) # AIC's
    vv=which(dis == max(dis))[1] # which has more variance components BIGGER
    vv2=c(1:2)[which(c(1:2)!= vv)] # SMALLER
    LR = (lls[vv] - lls[vv2])
    r.stat= abs(-2*((LR))) # -2(LL1 - LL2)
    df=dis[vv]-dis[vv2]
    chichi=pchisq((r.stat), df, lower.tail=FALSE)
    if(chichi > 1e-5){
      chichi <- round(chichi,5)
    }
    chichi2=paste(as.character(chichi),signifo(chichi), sep=" ")
    ### construct the table
    cat("Likelihood ratio test for mixed models\n")
    cat("==============================================================\n")
    result=data.frame(Df=c(dis[vv],dis[vv2]), AIC=c(aics[vv],aics[vv2]),
                      BIC=c(bics[vv],bics[vv2]), loLik=c(lls[vv],lls[vv2]),
                      Chisq=c("",as.character(round(r.stat,5))),
                      ChiDf=c("",as.character(df)), PrChisq=c("",chichi2 ))
    rownames(result) <- c(mods[vv],mods[vv2])
    print(result)
    cat("==============================================================\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

    #}
  }
  return(result)
}
#### =========== ####
## PLOTING FUNCTION #
#### =========== ####
plot.mmer <- function(x, stnd=TRUE, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  # provisional model
  prov <- mmer(fixed=x$call$fixed,
               #random=x$call$random,
               rcov=x$call$rcov,
               data=x$data, returnParam = TRUE,#reshape.results=TRUE,
               naMethodY = x$call$naMethodY,
               naMethodX = x$call$naMethodX)
  Xlist <- list()
  for(o in 1:length(prov[[3]])){
    Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
  }
  Xm <- do.call(cbind,Xlist)
  # std vs residuals, QQplot (std vs teor quantiles), sqrt(std residuals) vs fitted, std res vs leverage = cook's distance
  traits <- ncol(x$fitted)
  layout(matrix(1:4,2,2))

  resp <- x$terms$response[[1]]
  # ff <- fitted(x)
  rr <- residuals.mmer(x)
  for(i in 1:traits){

    plot(rr[,paste0(resp[i],".fitted")],scale(rr[,paste0(resp[i],".residuals")]),pch=20, col=transp("cadetblue"), ylab="Std Residuals", xlab="Fitted values", main="Residual vs Fitted", bty="n", ...); grid()
    plot(rr[,paste0(resp[i],".fitted")],sqrt(abs(scale(rr[,paste0(resp[i],".residuals")]))),pch=20, col=transp("thistle4"), ylab="Sqrt Abs Std Residuals", xlab="Fitted values", main="Scale-Location", bty="n", ...); grid()
    qqnorm(scale(rr[,paste0(resp[i],".residuals")]), pch=20, col=transp("tomato1"), ylab="Std Residuals", bty="n",...); grid()
    hat <- Xm%*%solve(t(Xm)%*%x$Vi%*%Xm)%*%t(Xm)%*%x$Vi # leverage including variance from random effects H= X(X'V-X)X'V-
    plot(diag(hat), scale(rr[,paste0(resp[i],".residuals")]), pch=20, col=transp("blue"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  }
  #####################
  layout(matrix(1,1,1))
}

plot.mmec <- function(x, stnd=TRUE, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  layout(matrix(1:4,2,2))
  # ff <- fitted(x)
  rr <- residuals.mmec(x)
  # for(i in 1:traits){

    plot(rr,scale(rr),pch=20, col=transp("cadetblue"), ylab="Std Residuals", xlab="Fitted values", main="Residual vs Fitted", bty="n", ...); grid()
    plot(rr,sqrt(abs(scale(rr))),pch=20, col=transp("thistle4"), ylab="Sqrt Abs Std Residuals", xlab="Fitted values", main="Scale-Location", bty="n", ...); grid()

    qqnorm(scale(rr), pch=20, col=transp("tomato1"), ylab="Std Residuals", bty="n",...); grid()
    # hat <- Xm%*%solve(t(Xm)%*%x$Vi%*%Xm)%*%t(Xm)%*%x$Vi # leverage including variance from random effects H= X(X'V-X)X'V-
    hat = x$W %*% x$Ci %*% t(x$W)
    plot(diag(hat), scale(rr), pch=20, col=transp("blue"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  # }
  #####################
  layout(matrix(1,1,1))
}
##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]   Solving Mixed Model Equations in R (sommer) 4.3.4 (2024-04-01) []",sep="")),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]   ------------- Multivariate Linear Mixed Models --------------  []")),appendLF=TRUE)
    packageStartupMessage(paste0(blue("[]   Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("*")), bgRed(white(" "))),"                      []")),appendLF=TRUE)
    packageStartupMessage(blue("[]   Published: PLoS ONE 2016, 11(6):1-15                           []"),appendLF=TRUE)
    packageStartupMessage(blue("[]   Dedicated to the University of Chapingo and UW-Madison         []"),appendLF=TRUE)
    packageStartupMessage(blue("[]   Type 'vignette('v1.sommer.quick.start')' for a short tutorial  []"),appendLF=TRUE)
    packageStartupMessage(blue("[]   Type 'citation('sommer')' to know how to cite sommer           []"),appendLF=TRUE)
    packageStartupMessage(blue(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(blue("sommer is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(blue("Current source is available at https://github.com/covaruber/sommer"),appendLF=TRUE)
    packageStartupMessage(blue("If needed, install as: devtools::install_github('covaruber/sommer')"),appendLF=TRUE)

  }
  invisible()
}
