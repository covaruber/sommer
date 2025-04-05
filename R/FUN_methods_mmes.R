##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    packageStartupMessage(magenta(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]   Solving Mixed Model Equations in R (sommer) 4.4.1 (2025-04-04) []",sep="")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]   ------------- Multivariate Linear Mixed Models --------------  []")),appendLF=TRUE)
    packageStartupMessage(paste0(magenta("[]   Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen(white(" ")), bgWhite(magenta("*")), bgRed(white(" "))),"                      []")),appendLF=TRUE)
    packageStartupMessage(magenta("[]   Published: PLoS ONE 2016, 11(6):1-15                           []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]   Dedicated to the University of Chapingo and UW-Madison         []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]   Type 'vignette('sommer.qg')' for a short tutorial              []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]   Type 'citation('sommer')' to know how to cite sommer           []"),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(magenta("sommer is updated on CRAN every 3-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(magenta("Current source is available at https://github.com/covaruber/sommer"),appendLF=TRUE)
    packageStartupMessage(magenta("If needed, install as: devtools::install_github('covaruber/sommer')"),appendLF=TRUE)
    
  }
  invisible()
}

#### =========== ####
## SUMMARY FUNCTION mmes #
#### =========== ####

"summary.mmes" <- function(object, ...) {

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
  varcomp <- as.data.frame(cbind(mys2,sqrt(diag(object$theta_se))))
  varcomp[,3] <- varcomp[,1]/varcomp[,2]
  colnames(varcomp) <- c("VarComp","VarCompSE","Zratio")

  constraints <- unlist(lapply(object$thetaC, as.vector))
  constraints <- constraints[which(constraints != 0)]
  varcomp$Constraint <- replace.values(constraints, 1:3, c("Positive","Unconstr","Fixed"))

  output <- list(varcomp=varcomp, betas=coef, method=method,logo=LLAIC)
  attr(output, "class")<-c("summary.mmes", "list")
  return(output)
}

"print.summary.mmes"<-function (x, digits = max(3, getOption("digits") - 3),  ...){

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
  cat(paste(rlh," sommer 4.4 ",rlh, "\n", collapse = ""))
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

"fitted.mmes" <- function(object,...){

  ff <- object$W %*% object$bu

  return(ff)
}

"print.fitted.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
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

"residuals.mmes" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  ff <- fitted.mmes(object)
  e <- object$y - ff
  return(e)
}

"print.residuals.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}
#### =========== ######
## RANEF FUNCTION #
#### =========== ######

"randef" <- function(object) {
  output<- object$u
  return(output)
}

#### =========== ####
## COEF FUNCTION ####
#### =========== ####

"coef.mmes" <- function(object, ...){
  object$b
}

"print.coef.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

#### =========== ####
## ANOVA FUNCTION ###
#### =========== ####

anova.mmes <- function(object, object2=NULL, ...) {
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

plot.mmes <- function(x, stnd=TRUE, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  layout(matrix(1:4,2,2))
  # ff <- fitted(x)
  rr <- residuals.mmes(x)
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

