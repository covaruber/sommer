##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "3.5.0"))
    stop("This package requires R 3.5.0 or later")
  if(interactive()) {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(magenta(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]  Solving Mixed Model Equations in R (sommer) ", desc$Version," (", desc$Date, ")  []",sep="")),appendLF=TRUE)
    packageStartupMessage(magenta(paste("[]  ------------- Multivariate Linear Mixed Models --------------   []")),appendLF=TRUE)
    packageStartupMessage(paste0(magenta("[]  Author: Giovanny Covarrubias-Pazaran",paste0(bgGreen
                                                                                           (white(" ")), bgWhite(magenta("M")), bgRed(white(" ")),"  ", bgRed(bold(yellow(" (") )),bgRed(bold(white("W"))), bgRed(bold(yellow(") "))) ) ,"                []")),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Published: PLoS ONE 2016, 11(6):1-15                            []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Dedicated to the University of Chapingo and UW-Madison          []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Type 'vignette('sommer.qg')' for a short tutorial               []"),appendLF=TRUE)
    packageStartupMessage(magenta("[]  Type 'citation('sommer')' to know how to cite sommer            []"),appendLF=TRUE)
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
  if(object$CiMode == 0){
    s2.beta <- rep(NA, length(1:nX))
  }else if(object$CiMode == 1){
    s2.beta <- rep(NA, length(1:nX))
  }else if(object$CiMode == 2){
    VarBeta <- object$Ci[1:nX,1:nX]
    s2.beta <- diag(as.matrix(VarBeta))
  }else{}
  
  coef$Std.Error <- sqrt(abs(s2.beta))
  coef$t.value <- coef$Estimate/coef$Std.Error
  
  varcomp <- object$covParNative

  # lapply(object$covStruct, function(x){x$free})
  # constraints <- unlist(lapply(object$thetaC, as.vector))
  # constraints <- constraints[which(constraints != 0)]
  # varcomp$Constraint <- replace.values(constraints, 1:3, c("Positive","Unconstr","Fixed"))

  output <- list(varcomp=varcomp, betas=coef, method=method,logo=LLAIC,
                 REML=if(is.null(object$REML)) TRUE else object$REML)
  attr(output, "class")<-c("summary.mmes", "list")
  return(output)
}

"print.summary.mmes"<-function (x, digits = max(3, getOption("digits") - 3),  ...){

  desc <- utils::packageDescription("sommer")
  nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$varcomp)),1,nchar))),na.rm = TRUE)

  # x$estimate <- round(x$estimate, digits = digits)
  # x$StdError    <- round(x$StdError   , digits = digits)
  # x$Zratio <- round(x$Zratio, digits = digits)
  # 
  # nmaxchar0 <- max(apply(x,1,function(y){
  #   nchar(paste(unlist(y), collapse = ""))
  # }) )
  
  if(nmaxchar0 < 26){
    nmaxchar0 <- 26
  } # + 26 spaces we have nmaxchar0+26  spaces to put the title

  nmaxchar <- nmaxchar0+44 ## add spaces from the 3 columns
  nmaxchar2 <- nmaxchar0+18
  nmaxchar3 <- nmaxchar0+34-46 #round(nmaxchar0/2)
  rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
  rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
  digits = max(3, getOption("digits") - 3)
  ################################################
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by ",
            if(isTRUE(x$REML)) "REML" else "ML", rlt,"\n", collapse = ""))
  cat(paste(rlh," sommer ",desc$Version,rlh, "\n", collapse = ""))
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

#### =========== ####################
## FACTOR-ANALYTIC LOADINGS/SCORES ##
#### =========== ####################

# Locate the single fam()/rrcm()-shaped term in a fitted mmes object and
# return its compiled factor descriptor together with the natural-scale
# parameter slice needed to rebuild loadings/specific variances.
#
# object$covStruct itself carries no names, but it is built in the same
# order as object$theta (random terms first, residual term last), which is
# named. Term names are therefore resolved through object$theta.
.mmes_fa_term <- function(object, term=NULL){

  if(!inherits(object, "mmes")){
    stop("object must be a fitted mmes model.", call.=FALSE)
  }

  structNames <- names(object$theta)
  if(is.null(structNames) || length(structNames) != length(object$covStruct)){
    stop("Internal mismatch between object$theta and object$covStruct.", call.=FALSE)
  }

  candidates <- character()
  for(i in seq_along(object$covStruct)){
    fs <- object$covStruct[[i]]$factors
    if(length(fs) == 1L && isTRUE(fs[[1]]$model %in% c("fa","rr"))){
      candidates <- c(candidates, structNames[i])
    }
  }

  if(is.null(term)){
    if(length(candidates) == 0L){
      stop("No fam()/rrcm() covariance term was found in this model.", call.=FALSE)
    }
    if(length(candidates) > 1L){
      stop(
        paste0(
          "Multiple fam()/rrcm() terms were found; please specify term as one of: ",
          paste(candidates, collapse=", ")
        ),
        call.=FALSE
      )
    }
    term <- candidates[1]
  }

  pos <- match(term, structNames)
  if(is.na(pos)){
    stop(
      paste0(
        "term '", term, "' was not found. Available terms: ",
        paste(structNames, collapse=", ")
      ),
      call.=FALSE
    )
  }

  factors <- object$covStruct[[pos]]$factors
  if(length(factors) != 1L || !isTRUE(factors[[1]]$model %in% c("fa","rr"))){
    stop(
      paste0(
        "term '", term, "' is not a single fam()/rrcm() covariance-shaping factor. ",
        "Terms combining fam()/rrcm() with additional shaping factors are not yet supported."
      ),
      call.=FALSE
    )
  }

  list(term=term, pos=pos, factor=factors[[1]])
}

# Reconstruct the normalized loadings (Lambda) and specific variances (Psi)
# of a fam()/rrcm() term such that sigma2*(Lambda %*% t(Lambda) + diag(Psi))
# reproduces object$theta[[term]] exactly (up to floating-point roundoff).
"loadings_mmes" <- function(object, term=NULL, varianceScale=TRUE, rotation=TRUE){

  located <- .mmes_fa_term(object, term)
  term <- located$term
  f <- located$factor

  q <- f$dim
  k <- f$order
  levels <- f$levels

  covPar <- object$covPar[[located$pos]]
  natural <- covPar[f$par_start:f$par_end]

  if(f$model == "fa"){
    nload <- f$fa_nload
    rows <- f$fa_row
    cols <- f$fa_col
  }else{
    nload <- f$rr_nload
    rows <- f$rr_row
    cols <- f$rr_col
  }

  loadingsRaw <- matrix(0, q, k)
  for(a in seq_len(nload)){
    loadingsRaw[rows[a], cols[a]] <- natural[a]
  }

  specificRaw <- rep(1, q)
  if(f$model == "fa" && q > 1L){
    specificRaw[-1] <- natural[nload + seq_len(q-1L)]
  }

  scale <- loadingsRaw[1,1]^2 + specificRaw[1]

  loadings <- loadingsRaw / sqrt(scale)
  specific <- specificRaw / scale

  dimnames(loadings) <- list(levels, paste0("F", seq_len(k)))
  names(specific) <- levels

  sigma2 <- unname(covPar[1])
  
  if(varianceScale){
    loadings <- loadings * sqrt(sigma2)
    specific <- specific * sqrt(sigma2)
  }
  
  if(rotation){
    V <- svd(loadings)$v
    loadings <- -loadings %*% V
    dimnames(loadings) <- list(levels, paste0("F", seq_len(k)))
  }
  
  list(
    loadings=loadings,
    specific=specific,
    sigma2=sigma2,
    model=f$model,
    term=term
  )
}

# Extract descriptor-defined covariance parameters on their native scale.
.covparams_mmes <- function(object, term=NULL){
  if(!inherits(object, "mmes")){
    stop("object must inherit from class 'mmes'.", call.=FALSE)
  }
  if(is.null(object$covStruct) || is.null(object$covPar)){
    stop("The fitted object does not contain covariance descriptors and parameters.",
         call.=FALSE)
  }

  termNames <- names(object$covStruct)
  if(is.null(termNames) || any(!nzchar(termNames))){
    termNames <- names(object$covPar)
  }
  if(is.null(termNames) || length(termNames) != length(object$covStruct)){
    termNames <- paste0("structure", seq_along(object$covStruct))
  }

  selected <- seq_along(object$covStruct)
  if(!is.null(term)){
    if(is.numeric(term)){
      selected <- as.integer(term)
      if(anyNA(selected) || any(selected < 1L | selected > length(termNames))){
        stop("Numeric term indices are outside the fitted covariance structures.",
             call.=FALSE)
      }
    }else{
      selected <- match(as.character(term), termNames)
      if(anyNA(selected)){
        stop("Unknown covariance term: ",
             paste(as.character(term)[is.na(selected)], collapse=", "),
             call.=FALSE)
      }
    }
  }

  rows <- list()
  outputIndex <- 0L

  for(i in selected){
    descriptor <- object$covStruct[[i]]
    current <- as.numeric(object$covPar[[i]])
    if(!length(current)) next
    scale <- current[1L]
    factors <- descriptor$factors

    if(!length(factors)){
      outputIndex <- outputIndex + 1L
      rows[[outputIndex]] <- data.frame(
        term=termNames[i], factor="sigma2", parameter="sigma2",
        estimate=scale, stringsAsFactors=FALSE
      )
      next
    }

    scaleAbsorbed <- FALSE
    for(j in seq_along(factors)){
      factor <- factors[[j]]
      start <- as.integer(factor$par_start)
      end <- as.integer(factor$par_end)
      factorPar <- if(end >= start) current[start:end] else numeric()
      reporter <- factor$native_report
      if(is.null(reporter) || !identical(reporter$backend, "R") ||
         !is.function(reporter$fun)){
        stop("Covariance factor ", j, " in term '", termNames[i],
             "' has no valid native reporting callback.", call.=FALSE)
      }

      absorbScale <- !scaleAbsorbed
      values <- reporter$fun(
        scale=scale,
        par=factorPar,
        factor=factor,
        absorb_scale=absorbScale
      )
      values <- unlist(values, use.names=TRUE)
      if(!length(values)) next
      if(is.null(names(values)) || any(!nzchar(names(values)))){
        stop("Native reporting callbacks must return named values.", call.=FALSE)
      }
      if(any(!is.finite(values))){
        stop("Native reporting callback returned non-finite values for term '",
             termNames[i], "'.", call.=FALSE)
      }

      outputIndex <- outputIndex + 1L
      rows[[outputIndex]] <- data.frame(
        term=termNames[i],
        factor=if(!is.null(factor$model) && nzchar(factor$model)) factor$model else
          paste0("factor", j),
        parameter=names(values),
        estimate=as.numeric(values),
        stringsAsFactors=FALSE
      )
      scaleAbsorbed <- scaleAbsorbed || absorbScale
    }

    if(!scaleAbsorbed){
      outputIndex <- outputIndex + 1L
      rows[[outputIndex]] <- data.frame(
        term=termNames[i], factor="sigma2", parameter="sigma2",
        estimate=scale, stringsAsFactors=FALSE
      )
    }
  }

  if(!length(rows)){
    return(data.frame(term=character(), factor=character(),
                      parameter=character(), estimate=numeric()))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

"covparams_mmes" <- function(object, term=NULL){
  .covparams_mmes(object, term)
}

.covparams_mmes_se <- function(object, term=NULL, rel_step=1e-6){
  if(length(rel_step) != 1L || !is.finite(rel_step) || rel_step <= 0){
    stop("rel_step must be one positive finite value.", call.=FALSE)
  }
  if(is.null(object$theta_se) || !is.matrix(object$theta_se)){
    stop("The fitted object does not contain a covariance-parameter uncertainty matrix.",
         call.=FALSE)
  }

  termNames <- names(object$covStruct)
  if(is.null(termNames) || any(!nzchar(termNames))){
    termNames <- names(object$covPar)
  }
  if(is.null(termNames) || length(termNames) != length(object$covStruct)){
    termNames <- paste0("structure", seq_along(object$covStruct))
  }

  selected <- seq_along(object$covStruct)
  if(!is.null(term)){
    if(is.numeric(term)){
      selected <- as.integer(term)
      if(anyNA(selected) || any(selected < 1L | selected > length(termNames))){
        stop("Numeric term indices are outside the fitted covariance structures.",
             call.=FALSE)
      }
    }else{
      selected <- match(as.character(term), termNames)
      if(anyNA(selected)){
        stop("Unknown covariance term: ",
             paste(as.character(term)[is.na(selected)], collapse=", "),
             call.=FALSE)
      }
    }
  }

  parameterCounts <- vapply(object$covPar, length, integer(1))
  totalParameters <- sum(parameterCounts)
  if(!all(dim(object$theta_se) == c(totalParameters, totalParameters))){
    stop("theta_se dimensions do not match the fitted covariance parameters.",
         call.=FALSE)
  }
  starts <- cumsum(c(1L, head(parameterCounts, -1L)))

  evaluateTerm <- function(index, parameterValues, template){
    candidate <- object
    candidate$covPar[[index]] <- parameterValues
    reported <- .covparams_mmes(candidate, index)
    if(nrow(reported) != nrow(template) ||
       !identical(reported$parameter, template$parameter) ||
       !identical(reported$factor, template$factor)){
      stop("Native reporting callback changed its output shape during numerical differentiation.",
           call.=FALSE)
    }
    reported$estimate
  }

  output <- vector("list", length(selected))
  for(outputIndex in seq_along(selected)){
    i <- selected[outputIndex]
    base <- .covparams_mmes(object, i)
    current <- as.numeric(object$covPar[[i]])
    nLocal <- length(current)
    jacobian <- matrix(0, nrow(base), nLocal)
    free <- as.logical(object$covStruct[[i]]$free)
    if(length(free) != nLocal){
      stop("Covariance descriptor free flags do not match covPar.", call.=FALSE)
    }

    for(k in which(free)){
      h <- rel_step * max(1, abs(current[k]))
      plus <- current
      minus <- current
      plus[k] <- plus[k] + h
      minus[k] <- minus[k] - h

      plusValue <- try(evaluateTerm(i, plus, base), silent=TRUE)
      minusValue <- try(evaluateTerm(i, minus, base), silent=TRUE)
      plusOK <- !inherits(plusValue, "try-error") && all(is.finite(plusValue))
      minusOK <- !inherits(minusValue, "try-error") && all(is.finite(minusValue))

      if(plusOK && minusOK){
        jacobian[,k] <- (plusValue - minusValue) / (2*h)
      }else if(plusOK){
        jacobian[,k] <- (plusValue - base$estimate) / h
      }else if(minusOK){
        jacobian[,k] <- (base$estimate - minusValue) / h
      }else{
        stop("Unable to numerically differentiate native covariance parameter '",
             base$parameter[1L], "' in term '", termNames[i], "'.", call.=FALSE)
      }
    }

    global <- starts[i] + seq_len(nLocal) - 1L
    covariance <- object$theta_se[global, global, drop=FALSE]
    covariance[!free,] <- 0
    covariance[,!free] <- 0
    nativeCovariance <- jacobian %*% covariance %*% t(jacobian)
    variance <- pmax(diag(nativeCovariance), 0)
    standardError <- sqrt(variance)
    zRatio <- rep(NA_real_, length(standardError))
    positiveSE <- standardError > 0
    zRatio[positiveSE] <- base$estimate[positiveSE] / standardError[positiveSE]

    base$StdError <- standardError
    base$Zratio <- zRatio
    attr(base, "vcov") <- nativeCovariance
    output[[outputIndex]] <- base
  }

  out <- do.call(rbind, output)
  rownames(out) <- NULL
  attr(out, "vcov") <- NULL
  out
}

"covparams_mmes_se" <- function(object, term=NULL, rel_step=1e-6){
  .covparams_mmes_se(object, term, rel_step)
}

# Predict per-level latent factor scores for a fam()/rrcm() term from its
# fitted loadings, covariance, and BLUPs. method="regression" (Thomson) uses
# the full fitted covariance; method="bartlett" uses only the specific
# (residual) variances and is the classic unbiased factor-score estimator.
"scores_mmes" <- function(object, term=NULL, method=c("regression","bartlett"),
                          varianceScale=TRUE, rotation=TRUE){

  method <- match.arg(method)
  fa <- loadings_mmes(object, term, varianceScale, rotation)
  term <- fa$term

  if(is.null(object$uList[[term]])){
    stop(
      paste0(
        "term '", term, "' has no BLUPs (it is a residual covariance structure); ",
        "scores_mmes() requires a random-effect fam()/rrcm() term."
      ),
      call.=FALSE
    )
  }

  L <- fa$loadings
  Sigma <- object$theta[[term]]
  U <- object$uList[[term]]

  if(!all(rownames(L) %in% colnames(U))){
    stop(
      paste0(
        "Internal mismatch between term '", term, "' loadings levels and BLUP levels."
      ),
      call.=FALSE
    )
  }
  U <- U[, rownames(L), drop=FALSE]

  if(method == "regression"){
    SigmaInv <- solve(Sigma)
    scores <- U %*% SigmaInv %*% L
  }else{
    PsiInv <- diag(1/(fa$sigma2 * fa$specific), nrow=length(fa$specific))
    scores <- U %*% PsiInv %*% L %*% solve(t(L) %*% PsiInv %*% L)
  }

  rownames(scores) <- rownames(U)
  colnames(scores) <- colnames(L)
  scores
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
    if(!is.null(object$REML) && !is.null(object2$REML) &&
       !identical(object$REML, object2$REML)){
      warning("Comparing a REML fit against a maximum-likelihood (REML=FALSE) fit is not a valid likelihood ratio test; refit both models with the same REML= setting.", call.=FALSE)
    }else if(isTRUE(object$REML) && isTRUE(object2$REML) &&
             !identical(deparse(object$args$fixed), deparse(object2$args$fixed))){
      warning("Both models were fit with REML=TRUE but specify different fixed effects. REML log-likelihoods are only comparable across models sharing the same fixed effects; refit both models with REML=FALSE for a valid likelihood ratio test on the fixed effects.", call.=FALSE)
    }
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

