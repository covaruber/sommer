
#### =========== ####
## SUMMARY FUNCTION mmes #
#### =========== ####

.mmes_require_covariance_information <- function(object){
  if(identical(object$covarianceInformation, "none") ||
     identical(object$control$output$covariance, "none") ||
     identical(object$engineControl$covarianceInformation, "none")){
    stop("Covariance information was omitted; refit with output=list(covariance='full') for this inference method.",
         call.=FALSE)
  }
  invisible(NULL)
}

.mmes_require_coefficient_matrix <- function(object){
  if(identical(object$retainCoefficientMatrix, FALSE) ||
     identical(object$control$output$coefficientMatrix, FALSE) ||
     identical(object$engineControl$retainCoefficientMatrix, FALSE)){
    stop("Coefficient matrix C was not retained; refit with output=list(coefficientMatrix=TRUE) for this method.",
         call.=FALSE)
  }
  invisible(NULL)
}

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
  laplaceUncertainty <- if(inherits(object, "mmes.glmm.laplace"))
    .mmes_laplace_uncertainty(object) else NULL
  if(inherits(object, "mmes.glmm.laplace")){
    LLAIC <- data.frame(logLik=object$llik,
      AIC=object$AIC, BIC=object$BIC, Method="Laplace",
      Converge=object$convergence, row.names="Value")
  }else if(inherits(object, "mmes.glmm")){
    LLAIC <- data.frame(Deviance=object$deviance,
                        Dispersion=if(length(object$dispersion) > 1L)
                          paste(names(object$dispersion), signif(object$dispersion, 4),
                                sep="=", collapse=", ") else object$dispersion,
                        Method="PQL", Converge=isTRUE(object$pqlConverged),
                        row.names="Value")
  }
  method="AI"
  coef <- data.frame(Estimate=object$b)

  ## se and t values for fixed effects
  nX <- length(object$b)
  if(inherits(object, "mmes.glmm.laplace")){
    coef$Std.Error <- laplaceUncertainty$betaSummary[, "Std. Error"]
    coef$t.value <- laplaceUncertainty$betaSummary[, "z value"]
    varcomp <- laplaceUncertainty$covParNative
    method <- "Laplace"
  }else if(identical(object$engine, "direct")){
    VarBeta <- object$VarBeta
    s2.beta <- rep(NA_real_, nX)
    if(!is.null(VarBeta)){
      varianceColumns <- if(is.null(object$fixedRank)) seq_len(nX) else object$fixedRank$retained
      s2.beta[varianceColumns] <- diag(as.matrix(VarBeta))
    }
  }else if(object$CiMode == 0){
    s2.beta <- rep(NA, length(1:nX))
  }else if(object$CiMode == 1){
    s2.beta <- rep(NA, length(1:nX))
  }else if(object$CiMode == 2){
    VarBeta <- object$Ci[1:nX,1:nX]
    s2.beta <- Matrix::diag(VarBeta)
  }else{}
  
  if(!inherits(object, "mmes.glmm.laplace")){
    coef$Std.Error <- sqrt(abs(s2.beta))
    coef$t.value <- coef$Estimate/coef$Std.Error
  }
  if(!is.null(object$fixedRank)){
    coef$Status <- "retained (conditional on zero constraints)"
    coef$Status[object$fixedRank$omitted] <- "aliased (constrained to zero)"
    coef$Std.Error[object$fixedRank$omitted] <- NA_real_
    coef$t.value[object$fixedRank$omitted] <- NA_real_
  }
  
if(!inherits(object, "mmes.glmm.laplace")){
  varcomp <- object$covParNative[,c("term", "parameter", "estimate", "StdError", "Zratio",
                                  intersect("status", names(object$covParNative)))]
}

  # lapply(object$covStruct, function(x){x$free})
  # constraints <- unlist(lapply(object$thetaC, as.vector))
  # constraints <- constraints[which(constraints != 0)]
  # varcomp$Constraint <- replace.values(constraints, 1:3, c("Positive","Unconstr","Fixed"))

  output <- list(varcomp=varcomp, betas=coef, method=method,logo=LLAIC,
                 remlOmissions=object$remlOmissions, fixedRank=object$fixedRank,
                 REML=if(is.null(object$REML)) TRUE else object$REML,
                 family=object$family,
                 familyParameters=if(is.null(laplaceUncertainty))
                   object$familyParameters else laplaceUncertainty$familyParameters,
                 zeroCoefficients=if(is.null(laplaceUncertainty))
                   object$zeroCoefficients else laplaceUncertainty$zeroCoefficients,
                 thresholds=if(is.null(laplaceUncertainty))
                   object$thresholds else laplaceUncertainty$thresholds,
                 laplaceDiagnostics=if(is.null(laplaceUncertainty))
                   object$laplaceDiagnostics else
                     utils::modifyList(object$laplaceDiagnostics, list(
                       curvatureValid=laplaceUncertainty$valid,
                       curvatureReason=laplaceUncertainty$reason,
                       uncertaintyValid=laplaceUncertainty$valid,
                       uncertaintyComputed=TRUE)),
                 pqlConverged=object$pqlConverged,
                 fitMethod=if(inherits(object, "mmes.glmm.laplace")) "Laplace" else NULL)
  if(identical(object$covarianceInformation, "none")){
    output$covarianceInformation <- "none"
  }
  if(identical(object$retainCoefficientMatrix, FALSE)){
    output$retainCoefficientMatrix <- FALSE
  }
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
  cat(paste("\n",rlt,
            if(identical(x$fitMethod, "Laplace"))
              "Generalized Linear Mixed Model fit by " else
              "Multivariate Linear Mixed Model fit by ",
            if(!is.null(x$fitMethod)) x$fitMethod else
            if(!is.null(x$pqlConverged)) "PQL" else
            if(isTRUE(x$REML)) "REML" else "ML", rlt,"\n", collapse = ""))
  cat(paste(rlh," sommer 4.4 ",rlh, "\n", collapse = ""))
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\n")
  cat("")
  print(x$logo)#, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nVariance-Covariance components:\n")
  if(!is.null(x$remlOmissions)){
    cat("Structurally unestimable REML components were omitted; their zeros are not variance estimates and their standard errors are unavailable.\n")
  }

  if(identical(x$covarianceInformation, "none")){
    cat("Covariance information omitted (covarianceInformation='none'); standard errors are unavailable.\n")
  }

  if(identical(x$retainCoefficientMatrix, FALSE)){
    cat("Coefficient matrix C was not retained; post-fit C-based inference requires a refit.\n")
  }
  print(x$varcomp, digits = digits)
  if(!is.null(x$familyParameters)){
    cat("\nFamily parameters:\n")
    print(x$familyParameters, digits=digits, row.names=FALSE)
  }
  if(!is.null(x$zeroCoefficients)){
    cat("\nZero-component coefficients (logit scale):\n")
    print(x$zeroCoefficients, digits=digits, row.names=FALSE)
  }
  if(!is.null(x$thresholds)){
    cat("\nOrdinal cutpoints (link scale):\n")
    print(x$thresholds, digits=digits, row.names=FALSE)
  }
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nFixed effects:\n")
  if(!is.null(x$fixedRank)){
    cat("Rank-reduced fit: aliased coefficients are constrained to zero, not estimated; retained coefficients are conditional on these constraints.\n")
  }
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

print.mmes.glmm.laplace <- function(x, ...){
  print(summary.mmes(x), ...)
  invisible(x)
}

vcov.mmes.glmm.laplace <- function(object, ...){
  if(length(list(...))){
    stop("vcov() for native Laplace fits does not accept additional arguments.",
         call.=FALSE)
  }
  uncertainty <- .mmes_laplace_uncertainty(object)
  if(is.null(uncertainty$fixedCovariance)){
    stop("Fixed-effect covariance is unavailable: ", uncertainty$reason, ".",
         call.=FALSE)
  }
  uncertainty$fixedCovariance
}

.anova_laplace_wald <- function(object){
  covariance <- vcov.mmes.glmm.laplace(object)
  termLabels <- attr(object$fixedTerms, "term.labels")
  assignment <- object$fixedAssignment
  testNames <- c("(Intercept)", termLabels)
  testRows <- lapply(seq_along(testNames), function(termPosition){
    columns <- if(termPosition == 1L) which(assignment == 0L) else
      which(assignment == termPosition - 1L)
    if(!length(columns)) return(NULL)
    estimate <- as.numeric(object$b[columns])
    termCovariance <- covariance[columns, columns, drop=FALSE]
    statistic <- tryCatch(drop(crossprod(estimate,
      solve(termCovariance, estimate))), error=function(condition){
        stop("Wald test for ", testNames[termPosition],
          " failed because its covariance submatrix is singular: ",
          conditionMessage(condition), call.=FALSE)
      })
    degreesFreedom <- length(columns)
    data.frame(Df=degreesFreedom, Chisq=statistic,
      `Pr(>Chisq)`=stats::pchisq(statistic, degreesFreedom,
        lower.tail=FALSE), check.names=FALSE,
      row.names=testNames[termPosition])
  })
  testRows <- Filter(Negate(is.null), testRows)
  if(!length(testRows)){
    stop("No estimable fixed-effect terms are available for Wald tests.",
         call.=FALSE)
  }
  result <- do.call(rbind, testRows)
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- paste(
    "Term-wise Wald chi-square tests using the Laplace fixed-effect covariance approximation.")
  result
}

.anova_laplace_compare <- function(smaller, larger, allowBoundary=FALSE){
  if(!isTRUE(smaller$convergence) || !isTRUE(larger$convergence)){
    stop("Likelihood comparisons require both native Laplace fits to converge.",
         call.=FALSE)
  }
  if(!allowBoundary && (isTRUE(smaller$laplaceDiagnostics$boundary) ||
     isTRUE(larger$laplaceDiagnostics$boundary))){
    stop("Likelihood comparisons are unavailable when either fit has a flagged boundary parameter; variance-component boundary tests require parametric bootstrap or other nonstandard calibration.",
         call.=FALSE)
  }
  familySignature <- function(fit) lapply(fit$familySpecifications,
    function(specification) list(family=specification$family,
      link=specification$link, parameter=specification$parameter,
      free=specification$free,
      zeroDesign=specification$zeroDesign, zeroOffset=specification$zeroOffset,
      nominalReference=specification$reference,
      nominalLevels=specification$levels,
      ordinalLevels=fit$ordinal$levels,
      thresholds=if(specification$family == "ordinal")
        if(specification$free) length(specification$thresholds) else specification$thresholds,
      fixedValue=if(specification$free) NULL else specification$logParameter))
  covarianceSignature <- function(fit) lapply(fit$covStruct, function(structure){
    factorDefinitions <- lapply(structure$factors, function(factor){
      factor$par <- NULL
      factor$free <- NULL
      factor
    })
    list(dim=structure$dim, par_names=structure$par_names,
      free=structure$free,
      fixedValues=structure$par[!structure$free],
      factors=factorDefinitions)
  })
  sameLikelihood <- identical(smaller$y, larger$y) &&
    identical(smaller$offset, larger$offset) &&
    identical(smaller$laplaceFamilyIndex, larger$laplaceFamilyIndex) &&
    identical(smaller$priorWeights, larger$priorWeights) &&
    identical(familySignature(smaller), familySignature(larger)) &&
    identical(deparse(smaller$randomFormula), deparse(larger$randomFormula)) &&
    identical(smaller$randomDesign, larger$randomDesign) &&
    identical(smaller$randomRelationships, larger$randomRelationships) &&
    identical(covarianceSignature(smaller), covarianceSignature(larger))
  if(!sameLikelihood){
    stop("Native Laplace likelihood comparisons require identical response rows, offsets, family/trial settings, random formula, and covariance structure; only nested fixed effects may differ.",
         call.=FALSE)
  }
  smallerTerms <- attr(smaller$fixedTerms, "term.labels")
  largerTerms <- attr(larger$fixedTerms, "term.labels")
  if(length(smallerTerms) >= length(largerTerms) ||
     !all(smallerTerms %in% largerTerms)){
    stop("Models must be strictly nested by fixed-effect terms, with object as the smaller model and object2 as the larger model.",
         call.=FALSE)
  }
  combinedDesign <- cbind(larger$fixedDesign, smaller$fixedDesign)
  if(Matrix::rankMatrix(combinedDesign)[1L] !=
     Matrix::rankMatrix(larger$fixedDesign)[1L]){
    stop("The smaller fixed-effect design is not nested in the larger model's column space.",
         call.=FALSE)
  }
  degreesFreedom <- ncol(larger$fixedDesign) - ncol(smaller$fixedDesign)
  likelihoodRatio <- 2 * (larger$llik - smaller$llik)
  if(likelihoodRatio < -1e-7){
    stop("The larger fixed-effect model has a lower fitted likelihood; refit or check optimizer convergence.",
         call.=FALSE)
  }
  likelihoodRatio <- max(0, likelihoodRatio)
  modelLabels <- vapply(list(smaller$fixedFormula, larger$fixedFormula),
    function(formula) paste(deparse(formula), collapse=""), character(1))
  result <- data.frame(Model=modelLabels,
    Df=c(ncol(smaller$fixedDesign), ncol(larger$fixedDesign)),
    AIC=c(smaller$AIC, larger$AIC), BIC=c(smaller$BIC, larger$BIC),
    logLik=c(smaller$llik, larger$llik),
    Chisq=c(NA_real_, likelihoodRatio),
    `Chi Df`=c(NA_real_, degreesFreedom),
    `Pr(>Chisq)`=c(NA_real_, stats::pchisq(likelihoodRatio,
      degreesFreedom, lower.tail=FALSE)), check.names=FALSE,
    row.names=c("smaller", "larger"))
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- paste(
    "Likelihood-ratio test for nested fixed effects in native Laplace models.")
  result
}

"fitted.mmes" <- function(object, type=c("response", "link"), ...){

  if(inherits(object, "mmes.glmm.laplace")){
    if(.laplace_is_ordinal(object)){
      type <- if(missing(type)) "response" else
        match.arg(type, c("response", "probabilities", "mean", "link"))
      return(.laplace_ordinal_output(object, object$linear.predictors, type))
    }
    type <- match.arg(type)
    return(if(type == "link") object$linear.predictors else object$fitted.values)
  }
  if(inherits(object, "mmes.glmm")){
    type <- match.arg(type)
    if(type == "link") return(object$linear.predictors)
    return(object$fitted.values)
  }

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

"residuals.mmes" <- function(object,
                               type=c("response", "deviance", "pearson", "working",
                                 "quantile"), ...) {
  if(inherits(object, "mmes.glmm")){
    type <- match.arg(type)
    if(inherits(object, "mmes.glmm.laplace")){
      if(.laplace_is_ordinal(object)) return(.laplace_ordinal_residuals(object, type))
      if(.laplace_is_nominal(object)){
        if(!type %in% c("response", "pearson")){
          stop("Nominal Laplace residuals support type='response' and type='pearson'; deviance, working and quantile residuals are not defined by this implementation.",
            call.=FALSE)
        }
        observed <- stats::model.matrix(~object$y - 1L)
        colnames(observed) <- levels(object$y)
        residual <- observed - object$category.probabilities
        if(type == "response") return(residual)
        probability <- object$category.probabilities
        if(any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)){
          stop("Nominal Pearson residuals are undefined for numerically zero or one category probabilities.",
            call.=FALSE)
        }
        return(residual / sqrt(probability * (1 - probability)))
      }
      if(type == "quantile") return(.laplace_quantile_residuals(object))
      response <- as.numeric(object$y)
      mean <- object$fitted.values
      families <- if(is.null(object$laplaceFamilies)) list(object$family) else object$laplaceFamilies
      familyIndex <- if(is.null(object$laplaceFamilyIndex))
        rep(1L, length(response)) else object$laplaceFamilyIndex
      trials <- if(is.null(object$priorWeights)) rep(1, length(response)) else object$priorWeights
      dispersion <- if(is.null(object$dispersion)) rep(1, length(families)) else object$dispersion
      if(type == "response") return(response - mean)
      if(type == "pearson"){
        variance <- .pql_apply(families, familyIndex, "variance", mean) *
          dispersion[familyIndex] / trials
        if(any(!is.finite(variance)) || any(variance <= 0)){
          stop("Pearson residuals are undefined at zero conditional variances.",
            call.=FALSE)
        }
        return((response - mean) / sqrt(variance))
      }
      if(type == "working"){
        derivative <- .pql_apply(families, familyIndex, "mu.eta", object$linear.predictors)
        if(any(!is.finite(derivative)) || any(derivative <= 0)){
          stop("Working residuals are undefined at zero link derivatives.",
            call.=FALSE)
        }
        return((response - mean) / derivative)
      }
      if(type == "quantile"){
        stop("Randomized quantile residuals are currently implemented for native Laplace fits only.",
             call.=FALSE)
      }
      contribution <- .laplace_deviance(response, mean, families, familyIndex, trials, dispersion)
      return(sign(response - mean) * sqrt(pmax(contribution, 0)))
    }
    y <- as.numeric(object$y)
    mu <- object$fitted.values
    pw <- if(is.null(object$priorWeights)) rep(1, length(y)) else object$priorWeights
    families <- if(is.null(object$pqlFamilies)) list(object$family) else object$pqlFamilies
    index <- if(is.null(object$pqlFamilyIndex)) rep(1L, length(y)) else object$pqlFamilyIndex
    if(type == "response") return(y - mu)
    if(type == "working") return(as.numeric(object$workingResponse) -
                    object$linear.predictorsNoOffset)
    if(type == "pearson") return((y - mu) * sqrt(pw / .pql_apply(families, index, "variance", mu)))
    contribution <- .pql_apply(families, index, "dev.resids", y, mu, pw)
    return(sign(y - mu) * sqrt(pmax(contribution, 0)))
  }
  if(length(type) == 1L && identical(type, "quantile")){
    stop("Randomized quantile residuals are currently implemented for native Laplace fits only.",
         call.=FALSE)
  }
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

rotate_back_mmes <- function(object){
  if(!inherits(object, "mmes")){
    stop("object must inherit from class 'mmes'.", call.=FALSE)
  }
  if(is.null(object$rotation) || is.null(object$uListEngine)){
    stop("object was not fitted with rotation=TRUE.", call.=FALSE)
  }

  term <- object$rotation$term
  out <- object$rotation$vectors %*% object$uListEngine[[term]]
  rownames(out) <- object$rotation$levels
  colnames(out) <- colnames(object$uListEngine[[term]])
  out
}

#### =========== ####
## COEF FUNCTION ####
#### =========== ####

"coef.mmes" <- function(object, ...){
  if(inherits(object, "mmes.glmm.laplace")) return(object$b)
  object$b
}

"logLik.mmes" <- function(object, ...){
  if(!inherits(object, "mmes.glmm.laplace")){
    stop("A marginal log-likelihood is not available for PQL fits.",
         call.=FALSE)
  }
  structure(object$llik,
    df=object$npar,
    nobs=object$nobs, class="logLik")
}

"print.coef.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}

#### =========== ####################
## FACTOR-ANALYTIC LOADINGS/SCORES ##
#### =========== ####################

# Locate the single fam()/rrm()-shaped term in a fitted mmes object and
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
      stop("No fam()/rrm() covariance term was found in this model.", call.=FALSE)
    }
    if(length(candidates) > 1L){
      stop(
        paste0(
          "Multiple fam()/rrm() terms were found; please specify term as one of: ",
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
        "term '", term, "' is not a single fam()/rrm() covariance-shaping factor. ",
        "Terms combining fam()/rrm() with additional shaping factors are not yet supported."
      ),
      call.=FALSE
    )
  }

  list(term=term, pos=pos, factor=factors[[1]])
}

# Reconstruct the normalized loadings (Lambda) and specific variances (Psi)
# of a fam()/rrm() term such that sigma2*(Lambda %*% t(Lambda) + diag(Psi))
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

  specificRaw <- if(isTRUE(f$zero_specific)) rep(0, q) else rep(1, q)
  if(f$model == "fa" && !identical(f$variance, "homogeneous") && q > 1L){
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
    specific <- specific * sigma2
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

.covparams_term_names <- function(object){
  termNames <- names(object$covStruct)
  if(is.null(termNames) || any(!nzchar(termNames))){
    termNames <- names(object$covPar)
  }
  if(is.null(termNames) || length(termNames) != length(object$covStruct)){
    termNames <- paste0("structure", seq_along(object$covStruct))
  }
  termNames
}

.covparams_select_terms <- function(termNames, term=NULL){
  selected <- seq_along(termNames)
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
  selected
}

.covparams_term_report <- function(descriptor, current, termName){
  if(!length(current)) return(list())
  scale <- current[1L]
  factors <- descriptor$factors
  rows <- list()
  if(!length(factors)){
    return(list(list(term=termName, factor="sigma2", section=NA_character_,
                     parameter="sigma2", estimate=scale)))
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
      stop("Covariance factor ", j, " in term '", termName,
           "' has no valid native reporting callback.", call.=FALSE)
    }
    absorbScale <- !scaleAbsorbed
    values <- reporter$fun(scale=scale, par=factorPar, factor=factor,
                           absorb_scale=absorbScale)
    values <- unlist(values, use.names=TRUE)
    if(!length(values)) next
    if(is.null(names(values)) || any(!nzchar(names(values)))){
      stop("Native reporting callbacks must return named values.", call.=FALSE)
    }
    if(any(!is.finite(values))){
      stop("Native reporting callback returned non-finite values for term '",
           termName, "'.", call.=FALSE)
    }
    parameterNames <- names(values)
    section <- rep(NA_character_, length(values))
    if(!is.null(factor$section_levels) || isTRUE(factor$section_owner)){
      section <- sub("^(.*)\\[([^]]*)\\]$", "\\2", parameterNames)
      parameterNames <- sub("^(.*)\\[([^]]*)\\]$", "\\1", parameterNames)
    }
    factorName <- if(!is.null(factor$label)) factor$label else
      if(!is.null(factor$model) && nzchar(factor$model)) factor$model else
        paste0("factor", j)
    rows[[length(rows) + 1L]] <- list(
      term=rep(termName, length(values)), factor=rep(factorName, length(values)),
      section=section, parameter=parameterNames, estimate=as.numeric(values))
    scaleAbsorbed <- scaleAbsorbed || absorbScale
  }
  if(!scaleAbsorbed){
    rows[[length(rows) + 1L]] <- list(
      term=termName, factor="sigma2", section=NA_character_,
      parameter="sigma2", estimate=scale)
  }
  rows
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
  termNames <- .covparams_term_names(object)
  selected <- .covparams_select_terms(termNames, term)
  rows <- list()
  for(i in selected){
    reported <- .covparams_term_report(
      object$covStruct[[i]], as.numeric(object$covPar[[i]]), termNames[i])
    if(length(reported)) rows <- c(rows, reported)
  }
  if(!length(rows)){
    return(data.frame(term=character(), factor=character(), section=character(),
                      parameter=character(), estimate=numeric()))
  }
  columns <- lapply(names(rows[[1L]]), function(column){
    unlist(lapply(rows, `[[`, column), use.names=FALSE)
  })
  names(columns) <- names(rows[[1L]])
  out <- as.data.frame(columns, stringsAsFactors=FALSE)
  rownames(out) <- NULL
  .mmes_reml_report_status(object, out)
}

"covparams_mmes" <- function(object, term=NULL){
  .covparams_mmes(object, term)
}

.covparams_mmes_se <- function(object, term=NULL, rel_step=1e-6){
  if(inherits(object, "mmes.glmm.laplace")){
    if(length(rel_step) != 1L || !is.finite(rel_step) ||
       !identical(as.numeric(rel_step), 1e-6)){
      stop("rel_step applies to Gaussian covariance reporting only; native Laplace uncertainty uses the fitted joint numerical Hessian.",
        call.=FALSE)
    }
    uncertainty <- .mmes_laplace_uncertainty(object)
    availableTerms <- unique(object$covParNative$term)
    selected <- .covparams_select_terms(availableTerms, term)
    rows <- which(object$covParNative$term %in% availableTerms[selected])
    out <- uncertainty$covParNative[rows, , drop=FALSE]
    attr(out, "vcov") <- uncertainty$covParCovariance[rows, rows, drop=FALSE]
    return(out)
  }
  .mmes_require_covariance_information(object)
  if(length(rel_step) != 1L || !is.finite(rel_step) || rel_step <= 0){
    stop("rel_step must be one positive finite value.", call.=FALSE)
  }
  if(is.null(object$theta_se) || !is.matrix(object$theta_se)){
    stop("The fitted object does not contain a covariance-parameter uncertainty matrix.",
         call.=FALSE)
  }

  termNames <- .covparams_term_names(object)
  selected <- .covparams_select_terms(termNames, term)

  if(!length(selected)){
    out <- .covparams_mmes(object, term)
    out$StdError <- out$Zratio <- numeric()
    attr(out, "vcov") <- matrix(numeric(), 0L, 0L)
    return(out)
  }

  parameterCounts <- vapply(object$covPar, length, integer(1))
  totalParameters <- sum(parameterCounts)
  if(!all(dim(object$theta_se) == c(totalParameters, totalParameters))){
    stop("theta_se dimensions do not match the fitted covariance parameters.",
         call.=FALSE)
  }
  starts <- cumsum(c(1L, head(parameterCounts, -1L)))

  output <- vector("list", length(selected))
  jacobians <- vector("list", length(selected))
  globalIndices <- vector("list", length(selected))
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

    reportValues <- function(values){
        rows <- .covparams_term_report(object$covStruct[[i]], values, termNames[i])
        if(!length(rows)) return(list(estimate=numeric(), factor=character(),
                                      section=character(), parameter=character()))
        flatten <- function(column){
          unlist(lapply(rows, `[[`, column), use.names=FALSE)
        }
        list(estimate=flatten("estimate"), factor=flatten("factor"),
             section=flatten("section"), parameter=flatten("parameter"))
    }
    sameShape <- function(reported){
      length(reported$estimate) == nrow(base) &&
        identical(as.character(reported$factor), as.character(base$factor)) &&
        identical(as.character(reported$section), as.character(base$section)) &&
        identical(as.character(reported$parameter), as.character(base$parameter))
    }
    for(k in which(free)){
      h <- rel_step * max(1, abs(current[k]))
      plus <- current
      minus <- current
      plus[k] <- plus[k] + h
      minus[k] <- minus[k] - h
      plusValue <- try(reportValues(plus), silent=TRUE)
      minusValue <- try(reportValues(minus), silent=TRUE)
      for(value in list(plusValue, minusValue)){
        if(!inherits(value, "try-error") && !sameShape(value)){
          stop("Native reporting callback changed output shape in term '",
               termNames[i], "'.", call.=FALSE)
        }
      }
      plusOK <- !inherits(plusValue, "try-error") &&
        length(plusValue$estimate) == nrow(base) &&
        all(is.finite(plusValue$estimate))
      minusOK <- !inherits(minusValue, "try-error") &&
        length(minusValue$estimate) == nrow(base) &&
        all(is.finite(minusValue$estimate))
      if(plusOK && minusOK){
        jacobian[,k] <- (plusValue$estimate - minusValue$estimate) / (2*h)
      }else if(plusOK){
        jacobian[,k] <- (plusValue$estimate - base$estimate) / h
      }else if(minusOK){
        jacobian[,k] <- (base$estimate - minusValue$estimate) / h
      }else{
        stop("Unable to numerically differentiate native covariance parameter '",
             base$parameter[1L], "' in term '", termNames[i], "'.", call.=FALSE)
      }
    }

    global <- starts[i] + seq_len(nLocal) - 1L
    jacobians[[outputIndex]] <- jacobian
    globalIndices[[outputIndex]] <- global
    base$StdError <- rep(NA_real_, nrow(base))
    base$Zratio <- rep(NA_real_, nrow(base))
    output[[outputIndex]] <- base
  }

  allIndices <- unlist(globalIndices, use.names=FALSE)
  allJacobian <- matrix(0, sum(vapply(output, nrow, integer(1))),
                        length(allIndices))
  rowStart <- 1L
  columnStart <- 1L
  for(i in seq_along(output)){
    rows <- nrow(output[[i]])
    columns <- length(globalIndices[[i]])
    if(rows && columns){
      allJacobian[rowStart:(rowStart + rows - 1L),
                  columnStart:(columnStart + columns - 1L)] <- jacobians[[i]]
    }
    rowStart <- rowStart + rows
    columnStart <- columnStart + columns
  }
  freeAll <- unlist(lapply(selected, function(i){
    as.logical(object$covStruct[[i]]$free)
  }), use.names=FALSE)
  covariance <- object$theta_se[allIndices, allIndices, drop=FALSE]
  covariance[!freeAll,] <- 0
  covariance[,!freeAll] <- 0
  unavailable <- is.na(covariance)
  covariance[unavailable] <- 0
  nativeCovariance <- allJacobian %*% covariance %*% t(allJacobian)
  if(any(unavailable)){
    depends <- allJacobian != 0
    nativeCovariance[(depends %*% unavailable %*% t(depends)) > 0] <- NA_real_
  }
  standardErrors <- sqrt(pmax(diag(nativeCovariance), 0))
  estimates <- unlist(lapply(output, `[[`, "estimate"), use.names=FALSE)
  zRatios <- rep(NA_real_, length(standardErrors))
  positiveSE <- is.finite(standardErrors) & standardErrors > 0
  zRatios[positiveSE] <- estimates[positiveSE] / standardErrors[positiveSE]
  rowStart <- 1L
  for(i in seq_along(output)){
    rows <- nrow(output[[i]])
    indices <- rowStart + seq_len(rows) - 1L
    output[[i]]$StdError <- standardErrors[indices]
    output[[i]]$Zratio <- zRatios[indices]
    rowStart <- rowStart + rows
  }
  out <- do.call(rbind, output)
  rownames(out) <- NULL
  attr(out, "vcov") <- nativeCovariance
  
  .mmes_reml_report_status(object, out)
}

"covparams_mmes_se" <- function(object, term=NULL, rel_step=1e-6){
  .covparams_mmes_se(object, term, rel_step)
}

# Predict per-level latent factor scores for a fam()/rrm() term from its
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
        "scores_mmes() requires a random-effect fam()/rrm() term."
      ),
      call.=FALSE
    )
  }

  L <- fa$loadings
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

  if(all(fa$specific == 0)){
    if(method == "bartlett"){
      stop("Bartlett scores are undefined for rrm(): specific variances are fixed at zero. Use method='regression'.",
           call.=FALSE)
    }
    scores <- U %*% t(solve(crossprod(L), t(L)))
  }else if(method == "regression"){
    # Woodbury contraction needs only a rank-by-rank solve, including when
    # the fit deliberately retains no dense trait covariance.
    specific <- fa$specific * if(varianceScale) 1 else fa$sigma2
    weightedL <- sweep(L, 1L, specific, `/`)
    core <- diag(ncol(L)) + crossprod(L, weightedL) *
      if(varianceScale) 1 else fa$sigma2
    scores <- (U %*% weightedL) %*% solve(core)
  }else{
    weights <- 1 / (fa$sigma2 * fa$specific)
    weightedL <- sweep(L, 1L, weights, `*`)
    scores <- U %*% weightedL %*% solve(crossprod(L, weightedL))
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
  if(inherits(object, "mmes.glmm.laplace")){
    if(length(list(...))){
      stop("Native Laplace anova() does not accept additional arguments.",
           call.=FALSE)
    }
    if(is.null(object2)) return(.anova_laplace_wald(object))
    if(!inherits(object2, "mmes.glmm.laplace")){
      stop("Both models must be native Laplace fits.", call.=FALSE)
    }
    return(.anova_laplace_compare(object, object2))
  }
  if(is.null(object2)){
    return(wald_mmes(object, ...))
  }
  .mmes_require_full_rank_inference(object, "Likelihood-ratio comparisons")
  .mmes_require_full_rank_inference(object2, "Likelihood-ratio comparisons")
  if(inherits(object, "mmes.glmm") || inherits(object2, "mmes.glmm")){
    stop("Likelihood ratio tests are not available for PQL fits: their log-likelihood is that of a working Gaussian model, not of the GLMM.", call.=FALSE)
  }
  {
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

.mmes_hat_diagonal <- function(W, Ci){
  rowSums((W %*% Ci) * W)
}

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
    if(is.null(x$Ci)){
      stop("plot.mmes() currently requires the Henderson engine's coefficient-matrix inverse (Ci); it is not available for the direct-inversion engine (control$formulation='direct').", call.=FALSE)
    }
    leverage <- .mmes_hat_diagonal(x$W, x$Ci)
    plot(leverage, scale(rr), pch=20, col=transp("blue"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  # }
  #####################
  layout(matrix(1,1,1))
}
