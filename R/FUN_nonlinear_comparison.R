.nonlinear_information_criterion <- function(object, others, k){
  fits <- c(list(object),others)
  if(!all(vapply(fits,inherits,logical(1),"mmes.nonlinear")))
    stop("Nonlinear information criteria require nonlinear mmes fits.",call.=FALSE)
  if(!is.null(k) && (length(k)!=1L || !is.finite(k) || k<0))
    stop("k must be one finite non-negative value.",call.=FALSE)
  value <- vapply(fits,function(fit)
    -2*as.numeric(stats::logLik(fit))+
      (if(is.null(k)) log(stats::nobs(fit)) else k)*attr(stats::logLik(fit),"df"),numeric(1))
  if(length(fits)==1L) return(unname(value))
  result <- data.frame(df=vapply(fits,function(fit) fit$npar,numeric(1)),value=value)
  names(result)[2L] <- if(is.null(k)) "BIC" else "AIC"
  result
}

.nonlinear_wald <- function(object,parm,null){
  if(is.null(parm) || is.null(null))
    stop("Nonlinear Wald tests require explicit parm and null on working coordinates; automatic nonlinear-term or boundary tests are unsupported.",
      call.=FALSE)
  available <- names(stats::coef(object))
  if(!is.character(parm) || !length(parm) || anyNA(parm) ||
     anyDuplicated(parm) || any(!parm %in% available) ||
     !is.numeric(null) || length(null)!=length(parm) || any(!is.finite(null)))
    stop("Select unique fixed/curve parameters and supply one finite working-coordinate null per parameter.",
      call.=FALSE)
  if(!is.null(names(null))){
    if(anyNA(names(null)) || anyDuplicated(names(null)) || !setequal(names(null),parm))
      stop("Named null values must match parm exactly.",call.=FALSE)
    null <- null[parm]
  }
  covariance <- stats::vcov(object,type="working")[parm,parm,drop=FALSE]
  difference <- object$optimizer$par[parm]-null
  statistic <- as.numeric(crossprod(difference,solve(covariance,difference)))
  result <- data.frame(Df=length(parm),Chisq=statistic,
    `Pr(>Chisq)`=stats::pchisq(statistic,length(parm),lower.tail=FALSE),
    check.names=FALSE)
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- paste("Approximate nonlinear Wald test on working coordinates:",
    paste(parm,collapse=", "),"; finite interior nulls only.")
  result
}

.nonlinear_compare <- function(smaller,larger){
  if(!inherits(larger,"mmes.nonlinear"))
    stop("Both models must be nonlinear mmes fits.",call.=FALSE)
  if(!identical(smaller$input,larger$input) ||
     !identical(smaller$nonlinear,larger$nonlinear) ||
     !identical(smaller$likelihoodMethod,larger$likelihoodMethod))
    stop("Nonlinear comparisons require identical curves, curve-random specifications, inputs and likelihood methods; only nested linear fixed effects may differ.",
      call.=FALSE)
  .nonlinear_uncertainty(smaller)
  .nonlinear_uncertainty(larger)
  proxy <- function(object){
    object <- .nonlinear_inference_proxy(object)
    object$laplaceParameterLayout <- object$nonlinearParameterLayout$covariance
    object$laplaceDiagnostics <- list(boundary=FALSE)
    object$AIC <- stats::AIC(object)
    object$BIC <- stats::BIC(object)
    object
  }
  result <- .anova_laplace_compare(proxy(smaller),proxy(larger))
  result$Df <- c(smaller$npar,larger$npar)
  attr(result,"heading") <- paste("Approximate likelihood-ratio test for nested linear",
    "fixed effects with the same nonlinear curve and random structure.")
  result
}
