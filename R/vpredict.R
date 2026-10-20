vpredict <- function(object, transform) {
  UseMethod("vpredict")
}

vpredict.mmes <- function(object, transform){
  .mmes_require_covariance_information(object)
  pframe <- as.list(object$covParNative$estimate)
  names(pframe) <- paste0("V", seq_along(pframe))
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe)
  gradient <- as.vector(attr(tvalue, "gradient"))
  active <- which(is.na(gradient) | gradient != 0)
  Vmat <- attr(object$covParNative, "vcov")
  se <- if(length(active)){
    sqrt(abs(drop(crossprod(gradient[active],
      Vmat[active, active, drop=FALSE] %*% gradient[active]))))
  }else 0
  out <- data.frame(Estimate=tvalue, SE=se)
  rownames(out) <- if(length(transform) == 3L) transform[[2L]] else ""
  class(out) <- c("vpredict.mmes", "data.frame")
  out
}