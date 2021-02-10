vpredict <- function(object, transform) {
  UseMethod("vpredict")
}

vpredict.mmer <- function (object, transform){
  
  # UseMethod("vpredict.mmer")
  
  # if(object$method %in% c("EMMA","EM")){
  #   stop("The pin function only works for 'NR' and 'AI' methods.",call. = FALSE)
  # }
  # 
  # if(object$method %in% c("MNR","MAI","MEMMA")){
  #   pframe <- as.list(object$sigma)#as.list(summary(object)[[3]][,1])
  # }else{
  #   pframe <- as.list(object$var.comp[,1])
  # }
  pframe <- as.list(summary(object)$varcomp[,1])
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  ## deriv creates a derivative for a simple expression
  # i.e. dx2x <- deriv(~ x^2, "x") ; dx2x
  ## eval evaluates a derivative expression for certain values specified in the deriv() expresion
  
  ## 1) get derivatives of the expression provided, 
  ## i.e. V1/(V1+V4) with respect to each of the variables existing, i.e. v1, v2, v3, v4
  ## 2) Evaluate the expression (derivatives) for values provided for each variable
  ## i.e. with the actual values of the variance components
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), 
                 pframe)
  X <- as.vector(attr(tvalue, "gradient")) # just make it a sinpe vector of derivatives
  tname <- if (length(transform) == 3) 
    transform[[2]]
  else ""
  n <- length(pframe) ## number of parameters available, i.e. V1,V2,V3,V4
  i <- rep(1:n, 1:n) ## repeat each parameter by its own
  j <- sequence(1:n) ## makes a sequence from 1 to the number provided, i.e. if sequence(1:2) = 1 1 2, because it makes the sequence for 1:1 and then 1:2
  k <- 1 + (i > j) # all where i <= j get a 1, all i > j get a 2
  Vmat <- object$sigmaSE
  toext <- upper.tri(Vmat)
  diag(toext) <- TRUE
  Vmat <- Vmat[which(toext,arr.ind = TRUE)] ## extract the upper triangular
  se <- sqrt(abs(sum(Vmat * X[i] * X[j] * k))) ## only X[i] and X[j] that match with the Vi's indicated are != than zero
  
  ### Vmat are the second derivatives
  ### X are the derivatives of the expression with respect to each parameter of interest
  ## X[i] * X[j] * k  multiplies once the var(var.i) and var(var.j) and twice the covar(covar.ij)
  ## then takes the sqrt(abs( sum[c(var(i),covar(i,j),var(j)] ))) , that's the SE
  ## Vmat * X[i] * X[j] * k --> 2nd derivatives * derivatives of i.e. h2 with respect to each term
  ## d''(x) * d'(x) * d'
  ## those var(vc.i) and covar(covar.ij) from the variance comp. come from the inverse if the second derivatives (Fisher's)
  toreturn2 <- data.frame(row.names = tname, Estimate = tvalue, SE = se)
  class(toreturn2) <- "vpredict.mmer"
  # attr(toreturn2, "class")<-c("vpredict.mmer", "data.frame")
  return(toreturn2)
}