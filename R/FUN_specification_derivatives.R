.spec_derivative_layout <- function(spec, data){
  names <- names(spec$parameters)
  variables <- names
  if(any(variables %in% spec$inputs))
    stop("Expression inputs conflict with parameter names.",call.=FALSE)
  .nonlinear_workspace(as.double(nrow(data))*
    (4*length(variables)^2+4*length(names)^2+length(variables)+length(names)),
    spec$control$execution$denseMemoryMB,"Specification derivative workspace")
  list(variables=variables,
    expressions=lapply(spec$expressions,function(term)
      stats::deriv3(term$expression,variables)))
}

.spec_predictor_derivatives <- function(parameters, mode, context){
  spec <- context$spec
  n <- context$n
  m <- length(spec$parameters)
  layout <- context$derivativeLayout
  coefficients <- parameters[seq_len(context$p)]
  splitMode <- .spec_split_mode(mode,context$widths)
  values <- list()
  first <- matrix(1,n,m)
  for(j in seq_len(m)){
    name <- names(spec$parameters)[j]
    declaration <- spec$parameters[[name]]
    eta <- as.numeric(context$design$parameters[[name]]$x %*%
      coefficients[declaration$positions])+context$design$parameters[[name]]$offset
    record <- declaration$randomIndex
    if(length(record))
      eta <- eta+as.numeric(context$randomDesign[[record]] %*% splitMode[[record]])
    values[[name]] <- if(declaration$link=="log") exp(eta) else eta
    if(declaration$link=="log") first[,j] <- values[[name]]
  }
  bindings <- new.env(parent=baseenv())
  for(name in spec$inputs) bindings[[name]] <- context$data[[name]]
  for(name in names(values)) bindings[[name]] <- values[[name]]
  value <- as.numeric(context$design$fixed$x %*%
    coefficients[spec$fixedPositions])+context$design$fixed$offset
  gradient <- matrix(0,n,m)
  hessian <- array(0,c(n,m,m))
  for(term in seq_along(spec$expressions)){
    expression <- eval(layout$expressions[[term]],bindings)
    sign <- spec$expressions[[term]]$sign
    value <- value+sign*as.numeric(expression)
    g <- attr(expression,"gradient")
    h <- attr(expression,"hessian")
    gradient <- gradient+sign*g
    hessian <- hessian+sign*h
  }
  jacobian <- matrix(0,n,context$q)
  starts <- cumsum(c(0L,head(context$widths,-1L)))
  columns <- lapply(seq_along(context$widths),function(i)
    starts[i]+seq_len(context$widths[i]))
  if(length(spec$additiveIndex)){
    record <- spec$additiveIndex
    value <- value+as.numeric(context$randomDesign[[record]] %*% splitMode[[record]])
    jacobian[,columns[[record]]] <- context$randomDesign[[record]]
  }
  for(j in seq_len(m)){
    record <- spec$parameters[[j]]$randomIndex
    if(length(record))
      jacobian[,columns[[record]]] <- jacobian[,columns[[record]],drop=FALSE]+
        context$randomDesign[[record]]*(gradient[,j]*first[,j])
  }
  list(eta=value,jacobian=jacobian,gradient=gradient,hessian=hessian,
    first=first,columns=columns)
}

.spec_conditional_hessian <- function(evaluation, family, precision, context){
  curvature <- family$blockCurvature
  if(is.null(curvature)) curvature <- Matrix::Diagonal(x=as.numeric(family$curvature))
  result <- precision+crossprod(evaluation$jacobian,curvature %*% evaluation$jacobian)
  parameters <- context$spec$parameters
  for(j in seq_along(parameters)){
    left <- parameters[[j]]$randomIndex
    if(!length(left)) next
    for(k in seq_along(parameters)){
      right <- parameters[[k]]$randomIndex
      if(!length(right)) next
      second <- evaluation$hessian[,j,k]*evaluation$first[,j]*evaluation$first[,k]
      if(j==k && parameters[[j]]$link=="log")
        second <- second+evaluation$gradient[,j]*evaluation$first[,j]
      result[evaluation$columns[[left]],evaluation$columns[[right]]] <-
        result[evaluation$columns[[left]],evaluation$columns[[right]],drop=FALSE]-
          crossprod(context$randomDesign[[left]],
            context$randomDesign[[right]]*(as.numeric(family$score)*second))
    }
  }
  as.matrix((result+t(result))/2)
}
