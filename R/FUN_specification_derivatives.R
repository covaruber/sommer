.spec_derivative_layout <- function(spec, data){
  names <- names(spec$parameters)
  graphs <- names(spec$predictors)
  placeholders <- if(length(graphs))
    stats::setNames(paste0(".specificationGraph",seq_along(graphs)),graphs) else character()
  replace <- function(expression){
    if(!is.call(expression)) return(expression)
    if(identical(expression[[1L]],as.name("graph_term")))
      return(as.name(placeholders[[as.character(expression[[2L]])]]))
    as.call(lapply(as.list(expression),replace))
  }
  variables <- c(names,unname(placeholders))
  if(any(variables %in% spec$inputs))
    stop("Expression inputs conflict with internal graph derivative coordinates.",call.=FALSE)
  .nonlinear_workspace(as.double(nrow(data))*
    (4*length(variables)^2+4*length(names)^2+length(variables)+length(names)),
    spec$control$execution$denseMemoryMB,"Specification derivative workspace")
  list(variables=variables,placeholders=placeholders,
    expressions=lapply(spec$expressions,function(term)
      stats::deriv3(replace(term$expression),variables)),
    graphs=lapply(spec$predictors,.nonlinear_graph_native_layout,data=data))
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
  graphDerivatives <- lapply(names(spec$predictors),function(name){
    graph <- spec$predictors[[name]]
    native <- layout$graphs[[name]]
    coordinates <- match(paste0(name,"::",graph$parameters),names(values))
    rows <- do.call(cbind,values[coordinates])
    k <- length(coordinates)
    directions <- matrix(0,k,3L*n)
    evaluate <- function(directions) sommer_nonlinear_graph_eval_cpp(native$nodes,
      native$output,native$inputs,unname(rows),native$transformed,directions,
      context$control$execution$denseMemoryMB)
    result <- evaluate(directions)
    bindings[[layout$placeholders[[name]]]] <- as.numeric(result$output)
    gradient <- matrix(0,n,m)
    gradient[,coordinates] <- result$gradient
    hessian <- array(0,c(n,m,m))
    for(j in seq_len(k)){
      directions[j,seq.int(1L,3L*n,by=3L)] <- 1
      directional <- evaluate(directions)
      hessian[,coordinates,coordinates[j]] <- directional$hessianVector
      directions[j,] <- 0
    }
    list(gradient=gradient,hessian=hessian)
  })
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
    if(!length(graphDerivatives)){
      gradient <- gradient+sign*g
      hessian <- hessian+sign*h
    }else for(i in seq_len(n)){
      projection <- rbind(diag(m),do.call(rbind,lapply(graphDerivatives,function(g)
        g$gradient[i,,drop=FALSE])))
      gradient[i,] <- gradient[i,]+sign*as.numeric(g[i,,drop=FALSE] %*% projection)
      second <- crossprod(projection,matrix(h[i,,],length(layout$variables)) %*% projection)
      for(j in seq_along(graphDerivatives))
        second <- second+g[i,m+j]*matrix(graphDerivatives[[j]]$hessian[i,,],m)
      hessian[i,,] <- hessian[i,,]+sign*second
    }
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
