.laplace_covariance_map <- function(structures, terms, vcc=NULL){
  table <- .vc_param_table(structures, terms)
  if(!length(structures)){
    if(!is.null(vcc))
      stop("vcc requires at least one random-effect covariance parameter.",
        call.=FALSE)
    return(list(structures=structures,table=table,
      map=list(T=matrix(numeric(),0L,0L),offset=numeric(),free=integer(),
        initial=numeric(),constrained=FALSE)))
  }
  if(is.null(vcc)){
    free <- which(table$free)
    tangent <- diag(nrow(table))[, free, drop=FALSE]
    offset <- unlist(lapply(structures, `[[`, "par"), use.names=FALSE)
    offset[free] <- 0
    initial <- unlist(lapply(structures, function(structure)
      structure$par[structure$free]), use.names=FALSE)
  }else{
    constraints <- .vc_constraint_map(structures, table, vcc)
    structures <- constraints$covStruct
    table <- .vc_param_table(structures, terms)
    table$group <- constraints$group
    tangent <- constraints$T
    offset <- constraints$o
    working <- unlist(lapply(structures, `[[`, "par"), use.names=FALSE)
    offset[!table$free] <- working[!table$free]
    initial <- if(ncol(tangent)) as.numeric(qr.solve(tangent, working-offset)) else numeric()
    if(any(!is.finite(working)) ||
       max(abs(as.numeric(offset + tangent %*% initial)-working)) > 1e-8){
      stop("vcc does not define finite consistent native Laplace starting coordinates.",
        call.=FALSE)
    }
  }
  names(initial) <- if(ncol(tangent)) vapply(seq_len(ncol(tangent)), function(column){
    row <- which(tangent[, column] != 0)[1L]
    paste(table$term[row], table$parameter[row], sep=":")
  }, character(1)) else character()
  list(structures=structures, table=table,
    map=list(T=tangent, offset=offset, free=which(table$free),
      initial=initial, constrained=!is.null(vcc)))
}

.laplace_expand_parameters <- function(parameters, layout){
  map <- layout$covarianceMap
  if(is.null(map) || !isTRUE(map$constrained)) return(parameters)
  fixed <- seq_len(layout$fixedCount)
  covariance <- layout$fixedCount + seq_len(layout$covarianceCount)
  family <- seq.int(layout$fixedCount + layout$covarianceCount + 1L,
    length(parameters) + 1L)
  family <- family[family <= length(parameters)]
  working <- as.numeric(map$offset + map$T %*% parameters[covariance])
  c(parameters[fixed], working[map$free], parameters[family])
}

.laplace_constrained_evaluate <- function(parameters, layout, evaluate){
  result <- evaluate(.laplace_expand_parameters(parameters, layout))
  map <- layout$covarianceMap
  if(!is.null(map) && isTRUE(map$constrained) && !is.null(result$gradient)){
    fixed <- seq_len(layout$fixedCount)
    covariance <- layout$fixedCount + seq_along(map$free)
    family <- seq.int(layout$fixedCount + length(map$free) + 1L,
      length(result$gradient) + 1L)
    family <- family[family <= length(result$gradient)]
    result$gradient <- c(result$gradient[fixed],
      as.numeric(crossprod(map$T[map$free, , drop=FALSE],
        result$gradient[covariance])), result$gradient[family])
  }
  result
}

.laplace_proposal_structures <- function(structures, parameters, layout){
  expanded <- .laplace_expand_parameters(parameters, layout)
  for(term in seq_along(structures)){
    positions <- layout$freePositions[[term]]
    if(length(positions)){
      structures[[term]]$par[positions] <- expanded[layout$fixedCount +
        layout$parameterStarts[term] + seq_along(positions)]
    }
  }
  structures
}

.laplace_covariance_signature <- function(object){
  map <- object$laplaceParameterLayout$covarianceMap
  if(is.null(map) || !isTRUE(map$constrained)) return(NULL)
  list(T=unname(map$T), offset=unname(map$offset), free=map$free)
}
