.nonlinear_joint_formula <- function(random, model, data){
  if(length(model$random)>1L || !is.null(model$joint)){
    formulas <- c(Filter(Negate(is.null),list(random)),unname(model$random))
    if(!is.null(model$joint)){
      formula <- model$joint$random
      reference <- intersect(all.vars(formula),names(data))
      if(!length(reference) || ".parameter" %in% names(data))
        stop("Joint curves require a data grouping variable and reserve .parameter.",
          call.=FALSE)
      replace <- function(value){
        if(identical(value,as.name(".parameter")))
          return(as.call(list(as.name(".nonlinear_parameter_axis"),
            as.name(reference[1L]),model$joint$parameters)))
        if(is.call(value)) return(as.call(lapply(as.list(value),replace)))
        value
      }
      formula[[2L]] <- replace(formula[[2L]])
      formulas <- c(formulas,list(formula))
    }
    env <- new.env(parent=environment(formulas[[1L]]))
    env$.nonlinear_parameter_axis <- .nonlinear_parameter_axis
    for(formula in formulas){
      for(variable in setdiff(all.vars(formula),c(names(data),"units"))){
        source <- environment(formula)
        if(exists(variable,source,inherits=TRUE)){
          value <- get(variable,source,inherits=TRUE)
          if(exists(variable,env,inherits=FALSE) &&
             !identical(get(variable,env),value))
            stop("Conflicting random-formula bindings for '",variable,"'.",call.=FALSE)
          assign(variable,value,env)
        }
      }
    }
    rhs <- Reduce(function(a,b) call("+",a,b),lapply(formulas,`[[`,2L))
    return(stats::as.formula(call("~",rhs),env=env))
  }
  if(!length(model$random)) return(random)
  curve <- model$random[[1L]]
  env <- new.env(parent=environment(curve))
  for(formula in Filter(Negate(is.null), list(random, curve))){
    source <- environment(formula)
    for(variable in setdiff(all.vars(formula), c(names(data), "units"))){
      if(exists(variable, source, inherits=TRUE)){
        value <- get(variable, source, inherits=TRUE)
        if(exists(variable, env, inherits=FALSE) &&
           !identical(get(variable, env), value)){
          stop("Conflicting random-formula bindings for '", variable, "'.", call.=FALSE)
        }
        assign(variable, value, env)
      }
    }
  }
  rhs <- if(is.null(random)) curve[[2L]] else call("+", random[[2L]], curve[[2L]])
  stats::as.formula(call("~", rhs), env=env)
}

.nonlinear_random_layout <- function(model, prepared){
  if(length(model$random)>1L || !is.null(model$joint))
    return(.nonlinear_multiple_layout(model,prepared))
  if(!length(model$random)) return(NULL)
  term <- length(prepared$rtermss)
  blocks <- which(prepared$Zind == term)
  design <- do.call(cbind, prepared$Z[blocks])
  entries <- Matrix::summary(as(as(design, "generalMatrix"), "CsparseMatrix"))
  if(nrow(entries) != nrow(design) || any(entries$x != 1) ||
     any(tabulate(entries$i, nrow(design)) != 1L)){
    stop("A random curve parameter requires one unit incidence per observation; slopes and overlapping loadings are not supported.",
      call.=FALSE)
  }
  widths <- vapply(prepared$Z, ncol, integer(1))
  first <- sum(widths[seq_len(min(blocks)-1L)])
  list(version=1L, parameter=names(model$random),
    coordinate=match(names(model$random), model$parameters),
    scale=if(model$positive[match(names(model$random), model$parameters)]) "log" else "identity",
    term=term, blocks=blocks, start=first, count=ncol(design),
    columns=colnames(design), design=design)
}

.nonlinear_random_parameters <- function(model, coordinates, layout, mode){
  if(is.null(layout)) return(NULL)
  positive <- model$positive
  if(identical(layout$version,2L)){
    values <- matrix(rep(as.numeric(coordinates),each=nrow(layout$groupMap)),
      nrow=nrow(layout$groupMap),dimnames=dimnames(layout$groupMap))
    selected <- layout$groupMap>0L
    values[selected] <- values[selected]+mode[layout$groupMap[selected]]
    values[,positive] <- exp(values[,positive,drop=FALSE])
    return(values)
  }
  values <- matrix(rep(as.numeric(coordinates), each=layout$count),
    nrow=layout$count, dimnames=list(layout$columns, model$parameters))
  values[, layout$coordinate] <- values[, layout$coordinate]+
    mode[layout$start+seq_len(layout$count)]
  values[, positive] <- exp(values[, positive, drop=FALSE])
  values
}

.nonlinear_parameter_axis <- function(reference, parameters){
  n <- if(is.null(dim(reference))) length(reference) else nrow(reference)
  Matrix::Matrix(1,nrow=n,ncol=length(parameters),sparse=TRUE,
    dimnames=list(NULL,parameters))
}

.nonlinear_term_labels <- function(model, terms){
  if(length(model$random)<=1L && is.null(model$joint)) return(terms)
  labels <- names(model$random)
  if(!is.null(model$joint))
    labels <- c(labels,paste(model$joint$parameters,collapse=","))
  selected <- length(terms)-length(labels)+seq_along(labels)
  if(!is.null(model$joint))
    terms[tail(selected,1L)] <- paste(deparse(model$joint$random[[2L]]),collapse="")
  terms[selected] <- paste0("curve[",labels,"]:",terms[selected])
  terms
}

.nonlinear_multiple_layout <- function(model,prepared){
  terms <- length(prepared$rtermss)
  number <- length(model$random)+as.integer(!is.null(model$joint))
  selected <- terms-number+seq_len(number)
  widths <- vapply(prepared$Z,ncol,integer(1))
  starts <- cumsum(c(0L,head(widths,-1L)))
  map <- matrix(0L,sum(widths),length(model$parameters))
  groups <- vector("list",length(model$parameters))
  names(groups) <- model$parameters
  for(index in seq_along(selected)){
    term <- selected[index]
    blocks <- which(prepared$Zind==term)
    design <- do.call(cbind,prepared$Z[blocks])
    columns <- unlist(lapply(blocks,function(b) starts[b]+seq_len(widths[b])))
    parameters <- if(index<=length(model$random)) names(model$random)[index] else
      model$joint$parameters
    factors <- prepared$covStruct[[term]]$factors
    labelAxes <- lapply(factors,function(f) as.character(f$levels))
    found <- integer()
    axis <- rep(1L,ncol(design))
    if(length(parameters)>1L){
      arguments <- as.list(model$joint$random[[2L]])[-1L]
      if(!is.null(names(arguments))) arguments <- arguments[!names(arguments) %in%
        c("Gu","sigma2","fixedSigma2","rotation","isFixed","verbose")]
      found <- which(vapply(arguments,function(value)
        ".parameter" %in% all.vars(value),logical(1)))
      if(length(found)!=1L || found>length(factors) ||
         !identical(as.character(factors[[found]]$levels),parameters))
        stop("Joint .parameter must retain the declared covariance-axis order.",
          call.=FALSE)
      dims <- c(vapply(factors,function(f) as.integer(f$dim),integer(1)),
        nrow(prepared$Ai[[term]]))
      axis <- rep(rep(seq_along(parameters),each=prod(dims[seq_along(dims)>found])),
        times=prod(dims[seq_along(dims)<found]))
    }
    if(length(found)) labelAxes <- labelAxes[-found]
    labelAxes <- c(labelAxes,list(rownames(prepared$Ai[[term]])))
    labels <- Reduce(function(a,b) as.vector(t(outer(a,b,paste,sep=":"))),labelAxes)
    if(anyDuplicated(labels))
      stop("Curve group labels must identify unique covariance coordinates.",call.=FALSE)
    for(j in seq_along(parameters)){
      local <- which(axis==j)
      incidence <- design[,local,drop=FALSE]
      entries <- Matrix::summary(as(as(incidence,"generalMatrix"),"CsparseMatrix"))
      if(nrow(entries)!=nrow(incidence) || any(entries$x!=1) ||
         any(tabulate(entries$i,nrow(incidence))!=1L))
        stop("Each random curve parameter requires one unit incidence per observation.",
          call.=FALSE)
      coordinate <- match(parameters[j],model$parameters)
      map[columns[local],coordinate] <- 1L
      groups[[coordinate]] <- stats::setNames(columns[local],labels)
    }
  }
  labels <- unique(unlist(lapply(groups,names),use.names=FALSE))
  groupMap <- matrix(0L,length(labels),length(model$parameters),
    dimnames=list(labels,model$parameters))
  columns <- character(nrow(map))
  for(j in seq_along(groups)){
    if(length(groups[[j]])){
      groupMap[names(groups[[j]]),j] <- groups[[j]]
      columns[groups[[j]]] <- paste0(model$parameters[j],":",names(groups[[j]]))
    }
  }
  first <- min(which(rowSums(map)>0))-1L
  list(version=2L,parameter=model$parameters[colSums(map)>0],
    coordinate=which(colSums(map)>0),start=first,count=nrow(map)-first,
    columns=columns[first+seq_len(nrow(map)-first)],
    coordinateMap=map,groupMap=groupMap)
}

.nonlinear_curve_bases <- function(object,design){
  layout <- object$nonlinearRandom
  k <- length(object$nonlinearCoordinates)
  lapply(seq_len(k),function(j){
    result <- design*0
    columns <- if(is.null(layout)) integer() else if(identical(layout$version,2L))
      which(layout$coordinateMap[,j]!=0) else
        if(j==layout$coordinate) layout$start+seq_len(layout$count) else integer()
    if(length(columns)) result[,columns] <- design[,columns,drop=FALSE]
    result
  })
}

.nonlinear_curve_value <- function(object,input,design,mode){
  bases <- .nonlinear_curve_bases(object,design)
  deviations <- do.call(cbind,lapply(bases,function(b) as.numeric(b%*%mode)))
  sommer_nonlinear_curve_joint_eval(object$nonlinear$curve,input,
    object$nonlinearCoordinates,deviations)
}

.nonlinear_has_reduced_rank <- function(structures){
  any(vapply(structures,function(structure)
    any(vapply(structure$factors,function(factor)
      isTRUE(factor$zero_specific) ||
        isTRUE(factor$representation$loadings$zero_specific),logical(1))),
    logical(1)))
}
