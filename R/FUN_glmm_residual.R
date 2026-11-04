.laplace_prepare_gaussian_residual <- function(fixed, data, family, rcov, W,
                                               weights, naMethodY, mmesArgs){
  structured <- !is.null(rcov) &&
    !identical(attr(stats::terms(rcov),"term.labels"),"units")
  if(!structured && is.null(W) && is.null(weights)) return(NULL)
  families <- if(inherits(family,"sommer_familym")) family$families else list(family)
  gaussian <- which(vapply(families,function(x)
    x$family == "gaussian" && x$link == "identity",logical(1)))
  if(!length(gaussian))
    stop("Structured rcov, W or weights require Gaussian observation blocks; non-Gaussian residual dependence is unsupported.",call.=FALSE)
  if(is.null(data) || !is.data.frame(data))
    stop("Gaussian observation residual blocks require a data frame.",call.=FALSE)
  frame <- stats::model.frame(fixed,data=data,na.action=stats::na.pass)
  response <- stats::model.response(frame)
  if(!is.numeric(response) || !is.null(dim(response)))
    stop("Gaussian observation blocks require a scalar numeric response in long format.",call.=FALSE)
  rows <- which(.laplace_family_index(family,data,nrow(data)) %in% gaussian)
  gaussianData <- droplevels(data[rows,,drop=FALSE])
  if(".gaussianResponse" %in% names(data))
    stop("Gaussian residual preparation reserves .gaussianResponse.",call.=FALSE)
  gaussianData$.gaussianResponse <- response[rows]
  args <- mmesArgs
  args$observationMask <- NULL
  prepared <- .mmes_prepare_model(.gaussianResponse~1,NULL,rcov,gaussianData,
    stats::gaussian(),naMethodY,args)
  keep <- rep(TRUE,nrow(data))
  keep[rows] <- prepared$obsInfo$included
  list(prepared=prepared,rows=rows,keep=keep,originalData=data)
}

.laplace_finish_gaussian_residual <- function(preparation, prepared, index,
                                               specifications, W, weights){
  if(is.null(preparation)) return(NULL)
  gaussian <- which(vapply(specifications,function(x) x$family=="gaussian",logical(1)))
  rows <- which(index %in% gaussian)
  if(!length(rows)) stop("No Gaussian observations remain for residual blocks.",call.=FALSE)
  original <- which(prepared$obsInfo$included)
  if(!is.null(prepared$categorical))
    original <- prepared$categorical$originalRows[prepared$categorical$rowMap]
  gaussianOriginal <- original[rows]
  available <- preparation$rows[preparation$prepared$obsInfo$included]
  select <- match(gaussianOriginal,available)
  if(anyNA(select)) stop("Gaussian residual row mapping is inconsistent.",call.=FALSE)
  descriptor <- preparation$prepared$covStruct[[1L]]
  block <- preparation$prepared$residualBlock[select]
  local <- preparation$prepared$residualIndex[select]
  originalCount <- nrow(preparation$originalData)
  if(!is.null(W)){
    if(length(dim(W))!=2L || nrow(W) != ncol(W))
      stop("Gaussian observation W must be square.",call.=FALSE)
    if(!isSymmetric(W) || any(!is.finite(W)))
      stop("Gaussian observation W must be finite and symmetric positive definite.",call.=FALSE)
    sources <- list(seq_len(originalCount),preparation$rows,
      sort(unique(original)),gaussianOriginal)
    sources <- unique(Filter(function(source) length(source)==nrow(W),sources))
    if(!length(sources))
      stop("W must match original or retained observations, or the original or retained Gaussian subset.",call.=FALSE)
    if(length(sources)>1L)
      stop("W dimensions ambiguously match different observation subsets; supply W for all original observations.",call.=FALSE)
    source <- sources[[1L]]
    nonGaussian <- which(!source %in% preparation$rows)
    entries <- Matrix::summary(to_sparse(W))
    outside <- entries$i %in% nonGaussian | entries$j %in% nonGaussian
    if(any(outside & (entries$i!=entries$j | abs(entries$x-1)>1e-12)) ||
       any(abs(diag(W)[nonGaussian]-1)>1e-12))
      stop("W must be identity on non-Gaussian observations and have no cross-family entries.",call.=FALSE)
    order <- match(gaussianOriginal,source)
    W <- W[order,order,drop=FALSE]
  }
  groups <- NULL
  if(!is.null(weights)){
    if(!inherits(weights,"formula") || length(weights)!=2L)
      stop("weights must be a one-sided formula describing independent blocks of W.",call.=FALSE)
    frame <- stats::model.frame(weights,preparation$originalData[gaussianOriginal,,drop=FALSE],
      na.action=stats::na.pass)
    if(!ncol(frame) || anyNA(frame))
      stop("The weights block formula must define nonmissing Gaussian row groups.",call.=FALSE)
    groups <- as.integer(do.call(interaction,c(frame,list(drop=TRUE,lex.order=TRUE))))
  }
  weight <- .mmes_prepare_weights(W,is.null(W),length(rows),rep(TRUE,length(rows)),
    length(rows),groups)
  list(rows=as.integer(rows),descriptor=descriptor,block=as.integer(block),
    index=as.integer(local),W=weight$W,weighted=weight$useH,weightBlocks=groups,
    originalRows=gaussianOriginal)
}

.laplace_add_residual_map <- function(mapping,model){
  if(is.null(model)) return(mapping)
  residual <- .laplace_covariance_map(list(model$descriptor),"residual")
  first <- mapping$map
  second <- residual$map
  mapping$structures <- c(mapping$structures,residual$structures)
  table <- mapping$table
  addition <- residual$table
  addition$structure <- addition$structure+length(mapping$structures)-1L
  if("group" %in% names(table)) addition$group <- NA_integer_
  mapping$table <- rbind(table,addition)
  mapping$table$index <- seq_len(nrow(mapping$table))
  mapping$map <- list(T=as.matrix(Matrix::bdiag(first$T,second$T)),
    offset=c(first$offset,second$offset),
    free=c(first$free,length(first$offset)+second$free),
    initial=c(first$initial,second$initial),constrained=TRUE)
  mapping
}

.laplace_covariance_structures <- function(object){
  c(object$covStruct,if(!is.null(object$laplaceResidual)) list(object$laplaceResidual$descriptor))
}

.laplace_covariance_terms <- function(object){
  c(object$randomMetadata,if(!is.null(object$laplaceResidual)) "residual")
}

.laplace_residual_matrix <- function(object,maxMemoryMB=NULL){
  if(!is.null(object$laplaceResidualCovariance)) return(object$laplaceResidualCovariance)
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  if(is.null(maxMemoryMB)) maxMemoryMB <- object$laplacePostPEVSettings$denseMemoryMB
  if(is.null(maxMemoryMB)) maxMemoryMB <- 1000
  sommer_glmm_residual_covariance(object$laplaceResidual,maxMemoryMB)
}

.laplace_residual_signature <- function(object){
  model <- object$laplaceResidual
  if(is.null(model)) return(NULL)
  descriptor <- model$descriptor
  list(rows=model$rows,block=model$block,index=model$index,W=model$W,
    free=descriptor$free,fixed=descriptor$par[!descriptor$free],
    factors=lapply(descriptor$factors,function(factor){
      factor$par <- NULL
      factor$native_report <- NULL
      factor
    }))
}
