.mmes_reml_omit <- function(X, Z, Zind, Ai, covStruct, terms, termNames){
  original <- list(Z=Z, covStruct=covStruct, terms=terms)
  mappings <- vector("list", length(terms))
  omitted <- list()
  fixedQR <- NULL
  inFixedSpace <- function(designs){
    if(!ncol(X)) return(all(vapply(designs, function(z) Matrix::nnzero(z) == 0L,
                                  logical(1))))
    if(is.null(fixedQR)) fixedQR <<- qr(X)
    for(z in designs){
      nonzero <- which(Matrix::colSums(abs(z)) > 0)
      if(!length(nonzero)) next
      for(first in seq.int(1L, length(nonzero), by=32L)){
        columns <- nonzero[seq.int(first, min(first + 31L, length(nonzero)))]
        panel <- as.matrix(z[, columns, drop=FALSE])
        residual <- qr.resid(fixedQR, panel)
        relative <- sqrt(colSums(residual^2) / colSums(panel^2))
        if(any(!is.finite(relative)) || any(relative > 1e-12)) return(FALSE)
      }
    }
    TRUE
  }
  activeTerms <- integer()
  newZ <- list()
  newIndex <- integer()
  newAi <- list()
  newCov <- list()
  newNames <- list()
  for(i in seq_along(terms)){
    cs <- covStruct[[i]]
    blocks <- which(Zind == i)
    keep <- seq_along(blocks)
    diagonal <- length(cs$factors) == 1L &&
      identical(cs$factors[[1L]]$model, "diag") &&
      is.null(cs$factors[[1L]]$section_levels) && all(cs$free)
    scalar <- isTRUE(cs$free[1L]) && !any(cs$free[-1L])
    if(diagonal){
      null <- vapply(Z[blocks], function(z) inFixedSpace(list(z)), logical(1))
      keep <- which(!null)
    }else if(scalar && inFixedSpace(Z[blocks])){
      keep <- integer()
    }
    removed <- setdiff(seq_along(blocks), keep)
    mappings[[i]] <- list(keep=keep, diagonal=diagonal,
                           dimension=cs$dim, width=ncol(Z[[blocks[1L]]]))
    if(length(removed)){
      parameters <- if(diagonal) paste0("variance[", cs$levels[removed], "]") else "sigma2"
      omitted[[length(omitted) + 1L]] <- data.frame(
        term=terms[i], parameter=parameters,
        reason="Random design lies in the fixed-effect space (REML)",
        stringsAsFactors=FALSE)
    }
    if(!length(keep)) next
    if(length(removed)){
      f <- cs$factors[[1L]]
      logVariances <- cs$par[1L] + c(0, cs$par[-1L])
      selected <- .compile_covfactor(list(type="diag", dim=length(keep),
        levels=f$levels[keep], par=logVariances[keep[-1L]] - logVariances[keep[1L]],
        free=rep(TRUE, length(keep) - 1L),
        par_names=if(length(keep) > 1L)
          paste0("variance_ratio[", f$levels[keep[-1L]], "]") else character()))
      selected$label <- f$label
      selected$par_start <- 2L
      selected$par_end <- length(keep)
      cs$factors <- list(selected)
      cs$par <- c(logVariances[keep[1L]], selected$par)
      cs$free <- rep(TRUE, length(keep))
      cs$par_names <- c("sigma2", if(length(keep) > 1L)
        paste(f$label, selected$par_names, sep=":"))
      cs$dim <- as.integer(length(keep))
      cs$levels <- cs$levels[keep]
    }
    activeTerms <- c(activeTerms, i)
    j <- length(activeTerms)
    mappings[[i]]$active <- j
    newZ <- c(newZ, Z[blocks[keep]])
    newIndex <- c(newIndex, rep(j, length(keep)))
    newAi[[j]] <- Ai[[i]]
    newCov[[j]] <- cs
    newNames[[j]] <- paste(paste(all.vars(stats::as.formula(paste0("~", terms[i]))),
                                collapse=":"), cs$par_names, sep=":")
  }
  if(!length(omitted)) return(NULL)
  residual <- length(covStruct)
  list(original=original, mappings=mappings, activeTerms=activeTerms,
       omitted=do.call(rbind, omitted), Z=newZ, Zind=newIndex, Ai=newAi,
       covStruct=c(newCov, covStruct[residual]),
       terms=terms[activeTerms],
       termNames=c(newNames, termNames[residual]))
}

.mmes_reml_public_factor <- function(factor, keep, levels){
  # The fitted parameters reference a retained level, not necessarily level one.
  factor$dim <- as.integer(length(levels))
  factor$levels <- levels
  factor$model <- "reml_diagonal"
  factor$evaluator <- list(backend="R", fun=local({
    selected <- keep
    dimension <- length(levels)
    function(par){
      values <- numeric(dimension)
      if(length(selected)) values[selected] <- c(1, exp(par))
      diag(values, dimension)
    }
  }))
  factor$derivative <- list(backend="numeric", rel_step=1e-6)
  factor$native_report <- list(backend="R", fun=local({
    selected <- keep
    labels <- levels
    function(scale, par, factor, absorb_scale=TRUE){
      values <- numeric(length(labels))
      if(length(selected)){
        values[selected] <- (if(absorb_scale) scale else 1) * c(1, par)
      }
      stats::setNames(values, paste0(if(absorb_scale) "variance" else
        "relative_variance", "[", labels, "]"))
    }
  }))
  factor
}

.mmes_reml_restore <- function(result, info, X){
  original <- info$original
  counts <- vapply(result$covPar, length, integer(1))
  starts <- cumsum(c(1L, counts))
  coefficients <- seq_len(ncol(X))
  publicRows <- seq_len(ncol(X))
  publicOffset <- ncol(X)
  parameterRows <- integer()
  publicParameters <- integer()
  parameterOffset <- 0L
  cs <- theta <- covPar <- working <- uList <- partitions <- vector("list", length(original$covStruct))
  pevs <- vector("list", length(original$terms))
  monitors <- vector("list", length(original$covStruct))
  for(i in seq_along(original$terms)){
    mapping <- info$mappings[[i]]
    old <- original$covStruct[[i]]
    dimension <- mapping$dimension
    width <- mapping$width
    ranges <- cbind(publicOffset + (seq_len(dimension) - 1L) * width + 1L,
                    publicOffset + seq_len(dimension) * width)
    partitions[[i]] <- ranges
    active <- mapping$active
    values <- matrix(0, width, dimension)
    pev <- matrix(0, width, dimension)
    covariance <- matrix(0, dimension, dimension)
    if(!is.null(active)){
      values[, mapping$keep] <- result$uList[[active]]
      if(length(result$uPevList[[active]])){
        pev[, mapping$keep] <- result$uPevList[[active]]
      }else pev <- matrix(numeric(), 0L, 0L)
      covariance[mapping$keep, mapping$keep] <- result$theta[[active]]
      source <- unlist(lapply(seq_len(nrow(result$partitions[[active]])), function(j){
        seq.int(result$partitions[[active]][j,1L], result$partitions[[active]][j,2L])
      }), use.names=FALSE)
      target <- unlist(lapply(mapping$keep, function(j) seq.int(ranges[j,1L], ranges[j,2L])),
                       use.names=FALSE)
      coefficients <- c(coefficients, source)
      publicRows <- c(publicRows, target)
      cs[[i]] <- info$covStruct[[active]]
      covPar[[i]] <- result$covPar[[active]]
      if(!is.null(result$covParWorking)) working[[i]] <- result$covParWorking[[active]]
      indices <- seq.int(starts[active], starts[active + 1L] - 1L)
      monitors[[i]] <- result$monitor[indices, , drop=FALSE]
      parameterRows <- c(parameterRows, indices)
      publicParameters <- c(publicParameters, parameterOffset + seq_along(indices))
      if(length(mapping$keep) < dimension){
        cs[[i]]$factors[[1L]] <- .mmes_reml_public_factor(
          cs[[i]]$factors[[1L]], mapping$keep, old$levels)
      }
    }else{
      cs[[i]] <- old
      cs[[i]]$par <- 0
      cs[[i]]$par_names <- "sigma2"
      cs[[i]]$free <- FALSE
      f <- .make_covfactor(dimension, old$levels,
        evaluator=list(backend="fixed", matrix=matrix(0, dimension, dimension)),
        native_report=function(scale, par, factor, absorb_scale=TRUE) c(sigma2=0),
        model="identity")
      if(mapping$diagonal){
        f <- .mmes_reml_public_factor(f, integer(), old$levels)
      }
      f$par_start <- 2L
      f$par_end <- 1L
      cs[[i]]$factors <- list(f)
      covPar[[i]] <- 1
      working[[i]] <- 0
      monitors[[i]] <- matrix(0, 1L, ncol(result$monitor))
    }
    cs[[i]]$dim <- old$dim
    cs[[i]]$levels <- old$levels
    if(length(mapping$keep) < dimension){
      cs[[i]]$reml_omitted <- setdiff(seq_len(dimension), mapping$keep)
    }
    theta[[i]] <- covariance
    uList[[i]] <- values
    pevs[[i]] <- pev
    parameterOffset <- parameterOffset + length(covPar[[i]])
    publicOffset <- publicOffset + dimension * width
  }
  residual <- length(original$covStruct)
  sourceResidual <- length(result$covPar)
  cs[[residual]] <- original$covStruct[[residual]]
  theta[[residual]] <- result$theta[[sourceResidual]]
  covPar[[residual]] <- result$covPar[[sourceResidual]]
  if(!is.null(result$covParWorking)) working[[residual]] <- result$covParWorking[[sourceResidual]]
  indices <- seq.int(starts[sourceResidual], starts[sourceResidual + 1L] - 1L)
  parameterRows <- c(parameterRows, indices)
  publicParameters <- c(publicParameters, parameterOffset + seq_along(indices))
  monitors[[residual]] <- result$monitor[indices, , drop=FALSE]
  parameterOffset <- parameterOffset + length(indices)
  map <- Matrix::sparseMatrix(i=publicRows, j=coefficients, x=1,
    dims=c(publicOffset, nrow(result$bu)))
  # Compose with the native latent-RR map; omitted effects have zero rows.
  result$coefficientMap <- if(is.null(result$coefficientMap)) map else
    map %*% result$coefficientMap
  if(!is.null(result$Ci) && nrow(result$Ci)) result$Ci <- map %*% result$Ci %*% t(map)
  result$bu <- as.matrix(map %*% result$bu)
  result$u <- result$bu[seq.int(ncol(X) + 1L, nrow(result$bu)), , drop=FALSE]
  result$uList <- uList[seq_along(original$terms)]
  result$uPevList <- pevs
  result$partitions <- partitions[seq_along(original$terms)]
  result$W <- do.call(cbind, c(list(X), original$Z))
  if(!is.null(result$Ci) && nrow(result$Ci)){
    diagonal <- Matrix::diag(result$Ci)
    result$uPevList <- lapply(result$partitions, function(ranges){
      matrix(diagonal[seq.int(ranges[1L,1L], ranges[nrow(ranges),2L])],
             nrow=ranges[1L,2L] - ranges[1L,1L] + 1L)
    })
  }
  parameterMap <- Matrix::sparseMatrix(i=publicParameters, j=parameterRows, x=1,
    dims=c(parameterOffset, sum(counts)))
  for(name in c("theta_se", "InfMat", "InfMatInv")){
    if(!is.null(result[[name]]) && nrow(result[[name]]) == sum(counts)){
      result[[name]] <- as.matrix(parameterMap %*% result[[name]] %*% t(parameterMap))
    }
  }
  for(name in c("percDelta", "monitorOriginalScale")){
    if(is.matrix(result[[name]]) && nrow(result[[name]]) == sum(counts)){
      result[[name]] <- as.matrix(parameterMap %*% result[[name]])
    }
  }
  if(length(result$dLu) == sum(counts)){
    result$dLu <- as.matrix(parameterMap %*% result$dLu)
  }
  if(!is.null(result$toBoundary) && ncol(result$toBoundary) == sum(counts)){
    result$toBoundary <- as.matrix(result$toBoundary %*% t(parameterMap))
  }
  result$theta <- theta
  result$covPar <- covPar
  if(!is.null(result$covParWorking)) result$covParWorking <- working
  result$monitor <- do.call(rbind, monitors)
  result$remlOmissions <- info$omitted
  result$remlActiveTerms <- c(info$activeTerms, residual)
  result$engineDiagnostics$remlOmissions <- info$omitted
  list(result=result, covStruct=cs)
}

.mmes_reml_report_status <- function(object, out){
  omitted <- object$remlOmissions
  if(is.null(omitted) || !nrow(out)) return(out)
  key <- function(x) paste(x$term, x$parameter, sep="\r")
  null <- key(out) %in% key(omitted)
  out$status <- ifelse(null, "unestimable (omitted)", "retained")
  if("StdError" %in% names(out)){
    out$StdError[null] <- NA_real_
    out$Zratio[null] <- NA_real_
  }
  covariance <- attr(out, "vcov")
  if(!is.null(covariance)){
    covariance[null,] <- NA_real_
    covariance[,null] <- NA_real_
    attr(out, "vcov") <- covariance
  }
  out
}
