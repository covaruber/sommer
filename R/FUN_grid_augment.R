# Missing-cell augmentation for gridded residuals (ASReml-style).
#
# Each missing cell of a residual block's observed-level rectangle gets one row
# with a zero response, zero random design and its own dummy fixed effect. The
# dummy absorbs the row exactly, so the REML likelihood, variance parameters,
# BLUEs/BLUPs and their covariances are unchanged, while every residual block
# becomes a complete Kronecker grid whose AR/ANTE precision stays banded.

.mmes_grid_factor_band <- function(factor){
  representation <- factor$representation
  kind <- if(is.null(representation$kind)) "dense" else representation$kind
  form <- if(is.null(representation$form)) "" else representation$form
  if(isTRUE(factor$structurally_diagonal) || identical(kind, "diagonal")) return(0L)
  if(identical(form, "ar1")) return(1L)
  if(identical(form, "stationary_ar")){
    order <- if(!is.null(representation$order)) representation$order else factor$order
    return(if(is.null(order)) -1L else as.integer(order))
  }
  if(identical(form, "modified_cholesky")){
    innovation <- representation$innovation
    if(is.null(innovation$rows) || is.null(innovation$columns)) return(-1L)
    return(as.integer(max(0L, innovation$rows - innovation$columns)))
  }
  -1L
}

.mmes_grid_augment <- function(mode, eligibility, covStruct, X, Z, yvar,
                               residualBlock, localIndex, sparseFactors){
  # Quiet skips are models that are not gridded AR/ANTE residuals at all.
  skip <- function(reason, quiet=FALSE){
    if(identical(mode, "force")){
      stop("representation$gridAugment='force' ", reason, ".", call.=FALSE)
    }
    list(applied=FALSE, reason=if(quiet) "" else reason)
  }
  if(identical(mode, "off")) return(list(applied=FALSE, reason=""))
  residual <- covStruct[[length(covStruct)]]
  factors <- residual$factors
  if(!length(factors)) return(skip("requires a factorized (vsm/dsumm) residual", TRUE))
  bands <- vapply(factors, .mmes_grid_factor_band, integer(1))
  if(any(bands < 0L)){
    return(skip("requires AR(p), ANTE or diagonal residual factors", TRUE))
  }
  if(!any(bands >= 1L)) return(skip("requires at least one AR(p)/ANTE residual factor", TRUE))
  if(nzchar(eligibility)) return(skip(eligibility))
  if(identical(sparseFactors, "off") && !identical(mode, "force")){
    return(list(applied=FALSE, reason="requires representation$sparseFactors other than 'off'"))
  }
  if(ncol(yvar) != 1L) return(skip("requires a single response column"))
  dims <- vapply(factors, function(f) as.integer(f$dim), integer(1))
  trails <- rev(cumprod(c(1, rev(dims))[seq_along(dims)]))
  local0 <- as.numeric(localIndex) - 1
  digits <- vapply(seq_along(dims), function(f) (local0 %/% trails[f]) %% dims[f],
                   numeric(length(local0)))
  if(!is.matrix(digits)) digits <- matrix(digits, nrow=length(local0))
  rowsByBlock <- split(seq_along(residualBlock), residualBlock)
  newBlock <- vector("list", length(rowsByBlock))
  newLocal <- vector("list", length(rowsByBlock))
  for(index in seq_along(rowsByBlock)){
    rows <- rowsByBlock[[index]]
    cells <- 0
    for(f in seq_along(dims)){
      levels <- sort(unique(digits[rows, f]))
      if(bands[f] > 1L && length(levels) > 1L && any(diff(levels) != 1)){
        return(skip("requires contiguous observed levels for AR(p)/ANTE factors of order above one"))
      }
      cells <- as.vector(outer(levels * trails[f], cells, "+"))
    }
    if(!.mmes_grid_fits(length(cells), length(rows)) && !identical(mode, "force")){
      return(skip("skips residual blocks whose observed-level grid is much larger than the data"))
    }
    missing <- setdiff(cells, local0[rows])
    newBlock[[index]] <- rep.int(as.integer(names(rowsByBlock)[index]), length(missing))
    newLocal[[index]] <- sort(missing)
  }
  newBlock <- unlist(newBlock, use.names=FALSE)
  newLocal <- unlist(newLocal, use.names=FALSE)
  m <- length(newLocal)
  observed <- nrow(X)
  if(!m) return(skip("found no missing grid cells", TRUE))
  if(m > observed && !identical(mode, "force")){
    return(skip("skips grids with fewer than half of the cells observed"))
  }
  dummy <- Matrix::sparseMatrix(i=seq_len(m), j=seq_len(m), x=1, dims=c(m, m))
  XEngine <- rbind(cbind(X, Matrix::Matrix(0, observed, m, sparse=TRUE)),
                   cbind(Matrix::Matrix(0, m, ncol(X), sparse=TRUE), dummy))
  XEngine <- as(as(XEngine, "generalMatrix"), "CsparseMatrix")
  ZEngine <- lapply(Z, function(z) rbind(z, Matrix::Matrix(0, m, ncol(z), sparse=TRUE)))
  yEngine <- rbind(yvar, Matrix::Matrix(0, m, 1L, sparse=TRUE))
  blockEngine <- c(as.integer(residualBlock), newBlock)
  localEngine <- c(as.integer(localIndex), as.integer(newLocal + 1))
  blockSizes <- tabulate(blockEngine)
  covStruct[[length(covStruct)]]$complete_residual_blocks <-
    length(blockSizes) > 0L && all(blockSizes == covStruct[[length(covStruct)]]$dim)
  attr(covStruct, "gridAugmentation") <- c(as.integer(observed), as.integer(ncol(X)))
  list(applied=TRUE, reason="", covStruct=covStruct, X=XEngine, Z=ZEngine,
       yvar=yEngine, W=Matrix::Diagonal(observed + m),
       residualBlock=blockEngine, localIndex=localEngine,
       info=list(observed=observed, fixed=ncol(X), cells=m,
                 blocks=length(rowsByBlock)))
}

# Mirrors ResidualPolicy::gridFits() in src/mme_paths.h.
.mmes_grid_fits <- function(rectangle, observed){
  rectangle <- as.numeric(rectangle)
  observed <- as.numeric(observed)
  rectangle <= 3 * observed + 1000 && rectangle - observed <= 4000
}

.mmes_grid_restore <- function(result, info){
  if(is.null(info)) return(result)
  p <- info$fixed
  m <- info$cells
  n <- info$observed
  engineEffects <- nrow(result$bu)
  publicIndices <- c(seq_len(p), seq.int(p + m + 1L, length.out=engineEffects - p - m))
  coefficientMap <- Matrix::sparseMatrix(i=seq_along(publicIndices), j=publicIndices,
    x=1, dims=c(length(publicIndices), engineEffects))
  result$coefficientMap <- if(is.null(result$coefficientMap)) coefficientMap else
    coefficientMap %*% result$coefficientMap
  if(!is.null(result$Ci) && nrow(result$Ci)){
    result$Ci <- coefficientMap %*% result$Ci %*% t(coefficientMap)
  }
  result$bu <- as.matrix(coefficientMap %*% result$bu)
  result$b <- result$bu[seq_len(p), , drop=FALSE]
  if(length(result$partitions)){
    result$partitions <- lapply(result$partitions, function(partition) partition - m)
  }
  if(!is.null(result$W) && nrow(result$W) == n + m){
    result$W <- result$W[seq_len(n), -(p + seq_len(m)), drop=FALSE]
  }
  if(length(result$yWorking) == n + m){
    result$yWorking <- if(is.matrix(result$yWorking))
      result$yWorking[seq_len(n), , drop=FALSE] else result$yWorking[seq_len(n)]
  }
  # The engine counted dummy columns/rows; information criteria use the data.
  finalLogLik <- result$llik[length(result$llik)]
  result$AIC <- -2 * finalLogLik + 2 * p
  result$BIC <- -2 * finalLogLik + log(n) * p
  if(is.list(result$engineDiagnostics)){
    result$engineDiagnostics$gridAugmentation <- list(active=TRUE,
      cells=m, observed=n, blocks=info$blocks)
  }
  result
}
