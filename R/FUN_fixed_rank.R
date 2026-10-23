.mmes_reduce_fixed <- function(fixedDesign){
  basis <- fixed_design_basis_cpp(
    as(as(as(fixedDesign, "dMatrix"), "generalMatrix"), "CsparseMatrix"))
  if(!length(basis$omitted)) return(list(X=fixedDesign, info=NULL))
  basis$columnNames <- colnames(fixedDesign)
  basis$assign <- attr(fixedDesign, "assign")
  basis$original <- fixedDesign
  reducedDesign <- fixedDesign[, basis$retained, drop=FALSE]
  attr(reducedDesign, "assign") <- basis$assign[basis$retained]
  attr(reducedDesign, "contrasts") <- attr(fixedDesign, "contrasts")
  list(X=reducedDesign, info=basis)
}

.mmes_restore_fixed <- function(result, basis){
  if(is.null(basis)) return(result)
  originalColumns <- ncol(basis$original)
  reducedColumns <- length(basis$retained)
  randomColumns <- nrow(result$bu) - reducedColumns
  publicIndices <- c(basis$retained, originalColumns + seq_len(randomColumns))
  coefficientMap <- Matrix::sparseMatrix(i=publicIndices,
    j=seq_along(publicIndices), x=1,
    dims=c(originalColumns + randomColumns, nrow(result$bu)))
  result$coefficientMap <- if(is.null(result$coefficientMap)) coefficientMap else
    coefficientMap %*% result$coefficientMap
  if(!is.null(result$Ci) && nrow(result$Ci)){
    result$Ci <- coefficientMap %*% result$Ci %*% t(coefficientMap)
  }
  result$bu <- as.matrix(coefficientMap %*% result$bu)
  result$b <- result$bu[seq_len(originalColumns), , drop=FALSE]
  if(length(result$partitions)){
    result$partitions <- lapply(result$partitions,
      function(partition) partition + originalColumns - reducedColumns)
  }
  randomDesign <- if(is.null(result$W)) NULL else
    result$W[, reducedColumns + seq_len(randomColumns), drop=FALSE]
  result$W <- if(randomColumns && !is.null(randomDesign))
    cbind(basis$original, randomDesign) else basis$original
  # Keep direct-engine VarBeta in reduced coordinates, avoiding a p-by-p expansion.
  basis$original <- NULL
  result$fixedRank <- basis
  result
}

.mmes_check_estimable_fixed <- function(object, contrastDesign){
  basis <- object$fixedRank
  if(is.null(basis)) return(invisible(NULL))
  fixedColumns <- length(basis$columnNames)
  if(ncol(contrastDesign) != nrow(object$bu)){
    stop("Prediction contrasts must have one column per original fixed/random coefficient.",
         call.=FALSE)
  }
  if(!nrow(contrastDesign)) return(invisible(NULL))
  omittedRelation <- basis$relation[, basis$omitted, drop=FALSE]
  invalidRows <- integer()
  for(firstRow in seq.int(1L, nrow(contrastDesign), by=32L)){
    panelRows <- seq.int(firstRow, min(firstRow + 31L, nrow(contrastDesign)))
    fixedPanel <- contrastDesign[panelRows, seq_len(fixedColumns), drop=FALSE]
    difference <- fixedPanel[, basis$omitted, drop=FALSE] -
      fixedPanel[, basis$retained, drop=FALSE] %*% omittedRelation
    discrepancy <- Matrix::rowSums(abs(difference))
    scale <- pmax(1, Matrix::rowSums(abs(fixedPanel)))
    invalidRows <- c(invalidRows,
      panelRows[!is.finite(discrepancy) |
                  discrepancy > max(1e-8, 100 * basis$tolerance) * scale])
  }
  if(length(invalidRows)){
    stop("Prediction/contrast rows ", paste(head(invalidRows, 10L), collapse=", "),
         if(length(invalidRows) > 10L) ", ..." else "",
         " are not estimable under the original rank-deficient fixed model. ",
         "Zero-constrained coefficients are not independently estimated effects. ",
         "Supply estimable contrasts (for example, observed-cell fitted means), ",
         "or refit an explicitly identified full-rank model.", call.=FALSE)
  }
  invisible(NULL)
}

.mmes_require_full_rank_inference <- function(object, method){
  if(!is.null(object$fixedRank)){
    stop(method, " is not available for automatically rank-reduced fixed effects. ",
         "Refit an explicitly identified full-rank model for this inference method; ",
         "estimable custom contrasts can be supplied to predict().",
         call.=FALSE)
  }
  invisible(NULL)
}
