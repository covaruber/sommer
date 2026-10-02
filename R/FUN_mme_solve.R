#### =========== ####
## KNOWN-COVARIANCE MIXED-MODEL EQUATIONS (mmes(solveOnly=TRUE))
#### =========== ####

# Same matrix as Matrix::sparse.model.matrix(fixed, mf), built from the
# distinct covariate rows: Matrix's interaction coding is superlinear in the
# number of records, while distinct rows are usually few (e.g. trait x herd).
.sparse_model_matrix_by_rows <- function(fixed, mf, contrasts, maxCells=2e7){
  fallback <- function() Matrix::sparse.model.matrix(fixed, data=mf, contrasts.arg=contrasts)
  rhs <- mf[-1L]
  if(!length(rhs) || anyNA(rhs) ||
     any(vapply(rhs, function(v) !is.null(dim(v)), logical(1)))){
    return(fallback())
  }
  codes <- lapply(rhs, function(v) if(is.factor(v)) as.integer(v) else match(v, unique(v)))
  key <- if(length(codes) == 1L) codes[[1L]] else do.call(paste, c(unname(codes), sep="\r"))
  first <- which(!duplicated(key))
  if(length(first) > maxCells / 10) return(fallback())
  mfu <- mf[first, , drop=FALSE]
  attr(mfu, "terms") <- attr(mf, "terms")
  Xu <- stats::model.matrix(attr(mf, "terms"), mfu, contrasts.arg=contrasts)
  if(as.double(nrow(Xu)) * ncol(Xu) > maxCells) return(fallback())
  rowOf <- match(key, key[first])
  X <- Matrix::t(Matrix::t(as(as(as(Xu, "dMatrix"), "generalMatrix"), "CsparseMatrix"))[, rowOf, drop=FALSE])
  X <- Matrix::drop0(X)
  dimnames(X) <- list(rownames(mf), colnames(Xu))
  attr(X, "assign") <- attr(Xu, "assign")
  attr(X, "contrasts") <- attr(Xu, "contrasts")
  X
}

.set_covstruct_par <- function(cs, work){
  cs$par <- as.numeric(work)
  cs$free[] <- FALSE
  for(k in seq_along(cs$factors)){
    s <- as.integer(cs$factors[[k]]$par_start)
    e <- as.integer(cs$factors[[k]]$par_end)
    if(length(s) && length(e) && e >= s) cs$factors[[k]]$par <- cs$par[s:e]
  }
  cs
}

.natural_vector_to_working <- function(cs, values, label){
  if(length(values) != length(cs$par) || any(!is.finite(values))){
    stop("covPar for term '", label, "' must contain ", length(cs$par),
         " finite values (", paste(cs$par_names, collapse=", "), ").", call.=FALSE)
  }
  vapply(seq_along(values), function(j){
    kind <- .vc_kind(cs, j)
    if(kind$kind %in% c("scale", "exp") && values[j] <= 0){
      stop("covPar value '", cs$par_names[j], "' of term '", label,
           "' must be positive.", call.=FALSE)
    }
    .vc_to_working(values[j], kind$kind, kind$lower, kind$upper)
  }, numeric(1))
}

# Working parameters reproducing a full covariance matrix K for terms whose
# only parameterized factor is usm() or dsm() (or that have no factors).
.matrix_to_working <- function(cs, K, label){
  K <- as.matrix(K)
  q <- as.integer(cs$dim)
  if(nrow(K) != q || ncol(K) != q || any(!is.finite(K)) || !isSymmetric(unname(K))){
    stop("The covariance matrix for term '", label, "' must be a finite symmetric ", q, " x ", q,
         " matrix.", call.=FALSE)
  }
  parameterized <- Filter(function(f) as.integer(f$par_end) >= as.integer(f$par_start), cs$factors)
  work <- NULL
  if(!length(parameterized) && q == 1L){
    work <- log(K[1, 1])
  }else if(length(cs$factors) == 1L && length(parameterized) == 1L && parameterized[[1L]]$dim == q){
    f <- parameterized[[1L]]
    if(!is.null(f$us_row)){
      ch <- tryCatch(chol(K), error=function(e) NULL)
      if(is.null(ch)) stop("The covariance matrix for term '", label, "' is not positive definite.", call.=FALSE)
      work <- c(log(K[1, 1]), .usm_par_from_K(K, f))
    }else if(identical(f$model, "diag")){
      d <- diag(K)
      if(any(d <= 0)) stop("The covariance matrix for term '", label, "' needs positive variances.", call.=FALSE)
      work <- c(log(d[1L]), log(d[-1L] / d[1L]))
    }
  }
  if(is.null(work) || length(work) != length(cs$par)){
    stop("A covariance matrix can be supplied in covPar only for ism(), dsm() and usm() terms; ",
         "give the natural parameters of term '", label, "' instead.", call.=FALSE)
  }
  if(max(abs(evaluate_covstruct_cpp(cs, work) - K)) > 1e-8 * max(1, abs(K))){
    stop("The covariance matrix for term '", label, "' is not representable by its structure ",
         "(e.g. nonzero covariances for dsm()).", call.=FALSE)
  }
  work
}

# Fixes every covariance parameter at known values. covPar may be a fitted
# mmes object; a list (named by term or in formula order) of natural-scale
# vectors as fit$covPar, or of covariance matrices for ism/dsm/usm terms; or a
# data frame with term, parameter and value columns (as the vcParams of an
# mmes(solveOnly=TRUE) result). NULL requires every parameter to be fixed in
# the formula already.
.mmes_known_covpar <- function(covStruct, labels, covPar){
  if(is.null(covPar)){
    free <- vapply(covStruct, function(cs) any(as.logical(cs$free)), logical(1))
    if(any(free)){
      stop("solveOnly=TRUE needs every covariance parameter. Supply covPar (e.g. a fitted mmes ",
           "object) or fix the parameters of: ", paste(labels[free], collapse=", "), call.=FALSE)
    }
    return(covStruct)
  }

  missingTerms <- function(found){
    if(anyNA(found)){
      stop("covPar does not contain the terms: ", paste(labels[is.na(found)], collapse=", "),
           ". Term labels must match the model formula.", call.=FALSE)
    }
  }

  if(inherits(covPar, "mmes")){
    work <- covPar$covParWorking
    source <- covPar$covStruct
    if(is.null(work) || is.null(source) || length(work) != length(source)){
      stop("The fitted mmes object does not carry working covariance parameters.", call.=FALSE)
    }
    found <- match(labels, names(source))
    missingTerms(found)
    for(i in seq_along(covStruct)){
      src <- source[[found[i]]]
      if(!identical(as.character(src$par_names), as.character(covStruct[[i]]$par_names)) ||
         !identical(as.integer(src$dim), as.integer(covStruct[[i]]$dim))){
        stop("The covariance structure of term '", labels[i],
             "' differs from the one in the fitted object.", call.=FALSE)
      }
      covStruct[[i]] <- .set_covstruct_par(covStruct[[i]], work[[found[i]]])
    }
    return(covStruct)
  }

  if(is.data.frame(covPar)){
    valueColumn <- intersect(c("value", "estimate", "start"), names(covPar))[1L]
    if(!all(c("term", "parameter") %in% names(covPar)) || is.na(valueColumn)){
      stop("A covPar data frame needs columns 'term', 'parameter' and 'value'.", call.=FALSE)
    }
    keys <- paste(covPar$term, covPar$parameter, sep="\r")
    if(anyDuplicated(keys)) stop("covPar contains duplicated term/parameter rows.", call.=FALSE)
    for(i in seq_along(covStruct)){
      cs <- covStruct[[i]]
      rows <- match(paste(labels[i], cs$par_names, sep="\r"), keys)
      if(anyNA(rows)){
        stop("covPar is missing parameters of term '", labels[i], "': ",
             paste(cs$par_names[is.na(rows)], collapse=", "), call.=FALSE)
      }
      values <- as.numeric(covPar[[valueColumn]][rows])
      covStruct[[i]] <- .set_covstruct_par(cs, .natural_vector_to_working(cs, values, labels[i]))
    }
    return(covStruct)
  }

  if(is.list(covPar)){
    found <- if(is.null(names(covPar)) && length(covPar) == length(covStruct)) seq_along(covStruct)
             else match(labels, names(covPar))
    missingTerms(found)
    for(i in seq_along(covStruct)){
      value <- covPar[[found[i]]]
      # fit$covPar stores natural vectors as one-column matrices.
      square <- (is.matrix(value) || inherits(value, "Matrix")) &&
        nrow(value) == ncol(value) && nrow(value) > 1L
      work <- if(square) .matrix_to_working(covStruct[[i]], value, labels[i])
              else .natural_vector_to_working(covStruct[[i]], as.numeric(value), labels[i])
      covStruct[[i]] <- .set_covstruct_par(covStruct[[i]], work)
    }
    return(covStruct)
  }

  stop("covPar must be a fitted mmes object, a named list or a data frame.", call.=FALSE)
}

.mmes_solve <- function(X, Z, Zind, Ai, yvar, W, useH, residualBlock, localIndex,
                        covStruct, labels, covPar, pcgTol, pcgMaxIters, verbose,
                        data, dataOriginal, obsInfo, partitionsX, call, args){
  if(ncol(yvar) != 1L){
    stop("solveOnly=TRUE requires one response column; use the long format with vsm(usm(trait), ...).",
         call.=FALSE)
  }
  covStruct <- .mmes_known_covpar(covStruct, labels, covPar)
  nRandom <- length(covStruct) - 1L
  residual <- covStruct[[length(covStruct)]]

  weights <- numeric(0)
  if(isTRUE(useH)){
    if(!Matrix::isDiagonal(W)){
      stop("solveOnly=TRUE supports diagonal weight matrices W only.", call.=FALSE)
    }
    weights <- as.numeric(Matrix::diag(W))
  }

  nFixed <- ncol(X)
  termStart <- integer(nRandom)
  termLevels <- integer(nRandom)
  lambda <- vector("list", nRandom)
  AiGeneral <- vector("list", nRandom)
  theta <- vector("list", length(covStruct))
  offset <- nFixed
  for(u in seq_len(nRandom)){
    cs <- covStruct[[u]]
    Zu <- Z[Zind == u]
    if(length(Zu) != cs$dim){
      stop("Internal error: term '", labels[u], "' has ", length(Zu),
           " design blocks for a covariance of dimension ", cs$dim, ".", call.=FALSE)
    }
    Sigma <- evaluate_covstruct_cpp(cs, cs$par)
    dimnames(Sigma) <- list(cs$levels, cs$levels)
    theta[[u]] <- Sigma
    ch <- tryCatch(chol(Sigma), error=function(e) NULL)
    if(is.null(ch)){
      stop("The covariance of term '", labels[u], "' is not positive definite at the supplied values.",
           call.=FALSE)
    }
    lambda[[u]] <- chol2inv(ch)
    termStart[u] <- offset
    termLevels[u] <- ncol(Zu[[1L]])
    AiGeneral[[u]] <- as(as(as(Ai[[u]], "dMatrix"), "generalMatrix"), "CsparseMatrix")
    offset <- offset + length(Zu) * termLevels[u]
  }
  sectioned <- any(vapply(residual$factors, function(f) !is.null(f$section_levels), logical(1)))
  if(!sectioned && residual$dim <= 500L){
    theta[[length(covStruct)]] <- evaluate_covstruct_cpp(residual, residual$par)
    dimnames(theta[[length(covStruct)]]) <- list(residual$levels, residual$levels)
  }
  names(theta) <- labels

  design <- to_sparse(do.call(cbind, c(list(X), Z)))
  y <- as.numeric(yvar[, 1L])

  message(crayon::blue("Engine selected: known-covariance PCG (solveOnly=TRUE)"))
  sol <- mme_pcg_solve(design, y, nFixed, termStart, termLevels, AiGeneral, lambda,
                       residual, residual$par, as.integer(residualBlock),
                       as.integer(localIndex), weights, pcgTol, as.integer(pcgMaxIters),
                       1000L, 10000L, numeric(0), verbose)
  if(!sol$converged){
    warning("PCG did not reach the requested tolerance (relative residual ",
            signif(sol$relres, 3), " after ", sol$iterations,
            " iterations); increase pcgMaxIters.", call.=FALSE)
  }

  effectNames <- c(colnames(X), unlist(lapply(Z, colnames), use.names=FALSE))
  bu <- matrix(sol$solution, ncol=1L, dimnames=list(effectNames, NULL))
  b <- bu[seq_len(nFixed), , drop=FALSE]
  u <- bu[-seq_len(nFixed), , drop=FALSE]
  uList <- vector("list", nRandom)
  partitions <- vector("list", nRandom)
  for(i in seq_len(nRandom)){
    q <- covStruct[[i]]$dim
    nl <- termLevels[i]
    first <- termStart[i] + (seq_len(q) - 1L) * nl + 1L
    partitions[[i]] <- cbind(first, first + nl - 1L)
    uList[[i]] <- matrix(sol$solution[termStart[i] + seq_len(q * nl)], nrow=nl, ncol=q,
                         dimnames=list(colnames(Z[Zind == i][[1L]]), covStruct[[i]]$levels))
  }
  names(uList) <- names(partitions) <- labels[seq_len(nRandom)]

  vcParams <- .vc_param_table(covStruct, labels)
  names(vcParams)[names(vcParams) == "start"] <- "value"
  names(covStruct) <- labels

  out <- list(b=b, u=u, bu=bu, uList=uList, partitions=partitions,
              partitionsX=partitionsX, theta=theta, covStruct=covStruct,
              vcParams=vcParams, fitted=sol$fitted, residuals=y - sol$fitted,
              y=yvar, data=data, dataOriginal=dataOriginal, obsInfo=obsInfo,
              pcg=sol[setdiff(names(sol), c("solution", "fitted"))],
              convergence=sol$converged, call=call, args=args)
  class(out) <- "mmesSolve"
  out
}

print.mmesSolve <- function(x, ...){
  cat("Known-covariance mixed-model solution (mmes, solveOnly=TRUE)\n")
  cat("Equations:", nrow(x$bu), " fixed:", nrow(x$b), " random:", nrow(x$u),
      " records:", length(x$fitted), "\n")
  cat("PCG", if(isTRUE(x$convergence)) "converged" else "did NOT converge",
      "in", x$pcg$iterations, "iterations; relative residual", signif(x$pcg$relres, 3), "\n")
  cat("Variance parameters used:\n")
  print(x$vcParams[, c("term", "parameter", "value")], row.names=FALSE)
  invisible(x)
}
