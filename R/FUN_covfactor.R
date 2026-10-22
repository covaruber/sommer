
# -------------------------------------------------------------------------
# Generic covariance-factor constructors for the Henderson-only interface.
# Each constructor supplies a dimensionless covariance SHAPE. vsm() owns the
# single overall variance scale, which avoids non-identifiability in arbitrary
# Kronecker products.
# -------------------------------------------------------------------------

.cov_dummy <- function(x, expr){
  if(is.matrix(x) || inherits(x, "Matrix")){
    dummy <- x
    if(is.null(colnames(dummy))) colnames(dummy) <- paste0("L", seq_len(ncol(dummy)))
  }else if(!is.character(x) && !is.factor(x)){
    dummy <- Matrix::Matrix(x, ncol=1, sparse=TRUE)
    colnames(dummy) <- expr
  }else{
    if(is.factor(x)){
      # Unused levels would add unidentifiable, data-free blocks to K (x) Gu.
      levs <- levels(x)[tabulate(as.integer(x), nlevels(x)) > 0L]
    }else{
      levs <- unique(na.omit(x))
    }
    if(length(levs) > 1L){
      xf <- factor(x, levels=levs)
      observed <- !is.na(xf)
      # Build the incidence matrix directly from (row, level) index pairs so
      # NA rows are simply omitted from the triplet list (zero row), instead
      # of allocating a full zero sparse matrix and using submatrix<-
      # assignment (dummy[observed,]<-...), which forces CHOLMOD/Matrix to
      # rebuild the whole sparse pattern and is O(n) per assigned column.
      dummy <- Matrix::sparseMatrix(i=which(observed), j=as.integer(xf[observed]),
                                    x=1, dims=c(length(x), length(levs)))
      colnames(dummy) <- levs
    }else{
      observed <- !is.na(x)
      dummy <- Matrix::sparseMatrix(i=which(observed), j=rep(1L, sum(observed)),
                                    x=1, dims=c(length(x), 1L))
      colnames(dummy) <- as.character(levs)
    }
  }
  if(!inherits(dummy, "dgCMatrix")){
    dummy <- as(as(as(dummy, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
  dummy
}


# -------------------------------------------------------------------------
# Low-level universal CovarianceFactor constructor.
#
# Package developers can use this directly for future covariance structures.
# In particular, evaluator=list(backend="R", fun=...) plus either an R
# derivative callback or derivative=list(backend="numeric") creates a new
# covariance model without any change to ai_mme_sp2().
# -------------------------------------------------------------------------
.covfactor_representation <- function(model, factor=list(),
                                     structurally_diagonal=FALSE){
  model <- if(length(model)) as.character(model)[1] else ""
  supplied <- factor$representation
  if(is.list(supplied) && is.character(supplied$kind) &&
     length(supplied$kind) == 1L){
    return(supplied)
  }

  kind <- if(isTRUE(structurally_diagonal) ||
             model %in% c("identity", "diag")){
    "diagonal"
  }else if(model %in% c("fa", "rr")){
    "lowrank_diagonal"
  }else if(model %in% c("csm")){
    "compound_symmetry"
  }else if(model == "ar1" && !isTRUE(factor$heterogeneous)){
    "banded_precision"
  }else if(model %in% c("arp", "ar1", "ar2", "ar3")){
    "autoregressive_precision"
  }else if(model == "ante"){
    "modified_cholesky_precision"
  }else if(model %in% c("sar", "car")){
    "spatial_precision"
  }else if(model == "ma"){
    "banded_covariance"
  }else if(model %in% c("us", "corg")){
    "dense_cholesky"
  }else{
    "dense"
  }

  out <- list(kind=kind)
  form <- switch(model,
    ar1=if(isTRUE(factor$heterogeneous)) "stationary_ar" else "ar1",
    arp="stationary_ar",
    ar2="stationary_ar",
    ar3="stationary_ar",
    ante="modified_cholesky",
    csm="compound_symmetry",
    sar="normalized_normal_equations",
    car="normalized_affine_precision",
    NULL
  )
  if(!is.null(form)) out$form <- form
  if(model == "ante"){
    out$innovation <- list(rows=factor$ante_row, columns=factor$ante_col,
                           coefficient_count=factor$ante_ncoef)
  }
  if(model %in% c("fa", "rr")){
    prefix <- paste0(model, "_")
    out$loadings <- list(rank=factor$order,
      count=factor[[paste0(prefix, "nload")]],
      rows=factor[[paste0(prefix, "row")]],
      columns=factor[[paste0(prefix, "col")]],
      diagonal=factor[[paste0(prefix, "diag")]],
      heterogeneous_specific=identical(model, "fa") &&
        !identical(factor$variance, "homogeneous"),
      zero_specific=isTRUE(factor$zero_specific))
  }
  if(model %in% c("sar", "car")){
    out$spatial <- list(weights=factor$W,
      lower=if(model == "sar") -factor$rho_bound else factor$rho_lower,
      upper=if(model == "sar") factor$rho_bound else factor$rho_upper,
      degree=if(model == "car") factor$degree else NULL)
  }
  if(kind %in% c("lowrank_diagonal", "autoregressive_precision",
                 "banded_covariance") && !is.null(factor$order)){
    out$order <- as.integer(factor$order)[1]
  }
  out
}

.make_covfactor <- function(dim, levels, par=numeric(), free=logical(),
                            par_names=character(), evaluator,
                            derivative=list(backend="numeric", rel_step=1e-6),
                            report=NULL, native_report=NULL, trust_cap=NULL,
                            structurally_diagonal=FALSE, model=NULL,
                            metadata=list()){
  dim <- as.integer(dim)[1]
  par <- as.numeric(par)
  free <- as.logical(free)
  par_names <- as.character(par_names)
  if(length(par) != length(free) || length(par) != length(par_names)){
    stop("CovarianceFactor par/free/par_names lengths differ.", call. = FALSE)
  }
  if(is.null(report)){
    report <- list(
      backend="builtin",
      transform=rep("identity", length(par)),
      lower=rep(NA_real_, length(par)),
      upper=rep(NA_real_, length(par))
    )
  }
  customNativeReport <- native_report
  if(is.null(trust_cap)) trust_cap <- rep(1.5, length(par))
  if(length(trust_cap) != length(par)){
    stop("CovarianceFactor trust_cap must have one value per parameter.",
         call. = FALSE)
  }
  out <- c(list(
    dim=dim,
    levels=as.character(levels),
    par=par,
    free=free,
    par_names=par_names,
    evaluator=evaluator,
    derivative=derivative,
    report=report,
    native_report=customNativeReport,
    trust_cap=as.numeric(trust_cap),
    structurally_diagonal=isTRUE(structurally_diagonal),
    representation=.covfactor_representation(
      model, metadata, structurally_diagonal
    ),
    descriptor_version=2L,
    model=model
  ), metadata)
  if(is.null(out$native_report)){
    out$native_report <- list(
      backend="R",
      fun=.covfactor_native_report(model, out)
    )
  }else if(is.function(out$native_report)){
    out$native_report <- list(backend="R", fun=out$native_report)
  }
  class(out) <- c("sommer_covfactor", "list")
  .validate_covfactor(out)
  out
}

.covfactor_native_report <- function(model, f){
  levels <- as.character(f$levels)

  if(model == "identity"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      if(absorb_scale) c(sigma2=scale) else numeric()
    })
  }

  if(model == "diag"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      multiplier <- if(absorb_scale) scale else 1
      prefix <- if(absorb_scale) "variance" else "relative_variance"
      stats::setNames(multiplier * c(1, par), paste0(prefix, "[", levels, "]"))
    })
  }

  if(model %in% c("ar1", "ar2", "ar3", "arp")){
    order <- if(model == "arp") as.integer(f$order) else
      as.integer(sub("ar", "", model))
    return(function(scale, par, factor, absorb_scale=TRUE){
      dependence <- par[seq_len(order)]
      dependenceNames <- if(order == 1L) "rho" else paste0("pacf[", seq_len(order), "]")
      names(dependence) <- dependenceNames
      if(length(par) == order){
        return(c(if(absorb_scale) c(variance=scale), dependence))
      }
      multiplier <- if(absorb_scale) scale else 1
      prefix <- if(absorb_scale) "variance" else "relative_variance"
      variances <- multiplier * c(1, par[-seq_len(order)])
      c(dependence,
        stats::setNames(variances, paste0(prefix, "[", levels, "]")))
    })
  }

  if(model == "csm"){
    heterogeneous <- identical(f$variance, "heterogeneous")
    return(function(scale, par, factor, absorb_scale=TRUE){
      if(!heterogeneous){
        return(c(if(absorb_scale) c(variance=scale), rho=par[1]))
      }
      multiplier <- if(absorb_scale) scale else 1
      prefix <- if(absorb_scale) "variance" else "relative_variance"
      c(rho=par[1], stats::setNames(multiplier * c(1, par[-1]),
                                    paste0(prefix, "[", levels, "]")))
    })
  }

  if(model == "us"){
    rows <- as.integer(f$us_row)
    cols <- as.integer(f$us_col)
    return(function(scale, par, factor, absorb_scale=TRUE){
      q <- length(levels)
      L <- matrix(0, q, q)
      L[1,1] <- 1
      for(k in seq_along(par)) L[rows[k], cols[k]] <- par[k]
      Sigma <- tcrossprod(L) * if(absorb_scale) scale else 1
      idx <- which(lower.tri(Sigma, diag=TRUE), arr.ind=TRUE)
      labels <- ifelse(idx[,1] == idx[,2],
                       paste0(if(absorb_scale) "variance" else "relative_variance",
                              "[", levels[idx[,1]], "]"),
                       paste0(if(absorb_scale) "covariance" else "relative_covariance",
                              "[", levels[idx[,1]], ",", levels[idx[,2]], "]"))
      stats::setNames(Sigma[idx], labels)
    })
  }

  if(model == "ma"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      c(if(absorb_scale) c(variance=scale),
        stats::setNames(par, factor$par_names))
    })
  }

  if(model == "corg"){
    rows <- as.integer(f$corg_row)
    cols <- as.integer(f$corg_col)
    return(function(scale, par, factor, absorb_scale=TRUE){
      q <- length(levels)
      A <- diag(q)
      for(k in seq_along(par)) A[rows[k], cols[k]] <- par[k]
      R <- stats::cov2cor(tcrossprod(A))
      idx <- which(lower.tri(R), arr.ind=TRUE)
      correlations <- stats::setNames(
        R[idx], paste0("correlation[", levels[idx[,1]], ",", levels[idx[,2]], "]"))
      c(if(absorb_scale) c(variance=scale), correlations)
    })
  }

  if(model == "fa"){
    rows <- as.integer(f$fa_row)
    cols <- as.integer(f$fa_col)
    nload <- as.integer(f$fa_nload)
    return(function(scale, par, factor, absorb_scale=TRUE){
      rawLoading <- par[seq_len(nload)]
      normalization <- rawLoading[1L]^2 + 1
      multiplier <- if(absorb_scale) sqrt(scale / normalization) else
        1 / sqrt(normalization)
      loading <- multiplier * rawLoading
      names(loading) <- paste0("loading[", levels[rows], ",F", cols, "]")
      specificRaw <- if(identical(factor$variance, "homogeneous")) rep(1, length(levels)) else
        c(1, par[-seq_len(nload)])
      specific <- (if(absorb_scale) scale else 1) * specificRaw / normalization
      if(identical(factor$variance, "homogeneous")){
        return(c(loading, common_specific_variance=specific[1L]))
      }
      names(specific) <- paste0(if(absorb_scale) "specific_variance" else
                                  "relative_specific_variance", "[", levels, "]")
      c(loading, specific)
    })
  }

  if(model == "rr"){
    rows <- as.integer(f$rr_row)
    cols <- as.integer(f$rr_col)
    return(function(scale, par, factor, absorb_scale=TRUE){
      remainder <- if(isTRUE(factor$zero_specific)) 0 else 1
      normalization <- par[1L]^2 + remainder
      multiplier <- if(absorb_scale) sqrt(scale / normalization) else
        1 / sqrt(normalization)
      loading <- multiplier * par
      names(loading) <- paste0("loading[", levels[rows], ",F", cols, "]")
      if(remainder > 0){
        return(c(loading,
          common_specific_variance=(if(absorb_scale) scale else 1) / normalization))
      }
      specific <- stats::setNames(rep(0, length(levels)),
        paste0(if(absorb_scale) "specific_variance" else
          "relative_specific_variance", "[", levels, "]"))
      c(loading, specific)
    })
  }

  if(model == "ante"){
    ncoef <- as.integer(f$ante_ncoef)
    return(function(scale, par, factor, absorb_scale=TRUE){
      coefficients <- if(ncoef) par[seq_len(ncoef)] else numeric()
      names(coefficients) <- factor$par_names[seq_len(ncoef)]
      ratios <- if(ncoef < length(par)) par[(ncoef+1L):length(par)] else numeric()
      multiplier <- if(absorb_scale) scale else 1
      prefix <- if(absorb_scale) "innovation_variance" else
        "relative_innovation_variance"
      innovations <- stats::setNames(multiplier * c(1, ratios),
                                     paste0(prefix, "[", levels, "]"))
      c(coefficients, innovations)
    })
  }

  if(model == "matern"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      c(if(absorb_scale) c(variance=scale),
        stats::setNames(as.numeric(par), factor$par_names))
    })
  }

  if(model == "metric"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      c(if(absorb_scale) c(variance=scale),
        stats::setNames(as.numeric(par), factor$par_names))
    })
  }

  if(model == "toeplitz"){
    return(function(scale, par, factor, absorb_scale=TRUE){
      c(if(absorb_scale) c(variance=scale),
        stats::setNames(par, paste0("pacf[", seq_along(par), "]")))
    })
  }

  if(model %in% c("sar", "car")){
    return(function(scale, par, factor, absorb_scale=TRUE){
      c(if(absorb_scale) c(variance=scale), rho=par[1])
    })
  }

  function(scale, par, factor, absorb_scale=TRUE){
    values <- stats::setNames(as.numeric(par), factor$par_names)
    if(absorb_scale) c(sigma2=scale, values) else values
  }
}

# -------------------------------------------------------------------------
# Universal CovarianceFactor compiler.
#
# Public covariance constructors create their convenient model-specific
# metadata in R, then this compiler converts it into the single interface
# consumed by vsm() and ai_mme_sp2():
#
#   evaluator  : how K(eta) is evaluated
#   derivative : how dK/deta_k is obtained
#   report     : transformation from working to reported coordinates
#   native_report: covariance-model interpretation on its natural scale
#   trust_cap  : per-parameter maximum proposal size in working coordinates
#
# ai_mme_sp2() never needs to know the statistical model label.  Built-in
# high-use structures use the native covariance backend; user-defined and new
# experimental structures can use R callbacks without changing the solver.
# -------------------------------------------------------------------------
.compile_covfactor <- function(f){
  if(inherits(f, "sommer_covfactor") &&
     !is.null(f$evaluator) && !is.null(f$derivative) &&
     !is.null(f$report) && !is.null(f$native_report) &&
     !is.null(f$trust_cap)){
    representation <- f$representation
    legacySpatial <- isTRUE(f$model %in% c("sar", "car")) &&
      identical(representation$kind, "dense") && is.null(representation$form)
    needsCapabilities <- is.null(representation) || legacySpatial ||
      (identical(representation$kind, "lowrank_diagonal") && is.null(representation$loadings)) ||
      (identical(representation$kind, "modified_cholesky_precision") && is.null(representation$innovation))
    if(needsCapabilities){
      metadata <- f
      metadata$representation <- NULL
      defaults <- .covfactor_representation(f$model, metadata, f$structurally_diagonal)
      if(legacySpatial) representation$kind <- defaults$kind
      f$representation <- if(is.null(representation)) defaults else utils::modifyList(defaults, representation)
      .validate_covfactor(f)
    }
    return(f)
  }

  if(is.null(f$type) || is.null(f$dim) || is.null(f$par) ||
     is.null(f$free) || is.null(f$par_names)){
    stop("Malformed covariance factor supplied to .compile_covfactor().",
         call. = FALSE)
  }

  model <- as.character(f$type)[1]
  q <- as.integer(f$dim)[1]
  p <- length(f$par)

  if(length(f$free) != p || length(f$par_names) != p){
    stop("Covariance factor par/free/par_names have inconsistent lengths.",
         call. = FALSE)
  }

  # Defaults: a native evaluator, numerical first derivative, identity
  # reporting, and a moderate working-coordinate trust cap.
  evaluator <- list(backend="native", op=model)
  derivative <- list(backend="numeric", rel_step=1e-6)
  transform <- rep("identity", p)
  lower <- rep(NA_real_, p)
  upper <- rep(NA_real_, p)
  trust <- rep(1.5, p)
  structurally_diagonal <- FALSE

  if(model == "identity"){
    evaluator <- list(backend="native", op="identity")
    derivative <- list(backend="none")
    trust <- numeric()
    structurally_diagonal <- TRUE

  }else if(model == "diag"){
    derivative <- list(backend="native", op="diag")
    structurally_diagonal <- TRUE
    transform[] <- "exp"
    trust[] <- 1.0

  }else if(model == "ar1"){
    derivative <- list(backend="native", op="ar1")
    transform[] <- "tanh"
    trust[] <- 1.0

  }else if(model == "us"){
    derivative <- list(backend="native", op="us")
    dd <- as.logical(f$us_diag)
    if(length(dd) != p) stop("Malformed US covariance metadata.", call. = FALSE)
    transform[dd] <- "exp"
    trust[dd] <- 1.0
    trust[!dd] <- 2.0

  }else if(model == "csm"){
    lo <- -1/(q-1)
    variance <- match.arg(f$variance, c("homogeneous", "heterogeneous"))
    if(variance == "homogeneous"){
      if(p != 1L) stop("Malformed homogeneous CSM covariance metadata.", call. = FALSE)
      evaluator$op <- "cor_uniform"
      transform[] <- "bounded_logit"
      lower[] <- lo
      upper[] <- 1
    }else{
      if(p != q) stop("Malformed heterogeneous CSM covariance metadata.", call. = FALSE)
      evaluator$op <- "corh"
      transform[1] <- "bounded_logit"
      lower[1] <- lo
      upper[1] <- 1
      if(p > 1) transform[2:p] <- "exp"
    }
    trust[] <- 1.0

  }else if(model == "arp"){
    transform[] <- "tanh"
    trust[] <- 1.0

  }else if(model == "ma"){
    trust[] <- 1.5

  }else if(model == "corg"){
    trust[] <- 1.5

  }else if(model == "fa"){
    nload <- as.integer(f$fa_nload)
    dd <- as.logical(f$fa_diag)
    if(length(dd) != nload || p < nload){
      stop("Malformed FA covariance metadata.", call. = FALSE)
    }
    if(nload){
      idx <- seq_len(nload)
      transform[idx[dd]] <- "exp"
      trust[idx[dd]] <- 1.0
      trust[idx[!dd]] <- 1.5
    }
    if(p > nload){
      transform[(nload+1L):p] <- "exp"
      trust[(nload+1L):p] <- 1.0
    }
    derivative <- list(backend="native", op="fa")

  }else if(model == "ante"){
    ncoef <- as.integer(f$ante_ncoef)
    if(ncoef < 0L || ncoef > p) stop("Malformed ANTE covariance metadata.", call. = FALSE)
    if(ncoef) trust[seq_len(ncoef)] <- 1.5
    if(p > ncoef){
      transform[(ncoef+1L):p] <- "exp"
      trust[(ncoef+1L):p] <- 1.0
    }

  }else if(model == "own"){
    if(!is.null(f$fixed_matrix)){
      evaluator <- list(backend="fixed", matrix=f$fixed_matrix)
      derivative <- list(backend="none")
      trust <- numeric()
      fm <- as.matrix(f$fixed_matrix)
      structurally_diagonal <- all(abs(fm[row(fm) != col(fm)]) <= 1e-14)
    }else{
      if(is.null(f$fun) || !is.function(f$fun)){
        stop("User covariance factor is missing fun(par).", call. = FALSE)
      }
      evaluator <- list(backend="R", fun=f$fun)
      if(!is.null(f$dfun)){
        derivative <- list(backend="R", fun=f$dfun)
      }else{
        derivative <- list(backend="numeric", rel_step=1e-6)
      }
      trust[] <- 1.5
    }

  }else if(model == "rr"){
    order <- as.integer(f$order)
    rows <- as.integer(f$rr_row)
    cols <- as.integer(f$rr_col)
    dd <- as.logical(f$rr_diag)
    nload <- as.integer(f$rr_nload)
    remainder <- if(isTRUE(f$zero_specific)) 0 else 1
    if(length(rows) != nload || length(cols) != nload ||
       length(dd) != nload || p != nload){
      stop("Malformed reduced-rank covariance metadata.", call. = FALSE)
    }

    rr_eval <- local({
      q0 <- q; k0 <- order; r0 <- rows; c0 <- cols; d0 <- dd; psi0 <- remainder
      function(par){
        L <- matrix(0, q0, k0)
        for(a in seq_along(par)){
          L[r0[a], c0[a]] <- if(d0[a]) exp(par[a]) else par[a]
        }
        M <- tcrossprod(L) + diag(psi0, q0)
        M / M[1,1]
      }
    })

    rr_d1 <- local({
      q0 <- q; k0 <- order; r0 <- rows; c0 <- cols; d0 <- dd; psi0 <- remainder
      function(par, k){
        L <- matrix(0, q0, k0)
        for(a in seq_along(par)){
          L[r0[a], c0[a]] <- if(d0[a]) exp(par[a]) else par[a]
        }
        dL <- matrix(0, q0, k0)
        dL[r0[k], c0[k]] <- if(d0[k]) exp(par[k]) else 1
        M <- tcrossprod(L) + diag(psi0, q0)
        D <- dL %*% t(L) + L %*% t(dL)
        s <- M[1,1]
        ds <- D[1,1]
        (D*s - M*ds)/(s*s)
      }
    })

    evaluator <- list(backend="R", fun=rr_eval)
    derivative <- list(backend="R", fun=rr_d1)
    transform[dd] <- "exp"
    trust[dd] <- 1.0
    trust[!dd] <- 1.5
    derivative <- list(backend="native", op="rr")

  }else{
    stop(paste0("Unknown covariance model label '", model,
                "'. New models should provide a sommer_covfactor descriptor ",
                "or use ownm()."), call. = FALSE)
  }

  f$model <- model                 # human-readable/debugging only
  f$type <- NULL                   # model label never crosses as solver dispatch
  f$evaluator <- evaluator
  f$derivative <- derivative
  f$representation <- .covfactor_representation(
    model, f, structurally_diagonal
  )
  if(f$representation$kind == "lowrank_diagonal"){
    f$precision <- list(backend="woodbury", kind=model)
  }else if(f$representation$kind == "compound_symmetry"){
    f$precision <- list(backend="woodbury", kind="compound_symmetry")
  }
  f$report <- list(
    backend="builtin",
    transform=as.character(transform),
    lower=as.numeric(lower),
    upper=as.numeric(upper)
  )
  if(is.function(f$native_report)){
    f$native_report <- list(backend="R", fun=f$native_report)
  }
  if(is.null(f$native_report)){
    f$native_report <- list(
      backend="R",
      fun=.covfactor_native_report(model, f)
    )
  }
  f$trust_cap <- as.numeric(trust)
  f$structurally_diagonal <- isTRUE(structurally_diagonal)
  f$descriptor_version <- 2L
  class(f) <- c("sommer_covfactor", "list")
  f
}

.mme_sparse_factor_kernels_auto <- function(factor){
  representation <- factor$representation
  dimension <- factor$dim
  if(!is.numeric(dimension) || length(dimension) != 1L ||
     !is.finite(dimension)){
    return(FALSE)
  }
  if(identical(representation$form, "modified_cholesky") &&
     dimension >= 64L && is.list(representation$innovation)){
    innovation <- representation$innovation
    rows <- innovation$rows
    columns <- innovation$columns
    count <- innovation$coefficient_count
    return(is.numeric(count) && length(count) == 1L && is.finite(count) &&
      count > 0L && count == length(rows) && count == length(columns) &&
      is.numeric(rows) && is.numeric(columns) &&
      all(is.finite(rows)) && all(is.finite(columns)) &&
      all(rows > columns) && max(rows - columns) <= 3L)
  }
  if(dimension >= 500L && is.character(representation$form) &&
     length(representation$form) == 1L &&
     representation$form %in% c("normalized_normal_equations",
                                "normalized_affine_precision") &&
     is.list(representation$spatial)){
    weights <- representation$spatial$weights
    return(is.matrix(weights) && identical(dim(weights), c(dimension, dimension)) &&
      all(is.finite(weights)) && mean(weights != 0) <= 0.01)
  }
  FALSE
}

.validate_covfactor <- function(f){
  required <- c("dim","levels","par","free","par_names","evaluator",
                "derivative","report","native_report","trust_cap",
                "structurally_diagonal","representation","descriptor_version")
  miss <- setdiff(required, names(f))
  if(length(miss)){
    stop(paste("Incomplete CovarianceFactor descriptor; missing:",
               paste(miss, collapse=", ")), call. = FALSE)
  }
  if(length(f$par) != length(f$free) || length(f$par) != length(f$par_names)){
    stop("CovarianceFactor par/free/par_names lengths differ.", call. = FALSE)
  }
  if(length(f$trust_cap) != length(f$par)){
    stop("CovarianceFactor trust_cap must have one value per working parameter.",
         call. = FALSE)
  }
  if(!is.list(f$representation) ||
     !is.character(f$representation$kind) ||
     length(f$representation$kind) != 1L ||
     !nzchar(f$representation$kind)){
    stop("CovarianceFactor representation must declare one family kind.",
         call. = FALSE)
  }
  innovation <- f$representation$innovation
  if(!is.null(innovation)){
    count <- innovation$coefficient_count
    rows <- innovation$rows
    columns <- innovation$columns
    if(length(count) != 1L || !is.finite(count) || count < 0 || count != as.integer(count) ||
       length(rows) != count || length(columns) != count ||
       anyNA(rows) || anyNA(columns) || any(rows != as.integer(rows)) ||
       any(columns != as.integer(columns)) || any(columns < 1L) ||
       any(rows > f$dim) || any(rows <= columns) ||
       anyDuplicated(data.frame(rows=rows, columns=columns))){
      stop("Invalid modified-Cholesky innovation capability.", call. = FALSE)
    }
    if(is.null(f$section_npar) && length(f$par) != count + f$dim - 1L){
      stop("Innovation capability has incompatible parameter count.", call. = FALSE)
    }
  }
  loadings <- f$representation$loadings
  if(!is.null(loadings)){
    count <- loadings$count
    rank <- loadings$rank
    rows <- loadings$rows
    columns <- loadings$columns
    if(length(count) != 1L || length(rank) != 1L || !is.finite(count) || !is.finite(rank) ||
       count < 1L || rank < 1L || count != as.integer(count) || rank != as.integer(rank) ||
       length(rows) != count || length(columns) != count || length(loadings$diagonal) != count ||
       anyNA(rows) || anyNA(columns) || anyNA(loadings$diagonal) ||
       any(rows != as.integer(rows)) || any(columns != as.integer(columns)) ||
       any(rows < 1L | rows > f$dim) || any(columns < 1L | columns > rank) ||
       anyDuplicated(data.frame(rows=rows, columns=columns)) ||
       !is.logical(loadings$diagonal) || !is.logical(loadings$heterogeneous_specific) ||
       length(loadings$heterogeneous_specific) != 1L || is.na(loadings$heterogeneous_specific)){
      stop("Invalid low-rank loading capability.", call. = FALSE)
    }
  }
  spatial <- f$representation$spatial
  if(!is.null(spatial)){
    if(!is.matrix(spatial$weights) || !is.numeric(spatial$weights) ||
       !identical(dim(spatial$weights), rep(as.integer(f$dim), 2L)) ||
       any(!is.finite(spatial$weights)) || length(spatial$lower) != 1L ||
       length(spatial$upper) != 1L || !is.finite(spatial$lower) ||
       !is.finite(spatial$upper) || spatial$lower >= spatial$upper){
      stop("Invalid normalized spatial precision capability.", call. = FALSE)
    }
    if(identical(f$representation$form, "normalized_affine_precision") &&
       (length(spatial$degree) != f$dim || any(!is.finite(spatial$degree)) || any(spatial$degree <= 0))){
      stop("Affine spatial precision requires positive finite degrees.", call. = FALSE)
    }
  }
  if(!is.list(f$native_report) || !identical(f$native_report$backend, "R") ||
     !is.function(f$native_report$fun)){
    stop("CovarianceFactor native_report must contain an R callback function.",
         call. = FALSE)
  }
  invisible(TRUE)
}
