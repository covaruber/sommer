

# -------------------------------------------------------------------------
# Descriptor-based covariance constructors for the Henderson/Kronecker engine.
# -------------------------------------------------------------------------

# Compound symmetry / uniform correlation:
#   K_ii = 1
#   K_ij = rho, i != j
# with rho mapped from an unconstrained working coordinate to the exact
# positive-definite interval (-1/(q-1), 1).
csm <- function(x, rho=0.10, fixed=FALSE,
                variance=c("homogeneous", "heterogeneous"), values=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)
  variance <- match.arg(variance)

  if(q < 2L){
    stop("csm() requires at least two levels.", call. = FALSE)
  }

  lo <- -1/(q-1)
  if(length(rho) != 1L || !is.finite(rho) || rho <= lo || rho >= 1){
    stop(sprintf("rho in csm() must lie strictly between %.6g and 1.", lo),
         call. = FALSE)
  }
  if(variance == "heterogeneous"){
    if(is.null(values)) values <- rep(1, q)
    if(length(values) != q || any(!is.finite(values)) || any(values <= 0)){
      stop("values in csm() must contain one positive finite variance per level.",
           call. = FALSE)
    }

    ratios <- values / values[1]
    eta_var <- log(ratios[-1])
    p <- (rho-lo)/(1-lo)
    eta_rho <- qlogis(p)
    par <- c(eta_rho, eta_var)

    if(is.null(fixed) || identical(fixed, FALSE)) fixed <- rep(FALSE, length(par))
    if(identical(fixed, TRUE)) fixed <- rep(TRUE, length(par))
    if(length(fixed) != length(par)){
      stop("fixed in heterogeneous csm() must have length q: rho plus q-1 variance ratios.",
           call. = FALSE)
    }

    return(list(
      Z=dummy,
      covFactor=.compile_covfactor(list(
        type="csm",
        variance="heterogeneous",
        dim=as.integer(q),
        levels=labs,
        par=par,
        free=!as.logical(fixed),
        par_names=c("rho", paste0("variance_ratio[", labs[-1], "]"))
      ))
    ))
  }

  if(length(fixed) != 1L){
    stop("fixed in homogeneous csm() must have length one.", call. = FALSE)
  }

  # inverse-logit map from (lo,1) to R
  p <- (rho-lo)/(1-lo)
  eta <- qlogis(p)

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
        type="csm",
        variance="homogeneous",
      dim=as.integer(q),
      levels=labs,
      par=c(eta_rho=eta),
      free=c(!isTRUE(fixed)),
      par_names="rho"
    ))
  )
}

# atm(): selected-level diagonal heterogeneity.
#
# The first level remains the product-scale reference (=1).  Among the
# remaining levels, only entries named in levs are estimable; all others are
# fixed at ratio 1.  If the first level is included in levs it remains the
# reference because vsm() owns the single overall scale.
atm <- function(x, levs, values=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  full_dummy <- .cov_dummy(x, expr)
  all_levels <- colnames(full_dummy)

  if(missing(levs) || is.null(levs)){
    stop("atm() requires levs: the levels to retain in this diagonal structure.",
         call. = FALSE)
  }

  if(is.numeric(levs)){
    if(any(levs < 1L | levs > ncol(full_dummy))){
      stop("Numeric levs in atm() are outside the available level range.",
           call. = FALSE)
    }
    levs <- all_levels[as.integer(levs)]
  }else{
    levs <- as.character(levs)
  }

  bad <- setdiff(levs, all_levels)
  if(length(bad)){
    stop(paste("Unknown levels in atm():", paste(bad, collapse=", ")),
         call. = FALSE)
  }

  # Preserve the user-requested level order. Observations belonging to levels
  # outside levs receive zeros in this factor design, matching the historical
  # purpose of atm(): fit the diagonal structure only for selected levels.
  dummy <- full_dummy[, levs, drop=FALSE]
  dummy <- to_sparse(dummy)
  q <- ncol(dummy)

  if(is.null(values)) values <- rep(1, q)
  if(length(values) != q || any(!is.finite(values)) || any(values <= 0)){
    stop("values in atm() must contain one positive finite value per selected level.",
         call. = FALSE)
  }

  ratios <- values / values[1]
  eta <- if(q > 1L) log(ratios[-1]) else numeric()

  if(is.null(fixed)) fixed <- rep(FALSE, length(eta))
  if(length(fixed) != length(eta)){
    stop("fixed in atm() must have length length(levs)-1.", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="diag",
      dim=as.integer(q),
      levels=colnames(dummy),
      par=eta,
      free=!as.logical(fixed),
      par_names=if(q > 1L)
        paste0("variance_ratio[", colnames(dummy)[-1], "]")
      else
        character(),
      selected_levels=levs,
      reference_level=colnames(dummy)[1]
    ))
  )
}

ism <- function(x){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="identity",
      dim=as.integer(q),
      levels=colnames(dummy),
      par=numeric(),
      free=logical(),
      par_names=character()
    ))
  )
}

dsm <- function(x, values=NULL, fixed=NULL, theta=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)

  if(!is.null(theta)){
    theta <- as.matrix(theta)
    if(!all(dim(theta) == c(q,q))){
      stop("theta supplied to dsm() must be q x q.", call. = FALSE)
    }
    if(any(abs(theta[row(theta) != col(theta)]) > 1e-10)){
      stop("theta supplied to dsm() must be diagonal.", call. = FALSE)
    }
    values <- diag(theta)
  }
  defaultStart <- is.null(values)
  if(is.null(values)) values <- rep(1, q)
  if(length(values) != q || any(!is.finite(values)) || any(values <= 0)){
    stop("values in dsm() must contain one positive finite value per level.", call. = FALSE)
  }

  # The first diagonal is the scale reference (=1). Remaining entries are
  # positive relative variances represented by log ratios.
  ratios <- values / values[1]
  eta <- if(q > 1L) log(ratios[-1]) else numeric()

  if(is.null(fixed)) fixed <- rep(FALSE, length(eta))
  if(length(fixed) != length(eta)){
    stop("fixed in dsm() must have length q-1.", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="diag",
      dim=as.integer(q),
      levels=colnames(dummy),
      par=eta,
      free=!as.logical(fixed),
      par_names=if(q > 1L) paste0("variance_ratio[", colnames(dummy)[-1], "]") else character(),
      default_start=defaultStart,
      ratio_index=seq_len(q - 1L)
    ))
  )
}

usm <- function(x, theta=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)

  if(is.null(theta)){
    K <- matrix(0.10, q, q)
    diag(K) <- 1
  }else{
    K <- as.matrix(theta)
    if(!all(dim(K) == c(q,q))){
      stop("theta supplied to usm() must be q x q.", call. = FALSE)
    }
    K <- (K + t(K))/2
    if(!is.finite(K[1,1]) || K[1,1] <= 0){
      stop("theta[1,1] supplied to usm() must be positive.", call. = FALSE)
    }
    K <- K / K[1,1]
  }

  ch <- try(chol(K), silent=TRUE)
  if(inherits(ch, "try-error")){
    stop("Initial theta supplied to usm() must be positive definite.", call. = FALSE)
  }
  L <- t(ch)
  # Numerical normalization guarantees L[1,1] == 1.
  L <- L / L[1,1]

  vals <- numeric()
  nm <- character()
  rows <- integer()
  cols <- integer()
  isdiag <- logical()

  for(i in seq_len(q)){
    for(j in seq_len(i)){
      if(i == 1L && j == 1L) next
      rows <- c(rows, i)
      cols <- c(cols, j)
      if(i == j){
        vals <- c(vals, log(L[i,i]))
        nm <- c(nm, paste0("chol_diag[", colnames(dummy)[i], "]"))
        isdiag <- c(isdiag, TRUE)
      }else{
        vals <- c(vals, L[i,j])
        nm <- c(nm, paste0("chol[", colnames(dummy)[i], ",", colnames(dummy)[j], "]"))
        isdiag <- c(isdiag, FALSE)
      }
    }
  }

  if(is.null(fixed)) fixed <- rep(FALSE, length(vals))
  if(length(fixed) != length(vals)){
    stop("fixed in usm() must have length q(q+1)/2 - 1.", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="us",
      dim=as.integer(q),
      levels=colnames(dummy),
      par=vals,
      free=!as.logical(fixed),
      par_names=nm,
      us_row=rows,
      us_col=cols,
      us_diag=isdiag,
      default_start=is.null(theta)
    ))
  )
}


# General positive-definite correlation matrix.
#
# We use an unconstrained unit-diagonal lower factor A and standardize
# A A' to correlation scale.  This spans the SPD correlation cone without
# requiring pairwise-correlation boundary repairs.
corgm <- function(x, theta=NULL, fixed=NULL,
                  variance=c("homogeneous", "heterogeneous"), values=NULL){
  variance <- match.arg(variance)
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)

  if(q < 2L){
    stop("corgm() requires at least two levels.", call. = FALSE)
  }

  rows <- integer()
  cols <- integer()
  for(i in 2:q){
    for(j in 1:(i-1L)){
      rows <- c(rows, i)
      cols <- c(cols, j)
    }
  }

  if(is.null(theta)){
    vals <- rep(0.05, length(rows))
  }else{
    K <- as.matrix(theta)
    if(!all(dim(K) == c(q,q))){
      stop("theta supplied to corgm() must be q x q.", call. = FALSE)
    }
    K <- (K+t(K))/2
    if(any(diag(K) <= 0) || any(!is.finite(K))){
      stop("theta supplied to corgm() must have positive finite diagonal.",
           call. = FALSE)
    }
    if(variance == "heterogeneous" && is.null(values)) values <- diag(K)
    K <- stats::cov2cor(K)
    ch <- try(chol(K), silent=TRUE)
    if(inherits(ch, "try-error")){
      stop("theta supplied to corgm() must be positive definite.", call. = FALSE)
    }
    A <- t(ch)
    A <- sweep(A, 1, diag(A), "/")
    vals <- mapply(function(i,j) A[i,j], rows, cols)
  }

  if(variance == "heterogeneous"){
    defaultStart <- is.null(values)
    if(is.null(values)) values <- rep(1, q)
    if(length(values) != q || any(!is.finite(values)) || any(values <= 0)){
      stop("values in corgm() must contain one positive variance per level.", call. = FALSE)
    }
    ratioPar <- log(values[-1] / values[1])
    nCor <- length(vals)
    if(is.null(fixed)) fixed <- rep(FALSE, nCor + q - 1L)
    if(length(fixed) != nCor + q - 1L){
      stop("fixed in heterogeneous corgm() must have length q(q-1)/2 + q-1 ",
           "(correlation coordinates, then variance ratios).", call. = FALSE)
    }
    corOf <- function(par){
      A <- diag(q)
      A[cbind(rows, cols)] <- par[seq_len(nCor)]
      S <- tcrossprod(A)
      s <- sqrt(diag(S))
      S / outer(s, s)
    }
    eval_fun <- function(par){
      d <- sqrt(c(1, exp(par[nCor + seq_len(q - 1L)])))
      d * t(d * corOf(par))
    }
    report_fun <- function(scale, par, factor, absorb_scale=TRUE){
      R <- corOf(par)
      multiplier <- if(absorb_scale) scale else 1
      prefix <- if(absorb_scale) "variance" else "relative_variance"
      c(stats::setNames(R[cbind(rows, cols)], paste0("correlation[", labs[rows], ",", labs[cols], "]")),
        stats::setNames(multiplier * c(1, par[nCor + seq_len(q - 1L)]),
                        paste0(prefix, "[", labs, "]")))
    }
    cf <- .make_covfactor(
      dim=q, levels=labs, par=c(as.numeric(vals), ratioPar),
      free=!as.logical(fixed),
      par_names=c(paste0("cor_chol[", labs[rows], ",", labs[cols], "]"),
                  paste0("variance_ratio[", labs[-1], "]")),
      evaluator=list(backend="R", fun=eval_fun),
      derivative=list(backend="numeric", rel_step=1e-6),
      report=list(backend="builtin",
                  transform=c(rep("identity", nCor), rep("exp", q - 1L)),
                  lower=rep(NA_real_, nCor + q - 1L), upper=rep(NA_real_, nCor + q - 1L)),
      native_report=report_fun,
      trust_cap=rep(1.0, nCor + q - 1L),
      structurally_diagonal=FALSE,
      model="corgh",
      metadata=list(ratio_index=nCor + seq_len(q - 1L), default_start=defaultStart)
    )
    return(list(Z=dummy, covFactor=cf))
  }

  if(is.null(fixed)) fixed <- rep(FALSE, length(vals))
  if(length(fixed) != length(vals)){
    stop("fixed in corgm() must have length q(q-1)/2.", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="corg",
      dim=as.integer(q),
      levels=labs,
      par=as.numeric(vals),
      free=!as.logical(fixed),
      par_names=paste0("cor_chol[", labs[rows], ",", labs[cols], "]"),
      corg_row=rows,
      corg_col=cols
    ))
  )
}


# User-defined covariance shape.
#
# Either:
#   ownm(x, K=<known SPD matrix>)
# or
#   ownm(x, fun=function(par) ..., par=..., dfun=function(par,k) ...)
#
# The wrapper automatically normalizes fun(par) by [1,1], so vsm() retains the
# unique overall sigma2. If dfun is omitted, ai_mme_sp2 uses a central numerical
# derivative of the normalized user function.
ownm <- function(x, K=NULL, fun=NULL, par=numeric(), fixed=NULL,
                 dfun=NULL, par_names=NULL, native_report=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)

  normalize_matrix <- function(M){
    M <- as.matrix(M)
    if(!all(dim(M) == c(q,q)) || any(!is.finite(M))){
      stop("ownm covariance function must return a finite q x q matrix.",
           call. = FALSE)
    }
    M <- (M+t(M))/2
    if(!is.finite(M[1,1]) || M[1,1] <= 0){
      stop("ownm covariance shape must have a positive [1,1] element.",
           call. = FALSE)
    }
    M / M[1,1]
  }

  if(!is.null(K)){
    if(!is.null(fun) || length(par)){
      stop("When K is supplied to ownm(), do not also supply fun or par.",
           call. = FALSE)
    }
    K0 <- normalize_matrix(K)
    if(any(eigen(K0, symmetric=TRUE, only.values=TRUE)$values <= 0)){
      stop("K supplied to ownm() must be positive definite.", call. = FALSE)
    }
    return(list(
      Z=dummy,
      covFactor=.compile_covfactor(list(
        type="own",
        dim=as.integer(q),
        levels=labs,
        par=numeric(),
        free=logical(),
        par_names=character(),
        fixed_matrix=K0
      ))
    ))
  }

  if(!is.function(fun)){
    stop("ownm() requires either K or a covariance function fun(par).",
         call. = FALSE)
  }

  par <- as.numeric(par)
  if(any(!is.finite(par))){
    stop("par in ownm() must be finite.", call. = FALSE)
  }

  raw_fun <- fun
  wrapped_fun <- function(p){
    M <- raw_fun(p)
    M <- as.matrix(M)
    M <- (M+t(M))/2
    M / M[1,1]
  }

  wrapped_dfun <- NULL
  if(!is.null(dfun)){
    if(!is.function(dfun)){
      stop("dfun in ownm() must be NULL or a function(par, k).", call. = FALSE)
    }
    raw_dfun <- dfun
    wrapped_dfun <- function(p, k){
      M <- as.matrix(raw_fun(p))
      M <- (M+t(M))/2
      D <- as.matrix(raw_dfun(p, k))
      D <- (D+t(D))/2
      s <- M[1,1]
      ds <- D[1,1]
      (D*s - M*ds)/(s*s)
    }
  }

  K0 <- normalize_matrix(raw_fun(par))
  if(any(eigen(K0, symmetric=TRUE, only.values=TRUE)$values <= 0)){
    stop("Initial ownm() covariance shape must be positive definite.",
         call. = FALSE)
  }

  if(is.null(fixed)) fixed <- rep(FALSE, length(par))
  if(length(fixed) != length(par)){
    stop("fixed in ownm() must have length(par).", call. = FALSE)
  }

  if(is.null(par_names)) par_names <- paste0("own[", seq_along(par), "]")
  if(length(par_names) != length(par)){
    stop("par_names in ownm() must have length(par).", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="own",
      dim=as.integer(q),
      levels=labs,
      par=par,
      free=!as.logical(fixed),
      par_names=as.character(par_names),
      fun=wrapped_fun,
      dfun=wrapped_dfun,
      native_report=native_report
    ))
  )
}
