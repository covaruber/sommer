
ar1m <- function(x, rho=0.30, fixed=FALSE,
                 variance=c("homogeneous", "heterogeneous"), values=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  variance <- match.arg(variance)
  if(q < 2L) stop("ar1m() requires at least two ordered levels.", call. = FALSE)
  levNum <- suppressWarnings(as.numeric(colnames(dummy)))
  if(!anyNA(levNum) && is.unsorted(levNum)){
    warning("ar1m(", expr, "): numeric levels are not in increasing order; AR1 lags follow levels(), ",
            "so reorder with factor(x, levels=sort(unique(as.numeric(as.character(x))))).",
            call. = FALSE)
  }
  if(length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1){
    stop("rho in ar1m() must be finite and strictly between -1 and 1.", call. = FALSE)
  }
  if(variance == "heterogeneous"){
    ar_fixed <- if(missing(fixed)) NULL else fixed
    return(.arp_m(dummy, order=1L, pacf=rho, fixed=ar_fixed, values=values,
                  heterogeneous=TRUE))
  }
  if(length(fixed) != 1L) stop("fixed in homogeneous ar1m() must have length one.", call. = FALSE)
  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="ar1",
      dim=as.integer(q),
      levels=colnames(dummy),
      par=c(atanh_rho=atanh(rho)),
      free=c(!isTRUE(fixed)),
      par_names="rho"
    ))
  )
}


# -------------------------------------------------------------------------
# Additional covariance-shape constructors.
#
# All structures are dimensionless shapes.  vsm() continues to own exactly
# one product-level sigma^2, preventing scale confounding in arbitrary
# Kronecker products.
# -------------------------------------------------------------------------

# Stable AR(p), p=2 or 3, parameterized through partial autocorrelations.
# PACF coordinates are tanh-transformed in C++, guaranteeing stationarity.
.ar_correlation_from_pacf <- function(pacf, q){
  order <- length(pacf)
  phi <- numeric()
  for(m in seq_len(order)){
    next_phi <- numeric(m)
    next_phi[m] <- pacf[m]
    if(m > 1L){
      next_phi[-m] <- phi - pacf[m] * rev(phi)
    }
    phi <- next_phi
  }

  A <- matrix(0, order, order)
  b <- numeric(order)
  for(k in seq_len(order)){
    A[k, k] <- A[k, k] + 1
    for(j in seq_len(order)){
      distance <- abs(k - j)
      if(distance == 0L){
        b[k] <- b[k] + phi[j]
      }else{
        A[k, distance] <- A[k, distance] - phi[j]
      }
    }
  }
  rho <- numeric(q)
  rho[1] <- 1
  rho[seq_len(order) + 1L] <- solve(A, b)
  if(q > order + 1L){
    for(h in (order + 1L):(q - 1L)){
      rho[h + 1L] <- sum(phi * rev(rho[(h - order + 1L):h]))
    }
  }
  stats::toeplitz(rho)
}

.arp_m <- function(x, order, pacf=NULL, fixed=NULL, values=NULL,
                   heterogeneous=FALSE){
  if(is.matrix(x) || inherits(x, "Matrix")){
    dummy <- x
  }else{
    expr <- as.character(substitute(x))
    dummy <- .cov_dummy(x, expr)
  }
  q <- ncol(dummy)

  order <- as.integer(order)
  if(length(order) != 1L || !order %in% c(1L,2L,3L)){
    stop("AR order must be 1, 2, or 3.", call. = FALSE)
  }
  if(q <= order){
    stop(sprintf("AR(%d) requires more than %d ordered levels.", order, order),
         call. = FALSE)
  }

  if(is.null(pacf)) pacf <- rep(0.20, order)
  if(length(pacf) != order || any(!is.finite(pacf)) || any(abs(pacf) >= 1)){
    stop("pacf must contain one finite value in (-1,1) per AR order.",
         call. = FALSE)
  }

  if(!heterogeneous){
    if(is.null(fixed)) fixed <- rep(FALSE, order)
    if(length(fixed) != order){
      stop("fixed must have length equal to the AR order.", call. = FALSE)
    }

    return(list(
      Z=dummy,
      covFactor=.compile_covfactor(list(
        type="arp",
        dim=as.integer(q),
        levels=colnames(dummy),
        order=order,
        par=atanh(pacf),
        free=!as.logical(fixed),
        par_names=paste0("pacf[", seq_len(order), "]")
      ))
    ))
  }

  if(is.null(values)) values <- rep(1, q)
  if(length(values) != q || any(!is.finite(values)) || any(values <= 0)){
    stop("values in heterogeneous AR() must contain one positive finite variance per level.",
         call. = FALSE)
  }
  par <- c(atanh(pacf), log((values / values[1])[-1]))
  if(is.null(fixed) || identical(fixed, FALSE)){
    fixed <- rep(FALSE, length(par))
  }
  if(identical(fixed, TRUE)){
    fixed <- rep(TRUE, length(par))
  }
  if(length(fixed) != length(par)){
    stop("fixed in heterogeneous AR() must have length order + q - 1.",
         call. = FALSE)
  }

  ar_eval <- local({
    q0 <- q
    order0 <- order
    function(eta){
      C <- .ar_correlation_from_pacf(tanh(eta[seq_len(order0)]), q0)
      variances <- c(1, exp(eta[order0 + seq_len(q0 - 1L)]))
      tcrossprod(sqrt(variances)) * C
    }
  })

  list(
    Z=dummy,
    covFactor=.make_covfactor(
      dim=as.integer(q),
      levels=colnames(dummy),
      par=par,
      free=!as.logical(fixed),
      par_names=c(paste0("pacf[", seq_len(order), "]"),
                  paste0("variance_ratio[", colnames(dummy)[-1], "]")),
      evaluator=list(backend="R", fun=ar_eval),
      derivative=list(backend="numeric", rel_step=1e-6),
      report=list(
        backend="builtin",
        transform=c(rep("tanh", order), rep("exp", q - 1L)),
        lower=rep(NA_real_, length(par)),
        upper=rep(NA_real_, length(par))
      ),
      trust_cap=rep(1, length(par)),
      model=paste0("ar", order),
      metadata=list(order=order, heterogeneous=TRUE)
    )
  )
}

ar2m <- function(x, pacf=c(0.20,0.10), fixed=NULL,
                 variance=c("homogeneous", "heterogeneous"), values=NULL){
  .arp_m(x, order=2L, pacf=pacf, fixed=fixed, values=values,
         heterogeneous=match.arg(variance) == "heterogeneous")
}

ar3m <- function(x, pacf=c(0.20,0.10,0.05), fixed=NULL,
                 variance=c("homogeneous", "heterogeneous"), values=NULL){
  .arp_m(x, order=3L, pacf=pacf, fixed=fixed, values=values,
         heterogeneous=match.arg(variance) == "heterogeneous")
}


# Stationary moving-average covariance of order 1 or 2.
#
# theta are the ordinary MA polynomial coefficients in
#   e_t + theta[1] e_{t-1} + theta[2] e_{t-2}.
# Any finite coefficients produce a valid covariance; invertibility is not
# required to define the covariance matrix.
mam <- function(x, order=1L, theta=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)

  order <- as.integer(order)
  if(length(order) != 1L || !order %in% c(1L,2L)){
    stop("mam() currently supports order=1 or order=2.", call. = FALSE)
  }
  if(q <= order){
    stop(sprintf("MA(%d) requires more than %d ordered levels.", order, order),
         call. = FALSE)
  }

  if(is.null(theta)) theta <- rep(0.15, order)
  if(length(theta) != order || any(!is.finite(theta))){
    stop("theta in mam() must contain one finite coefficient per MA order.",
         call. = FALSE)
  }

  if(is.null(fixed)) fixed <- rep(FALSE, order)
  if(length(fixed) != order){
    stop("fixed in mam() must have length equal to order.", call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="ma",
      dim=as.integer(q),
      levels=colnames(dummy),
      order=order,
      par=as.numeric(theta),
      free=!as.logical(fixed),
      par_names=paste0("ma[", seq_len(order), "]")
    ))
  )
}

ma1m <- function(x, theta=0.15, fixed=FALSE){
  mam(x, order=1L, theta=theta, fixed=fixed)
}

ma2m <- function(x, theta=c(0.15,0.05), fixed=NULL){
  mam(x, order=2L, theta=theta, fixed=fixed)
}


# Antedependence of order k using a modified-Cholesky representation
#
#   T y = e,  Cov(e)=D,
#   K = T^{-1} D T^{-T},
#
# with T unit lower triangular and regression coefficients only in the
# requested k subdiagonals. D[1,1]=1 is the scale reference.
antem <- function(x, order=1L, beta=NULL, innovations=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)

  order <- as.integer(order)
  if(length(order) != 1L || order < 1L || order >= q){
    stop("order in antem() must satisfy 1 <= order < number of levels.",
         call. = FALSE)
  }

  rows <- integer()
  cols <- integer()
  for(i in 2:q){
    j0 <- max(1L, i-order)
    for(j in j0:(i-1L)){
      rows <- c(rows, i)
      cols <- c(cols, j)
    }
  }

  ncoef <- length(rows)

  if(is.null(beta)) beta <- rep(0.10, ncoef)
  if(length(beta) != ncoef || any(!is.finite(beta))){
    stop("beta in antem() has the wrong length or contains non-finite values.",
         call. = FALSE)
  }

  if(is.null(innovations)) innovations <- rep(1, q)
  if(length(innovations) != q || any(!is.finite(innovations)) ||
     any(innovations <= 0)){
    stop("innovations in antem() must contain q positive finite values.",
         call. = FALSE)
  }

  ratios <- innovations / innovations[1]
  par <- c(as.numeric(beta), log(ratios[-1]))
  par_names <- c(
    paste0("ante[", labs[rows], "<-", labs[cols], "]"),
    paste0("innovation_ratio[", labs[-1], "]")
  )

  if(is.null(fixed)) fixed <- rep(FALSE, length(par))
  if(length(fixed) != length(par)){
    stop("fixed in antem() must match the antedependence parameter count.",
         call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="ante",
      dim=as.integer(q),
      levels=labs,
      order=order,
      par=par,
      free=!as.logical(fixed),
      par_names=par_names,
      ante_ncoef=as.integer(ncoef),
      ante_row=rows,
      ante_col=cols
    ))
  )
}

# -------------------------------------------------------------------------
# General positive-definite Toeplitz correlation covariance.
#
# A q x q Toeplitz correlation matrix is parameterized by q-1 reflection
# coefficients / partial autocorrelations in (-1,1).  This gives an
# unconstrained working parameterization through atanh(pacf) while preserving
# positive definiteness.
# -------------------------------------------------------------------------
toeplitzm <- function(x, pacf=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  if(q < 2L){
    stop("toeplitzm() requires at least two ordered levels.", call. = FALSE)
  }

  p <- q-1L
  if(is.null(pacf)) pacf <- rep(0.10, p)
  if(length(pacf) != p || any(!is.finite(pacf)) || any(abs(pacf) >= 1)){
    stop("pacf in toeplitzm() must contain q-1 finite values strictly inside (-1,1).",
         call. = FALSE)
  }
  if(is.null(fixed)) fixed <- rep(FALSE, p)
  if(length(fixed) != p){
    stop("fixed in toeplitzm() must have length q-1.", call. = FALSE)
  }

  eval_fun <- local({
    q0 <- q
    function(par){
      kappa <- tanh(par)
      phi <- numeric()
      for(m in seq_along(kappa)){
        nxt <- numeric(m)
        nxt[m] <- kappa[m]
        if(m > 1L){
          for(j in seq_len(m-1L)){
            nxt[j] <- phi[j] - kappa[m] * phi[m-j]
          }
        }
        phi <- nxt
      }

      ord <- length(phi)
      A <- matrix(0, ord, ord)
      b <- numeric(ord)
      for(kk in seq_len(ord)){
        A[kk,kk] <- A[kk,kk] + 1
        for(jj in seq_len(ord)){
          d <- abs(kk-jj)
          if(d == 0L) b[kk] <- b[kk] + phi[jj]
          else A[kk,d] <- A[kk,d] - phi[jj]
        }
      }
      rho <- c(1, as.numeric(solve(A,b)))
      stats::toeplitz(rho[seq_len(q0)])
    }
  })

  cf <- .make_covfactor(
    dim=q,
    levels=colnames(dummy),
    par=atanh(pacf),
    free=!as.logical(fixed),
    par_names=paste0("pacf[", seq_len(p), "]"),
    evaluator=list(backend="R", fun=eval_fun),
    derivative=list(backend="numeric", rel_step=1e-6),
    report=list(backend="builtin",
                transform=rep("tanh",p),
                lower=rep(NA_real_,p),
                upper=rep(NA_real_,p)),
    trust_cap=rep(1.0,p),
    structurally_diagonal=FALSE,
    model="toeplitz"
  )

  list(Z=dummy, covFactor=cf)
}
