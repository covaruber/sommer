
# =========================================================================
# Additional generic CovarianceFactor structures
#   maternm(), toeplitzm(), sar(), car()
#
# These constructors intentionally use the universal CovarianceFactor v2
# callback interface.  No ai_mme_sp2() change is required.
# =========================================================================

# Build an incidence matrix for unique spatial coordinates and return the
# coordinate matrix in exactly the same column order as the incidence matrix.
.spatial_cov_dummy <- function(x, expr){
  if(is.data.frame(x)) x <- as.matrix(x)

  if(is.vector(x) && !is.list(x)){
    x <- matrix(as.numeric(x), ncol=1)
  }else{
    x <- as.matrix(x)
  }

  if(!is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L || any(!is.finite(x))){
    stop("Spatial coordinates must be a finite numeric vector, matrix, or data frame.",
         call. = FALSE)
  }

  # unique.data.frame() preserves first occurrence, which keeps model ordering
  # deterministic and consistent with the observation design.
  ux <- unique(as.data.frame(x, check.names=FALSE))
  coords <- as.matrix(ux)
  storage.mode(coords) <- "double"

  key <- function(M){
    apply(M, 1L, function(z) paste(format(z, digits=17, scientific=FALSE,
                                           trim=TRUE), collapse="\r"))
  }
  lev_key <- key(coords)
  obs_key <- key(x)
  idx <- match(obs_key, lev_key)
  if(anyNA(idx)){
    stop("Internal spatial-coordinate matching failure.", call. = FALSE)
  }

  Z <- Matrix::sparseMatrix(i=seq_len(nrow(x)), j=idx, x=1,
                            dims=c(nrow(x), nrow(coords)))

  if(ncol(coords) == 1L){
    labs <- format(coords[,1], digits=12, trim=TRUE)
  }else{
    labs <- apply(coords, 1L, function(z) paste(format(z, digits=10,
                                                        trim=TRUE), collapse=":"))
  }
  colnames(Z) <- make.unique(as.character(labs))
  Z <- to_sparse(Z)

  list(Z=Z, coords=coords, levels=colnames(Z), expr=expr)
}

# Align a user-supplied spatial weights/adjacency matrix to the levels of x.
# Named matrices are reordered by level name.  Unnamed matrices are assumed
# already to follow the incidence-column order.
.align_spatial_weights <- function(W, levels, symmetric=FALSE,
                                   nonnegative=FALSE,
                                   zero_diagonal=FALSE,
                                   name="W"){
  W <- as.matrix(W)
  q <- length(levels)
  if(!all(dim(W) == c(q,q)) || any(!is.finite(W))){
    stop(sprintf("%s must be a finite %d x %d matrix.", name, q, q),
         call. = FALSE)
  }

  if(!is.null(rownames(W)) || !is.null(colnames(W))){
    if(is.null(rownames(W)) || is.null(colnames(W))){
      stop(sprintf("If %s has dimnames, both row and column names are required.", name),
           call. = FALSE)
    }
    miss <- union(setdiff(levels, rownames(W)), setdiff(levels, colnames(W)))
    if(length(miss)){
      stop(sprintf("Levels missing from %s: %s", name, paste(miss, collapse=", ")),
           call. = FALSE)
    }
    W <- W[levels, levels, drop=FALSE]
  }

  if(symmetric && max(abs(W-t(W))) > 1e-10){
    stop(sprintf("%s must be symmetric for this covariance structure.", name),
         call. = FALSE)
  }
  if(nonnegative && any(W < -1e-12)){
    stop(sprintf("%s must contain non-negative adjacency weights.", name),
         call. = FALSE)
  }
  if(zero_diagonal && any(abs(diag(W)) > 1e-12)){
    stop(sprintf("%s must have a zero diagonal.", name), call. = FALSE)
  }

  W
}

# -------------------------------------------------------------------------
# Matern correlation covariance shape.
#
# x can be a numeric coordinate vector or an n x d coordinate matrix/data
# frame. Repeated coordinates share one covariance-product level.
#
# K(d) = 2^(1-nu)/Gamma(nu) * z^nu * BesselK_nu(z),
# z = sqrt(2*nu) * d/range, with K(0)=1.
# -------------------------------------------------------------------------
maternm <- function(x, range=NULL, nu=0.5, fixed=c(FALSE,FALSE),
                    distance=NULL, anisotropy=c("none", "geometric"),
                    angle=pi/2, ratio=1.5){
  anisotropy <- match.arg(anisotropy)
  expr <- paste(deparse(substitute(x)), collapse="")
  sx <- .spatial_cov_dummy(x, expr)
  Z <- sx$Z
  coords <- sx$coords
  q <- ncol(Z)
  geometric <- anisotropy == "geometric"
  if(geometric && (ncol(coords) != 2L || !is.null(distance))){
    stop("Geometric anisotropy in maternm() needs two-dimensional coordinates and no distance matrix.",
         call. = FALSE)
  }

  if(q < 2L){
    stop("maternm() requires at least two distinct spatial locations.", call. = FALSE)
  }

  if(is.null(distance)){
    D <- as.matrix(stats::dist(coords))
  }else{
    D <- as.matrix(distance)
    if(!all(dim(D) == c(q,q)) || any(!is.finite(D)) ||
       max(abs(D-t(D))) > 1e-10 || any(D < -1e-12) ||
       any(abs(diag(D)) > 1e-10)){
      stop("distance in maternm() must be a finite symmetric q x q distance matrix with zero diagonal.",
           call. = FALSE)
    }
  }

  posd <- D[D > 0]
  if(!length(posd)){
    stop("maternm() requires at least one positive inter-location distance.",
         call. = FALSE)
  }
  if(is.null(range)) range <- stats::median(posd)

  if(length(range) != 1L || !is.finite(range) || range <= 0){
    stop("range in maternm() must be one positive finite value.", call. = FALSE)
  }
  if(length(nu) != 1L || !is.finite(nu) || nu <= 0){
    stop("nu in maternm() must be one positive finite value.", call. = FALSE)
  }
  nPar <- if(geometric) 4L else 2L
  if(length(fixed) == 1L) fixed <- rep(fixed, nPar)
  if(length(fixed) == 2L && geometric) fixed <- c(fixed, FALSE, FALSE)
  if(length(fixed) != nPar){
    stop("fixed in maternm() must have length 1 or ", nPar,
         if(geometric) " (range, nu, angle, ratio)." else " (range, nu).",
         call. = FALSE)
  }
  if(geometric){
    if(length(angle) != 1L || !is.finite(angle) || angle <= 0 || angle >= pi){
      stop("angle in maternm() must lie strictly between 0 and pi.", call. = FALSE)
    }
    if(length(ratio) != 1L || !is.finite(ratio) || ratio <= 0){
      stop("ratio in maternm() must be one positive value.", call. = FALSE)
    }
  }

  eval_fun <- local({
    D0 <- D
    C0 <- coords
    function(par){
      r <- exp(par[1])
      v <- exp(par[2])
      if(geometric) D0 <- .geometric_distance(C0, pi*stats::plogis(par[3]), exp(par[4]))$d
      z <- sqrt(2*v) * D0 / r
      K <- matrix(1, nrow(D0), ncol(D0))
      use <- z > 1e-10
      if(any(use)){
        zu <- z[use]
        # Scaled Bessel K avoids premature underflow; exp(-z) restores K_nu.
        bk <- besselK(zu, nu=v, expon.scaled=TRUE)
        logc <- (1-v)*log(2) - lgamma(v) + v*log(zu) + log(bk) - zu
        val <- exp(logc)
        val[!is.finite(val)] <- 0
        K[use] <- val
      }
      diag(K) <- 1
      K <- (K+t(K))/2
      K
    }
  })

  maternPar <- c(log_range=log(range), log_nu=log(nu))
  maternNames <- c("range","nu")
  maternTransform <- c("exp","exp")
  maternLower <- c(NA_real_,NA_real_)
  maternUpper <- c(NA_real_,NA_real_)
  if(geometric){
    maternPar <- c(maternPar, eta_angle=stats::qlogis(angle/pi), log_ratio=log(ratio))
    maternNames <- c(maternNames, "angle", "ratio")
    maternTransform <- c(maternTransform, "bounded_logit", "exp")
    maternLower <- c(maternLower, 0, NA_real_)
    maternUpper <- c(maternUpper, pi, NA_real_)
  }

  cf <- .make_covfactor(
    dim=q,
    levels=sx$levels,
    par=maternPar,
    free=!as.logical(fixed),
    par_names=maternNames,
    evaluator=list(backend="R", fun=eval_fun),
    derivative=list(backend="numeric", rel_step=1e-6),
    report=list(backend="builtin",
                transform=maternTransform,
                lower=maternLower,
                upper=maternUpper),
    trust_cap=c(1.0,0.75,rep(1.0, nPar-2L)),
    structurally_diagonal=FALSE,
    model="matern",
    metadata=list(distance=D, coordinates=coords)
  )

  list(Z=Z, covFactor=cf)
}

# Rotated/stretched pairwise distances for geometric anisotropy; also returns
# the pieces needed for analytic derivatives with respect to angle and ratio.
.geometric_distance <- function(coords, angle, ratio){
  dx <- outer(coords[,1], coords[,1], "-")
  dy <- outer(coords[,2], coords[,2], "-")
  u <- dx*cos(angle) + dy*sin(angle)
  v <- -dx*sin(angle) + dy*cos(angle)
  d <- sqrt(u^2 + (v/ratio)^2)
  list(d=d, u=u, v=v)
}

.metric_shape <- function(model){
  switch(model,
    exponential=list(g=function(h) exp(-h), dg=function(h) -exp(-h)),
    gaussian=list(g=function(h) exp(-h^2), dg=function(h) -2*h*exp(-h^2)),
    spherical=list(g=function(h) ifelse(h < 1, 1 - 1.5*h + 0.5*h^3, 0),
                   dg=function(h) ifelse(h < 1, -1.5 + 1.5*h^2, 0)),
    circular=list(g=function(h){
                    hh <- pmin(h, 1)
                    ifelse(h < 1, (2/pi)*(acos(hh) - hh*sqrt(1 - hh^2)), 0)
                  },
                  dg=function(h) ifelse(h < 1, -(4/pi)*sqrt(1 - pmin(h, 1)^2), 0)))
}

# -------------------------------------------------------------------------
# Distance-based (metric) correlation structures: exponential, power,
# gaussian, spherical and circular, isotropic or anisotropic.
# -------------------------------------------------------------------------
metricm <- function(x, model=c("exponential", "power", "gaussian", "spherical", "circular"),
                    range=NULL, rho=NULL, anisotropy=c("none", "product", "geometric"),
                    angle=pi/2, ratio=1.5, metric=c("euclidean", "manhattan"),
                    distance=NULL, fixed=NULL){
  model <- match.arg(model)
  anisotropy <- match.arg(anisotropy)
  metric <- match.arg(metric)
  expr <- paste(deparse(substitute(x)), collapse="")
  sx <- .spatial_cov_dummy(x, expr)
  coords <- sx$coords
  q <- ncol(sx$Z)
  nd <- ncol(coords)
  if(q < 2L) stop("metricm() requires at least two distinct locations.", call. = FALSE)
  if(anisotropy != "none" && nd != 2L){
    stop("Anisotropic metricm() models require two-dimensional coordinates.", call. = FALSE)
  }
  if(anisotropy == "product" && !(model %in% c("exponential", "power", "gaussian"))){
    stop("Product anisotropy is available for the exponential, power and gaussian models.",
         call. = FALSE)
  }
  if(model == "spherical" && nd > 3L){
    stop("The spherical model is only valid in up to three dimensions.", call. = FALSE)
  }
  if(model == "circular" && nd > 2L){
    stop("The circular model is only valid in up to two dimensions.", call. = FALSE)
  }
  if(metric == "manhattan" && (anisotropy != "none" || !(model %in% c("exponential", "power")))){
    stop("The manhattan metric is only valid for isotropic exponential and power models.",
         call. = FALSE)
  }
  if(!is.null(distance) && anisotropy != "none"){
    stop("A user distance matrix can only be used with anisotropy='none'.", call. = FALSE)
  }

  if(!is.null(distance)){
    D <- as.matrix(distance)
    if(!all(dim(D) == c(q,q)) || any(!is.finite(D)) || max(abs(D-t(D))) > 1e-10 ||
       any(D < -1e-12) || any(abs(diag(D)) > 1e-10)){
      stop("distance in metricm() must be a finite symmetric q x q distance matrix with zero diagonal.",
           call. = FALSE)
    }
  }else{
    D <- as.matrix(stats::dist(coords, method=metric))
  }
  posd <- D[D > 0]
  if(!length(posd)) stop("metricm() requires a positive inter-location distance.", call. = FALSE)
  med <- stats::median(posd)
  isPower <- model == "power"

  # one parameter per axis for product anisotropy, otherwise one
  nAxis <- if(anisotropy == "product") 2L else 1L
  axisDist <- if(anisotropy == "product"){
    list(abs(outer(coords[,1], coords[,1], "-")), abs(outer(coords[,2], coords[,2], "-")))
  }else NULL
  if(isPower){
    if(is.null(rho)) rho <- rep(exp(-1/med), nAxis)
    if(length(rho) == 1L) rho <- rep(rho, nAxis)
    if(length(rho) != nAxis || any(!is.finite(rho)) || any(rho <= 0 | rho >= 1)){
      stop("rho in metricm() must lie strictly between 0 and 1 (one value per axis).",
           call. = FALSE)
    }
    base <- stats::qlogis(rho)
    baseNames <- if(nAxis == 2L) c("rho_x", "rho_y") else "rho"
    baseReport <- rep("bounded_logit", nAxis)
  }else{
    if(is.null(range)) range <- rep(med, nAxis)
    if(length(range) == 1L) range <- rep(range, nAxis)
    if(length(range) != nAxis || any(!is.finite(range)) || any(range <= 0)){
      stop("range in metricm() must be positive (one value per axis).", call. = FALSE)
    }
    base <- log(range)
    baseNames <- if(nAxis == 2L) c("range_x", "range_y") else "range"
    baseReport <- rep("exp", nAxis)
  }
  par <- base
  parNames <- baseNames
  transform <- baseReport
  lower <- ifelse(baseReport == "bounded_logit", 0, NA_real_)
  upper <- ifelse(baseReport == "bounded_logit", 1, NA_real_)
  if(anisotropy == "geometric"){
    if(length(angle) != 1L || !is.finite(angle) || angle <= 0 || angle >= pi){
      stop("angle in metricm() must lie strictly between 0 and pi.", call. = FALSE)
    }
    if(length(ratio) != 1L || !is.finite(ratio) || ratio <= 0){
      stop("ratio in metricm() must be one positive value.", call. = FALSE)
    }
    par <- c(par, stats::qlogis(angle/pi), log(ratio))
    parNames <- c(parNames, "angle", "ratio")
    transform <- c(transform, "bounded_logit", "exp")
    lower <- c(lower, 0, NA_real_)
    upper <- c(upper, pi, NA_real_)
  }
  if(is.null(fixed)) fixed <- rep(FALSE, length(par))
  if(length(fixed) == 1L) fixed <- rep(fixed, length(par))
  if(length(fixed) != length(par)){
    stop("fixed in metricm() must have one value per parameter: ",
         paste(parNames, collapse=", "), call. = FALSE)
  }
  shape <- if(isPower) NULL else .metric_shape(model)

  # K and dK/dpar for the current working parameters.
  kernel <- function(par, k=0L){
    if(anisotropy == "product"){
      if(isPower){
        r <- stats::plogis(par[1:2])
        K <- r[1]^axisDist[[1]] * r[2]^axisDist[[2]]
        if(k == 0L) return(K)
        return(K * axisDist[[k]] * (1 - r[k]))
      }
      ph <- exp(par[1:2])
      h <- lapply(1:2, function(a) axisDist[[a]]/ph[a])
      if(model == "exponential"){
        K <- exp(-h[[1]] - h[[2]])
        if(k == 0L) return(K)
        return(K * h[[k]])
      }
      K <- exp(-h[[1]]^2 - h[[2]]^2)
      if(k == 0L) return(K)
      return(K * 2 * h[[k]]^2)
    }
    if(anisotropy == "geometric"){
      a <- pi*stats::plogis(par[2]); lam <- exp(par[3])
      gd <- .geometric_distance(coords, a, lam)
      d <- gd$d
    }else{
      d <- D
    }
    if(isPower){
      r <- stats::plogis(par[1])
      K <- r^d
      dKdd <- K*log(r)
      dK1 <- K*d*(1 - r)
    }else{
      ph <- exp(par[1])
      h <- d/ph
      K <- shape$g(h)
      dKdd <- shape$dg(h)/ph
      dK1 <- -shape$dg(h)*h
    }
    if(k == 0L){ diag(K) <- 1; return(K) }
    if(k == 1L){ diag(dK1) <- 0; return(dK1) }
    safe <- d > 0
    dd <- matrix(0, q, q)
    if(k == 2L){
      sa <- stats::plogis(par[2])
      dd[safe] <- (gd$u*gd$v*(1 - 1/lam^2))[safe]/d[safe] * pi*sa*(1 - sa)
    }else{
      dd[safe] <- -(gd$v^2/lam^2)[safe]/d[safe]
    }
    dKdd * dd
  }

  K0 <- kernel(par)
  if(model == "gaussian" && rcond(K0) < 1e-10){
    warning("The starting gaussian correlation matrix is nearly singular; ",
            "keep a residual (nugget) term or use a smaller starting range.", call. = FALSE)
  }

  cf <- .make_covfactor(
    dim=q,
    levels=sx$levels,
    par=stats::setNames(par, parNames),
    free=!as.logical(fixed),
    par_names=parNames,
    evaluator=list(backend="R", fun=function(par) kernel(par)),
    derivative=list(backend="R", fun=function(par, k) kernel(par, as.integer(k))),
    report=list(backend="builtin", transform=transform, lower=lower, upper=upper),
    trust_cap=rep(1.0, length(par)),
    structurally_diagonal=FALSE,
    model="metric",
    metadata=list(distance=D, coordinates=coords, metric_model=model,
                  anisotropy=anisotropy)
  )

  list(Z=sx$Z, covFactor=cf)
}

# -------------------------------------------------------------------------
# Simultaneous autoregressive (SAR) covariance.
#
#   B = I - rho W
#   M = B^{-1} B^{-T}
#   K = M / M[1,1]
#
# rho is restricted to (-1/r(W), 1/r(W)), where r(W) is the spectral radius.
# This conservative interval guarantees B is nonsingular for any real square W.
# -------------------------------------------------------------------------
sar <- function(x, W, rho=0.10, fixed=FALSE){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  if(q < 2L) stop("sar() requires at least two spatial levels.", call. = FALSE)

  W <- .align_spatial_weights(W, colnames(dummy), symmetric=FALSE,
                              nonnegative=FALSE, zero_diagonal=FALSE,
                              name="W")
  ev <- eigen(W, only.values=TRUE)$values
  radius <- max(Mod(ev))
  if(!is.finite(radius) || radius <= 0){
    stop("W in sar() must have a positive finite spectral radius.", call. = FALSE)
  }
  bound <- 1/radius
  if(length(rho) != 1L || !is.finite(rho) || abs(rho) >= bound){
    stop(sprintf("rho in sar() must lie strictly between %.6g and %.6g.",
                 -bound, bound), call. = FALSE)
  }
  if(length(fixed) != 1L) stop("fixed in sar() must have length one.", call. = FALSE)

  # Symmetric bounded-logit map (-bound,bound) -> R.
  p0 <- (rho + bound)/(2*bound)
  eta <- qlogis(p0)

  eval_fun <- local({
    W0 <- W; b0 <- bound; q0 <- q
    function(par){
      s <- plogis(par[1])
      r <- -b0 + 2*b0*s
      A <- solve(diag(q0) - r*W0)
      M <- tcrossprod(A)
      M / M[1,1]
    }
  })

  d1_fun <- local({
    W0 <- W; b0 <- bound; q0 <- q
    function(par, k){
      if(k != 1L) stop("Invalid SAR derivative index.", call. = FALSE)
      s <- plogis(par[1])
      r <- -b0 + 2*b0*s
      dr <- 2*b0*s*(1-s)
      A <- solve(diag(q0) - r*W0)
      dA <- dr * A %*% W0 %*% A
      M <- tcrossprod(A)
      D <- dA %*% t(A) + A %*% t(dA)
      sc <- M[1,1]; dsc <- D[1,1]
      (D*sc - M*dsc)/(sc*sc)
    }
  })

  cf <- .make_covfactor(
    dim=q,
    levels=colnames(dummy),
    par=c(eta_rho=eta),
    free=c(!isTRUE(fixed)),
    par_names="rho",
    evaluator=list(backend="R", fun=eval_fun),
    derivative=list(backend="R", fun=d1_fun),
    report=list(backend="builtin",
                transform="bounded_logit",
                lower=-bound,
                upper=bound),
    trust_cap=1.0,
    structurally_diagonal=FALSE,
    model="sar",
    metadata=list(W=W, rho_bound=bound)
  )

  list(Z=dummy, covFactor=cf)
}

# -------------------------------------------------------------------------
# Proper conditional autoregressive (CAR) covariance.
#
# For a symmetric nonnegative adjacency matrix W with zero diagonal:
#   D = diag(rowSums(W))
#   Q = D - rho W
#   M = Q^{-1}
#   K = M / M[1,1]
#
# Using S = D^{-1/2} W D^{-1/2}, rho is mapped to the exact open interval
# for which I-rho*S (and hence Q) is positive definite.
# -------------------------------------------------------------------------
car <- function(x, W, rho=0.10, fixed=FALSE){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  if(q < 2L) stop("car() requires at least two spatial levels.", call. = FALSE)

  W <- .align_spatial_weights(W, colnames(dummy), symmetric=TRUE,
                              nonnegative=TRUE, zero_diagonal=TRUE,
                              name="W")
  rs <- rowSums(W)
  if(any(!is.finite(rs)) || any(rs <= 0)){
    stop("Every level in car() must have positive total adjacency weight (no isolated levels).",
         call. = FALSE)
  }
  S <- diag(1/sqrt(rs)) %*% W %*% diag(1/sqrt(rs))
  lam <- eigen(S, symmetric=TRUE, only.values=TRUE)$values
  lam_pos <- lam[lam > 1e-12]
  lam_neg <- lam[lam < -1e-12]
  upper <- if(length(lam_pos)) 1/max(lam_pos) else Inf
  lower <- if(length(lam_neg)) 1/min(lam_neg) else -Inf

  # A finite two-sided transform is desirable for optimizer/reporting.  For
  # ordinary undirected adjacency matrices S has both signs; if not, use a
  # conservative symmetric finite interval based on the spectral radius.
  if(!is.finite(lower) || !is.finite(upper)){
    rad <- max(abs(lam))
    if(!is.finite(rad) || rad <= 0){
      stop("Unable to determine a valid CAR dependence interval.", call. = FALSE)
    }
    lower <- -1/rad
    upper <-  1/rad
  }

  eps <- 1e-12 * max(1, abs(lower), abs(upper))
  if(length(rho) != 1L || !is.finite(rho) ||
     rho <= lower+eps || rho >= upper-eps){
    stop(sprintf("rho in car() must lie strictly between %.6g and %.6g.",
                 lower, upper), call. = FALSE)
  }
  if(length(fixed) != 1L) stop("fixed in car() must have length one.", call. = FALSE)

  p0 <- (rho-lower)/(upper-lower)
  eta <- qlogis(p0)
  Dmat <- diag(rs)

  eval_fun <- local({
    W0 <- W; D0 <- Dmat; lo0 <- lower; hi0 <- upper
    function(par){
      s <- plogis(par[1])
      r <- lo0 + (hi0-lo0)*s
      M <- solve(D0 - r*W0)
      M <- (M+t(M))/2
      M / M[1,1]
    }
  })

  d1_fun <- local({
    W0 <- W; D0 <- Dmat; lo0 <- lower; hi0 <- upper
    function(par, k){
      if(k != 1L) stop("Invalid CAR derivative index.", call. = FALSE)
      s <- plogis(par[1])
      r <- lo0 + (hi0-lo0)*s
      dr <- (hi0-lo0)*s*(1-s)
      M <- solve(D0 - r*W0)
      D <- dr * M %*% W0 %*% M
      M <- (M+t(M))/2
      D <- (D+t(D))/2
      sc <- M[1,1]; dsc <- D[1,1]
      (D*sc - M*dsc)/(sc*sc)
    }
  })

  cf <- .make_covfactor(
    dim=q,
    levels=colnames(dummy),
    par=c(eta_rho=eta),
    free=c(!isTRUE(fixed)),
    par_names="rho",
    evaluator=list(backend="R", fun=eval_fun),
    derivative=list(backend="R", fun=d1_fun),
    report=list(backend="builtin",
                transform="bounded_logit",
                lower=lower,
                upper=upper),
    trust_cap=1.0,
    structurally_diagonal=FALSE,
    model="car",
    metadata=list(W=W, degree=rs,
                  rho_lower=lower, rho_upper=upper)
  )

  list(Z=dummy, covFactor=cf)
}
