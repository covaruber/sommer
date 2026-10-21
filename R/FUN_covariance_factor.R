

# Factor analytic covariance shape of rank k:
#   M = Lambda Lambda' + Psi,  Psi diagonal positive,
# then K = M / M[1,1].
#
# Lambda is lower-triangular in its leading k x k block for rotational
# identification. Leading loading diagonals are positive (log parameterized).
fam <- function(x, k=1L, loadings=NULL, specific=NULL, fixed=NULL,
                variance=c("heterogeneous", "homogeneous"), shared=NULL){
  variance <- match.arg(variance)
  if(!is.null(shared) && (!is.character(shared) || length(shared) != 1L ||
     is.na(shared) || !nzchar(shared))){
    stop("shared in fam() must be NULL or one non-empty loading-group name.", call.=FALSE)
  }
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)

  k <- as.integer(k)
  if(length(k) != 1L || k < 1L || k >= q){
    stop("k in fam() must satisfy 1 <= k < number of levels.", call. = FALSE)
  }

  rows <- integer()
  cols <- integer()
  diag_loading <- logical()

  for(j in seq_len(k)){
    for(i in j:q){
      rows <- c(rows, i)
      cols <- c(cols, j)
      diag_loading <- c(diag_loading, i == j)
    }
  }

  nload <- length(rows)

  if(is.null(loadings)){
    L <- matrix(0, q, k)
    for(j in seq_len(k)){
      L[j,j] <- 0.60
      if(j < q) L[(j+1L):q,j] <- 0.10
    }
  }else{
    L <- as.matrix(loadings)
    if(!all(dim(L) == c(q,k)) || any(!is.finite(L))){
      stop("loadings in fam() must be a finite q x k matrix.", call. = FALSE)
    }
    for(j in seq_len(k)){
      if(j > 1L) L[seq_len(j-1L),j] <- 0
      if(L[j,j] <= 0){
        stop("Leading diagonal loadings in fam() must be positive.", call. = FALSE)
      }
    }
  }

  if(is.null(specific)) specific <- rep(0.50, q)
  if(variance == "homogeneous" && length(specific) == 1L){
    specific <- rep(specific, q)
  }
  if(length(specific) != q || any(!is.finite(specific)) || any(specific <= 0)){
    stop("specific in fam() must contain q positive finite specific variances (or one value when variance='homogeneous').",
         call. = FALSE)
  }
  if(variance == "homogeneous" && any(specific != specific[1L])){
    stop("specific in homogeneous fam() must be a single value or q equal values.",
         call. = FALSE)
  }

  # Remove the otherwise redundant internal FA scale.  The first specific
  # variance is the reference (=1); loadings are expressed in its SD units.
  specific_scale <- specific[1]
  specific <- specific / specific_scale
  L <- L / sqrt(specific_scale)

  load_par <- numeric(nload)
  load_names <- character(nload)

  for(a in seq_len(nload)){
    i <- rows[a]; j <- cols[a]
    if(diag_loading[a]){
      load_par[a] <- log(L[i,j])
      load_names[a] <- paste0("loading_diag[", labs[i], ",F", j, "]")
    }else{
      load_par[a] <- L[i,j]
      load_names[a] <- paste0("loading[", labs[i], ",F", j, "]")
    }
  }

  par <- if(variance == "heterogeneous") c(load_par, log(specific[-1])) else load_par
  par_names <- if(variance == "heterogeneous"){
    c(load_names, paste0("specific_ratio[", labs[-1], "]"))
  }else load_names

  if(is.null(fixed)) fixed <- rep(FALSE, length(par))
  if(length(fixed) != length(par)){
    stop("fixed in fam() must match the number of loading plus specific-variance parameters.",
         call. = FALSE)
  }

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="fa",
      dim=as.integer(q),
      levels=labs,
      order=k,
      par=par,
      free=!as.logical(fixed),
      par_names=par_names,
      fa_nload=as.integer(nload),
      fa_row=rows,
      fa_col=cols,
      fa_diag=diag_loading,
      variance=variance,
      shared=shared,
      fa_specific_reference=labs[1]
    ))
  )
}



# Pure reduced-rank covariance: K = Lambda Lambda', with Lambda[1,1]=1.
# The reference loading removes the internal scale already supplied by vsm().
rrm <- function(x, k=1L, loadings=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)

  if(length(k) != 1L || !is.finite(k) || k != as.integer(k) || k < 1L || k >= q){
    stop("k in rrm() must satisfy 1 <= k < number of levels.", call. = FALSE)
  }
  k <- as.integer(k)

  rows <- integer()
  cols <- integer()
  diag_loading <- logical()

  for(j in seq_len(k)){
    for(i in j:q){
      rows <- c(rows, i)
      cols <- c(cols, j)
      diag_loading <- c(diag_loading, i == j)
    }
  }

  nload <- length(rows)

  if(is.null(loadings)){
    L <- matrix(0, q, k)
    for(j in seq_len(k)){
      L[j,j] <- 0.60
      if(j < q) L[(j+1L):q,j] <- 0.10
    }
  }else{
    L <- as.matrix(loadings)
    if(!all(dim(L) == c(q,k)) || any(!is.finite(L))){
      stop("loadings in rrm() must be a finite q x k matrix.", call. = FALSE)
    }
    for(j in seq_len(k)){
      if(j > 1L) L[seq_len(j-1L),j] <- 0
      if(L[j,j] <= 0){
        stop("Leading diagonal loadings in rrm() must be positive.", call. = FALSE)
      }
    }
  }

  L <- L / L[1L, 1L]
  vals <- numeric(nload)
  nm <- character(nload)

  for(a in seq_len(nload)){
    i <- rows[a]; j <- cols[a]
    if(diag_loading[a]){
      vals[a] <- log(L[i,j])
      nm[a] <- paste0("loading_diag[", labs[i], ",F", j, "]")
    }else{
      vals[a] <- L[i,j]
      nm[a] <- paste0("loading[", labs[i], ",F", j, "]")
    }
  }

  if(is.null(fixed)) fixed <- rep(FALSE, length(vals))
  if(length(fixed) != length(vals)){
    stop("fixed in rrm() must match the number of loading parameters.", call. = FALSE)
  }
  if(!is.logical(fixed) || anyNA(fixed)){
    stop("fixed in rrm() must be a logical vector without missing values.", call. = FALSE)
  }
  fixed[1L] <- TRUE

  list(
    Z=dummy,
    covFactor=.compile_covfactor(list(
      type="rr",
      dim=as.integer(q),
      levels=labs,
      order=k,
      par=vals,
      free=!as.logical(fixed),
      par_names=nm,
      rr_nload=as.integer(nload),
      rr_row=rows,
      rr_col=cols,
      rr_diag=diag_loading,
      zero_specific=TRUE
    ))
  )
}
