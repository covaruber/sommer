#### =========== ####
## MULTIVARIATE (LONG-FORMAT) HELPERS
#### =========== ####

stackTraits <- function(data, traits, keep=NULL, trait="trait", value="value",
                        record="record"){
  data <- as.data.frame(data)
  traits <- as.character(traits)
  missing <- setdiff(traits, names(data))
  if(length(missing)) stop("Traits not found in data: ", paste(missing, collapse=", "), call.=FALSE)
  if(length(traits) < 1L) stop("At least one trait is required.", call.=FALSE)
  if(is.null(keep)) keep <- setdiff(names(data), traits)
  clash <- intersect(c(trait, value, record), keep)
  if(length(clash)){
    stop("Column names already used in data: ", paste(clash, collapse=", "),
         "; choose other trait/value/record names.", call.=FALSE)
  }
  n <- nrow(data)
  rec <- if(!is.null(rownames(data)) && !anyDuplicated(rownames(data)) &&
            !identical(rownames(data), as.character(seq_len(n)))) rownames(data) else
              paste0("r", seq_len(n))
  out <- data[rep(seq_len(n), times=length(traits)), keep, drop=FALSE]
  out[[record]] <- factor(rep(rec, times=length(traits)), levels=rec)
  out[[trait]] <- factor(rep(traits, each=n), levels=traits)
  out[[value]] <- unlist(lapply(traits, function(t) data[[t]]), use.names=FALSE)
  rownames(out) <- NULL
  out
}

covmatrix_mmes <- function(object, term, se=TRUE, max.dim=500L, rel_step=1e-5){
  if(!inherits(object, "mmes")) stop("object must inherit from class 'mmes'.", call.=FALSE)
  terms <- names(object$covStruct)
  if(is.numeric(term)) term <- terms[term]
  i <- match(term, terms)
  if(is.na(i)) stop("Unknown term. Available: ", paste(terms, collapse=", "), call.=FALSE)
  cs <- object$covStruct[[i]]
  if(cs$dim > max.dim){
    stop("The covariance of this term has dimension ", cs$dim, " (> max.dim).", call.=FALSE)
  }
  if(any(vapply(cs$factors, function(f) !is.null(f$section_levels), logical(1)))){
    stop("covmatrix_mmes() does not support dsumm() sections.", call.=FALSE)
  }
  natural <- as.numeric(object$covPar[[i]])
  evalAt <- function(v){
    K <- evaluate_covstruct_cpp(cs, .natural_to_working(cs, v))
    dimnames(K) <- list(cs$levels, cs$levels)
    K
  }
  corOf <- function(K){
    s <- sqrt(diag(K))
    R <- K / outer(s, s)
    diag(R) <- 1
    R
  }
  V <- evalAt(natural)
  out <- list(covariance=V, correlation=corOf(V))
  if(se){
    counts <- vapply(object$covPar, length, integer(1))
    global <- sum(counts[seq_len(i - 1L)]) + seq_len(counts[i])
    W <- as.matrix(object$theta_se)[global, global, drop=FALSE]
    free <- as.logical(cs$free)
    W[!free, ] <- 0
    W[, !free] <- 0
    q <- nrow(V)
    Jc <- matrix(0, q*q, length(natural))
    Jr <- Jc
    for(k in which(free)){
      h <- rel_step * max(abs(natural[k]), 1e-3)
      up <- natural; dn <- natural
      up[k] <- up[k] + h; dn[k] <- dn[k] - h
      Vu <- evalAt(up); Vd <- evalAt(dn)
      Jc[, k] <- as.vector(Vu - Vd) / (2*h)
      Jr[, k] <- as.vector(corOf(Vu) - corOf(Vd)) / (2*h)
    }
    toMat <- function(J){
      M <- matrix(sqrt(pmax(diag(J %*% W %*% t(J)), 0)), q, q)
      dimnames(M) <- dimnames(V)
      M
    }
    out$covariance.se <- toMat(Jc)
    out$correlation.se <- toMat(Jr)
    diag(out$correlation.se) <- 0
  }
  out
}

# Observation -> level of factor fidx for product coordinates (1-based).
.factor_digit <- function(localIndex, factors, fidx){
  dims <- vapply(factors, function(f) as.integer(f$dim), integer(1))
  trailing <- rev(cumprod(rev(c(dims[-1L], 1L))))
  ((localIndex - 1L) %/% trailing[fidx]) %% dims[fidx] + 1L
}

.usm_par_from_K <- function(K, f){
  L <- t(chol(K))
  L <- L / L[1, 1]
  ifelse(f$us_diag, log(L[cbind(f$us_row, f$us_row)]), L[cbind(f$us_row, f$us_col)])
}

# Heterogeneous starting shapes for default usm()/dsm() factors from the
# mean squares of fixed-effect residuals within each level of the factor
# (mean squares, unlike var(), are invariant to eigen rotations).
.mv_start_shapes <- function(cs, localIndex, r){
  if(!length(cs$factors) || length(localIndex) != length(r)) return(cs)
  for(fidx in seq_along(cs$factors)){
    f <- cs$factors[[fidx]]
     if(!isTRUE(f$default_start) || !(f$model %in% c("us", "diag", "corgh")) || f$dim < 2L ||
       !is.null(f$section_levels)) next
    lev <- .factor_digit(localIndex, cs$factors, fidx)
    ok <- !is.na(lev)
    v <- vapply(seq_len(f$dim), function(l){
      x <- r[ok & lev == l]
      if(length(x) < 3L) NA_real_ else mean(x^2)
    }, numeric(1))
    if(any(!is.finite(v)) || any(v <= 0)) next
    ratio <- v / v[1L]
    s <- f$par_start:f$par_end
    if(f$model == "diag"){
      new <- log(ratio[-1L])
    }else if(f$model == "corgh"){
      new <- cs$par[s]
      new[f$ratio_index] <- log(ratio[-1L])
    }else{
      C <- matrix(0.10, f$dim, f$dim)
      diag(C) <- 1
      K <- sqrt(ratio) * t(sqrt(ratio) * C)
      new <- .usm_par_from_K(K, f)
    }
    keepFixed <- !as.logical(cs$free[s])
    new[keepFixed] <- cs$par[s][keepFixed]
    cs$par[s] <- new
    cs$factors[[fidx]]$par <- new
    if(isTRUE(cs$sigma2_is_default) && isTRUE(cs$free[1L])){
      cs$par[1L] <- cs$par[1L] - log(mean(ratio))
    }
  }
  cs
}
