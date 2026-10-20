#### =========== ####
## MARKER EFFECTS FROM A FITTED RELATIONSHIP-BASED (GBLUP) TERM
#### =========== ####

# Marker transform used by A.mat(): imputation, MAF and monomorphic filters,
# centring by 2p and the VanRaden scale 2*sum(p*q).
.amat_marker_transform <- function(M, min.MAF){
  M <- .relationship_impute(M)
  p <- colMeans(M + 1) / 2
  maf <- pmin(p, 1 - p)
  keep <- maf > min.MAF
  keep[keep] <- apply(M[, keep, drop=FALSE], 2, stats::var) > 0
  pk <- p[keep]
  Z <- sweep(M[, keep, drop=FALSE] + 1, 2, 2 * pk)
  list(Z=Z, keep=keep, k=2 * sum(pk * (1 - pk)))
}

meffects_mmes <- function(object, term, M, method=c("A.mat", "custom"),
                          min.MAF=0, Z=NULL, scale=NULL, blend=0, Gu=NULL,
                          se=FALSE, chunk=1000){
  if(!inherits(object, "mmes")) stop("object must inherit from class 'mmes'.", call.=FALSE)
  method <- match.arg(method)
  terms <- names(object$uList)
  if(is.numeric(term)) term <- terms[term]
  ti <- match(term, terms)
  if(is.na(ti)) stop("Unknown random term. Available: ", paste(terms, collapse=", "), call.=FALSE)
  if(length(blend) != 1L || !is.finite(blend) || blend < 0 || blend >= 1){
    stop("blend must be one value in [0, 1).", call.=FALSE)
  }
  chunk <- .mme_positive(chunk, "chunk", integer=TRUE)
  U <- as.matrix(object$uList[[ti]])
  levels <- rownames(U)

  # G^{-1} exactly as used in the fit
  if(!is.null(Gu)){
    if(.is_giv(Gu)) Gu <- .giv_to_precision(Gu)
    Ginv <- if(isTRUE(attr(Gu, "inverse"))) Gu else solve(as.matrix(Gu))
  }else if(!is.null(object$rotation) &&
           as.character(object$rotation$term) %in% c(as.character(ti), term)){
    V <- object$rotation$vectors
    Ginv <- V %*% (object$rotation$precision * t(V))
    dimnames(Ginv) <- list(object$rotation$levels, object$rotation$levels)
  }else{
    Ginv <- .mmes_random_precision(object)[[ti]]
    if(is.null(rownames(Ginv))) dimnames(Ginv) <- list(levels, levels)
  }
  if(!all(levels %in% rownames(Ginv))){
    stop("The relationship matrix does not contain all levels of the term.", call.=FALSE)
  }
  Ginv <- Ginv[levels, levels, drop=FALSE]

  if(method == "A.mat"){
    if(missing(M) || is.null(rownames(M))) stop("M must be a marker matrix with rownames.", call.=FALSE)
    M <- as.matrix(M)
    absent <- setdiff(levels, rownames(M))
    if(length(absent)){
      stop("Individuals of the term missing from M: ",
           paste(utils::head(absent, 10), collapse=", "), call.=FALSE)
    }
    extra <- setdiff(rownames(M), levels)
    if(length(extra)) message(length(extra), " rows of M are not levels of the term and are ignored.")
    # allele frequencies are computed on the rows used to build G (all rows of M)
    tr <- .amat_marker_transform(M, min.MAF)
    Zm <- tr$Z[levels, , drop=FALSE]
    s <- (1 - blend) / tr$k
    keep <- tr$keep
    markerNames <- if(is.null(colnames(M))) paste0("m", seq_len(ncol(M))) else colnames(M)
  }else{
    if(is.null(Z) || is.null(scale)) stop("method='custom' needs Z and scale.", call.=FALSE)
    Z <- as.matrix(Z)
    if(is.null(rownames(Z)) || !all(levels %in% rownames(Z))){
      stop("Z must have rownames containing every level of the term.", call.=FALSE)
    }
    Zm <- Z[levels, , drop=FALSE]
    s <- (1 - blend) * scale
    keep <- rep(TRUE, ncol(Z))
    markerNames <- if(is.null(colnames(Z))) paste0("m", seq_len(ncol(Z))) else colnames(Z)
  }

  effects <- s * crossprod(Zm, Matrix::crossprod(Ginv, U))
  coords <- colnames(U)
  if(is.null(coords)) coords <- paste0("c", seq_len(ncol(U)))

  out <- data.frame(marker=rep(markerNames, times=ncol(U)),
                    coordinate=rep(coords, each=length(markerNames)),
                    effect=NA_real_, stringsAsFactors=FALSE)
  rows <- function(c) (c - 1L) * length(markerNames) + which(keep)
  for(c in seq_len(ncol(U))) out$effect[rows(c)] <- effects[, c]

  if(se){
    .mmes_require_coefficient_matrix(object)
    if(is.null(object$C) || nrow(object$C) == 0){
      stop("se=TRUE requires the coefficient matrix C in the fitted object.", call.=FALSE)
    }
    theta <- as.matrix(object$theta[[ti]])
    part <- object$partitions[[ti]]
    nEff <- nrow(object$bu)
    solveVcov <- .mmes_vcov_solver(object)
    out$se <- NA_real_
    for(start in seq.int(1L, max(1L, ncol(Zm)), by=chunk)){
      if(!ncol(Zm)) break
      idx <- seq.int(start, min(ncol(Zm), start + chunk - 1))
      GiZ <- Ginv %*% Zm[, idx, drop=FALSE]
      gPart <- s^2 * colSums(Zm[, idx, drop=FALSE] * GiZ)
      H <- s * GiZ
      for(c in seq_len(ncol(U))){
        cols <- part[c, 1]:part[c, 2]
        D <- Matrix::sparseMatrix(i=rep(seq_along(idx), each=length(cols)),
                                  j=rep(cols, times=length(idx)),
                                  x=as.vector(H),
                                  dims=c(length(idx), nEff))
        pev <- diag(solveVcov(D))
        v <- theta[c, c] * gPart - pev
        out$se[rows(c)[idx]] <- sqrt(pmax(v, 0))
      }
    }
    out$z <- out$effect / out$se
    out$p.value <- 2 * stats::pnorm(-abs(out$z))
  }
  if(ncol(U) == 1L) out$coordinate <- NULL
  attr(out, "scale") <- s
  attr(out, "markersKept") <- markerNames[keep]
  out
}
