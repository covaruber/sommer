vsm <- function(..., Gu=NULL, sigma2=NULL, fixedSigma2=FALSE,
                rotation=FALSE, isFixed=FALSE, verbose=TRUE){

  init <- list(...)
  expr_names <- as.character(substitute(list(...)))[-1L]

  if(length(init) < 1L){
    stop("vsm() requires at least one structured term.", call. = FALSE)
  }

  if(!all(vapply(init, is.list, logical(1)))){
    stop(
      paste0(
        "Every term supplied to vsm() must be wrapped in a covariance ",
        "constructor returning a CovarianceFactor descriptor, such as ",
        "ism(), dsm(), usm(), ar1m(), csm(), rrm(), fam(), ",
        "maternm(), toeplitzm(), sar(), car(), or ownm()."
      ),
      call. = FALSE
    )
  }

  # sigma2=NULL (the default) is a marker meaning "let mmes() replace this
  # with a data-driven starting value"; an explicit user value is always
  # honored as-is and never overridden.
  sigma2IsDefault <- is.null(sigma2)
  if(sigma2IsDefault) sigma2 <- 0.15

  if(!is.finite(sigma2) || length(sigma2) != 1L || sigma2 <= 0){
    stop("sigma2 in vsm() must be one positive finite value.",
         call. = FALSE)
  }

  if(length(fixedSigma2) != 1L){
    stop("fixedSigma2 must have length one.", call. = FALSE)
  }

  if(length(rotation) != 1L || !is.logical(rotation) || is.na(rotation)){
    stop("rotation must be a single TRUE/FALSE value.", call. = FALSE)
  }


  # ======================================================================
  # Main-effect term and covariance-shaping factors
  # ======================================================================

  # Last term is the main-effect incidence. Every preceding term is a
  # covariance-shaping factor. This makes Kronecker depth unlimited.
  main <- init[[length(init)]]

  if(is.null(main$Z)){
    stop("The last vsm() term must provide a Z matrix.",
         call. = FALSE)
  }

  factor_terms <-
    if(length(init) > 1L) init[-length(init)] else list()


  main_vars <-
    all.vars(
      as.formula(
        paste0("~", expr_names[length(expr_names)])
      )
    )

  all_vars <- unique(
    unlist(
      lapply(
        expr_names,
        function(z){
          all.vars(as.formula(paste0("~", z)))
        }
      )
    )
  )

  is.residual <-
    "units" %in% main_vars ||
    "units" %in% all_vars

  n <- nrow(main$Z)


  # ======================================================================
  # Normalize main-effect design
  # ======================================================================

  mainZ <- to_sparse(main$Z)


  # ======================================================================
  # Relationship / known precision matrix
  #
  # IMPORTANT:
  # Matrix coercion and subsetting may drop arbitrary attributes.
  # Therefore we record the precision status BEFORE doing any operation
  # on Gu and explicitly restore it afterwards.
  # ======================================================================

  Gu_is_inverse <- FALSE
  Gu_was_supplied <- !is.null(Gu)

  if(.is_giv(Gu)) Gu <- .giv_to_precision(Gu)

  if(!is.null(Gu)){

    Gu_is_inverse <- isTRUE(attr(Gu, "inverse"))

    if(!Gu_is_inverse){
      stop(
        "Gu must have attr(Gu, 'inverse')=TRUE for the Henderson solver.",
        call. = FALSE
      )
    }

    # Conversion can drop custom attributes.
    Gu <- to_precision_sparse(Gu)

    # Restore immediately so the invariant is maintained throughout vsm().
    attr(Gu, "inverse") <- TRUE

    if(is.null(colnames(Gu)) || is.null(rownames(Gu))){
      stop(
        "Gu must have row and column names matching the main-effect levels.",
        call. = FALSE
      )
    }

    if(nrow(Gu) != ncol(Gu)){
      stop("Gu must be a square precision matrix.",
           call. = FALSE)
    }

    if(!identical(rownames(Gu), colnames(Gu))){
      stop(
        "Gu must have identical row and column level names in the same order.",
        call. = FALSE
      )
    }

    miss <- setdiff(colnames(mainZ), colnames(Gu))

    if(length(miss)){
      stop(
        paste(
          "Levels missing from Gu:",
          paste(miss, collapse=", ")
        ),
        call. = FALSE
      )
    }

    # Preserve historical ability to predict levels contained in Gu that
    # are not observed in the data by adding zero incidence columns.
    extra <- setdiff(colnames(Gu), colnames(mainZ))

    if(length(extra)){

      if(verbose){
        cat(
          "Adding", length(extra), "additional Gu levels to the main-effect model matrix:",
          paste(utils::head(extra, 10L), collapse=", "),
          if(length(extra) > 10L) "..." else "",
          "\n"
        )
      }

      add <- Matrix::Matrix(
        0,
        nrow=nrow(mainZ),
        ncol=length(extra),
        sparse=TRUE
      )

      colnames(add) <- extra

      mainZ <- cbind(mainZ, add)
      mainZ <- to_sparse(mainZ)
    }
  }


  # ======================================================================
  # Row-wise Kronecker / Khatri-Rao product
  #
  # Ordering agrees with kronecker(K1,K2,...):
  # earlier factors are slow indices and later factors are fast indices.
  # ======================================================================

  row_kron <- function(A, B){

    A <- to_sparse(A)
    B <- to_sparse(B)

    if(nrow(A) != nrow(B)){
      stop(
        paste0(
          "All covariance factors inside vsm() must have ",
          "the same number of rows."
        ),
        call. = FALSE
      )
    }

    sa <- Matrix::summary(A)
    sb <- Matrix::summary(B)
    m <- merge(
      data.frame(i=sa$i, ja=sa$j, xa=sa$x),
      data.frame(i=sb$i, jb=sb$j, xb=sb$x),
      by="i", sort=FALSE
    )
    nb <- ncol(B)
    ans <- Matrix::sparseMatrix(
      i=m$i, j=(m$ja - 1L) * nb + m$jb, x=m$xa * m$xb,
      dims=c(nrow(A), ncol(A) * nb)
    )
    colnames(ans) <- as.vector(t(outer(colnames(A), colnames(B), paste, sep=":")))

    to_sparse(ans)
  }


  # ======================================================================
  # Compile covariance-shaping factors
  # ======================================================================

  if(length(factor_terms) == 0L){

    Z0 <- Matrix::Matrix(
      1,
      nrow=n,
      ncol=1,
      sparse=TRUE
    )

    colnames(Z0) <- "1"

    factors <- list()

  }else{

    if(any(
      vapply(
        factor_terms,
        function(z) is.null(z$covFactor),
        logical(1)
      )
    )){
      stop(
        paste0(
          "Every covariance-shaping term in vsm() must ",
          "supply a covFactor descriptor."
        ),
        call. = FALSE
      )
    }

    Z0 <- to_sparse(factor_terms[[1]]$Z)

    if(length(factor_terms) > 1L){

      for(i in 2:length(factor_terms)){
        Z0 <- row_kron(
          Z0,
          factor_terms[[i]]$Z
        )
      }
    }

    factors <- lapply(
      factor_terms,
      function(z){

        f <- .compile_covfactor(z$covFactor)

        .validate_covfactor(f)

        f
      }
    )
  }


  # ======================================================================
  # Random-effect design
  #
  # Build one Z block for each covariance-product coordinate.
  #
  # Residual structures only need the product-coordinate layout and do not
  # materialize an observation-identity random-effect design.
  # ======================================================================

  if(is.residual){

    Z <- list()

  }else{

    Z <- vector(
      "list",
      ncol(Z0)
    )

    for(j in seq_len(ncol(Z0))){

      # Row scaling keeps O(nnz); an outer-product mask would be dense n x levels.
      Z[[j]] <-
        to_sparse(
          Matrix::drop0(
            Matrix::Diagonal(x=as.numeric(Z0[,j])) %*% mainZ
          )
        )

      colnames(Z[[j]]) <-
        colnames(mainZ)
    }
  }


  # ======================================================================
  # Final Gu normalization
  #
  # This section establishes an explicit invariant:
  #
  # Every vsm() object crossing the R/C++ or vsm()/covm() boundary has:
  #
  #   1. Gu represented as dgCMatrix
  #   2. attr(Gu,"inverse") == TRUE
  #
  # Matrix operations above are never trusted to preserve the attribute.
  # ======================================================================

  if(is.residual){

    Gu <- to_sparse(
      Matrix::Matrix(
        1,
        1,
        1,
        sparse=TRUE
      )
    )

    attr(Gu, "inverse") <- TRUE

  }else if(is.null(Gu)){

    Gu <- to_sparse(
      Matrix::Diagonal(
        n=ncol(mainZ),
        x=1
      )
    )

    colnames(Gu) <-
      rownames(Gu) <-
      colnames(mainZ)

    attr(Gu, "inverse") <- TRUE

  }else{

    # Subsetting is required to put Gu into exactly the same ordering as
    # mainZ. Subsetting/coercion may drop custom attributes, so restore the
    # precision marker AFTER the operation.
    if(!identical(colnames(Gu), colnames(mainZ))){
      guOrder <- match(colnames(mainZ), colnames(Gu))
      Gu <- Gu[guOrder, guOrder, drop=FALSE]
    }

    Gu <- to_precision_sparse(Gu)

    attr(Gu, "inverse") <- TRUE
  }


  # Final representation validation.
  if(!inherits(Gu, c("dgCMatrix", "dsCMatrix"))){
    stop(
      paste0(
        "Internal vsm() error: Gu was not normalized ",
        "to a sparse precision matrix."
      ),
      call. = FALSE
    )
  }

  if(!isTRUE(attr(Gu, "inverse"))){
    stop(
      paste0(
        "Internal vsm() error: Gu lost its ",
        "inverse/precision attribute."
      ),
      call. = FALSE
    )
  }

  rotationInfo <- NULL
  GuRot <- NULL

  if(isTRUE(rotation)){
    if(is.residual){
      stop("rotation=TRUE is only available for random-effect vsm() terms.",
           call. = FALSE)
    }
    if(!Gu_was_supplied){
      stop("rotation=TRUE requires a supplied Gu precision matrix.",
           call. = FALSE)
    }
    if(!isSymmetric(Gu, tol=1e-10)){
      stop("rotation=TRUE requires a symmetric Gu precision matrix.",
           call. = FALSE)
    }

    eig <- eigen(as.matrix(Gu), symmetric=TRUE)
    if(any(!is.finite(eig$values)) || any(eig$values <= 0)){
      stop("rotation=TRUE requires a positive-definite Gu precision matrix.",
           call. = FALSE)
    }

    U <- eig$vectors
    for(j in seq_len(ncol(U))){
      pivot <- which.max(abs(U[,j]))
      if(U[pivot,j] < 0) U[,j] <- -U[,j]
    }
    modeNames <- paste0(".eigen", seq_along(eig$values))
    rownames(U) <- colnames(Gu)
    colnames(U) <- modeNames

    GuRot <- to_sparse(Matrix::Diagonal(x=eig$values))
    rownames(GuRot) <- colnames(GuRot) <- modeNames
    attr(GuRot, "inverse") <- TRUE

    rotationInfo <- list(
      vectors=U,
      precision=eig$values,
      covariance=1/eig$values,
      levels=colnames(Gu),
      modes=modeNames
    )
  }


  # ======================================================================
  # Flatten covariance-product descriptor
  #
  # All optimizer coordinates are unconstrained working parameters:
  #
  #   log(sigma2)
  #   atanh(rho)
  #   log variance ratios
  #   normalized-Cholesky coordinates
  #   etc.
  # ======================================================================

  par <- c(
    log_sigma2=log(sigma2)
  )

  free <- c(
    !isTRUE(fixedSigma2)
  )

  par_names <- "sigma2"


  if(length(factors)){

    for(i in seq_along(factors)){

      f <- factors[[i]]

      f$par_start <-
        length(par) + 1L

      if(length(f$par)){

        par <- c(
          par,
          f$par
        )

        free <- c(
          free,
          f$free
        )

        prefix <-
          if(length(expr_names) >= i){
            expr_names[i]
          }else{
            paste0("factor", i)
          }

        par_names <- c(
          par_names,
          paste(
            prefix,
            f$par_names,
            sep=":"
          )
        )
      }

      f$par_end <- length(par)
      f$label <- if(length(expr_names) >= i) expr_names[i] else paste0("factor", i)

      factors[[i]] <- f
    }
  }


  product_dim <-
    if(length(factors)){
      prod(
        vapply(
          factors,
          `[[`,
          numeric(1),
          "dim"
        )
      )
    }else{
      1L
    }


  if(product_dim != ncol(Z0)){
    stop(
      paste0(
        "Internal vsm() error: covariance product dimension ",
        "does not match the combined design."
      ),
      call. = FALSE
    )
  }


  covStruct <- list(

    type="kron",

    par=as.numeric(par),

    free=as.logical(free),

    par_names=par_names,

    factors=factors,

    dim=as.integer(product_dim),

    levels=colnames(Z0),

    scale_index=1L,

    descriptor_version=2L,

    factor_interface="CovarianceFactor",

    parameterization="working",

    main_levels=colnames(mainZ),

    # Lets mmes() know it may replace par[1] (log_sigma2) with a
    # data-driven starting value; never set when the user supplied sigma2.
    sigma2_is_default=sigma2IsDefault
  )

  names(covStruct$par) <-
    par_names


  # ======================================================================
  # Residual product-coordinate mapping
  #
  # Residual covariance factors must define exactly one covariance-product
  # coordinate for every observation.
  # ======================================================================

  residualLocalIndex <- NULL

  if(is.residual){

    ss <- Matrix::summary(Z0)

    if(anyDuplicated(ss$i) || any(abs(ss$x - 1) > 1e-12)){
      stop(
        paste0(
          "Residual covariance factors must define exactly ",
          "one covariance-product level per observation."
        ),
        call. = FALSE
      )
    }

    residualLocalIndex <- rep(NA_integer_, nrow(Z0))
    residualLocalIndex[ss$i] <- as.integer(ss$j)
  }


  # ======================================================================
  # Return object
  # ======================================================================

  output <- list(

    Z=Z,

    Gu=Gu,

    GuRot=GuRot,

    rotation=rotationInfo,

    covStruct=covStruct,

    residualLocalIndex=residualLocalIndex,

    productDesign=Z0,

    partitionsR=NULL
  )


  if(isFixed){

    return(
      as.matrix(
        do.call(
          cbind,
          Z
        )
      )
    )
  }


  output
}

to_sparse <- function(z){
  if(!inherits(z, "dgCMatrix")){
    z <- as(as(z, "generalMatrix"), "CsparseMatrix") # methods::as(z, "dgCMatrix")
  }
  z
}

# Symmetric precision matrices keep one triangle (dsCMatrix); C++ expands them.
to_precision_sparse <- function(z){
  if(inherits(z, "dsCMatrix")) return(z)
  if(!inherits(z, "Matrix")) z <- Matrix::Matrix(z, sparse=TRUE)
  if(Matrix::isSymmetric(z)){
    z <- as(Matrix::forceSymmetric(as(z, "CsparseMatrix")), "CsparseMatrix")
    if(inherits(z, "dsCMatrix")) return(z)
  }
  to_sparse(z)
}

# 3-column (row, column, value) precision input, e.g. ASReml ainverse()/nadiv
# listAinv. A 3x3 full matrix with dimnames and no rowNames attr is not giv.
.is_giv <- function(x){
  if(is.null(x) || inherits(x, "Matrix")) return(FALSE)
  if(!(is.matrix(x) || is.data.frame(x)) || ncol(x) != 3L) return(FALSE)
  if(!is.null(attr(x, "rowNames")) || is.data.frame(x)) return(TRUE)
  !(nrow(x) == 3L && !is.null(rownames(x)) && !is.null(colnames(x)))
}

.giv_to_precision <- function(x){
  levs <- attr(x, "rowNames")
  if(is.null(levs)){
    stop("Gu supplied in 3-column (row, column, value) format must carry a ",
         "'rowNames' attribute with the level names.", call. = FALSE)
  }
  levs <- as.character(levs)
  if(anyNA(levs) || anyDuplicated(levs)){
    stop("attr(Gu, 'rowNames') must contain unique, non-missing level names.",
         call. = FALSE)
  }
  if(isFALSE(attr(x, "inverse")) || isFALSE(attr(x, "INVERSE"))){
    stop("Gu in 3-column format must be a precision (inverse) matrix.",
         call. = FALSE)
  }
  nl <- length(levs)
  ri <- x[, 1L]; ci <- x[, 2L]; v <- x[, 3L]
  if(!is.numeric(ri) || !is.numeric(ci) || !is.numeric(v)){
    stop("All three columns of a 3-column Gu must be numeric.", call. = FALSE)
  }
  if(anyNA(ri) || anyNA(ci) || any(ri != round(ri)) || any(ci != round(ci)) ||
     any(ri < 1) || any(ci < 1) || any(ri > nl) || any(ci > nl)){
    stop("Row/Column indices of a 3-column Gu must be integers between 1 and ",
         "length(attr(Gu, 'rowNames')) = ", nl, ".", call. = FALSE)
  }
  if(any(!is.finite(v))){
    stop("Values of a 3-column Gu must be finite.", call. = FALSE)
  }
  # Canonicalize to the upper triangle; lower-only, upper-only or both are fine.
  i <- as.integer(pmin(ri, ci)); j <- as.integer(pmax(ri, ci))
  key <- (as.numeric(j) - 1) * nl + i
  dup <- duplicated(key)
  if(any(dup)){
    first <- v[match(key[dup], key)]
    if(any(abs(v[dup] - first) > 1e-8 * pmax(1, abs(first)))){
      stop("3-column Gu contains conflicting values for the same cell ",
           "(e.g. (i,j) and (j,i) differ); the matrix must be symmetric.",
           call. = FALSE)
    }
    i <- i[!dup]; j <- j[!dup]; v <- v[!dup]
  }
  missDiag <- setdiff(seq_len(nl), i[i == j])
  if(length(missDiag)){
    stop("3-column Gu is missing diagonal values for ", length(missDiag),
         " level(s), e.g.: ", paste(head(levs[missDiag], 5), collapse=", "),
         call. = FALSE)
  }
  out <- Matrix::sparseMatrix(i=i, j=j, x=as.numeric(v), dims=c(nl, nl),
                              dimnames=list(levs, levs), symmetric=TRUE)
  attr(out, "inverse") <- TRUE
  out
}

# Combine two vsm() random-effect structures into one correlated random structure.
#
# The returned object uses exactly the same CovarianceFactor-v2 contract as vsm().
# ai_mme_sp2() therefore does not require any covm-specific code.
#
# Current implementation intentionally requires SIMPLE vsm() structures
# (one covariance-product coordinate on each side).  Thus the two effects may
# have different incidence matrices, but they must act on the same coefficient
# space and use the same Gu precision matrix.
#
# Covariance model:
#
#   Var([u1',u2']') = sigma2 * K_effect %x% A
#
# where K_effect is a normalized 2 x 2 unstructured covariance shape and
# sigma2 is the variance of the first effect.  The normalized-Cholesky
# parameterization guarantees positive definiteness.
#
covm <- function(ran1, ran2, thetaC=NULL, theta=NULL,
                 fixed=NULL, fixedSigma2=FALSE,
                 labels=c("ran1","ran2"), tol=1e-10){

  # ======================================================================
  # Validate vsm() objects
  # ======================================================================

  check_vsm <- function(x, nm){

    if(
      !is.list(x) ||
      is.null(x$Z) ||
      is.null(x$Gu) ||
      is.null(x$covStruct)
    ){
      stop(
        nm,
        " must be the result of vsm().",
        call.=FALSE
      )
    }

    cs <- x$covStruct

    if(
      !identical(cs$type, "kron") ||
      is.null(cs$descriptor_version) ||
      cs$descriptor_version < 2L
    ){
      stop(
        nm,
        " must use the CovarianceFactor-v2 vsm() interface.",
        call.=FALSE
      )
    }

    # Current implementation intentionally handles simple random effects.
    if(
      length(x$Z) != 1L ||
      cs$dim != 1L ||
      length(cs$factors) != 0L
    ){
      stop(
        paste0(
          "The current covm() implementation combines simple vsm() ",
          "random effects only. Each side must have one covariance-product ",
          "coordinate, e.g. vsm(ism(effect), Gu=Ai). Shared structured ",
          "factors can be added later without changing ai_mme_sp2()."
        ),
        call.=FALSE
      )
    }

    if(!isTRUE(attr(x$Gu, "inverse"))){
      stop(
        nm,
        "$Gu must be an inverse/precision matrix with ",
        "attr(Gu,'inverse')=TRUE.",
        call.=FALSE
      )
    }

    invisible(TRUE)
  }


  check_vsm(ran1, "ran1")
  check_vsm(ran2, "ran2")


  # ======================================================================
  # Effect labels
  # ======================================================================

  if(
    length(labels) != 2L ||
    anyNA(labels) ||
    any(!nzchar(labels)) ||
    anyDuplicated(labels)
  ){
    stop(
      "labels must contain two different non-empty names.",
      call.=FALSE
    )
  }

  labels <- as.character(labels)


  # ======================================================================
  # Capture precision status BEFORE Matrix operations
  # ======================================================================

  inv1 <- isTRUE(
    attr(ran1$Gu, "inverse")
  )

  inv2 <- isTRUE(
    attr(ran2$Gu, "inverse")
  )

  if(!inv1 || !inv2){
    stop(
      paste0(
        "Both Gu matrices must be inverse/precision matrices ",
        "for the Henderson solver."
      ),
      call.=FALSE
    )
  }


  # ======================================================================
  # Normalize matrices
  #
  # to_sparse() is not assumed to preserve arbitrary R attributes.
  # ======================================================================

  Z1 <- to_sparse(
    ran1$Z[[1L]]
  )

  Z2 <- to_sparse(
    ran2$Z[[1L]]
  )

  G1 <- to_precision_sparse(
    ran1$Gu
  )

  G2 <- to_precision_sparse(
    ran2$Gu
  )

  # Explicitly restore precision metadata after coercion.
  attr(G1, "inverse") <- TRUE
  attr(G2, "inverse") <- TRUE


  # ======================================================================
  # Incidence compatibility
  # ======================================================================

  if(nrow(Z1) != nrow(Z2)){
    stop(
      paste0(
        "The two random-effect incidence matrices must ",
        "have the same number of observations."
      ),
      call.=FALSE
    )
  }

  if(ncol(Z1) != ncol(Z2)){
    stop(
      paste0(
        "The two random effects must have the same ",
        "coefficient dimension."
      ),
      call.=FALSE
    )
  }


  # ======================================================================
  # Main-effect level compatibility
  # ======================================================================

  lev1 <- ran1$covStruct$main_levels
  lev2 <- ran2$covStruct$main_levels

  if(is.null(lev1)){
    lev1 <- colnames(Z1)
  }

  if(is.null(lev2)){
    lev2 <- colnames(Z2)
  }

  lev1 <- as.character(lev1)
  lev2 <- as.character(lev2)


  if(!identical(lev1, lev2)){
    stop(
      paste0(
        "The two random effects must have identical ",
        "main-effect levels in the same order."
      ),
      call.=FALSE
    )
  }


  # ======================================================================
  # Gu compatibility
  # ======================================================================

  if(
    !all(dim(G1) == dim(G2)) ||
    !identical(rownames(G1), rownames(G2)) ||
    !identical(colnames(G1), colnames(G2))
  ){
    stop(
      "ran1 and ran2 must use the same Gu coefficient space.",
      call.=FALSE
    )
  }


  DG <- Matrix::drop0(
    G1 - G2,
    tol=tol
  )

  if(length(DG@x)){
    stop(
      "ran1 and ran2 must use the same Gu precision matrix.",
      call.=FALSE
    )
  }


  # ======================================================================
  # Initial natural-scale 2 x 2 covariance
  #
  # Preserve the historical covm() starting covariance.
  # ======================================================================

  if(is.null(theta)){

    theta <-
      diag(2) * 0.15 +
      matrix(
        0.015,
        2,
        2
      )

  }else{

    theta <- as.matrix(theta)
  }


  if(
    !all(dim(theta) == c(2L,2L)) ||
    any(!is.finite(theta))
  ){
    stop(
      "theta must be a finite 2 x 2 covariance matrix.",
      call.=FALSE
    )
  }


  theta <-
    (theta + t(theta)) / 2


  ev <- eigen(
    theta,
    symmetric=TRUE,
    only.values=TRUE
  )$values


  if(min(ev) <= tol){
    stop(
      "theta supplied to covm() must be positive definite.",
      call.=FALSE
    )
  }


  # ======================================================================
  # Convert natural covariance to the standard vsm() representation
  #
  #     Sigma = sigma2 * K_effect %x% A
  #
  # K_effect is normalized so K[1,1] = 1.
  # ======================================================================

  sigma2 <- theta[1L,1L]

  if(
    !is.finite(sigma2) ||
    sigma2 <= 0
  ){
    stop(
      "theta[1,1] must be positive.",
      call.=FALSE
    )
  }


  K <- theta / sigma2


  ch <- try(
    chol(K),
    silent=TRUE
  )

  if(inherits(ch, "try-error")){
    stop(
      "theta supplied to covm() must be positive definite.",
      call.=FALSE
    )
  }


  L <- t(ch)

  # Numerical normalization guarantees L11 = 1.
  L <- L / L[1L,1L]


  # ======================================================================
  # Normalized-Cholesky coordinates
  #
  # For q=2:
  #
  #     L = [ 1       0       ]
  #         [ L21   exp(eta22) ]
  #
  # Thus there are two factor coordinates:
  #
  #     L21
  #     log(L22)
  # ======================================================================

  us_par <- c(
    L[2L,1L],
    log(L[2L,2L])
  )


  us_names <- c(

    paste0(
      "chol[",
      labels[2L],
      ",",
      labels[1L],
      "]"
    ),

    paste0(
      "chol_diag[",
      labels[2L],
      "]"
    )
  )


  # ======================================================================
  # Fixed parameters
  # ======================================================================

  if(is.null(fixed)){
    fixed <- c(
      FALSE,
      FALSE
    )
  }


  if(
    length(fixed) != 2L ||
    anyNA(fixed)
  ){
    stop(
      paste0(
        "fixed must be a logical vector of length 2 for ",
        "the normalized-Cholesky coordinates."
      ),
      call.=FALSE
    )
  }


  fixed <- as.logical(fixed)


  # ======================================================================
  # Legacy thetaC
  #
  # Cell-wise thetaC constraints cannot generally be mapped exactly onto
  # normalized-Cholesky coordinates.
  # ======================================================================

  if(!is.null(thetaC)){

    stop(
      paste0(
        "thetaC is not supported by the CovarianceFactor-v2 covm() ",
        "parameterization. Use fixedSigma2 to fix the first-effect ",
        "variance and fixed to fix the two normalized-Cholesky ",
        "coordinates. Cell-wise thetaC codes cannot in general be ",
        "translated exactly to Cholesky-coordinate constraints."
      ),
      call.=FALSE
    )
  }


  # ======================================================================
  # Construct generic effect covariance factor
  #
  # This is exactly the same generic 'us' factor consumed by the covariance
  # engine. ai_mme_sp2() therefore has no knowledge of covm().
  # ======================================================================

  effectFactor <- .compile_covfactor(
    list(

      type="us",

      dim=2L,

      levels=labels,

      par=us_par,

      free=!fixed,

      par_names=us_names,

      us_row=c(
        2L,
        2L
      ),

      us_col=c(
        1L,
        2L
      ),

      us_diag=c(
        FALSE,
        TRUE
      )
    )
  )


  .validate_covfactor(
    effectFactor
  )


  # ======================================================================
  # Flatten into standard vsm() CovarianceFactor-v2 descriptor
  # ======================================================================

  par <- c(
    log_sigma2=log(sigma2),
    effectFactor$par
  )


  free <- c(
    !isTRUE(fixedSigma2),
    effectFactor$free
  )


  par_names <- c(
    "sigma2",
    effectFactor$par_names
  )


  effectFactor$par_start <- 2L
  effectFactor$par_end <- length(par)


  covStruct <- list(

    type="kron",

    par=as.numeric(par),

    free=as.logical(free),

    par_names=par_names,

    factors=list(
      effectFactor
    ),

    dim=2L,

    levels=labels,

    scale_index=1L,

    descriptor_version=2L,

    factor_interface="CovarianceFactor",

    parameterization="working",

    main_levels=lev1
  )


  names(covStruct$par) <-
    par_names


  # ======================================================================
  # Final Gu
  #
  # K_effect %x% A uses one shared precision matrix A.
  # Explicitly preserve precision metadata in the returned object.
  # ======================================================================

  Gu <- G1

  attr(Gu, "inverse") <- TRUE


  if(!inherits(Gu, c("dgCMatrix", "dsCMatrix"))){
    stop(
      paste0(
        "Internal covm() error: Gu was not normalized ",
        "to a sparse precision matrix."
      ),
      call.=FALSE
    )
  }


  if(!isTRUE(attr(Gu, "inverse"))){
    stop(
      paste0(
        "Internal covm() error: Gu lost its ",
        "inverse/precision attribute."
      ),
      call.=FALSE
    )
  }


  # ======================================================================
  # Return standard random-structure object
  #
  # Ordering:
  #
  #     effect coordinate first
  #     coefficient level second
  #
  # giving
  #
  #     K_effect %x% A
  #
  # Z1 and Z2 are therefore the two covariance-product coordinates.
  # ======================================================================

  list(

    Z=list(
      Z1,
      Z2
    ),

    Gu=Gu,

    covStruct=covStruct,

    residualLocalIndex=NULL,

    productDesign=NULL,

    partitionsR=NULL,

    covm=TRUE,

    covm_labels=labels
  )
}

# -------------------------------------------------------------------------
# Covariance structure spanning several random terms (ASReml str()):
#   Var(u_1,...,u_T) = sigma2 * K_terms %x% K_inner %x% A
# Terms must share the coefficient space, Gu and inner covariance factors.
# -------------------------------------------------------------------------
strm <- function(..., cov=usm, Gu=NULL, labels=NULL, sigma2=NULL,
                 fixedSigma2=FALSE){
  terms <- list(...)
  exprs <- vapply(as.list(substitute(list(...)))[-1L],
                  function(e) paste(deparse(e), collapse=""), character(1))
  nT <- length(terms)
  if(nT < 2L) stop("strm() needs at least two vsm() terms.", call.=FALSE)
  for(t in seq_len(nT)){
    x <- terms[[t]]
    if(!is.list(x) || is.null(x$covStruct) || !identical(x$covStruct$type, "kron") ||
       is.null(x$covStruct$descriptor_version) || x$covStruct$descriptor_version < 2L){
      stop("Every argument of strm() must be a vsm() term.", call.=FALSE)
    }
    if(!length(x$Z)) stop("strm() terms must be random effects, not residual terms.", call.=FALSE)
    if(!is.null(x$rotation)) stop("rotation=TRUE is not available inside strm().", call.=FALSE)
    if(isFALSE(x$covStruct$sigma2_is_default) || !isTRUE(x$covStruct$free[1L])){
      warning("sigma2/fixedSigma2 of individual strm() terms are ignored; use those of strm().",
              call.=FALSE)
    }
  }
  if(is.null(labels)){
    labels <- if(!is.null(names(terms)) && all(nzchar(names(terms)))) names(terms) else
      paste0("t", seq_len(nT))
  }
  labels <- as.character(labels)
  if(length(labels) != nT || anyDuplicated(labels) || any(!nzchar(labels))){
    stop("labels must contain one distinct non-empty name per term.", call.=FALSE)
  }

  # inner covariance factors must coincide across terms
  inner <- terms[[1L]]$covStruct$factors
  signature <- function(fs) vapply(fs, function(f) paste(f$model, f$dim,
                                   paste(f$levels, collapse="\r")), character(1))
  for(t in seq_len(nT)[-1L]){
    if(!identical(signature(terms[[t]]$covStruct$factors), signature(inner))){
      stop("All strm() terms must use the same covariance factors (same constructors and levels).",
           call.=FALSE)
    }
  }
  qInner <- length(terms[[1L]]$Z)

  isIdentity <- function(G){
    Matrix::isDiagonal(G) && all(abs(Matrix::diag(G) - 1) < 1e-12)
  }
  if(!is.null(Gu)){
    if(.is_giv(Gu)) Gu <- .giv_to_precision(Gu)
    if(!isTRUE(attr(Gu, "inverse"))){
      stop("Gu must have attr(Gu, 'inverse')=TRUE for the Henderson solver.", call.=FALSE)
    }
    Gu <- to_precision_sparse(Gu)
    if(is.null(rownames(Gu)) || !identical(rownames(Gu), colnames(Gu))){
      stop("Gu must have identical row and column names.", call.=FALSE)
    }
    common <- rownames(Gu)
  }else{
    termGu <- lapply(terms, `[[`, "Gu")
    known <- which(!vapply(termGu, isIdentity, logical(1)))
    if(length(known)){
      Gu <- termGu[[known[1L]]]
      common <- rownames(Gu)
      for(k in known[-1L]){
        Gk <- termGu[[k]]
        if(!identical(rownames(Gk), common) ||
           length(Matrix::drop0(Gk - Gu, tol=1e-10)@x)){
          stop("strm() terms with a Gu must share the same Gu; supply it once through strm(Gu=).",
               call.=FALSE)
        }
      }
    }else{
      common <- unique(unlist(lapply(terms, function(x) x$covStruct$main_levels)))
      Gu <- to_sparse(Matrix::Diagonal(n=length(common), x=1))
      rownames(Gu) <- colnames(Gu) <- common
    }
  }
  missing <- setdiff(unlist(lapply(terms, function(x) x$covStruct$main_levels)), common)
  if(length(missing)){
    stop("Levels missing from Gu: ", paste(utils::head(missing, 10), collapse=", "), call.=FALSE)
  }
  Gu <- to_precision_sparse(Gu[common, common, drop=FALSE])
  attr(Gu, "inverse") <- TRUE

  align <- function(Zt){
    s <- Matrix::summary(to_sparse(Zt))
    Zc <- Matrix::sparseMatrix(i=s$i, j=match(colnames(Zt)[s$j], common), x=s$x,
                               dims=c(nrow(Zt), length(common)))
    colnames(Zc) <- common
    to_sparse(Zc)
  }
  Z <- unlist(lapply(terms, function(x) lapply(x$Z, align)), recursive=FALSE)
  n <- vapply(Z, nrow, integer(1))
  if(length(unique(n)) != 1L){
    stop("All strm() terms must have the same number of observations.", call.=FALSE)
  }

  if(!is.function(cov)) stop("cov must be a covariance constructor such as usm or function(x) fam(x, 1).",
                             call.=FALSE)
  termFactor <- cov(factor(labels, levels=labels))$covFactor
  if(is.null(termFactor)) stop("cov must return a CovarianceFactor descriptor.", call.=FALSE)
  if(termFactor$dim != nT) stop("The term covariance factor must have one level per term.", call.=FALSE)
  factors <- c(list(termFactor), inner)

  sigma2IsDefault <- is.null(sigma2)
  if(sigma2IsDefault) sigma2 <- 0.15
  par <- c(log_sigma2=log(sigma2))
  free <- !isTRUE(fixedSigma2)
  parNames <- "sigma2"
  factorLabels <- c(paste0("strm(", paste(labels, collapse=","), ")"),
                    vapply(inner, function(f) if(!is.null(f$label)) f$label else f$model,
                           character(1)))
  for(i in seq_along(factors)){
    f <- .compile_covfactor(factors[[i]])
    .validate_covfactor(f)
    f$par_start <- length(par) + 1L
    if(length(f$par)){
      par <- c(par, f$par)
      free <- c(free, f$free)
      parNames <- c(parNames, paste(factorLabels[i], f$par_names, sep=":"))
    }
    f$par_end <- length(par)
    f$label <- factorLabels[i]
    factors[[i]] <- f
  }
  innerLevels <- terms[[1L]]$covStruct$levels
  coordLevels <- if(qInner == 1L) labels else
    as.vector(t(outer(labels, innerLevels, paste, sep=":")))

  covStruct <- list(
    type="kron", par=stats::setNames(as.numeric(par), parNames), free=as.logical(free),
    par_names=parNames, factors=factors, dim=as.integer(nT*qInner),
    levels=coordLevels, scale_index=1L, descriptor_version=2L,
    factor_interface="CovarianceFactor", parameterization="working",
    main_levels=common, sigma2_is_default=sigma2IsDefault
  )

  list(Z=Z, Gu=Gu, GuRot=NULL, rotation=NULL, covStruct=covStruct,
       residualLocalIndex=NULL, productDesign=NULL, partitionsR=NULL,
       strm_labels=labels, strm_terms=exprs)
}

###############
## VS structures for mmec

H <- function(timevar=NULL, idvar=NULL, response=NULL, Gu=NULL){
  if(is.null(timevar) ){stop("Please provide the timevar argument.", call. = FALSE)}
  if(is.null(idvar) ){stop("Please provide the idvar argument.", call. = FALSE)}
  if(is.null(response) ){stop("Please provide the response argument.", call. = FALSE)}

  dtx <- data.frame(timevar=timevar, idvar=idvar, v.names=response)
  dtx2 <- aggregate(v.names~timevar+idvar, data=dtx, FUN=mean, na.rm=TRUE)
  wide <- reshape(dtx2, direction = "wide", idvar = "idvar",
                  timevar = "timevar", v.names = "v.names", sep= "_")
  rowNamesWide <-  wide[,1]
  rownames(wide) <- rowNamesWide
  wide <- wide[,-1]
  # if user doesn't provide the a Gu we impute simply and use the correlation matrix as a Gu
  if(is.null(Gu)){
    X <- apply(wide, 2, imputev)
    Gu <- cor(t(X))
  }else{
    Gu = cov2cor(Gu)
  }
  # impute missing data using a relationship matrix
  if(is.null(rownames(Gu))){stop("Gu needs to have row names.", call. = FALSE)}
  if(is.null(colnames(Gu))){stop("Gu needs to have column names.", call. = FALSE)}
  for(iEnv in 1:ncol(wide)){ # iEnv=1
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    imputationVector <- as.numeric(Gu[as.character(rowNamesWide),as.character(rowNamesWide[withData])] %*% as.matrix(wide[withData,iEnv]))
    wide[,iEnv] <- imputationVector  # wide[withoutData,iEnv] <- imputationVector[withoutData]
    # scaleFactor=imputationVector[withData[1]] / wide[withData[1],iEnv]
  }
  colnames(wide) <- gsub("v.names_","", colnames(wide))
  return(wide)
}

# Direct sum of a residual vsm() term over the levels of `by`:
#   R = blockdiag_t sigma2_t * kron_f K_f(theta_{f,t})
# Levels outside `levels` receive an iid residual sigma2_t * I.
dsumm <- function(x, by, levels=NULL){
  byExpr <- deparse(substitute(by))

  if(!is.list(x) || is.null(x$covStruct) || !identical(x$covStruct$type, "kron") ||
     is.null(x$residualLocalIndex)){
    stop("dsumm() expects a residual vsm() term, e.g. ",
         "dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial).", call. = FALSE)
  }
  inner <- x$covStruct
  if(any(vapply(inner$factors, function(f) !is.null(f$section_levels), logical(1)))){
    stop("dsumm() terms cannot be nested.", call. = FALSE)
  }

  n <- length(x$residualLocalIndex)
  if(length(by) != n){
    stop("by in dsumm() must have one value per observation.", call. = FALSE)
  }
  byFactor <- if(is.factor(by)) droplevels(by) else factor(by)
  sections <- levels(byFactor)
  nSections <- length(sections)
  if(nSections < 1L){
    stop("by in dsumm() has no non-missing levels.", call. = FALSE)
  }

  listed <- if(is.null(levels)) sections else unique(as.character(levels))
  unknown <- setdiff(listed, sections)
  if(length(unknown)){
    stop("Unknown levels in dsumm(): ", paste(unknown, collapse=", "), call. = FALSE)
  }
  if(!length(listed)){
    stop("levels in dsumm() must select at least one level of by.", call. = FALSE)
  }
  nListed <- length(listed)
  sectionIndex <- match(sections, listed) - 1L
  sectionIndex[is.na(sectionIndex)] <- -1L

  qInner <- as.integer(inner$dim)
  byIndex <- as.integer(byFactor)
  local <- rep(NA_integer_, n)
  isListed <- !is.na(byIndex) & sectionIndex[pmax(byIndex, 1L)] >= 0L
  local[isListed] <- (byIndex[isListed] - 1L) * qInner + x$residualLocalIndex[isListed]
  for(t in which(sectionIndex < 0L)){
    rows <- which(byIndex == t)
    if(length(rows) > qInner){
      stop("Level '", sections[t], "' of dsumm(by=) has more observations than ",
           "coordinates in the inner vsm() term.", call. = FALSE)
    }
    # iid sections need only distinct coordinates, not their spatial labels.
    local[rows] <- (t - 1L) * qInner + seq_along(rows)
  }

  owner <- .compile_covfactor(list(
    type="diag",
    dim=nSections,
    levels=sections,
    par=rep(0, nSections - 1L),
    free=rep(TRUE, nSections - 1L),
    par_names=if(nSections > 1L) paste0("variance_ratio[", sections[-1L], "]") else character()
  ))
  owner$label <- byExpr
  owner$section_owner <- TRUE

  par <- inner$par[1L]
  free <- inner$free[1L]
  parNames <- "sigma2"
  append_factor <- function(f){
    f$par_start <- length(par) + 1L
    par <<- c(par, f$par)
    free <<- c(free, f$free)
    if(length(f$par)) parNames <<- c(parNames, paste(f$label, f$par_names, sep=":"))
    f$par_end <- length(par)
    f
  }

  factors <- list(append_factor(owner))
  for(f in inner$factors){
    if(identical(f$model, "identity")){
      factors[[length(factors) + 1L]] <- append_factor(f)
      next
    }
    localIdx <- if(f$par_end >= f$par_start) seq.int(f$par_start, f$par_end) else integer()
    nPar <- length(localIdx)
    reportOne <- f$native_report$fun
    template <- f

    g <- f
    g$par <- rep(inner$par[localIdx], nListed)
    g$free <- rep(inner$free[localIdx], nListed)
    g$par_names <- paste0(rep(f$par_names, nListed), "[", rep(listed, each=nPar), "]")
    g$trust_cap <- rep(f$trust_cap, nListed)
    g$report$transform <- rep(f$report$transform, nListed)
    g$report$lower <- rep(f$report$lower, nListed)
    g$report$upper <- rep(f$report$upper, nListed)
    g$section_var <- byExpr
    g$section_levels <- sections
    g$section_index <- sectionIndex
    g$section_npar <- nPar
    g$section_factor <- 1L
    g$native_report <- list(backend="R", fun=function(scale, par, factor, absorb_scale=TRUE){
      values <- lapply(seq_len(nListed), function(s){
        v <- reportOne(scale=scale, par=par[(s - 1L) * nPar + seq_len(nPar)],
                       factor=template, absorb_scale=absorb_scale)
        stats::setNames(v, paste0(names(v), "[", listed[s], "]"))
      })
      unlist(values)
    })
    .validate_covfactor(g)
    factors[[length(factors) + 1L]] <- append_factor(g)
  }

  covStruct <- inner
  covStruct$par <- stats::setNames(as.numeric(par), parNames)
  covStruct$free <- as.logical(free)
  covStruct$par_names <- parNames
  covStruct$factors <- factors
  covStruct$dim <- as.integer(nSections * qInner)
  covStruct$levels <- paste(rep(sections, each=qInner), rep(inner$levels, nSections), sep=":")

  observed <- which(!is.na(local))
  productDesign <- Matrix::sparseMatrix(
    i=observed, j=local[observed], x=1, dims=c(n, nSections * qInner)
  )

  list(
    Z=list(),
    Gu=x$Gu,
    GuRot=NULL,
    rotation=NULL,
    covStruct=covStruct,
    residualLocalIndex=local,
    residualSelfMasked=TRUE,
    productDesign=productDesign,
    partitionsR=NULL
  )
}
