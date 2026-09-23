vsm <- function(..., Gu=NULL, sigma2=NULL, fixedSigma2=FALSE,
                isFixed=FALSE, verbose=TRUE){
  
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
        "ism(), dsm(), usm(), ar1m(), csm(), rrcm(), fam(), ",
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
  
  if(!is.null(Gu)){
    
    Gu_is_inverse <- isTRUE(attr(Gu, "inverse"))
    
    if(!Gu_is_inverse){
      stop(
        "Gu must have attr(Gu, 'inverse')=TRUE for the Henderson solver.",
        call. = FALSE
      )
    }
    
    # Conversion can drop custom attributes.
    Gu <- to_sparse(Gu)
    
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
          "Adding additional Gu levels to the main-effect model matrix:",
          paste(extra, collapse=", "),
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
    
    out <- vector(
      "list",
      ncol(A) * ncol(B)
    )
    
    nm <- character(length(out))
    
    cc <- 1L
    
    for(i in seq_len(ncol(A))){
      for(j in seq_len(ncol(B))){
        
        out[[cc]] <-
          A[,i,drop=FALSE] *
          B[,j,drop=FALSE]
        
        nm[cc] <-
          paste(
            colnames(A)[i],
            colnames(B)[j],
            sep=":"
          )
        
        cc <- cc + 1L
      }
    }
    
    ans <- do.call(cbind, out)
    
    colnames(ans) <- nm
    
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
      
      mask <-
        Z0[,j,drop=FALSE] %*%
        Matrix::Matrix(
          1,
          1,
          ncol(mainZ)
        )
      
      Z[[j]] <-
        to_sparse(
          mainZ * mask
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
    Gu <- Gu[
      colnames(mainZ),
      colnames(mainZ),
      drop=FALSE
    ]
    
    Gu <- to_sparse(Gu)
    
    attr(Gu, "inverse") <- TRUE
  }
  
  
  # Final representation validation.
  if(!inherits(Gu, "dgCMatrix")){
    stop(
      paste0(
        "Internal vsm() error: Gu was not normalized ",
        "to dgCMatrix."
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
      
      f <- .compile_covfactor(
        factors[[i]]
      )
      
      .validate_covfactor(f)
      
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
    
    byrow <- split(
      seq_len(nrow(ss)),
      ss$i
    )
    
    residualLocalIndex <-
      rep(
        NA_integer_,
        nrow(Z0)
      )
    
    for(rr in seq_len(nrow(Z0))){
      
      hits <-
        byrow[[as.character(rr)]]
      
      if(is.null(hits)){
        next
      }
      
      if(
        length(hits) != 1L ||
        abs(ss$x[hits] - 1) > 1e-12
      ){
        stop(
          paste0(
            "Residual covariance factors must define exactly ",
            "one covariance-product level per observation."
          ),
          call. = FALSE
        )
      }
      
      residualLocalIndex[rr] <-
        ss$j[hits]
    }
  }
  
  
  # ======================================================================
  # Return object
  # ======================================================================
  
  output <- list(
    
    Z=Z,
    
    Gu=Gu,
    
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

## small matrix constructors
unsm <- function(x, reps=NULL){
  mm <- matrix(1,x,x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 2
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
}

unsm2 <- function (x, reps = NULL) {
  mm <- matrix(1, x, x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 0
  if (!is.null(reps)) {
    return(rep(list(mm), reps))
  }
  else {
    return(mm)
  }
}

fixm <- function(x, reps=NULL){
  mm <- matrix(3,x,x)
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
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
  
  G1 <- to_sparse(
    ran1$Gu
  )
  
  G2 <- to_sparse(
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
  
  Gu <- to_sparse(Gu)
  
  attr(Gu, "inverse") <- TRUE
  
  
  if(!inherits(Gu, "dgCMatrix")){
    stop(
      paste0(
        "Internal covm() error: Gu was not normalized ",
        "to dgCMatrix."
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


replace.values <- function(Values,Search,Replace){
  dd0 <- data.frame(Values)
  vv <- which(Values%in%Search)
  dd <- data.frame(Search,Replace)
  rownames(dd) <- Search
  dd0[vv,"Values"] <- as.character(dd[Values[vv],"Replace"])
  return(dd0[,1])
}

myformula <- function(x){
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  yuyuf <- strsplit(as.character(x[3]), split = "[+]")[[1]]
  termss <- apply(data.frame(yuyuf),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  newtermss <- apply(data.frame(yuyuf),1,function(y){
    newy <- expi(y)
    if(length(newy) > 0){
      newy <- gsub(",.*","",newy)
    }else{newy <- y}
    return(newy)
  })
  resp <- strsplit(as.character(x[2]), split = "[+]")[[1]]
  newx <- paste(resp, "~",paste(newtermss,collapse = "+"))
  return(newx)
}

##############
## na.methods

subdata <- function(data,fixed,na.method.Y=NULL,na.method.X=NULL){
  
  # silently change all columns that are defined as character into factors
  # columnTypes <- unlist(lapply(data, class))
  # columnTypesC <- which(columnTypes == "character")
  # if(length(columnTypesC) > 0){ # if there's character types change them to factor
  #   for(cti in 1:length(columnTypesC)){
  #     data[,cti] <- as.factor(data[,cti])
  #   }
  # }
  ####
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (is(mfna, "try-error") ) { # class(mfna) == "try-error"
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, parent.frame())
  yvar <- as.matrix(model.response(mfna))
  nt <- ncol(yvar)
  good <- 1:nrow(data)
  if(nt==1){colnames(yvar) <- response}
  if(na.method.Y=="include"){
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      use <- touse[i]
      # print(iname)
      data[,use] <- imputev(data[,use])
    }
  }else if(na.method.Y=="include2"){
    tlist <- list()
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      # print(touse[i])
      use <- touse[i]
      vivi <- as.vector(data[,use])
      tlist[[i]] <- which(!is.na(vivi))
    }
    # print(tlist)
    good <- sort(unique(unlist(tlist)))
    data <- data[good,]
    for(i in 1:length(touse)){
      use <- touse[i]
      # print(iname)
      data[,use] <- imputev(data[,use])
    }
  }else if(na.method.Y=="exclude"){
    tlist <- list()
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      # print(touse[i])
      use <- touse[i]
      vivi <- as.vector(data[,use])
      tlist[[i]] <- which(!is.na(vivi))
    }
    # print(tlist)
    if(length(tlist)==1){ #only one trait
      good <- tlist[[1]]
    }else{#more than one trait
      good <- Reduce(intersect,tlist)
    }
    data <- data[good,]
  }else{stop("na.method.Y not recognized")}
  data <- data.frame(data)
  
  ##########
  ## na.method x
  yuyu <- strsplit(as.character(fixed[3]), split = "[+]")[[1]]
  xtermss <- apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  xtermss2 <- apply(data.frame(xtermss),1,function(x){gsub(",.*","",expi2(x))})
  xtermss2[which(xtermss2 == "")] <- xtermss[which(xtermss2 == "")]
  
  xtermss2 <- intersect(colnames(data),xtermss2) # only focus on the terms that are in teh dataset so we can skip overlay and weird vs structures
  
  # print(xtermss2)
  if(length(xtermss2) > 0){
    mycl <- as.vector(unlist(lapply(data.frame(data[,xtermss2]),class)))
    
    if(na.method.X=="include"){
      touse <- xtermss2
      for(i in 1:length(touse)){
        use <- touse[i]
        usecl <- mycl[i]
        if(usecl == "factor"){data[,use] <- as.factor(imputev(data[,use]))}else{data[,use] <- imputev(data[,use])}
      }
    }else if(na.method.X=="exclude"){
      tlist <- list()
      touse <- xtermss2
      for(i in 1:length(touse)){
        # print(touse[i])
        use <- touse[i]
        vivi <- as.vector(data[,use])
        tlist[[i]] <- which(!is.na(vivi))
      }
      # print(tlist)
      if(length(tlist)==1){ #only one trait
        good <- tlist[[1]]
      }else{#more than one trait
        good <- Reduce(intersect,tlist)
      }
      data <- data[good,]
    }else{stop("na.method.Y not recognized")}
    data <- data.frame(data)
  }
  
  return(list(datar=data,good=good))
  
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
      levs <- levels(x)
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
.make_covfactor <- function(dim, levels, par=numeric(), free=logical(),
                            par_names=character(), evaluator,
                            derivative=list(backend="numeric", rel_step=1e-6),
                            report=NULL, trust_cap=NULL,
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
    trust_cap=as.numeric(trust_cap),
    structurally_diagonal=isTRUE(structurally_diagonal),
    descriptor_version=2L,
    model=model
  ), metadata)
  class(out) <- c("sommer_covfactor", "list")
  .validate_covfactor(out)
  out
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
#   trust_cap  : per-parameter maximum proposal size in working coordinates
#
# ai_mme_sp2() never needs to know the statistical model label.  Built-in
# high-use structures use the native covariance backend; user-defined and new
# experimental structures can use R callbacks without changing the solver.
# -------------------------------------------------------------------------
.compile_covfactor <- function(f){
  if(inherits(f, "sommer_covfactor") &&
     !is.null(f$evaluator) && !is.null(f$derivative) &&
     !is.null(f$report) && !is.null(f$trust_cap)){
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
    # Reduced-rank is intentionally compiled to the generic callback backend.
    # This is the proof-of-concept for future covariance structures: adding
    # rrcm() requires no native C++ evaluator or solver-specific branch.
    order <- as.integer(f$order)
    rows <- as.integer(f$rr_row)
    cols <- as.integer(f$rr_col)
    dd <- as.logical(f$rr_diag)
    nload <- as.integer(f$rr_nload)
    if(length(rows) != nload || length(cols) != nload ||
       length(dd) != nload || p != nload){
      stop("Malformed reduced-rank covariance metadata.", call. = FALSE)
    }

    rr_eval <- local({
      q0 <- q; k0 <- order; r0 <- rows; c0 <- cols; d0 <- dd
      function(par){
        L <- matrix(0, q0, k0)
        for(a in seq_along(par)){
          L[r0[a], c0[a]] <- if(d0[a]) exp(par[a]) else par[a]
        }
        M <- tcrossprod(L) + diag(q0)
        M / M[1,1]
      }
    })

    rr_d1 <- local({
      q0 <- q; k0 <- order; r0 <- rows; c0 <- cols; d0 <- dd
      function(par, k){
        L <- matrix(0, q0, k0)
        for(a in seq_along(par)){
          L[r0[a], c0[a]] <- if(d0[a]) exp(par[a]) else par[a]
        }
        dL <- matrix(0, q0, k0)
        dL[r0[k], c0[k]] <- if(d0[k]) exp(par[k]) else 1
        M <- tcrossprod(L) + diag(q0)
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

  }else{
    stop(paste0("Unknown covariance model label '", model,
                "'. New models should provide a sommer_covfactor descriptor ",
                "or use ownm()."), call. = FALSE)
  }

  f$model <- model                 # human-readable/debugging only
  f$type <- NULL                   # model label never crosses as solver dispatch
  f$evaluator <- evaluator
  f$derivative <- derivative
  f$report <- list(
    backend="builtin",
    transform=as.character(transform),
    lower=as.numeric(lower),
    upper=as.numeric(upper)
  )
  f$trust_cap <- as.numeric(trust)
  f$structurally_diagonal <- isTRUE(structurally_diagonal)
  f$descriptor_version <- 2L
  class(f) <- c("sommer_covfactor", "list")
  f
}

.validate_covfactor <- function(f){
  required <- c("dim","levels","par","free","par_names","evaluator",
                "derivative","report","trust_cap","structurally_diagonal",
                "descriptor_version")
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
  invisible(TRUE)
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

ar1m <- function(x, rho=0.30, fixed=FALSE,
                 variance=c("homogeneous", "heterogeneous"), values=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  variance <- match.arg(variance)
  if(q < 2L) stop("ar1m() requires at least two ordered levels.", call. = FALSE)
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
      par_names=if(q > 1L) paste0("variance_ratio[", colnames(dummy)[-1], "]") else character()
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
      us_diag=isdiag
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
  rho[abs(row(matrix(0, q, q)) - col(matrix(0, q, q))) + 1L]
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
      model=paste0("ar", order)
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


# General positive-definite correlation matrix.
#
# We use an unconstrained unit-diagonal lower factor A and standardize
# A A' to correlation scale.  This spans the SPD correlation cone without
# requiring pairwise-correlation boundary repairs.
corgm <- function(x, theta=NULL, fixed=NULL){
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
    K <- stats::cov2cor(K)
    ch <- try(chol(K), silent=TRUE)
    if(inherits(ch, "try-error")){
      stop("theta supplied to corgm() must be positive definite.", call. = FALSE)
    }
    A <- t(ch)
    A <- sweep(A, 1, diag(A), "/")
    vals <- mapply(function(i,j) A[i,j], rows, cols)
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


# Factor analytic covariance shape of rank k:
#   M = Lambda Lambda' + Psi,  Psi diagonal positive,
# then K = M / M[1,1].
#
# Lambda is lower-triangular in its leading k x k block for rotational
# identification. Leading loading diagonals are positive (log parameterized).
fam <- function(x, k=1L, loadings=NULL, specific=NULL, fixed=NULL){
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
  if(length(specific) != q || any(!is.finite(specific)) || any(specific <= 0)){
    stop("specific in fam() must contain q positive finite specific variances.",
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
  
  par <- c(load_par, log(specific[-1]))
  par_names <- c(load_names, paste0("specific_ratio[", labs[-1], "]"))
  
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
      fa_specific_reference=labs[1]
    ))
  )
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
                 dfun=NULL, par_names=NULL){
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
      dfun=wrapped_dfun
    ))
  )
}



# Reduced-rank covariance approximation of rank k:
#   M = Lambda Lambda' + I,
# followed by K = M / M[1,1].
#
# The rank-k component Lambda Lambda' captures the dominant covariance pattern,
# while the identity term supplies a common isotropic remainder.  This keeps the
# covariance strictly positive definite, which is required by the current
# Henderson precision-based implementation.  The leading k x k block of Lambda
# is lower triangular for rotational identification and its diagonal entries are
# positive through a log parameterization.
rrcm <- function(x, k=1L, loadings=NULL, fixed=NULL){
  expr <- as.character(substitute(x))
  dummy <- .cov_dummy(x, expr)
  q <- ncol(dummy)
  labs <- colnames(dummy)
  
  k <- as.integer(k)
  if(length(k) != 1L || k < 1L || k >= q){
    stop("k in rrcm() must satisfy 1 <= k < number of levels.", call. = FALSE)
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
      stop("loadings in rrcm() must be a finite q x k matrix.", call. = FALSE)
    }
    for(j in seq_len(k)){
      if(j > 1L) L[seq_len(j-1L),j] <- 0
      if(L[j,j] <= 0){
        stop("Leading diagonal loadings in rrcm() must be positive.", call. = FALSE)
      }
    }
  }
  
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
    stop("fixed in rrcm() must match the number of loading parameters.", call. = FALSE)
  }
  
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
      rr_diag=diag_loading
    ))
  )
}

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
                    distance=NULL){
  expr <- paste(deparse(substitute(x)), collapse="")
  sx <- .spatial_cov_dummy(x, expr)
  Z <- sx$Z
  coords <- sx$coords
  q <- ncol(Z)

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
  if(length(fixed) == 1L) fixed <- rep(fixed, 2L)
  if(length(fixed) != 2L){
    stop("fixed in maternm() must have length 1 or 2 (range, nu).",
         call. = FALSE)
  }

  eval_fun <- local({
    D0 <- D
    function(par){
      r <- exp(par[1])
      v <- exp(par[2])
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

  cf <- .make_covfactor(
    dim=q,
    levels=sx$levels,
    par=c(log_range=log(range), log_nu=log(nu)),
    free=!as.logical(fixed),
    par_names=c("range","nu"),
    evaluator=list(backend="R", fun=eval_fun),
    derivative=list(backend="numeric", rel_step=1e-6),
    report=list(backend="builtin",
                transform=c("exp","exp"),
                lower=c(NA_real_,NA_real_),
                upper=c(NA_real_,NA_real_)),
    trust_cap=c(1.0,0.75),
    structurally_diagonal=FALSE,
    model="matern",
    metadata=list(distance=D, coordinates=coords)
  )

  list(Z=Z, covFactor=cf)
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
