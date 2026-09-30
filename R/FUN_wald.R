#### =========== ####
## WALD TESTS AND DENOMINATOR DEGREES OF FREEDOM
#### =========== ####

.wald_fixed_vcov <- function(object){
  p <- nrow(object$b)
  if(identical(object$engine, "direct")){
    if(is.null(object$VarBeta)) stop("The fitted object does not contain VarBeta.", call.=FALSE)
    return(as.matrix(object$VarBeta))
  }
  D <- cbind(Matrix::Diagonal(p),
             Matrix::Matrix(0, p, nrow(object$bu) - p, sparse=TRUE))
  D <- as(as(D, "generalMatrix"), "CsparseMatrix")
  as.matrix(predict_mmes_vcov_cpp(object, .mmes_engine_contrast(object, D)))
}

.wald_term_columns <- function(object){
  cols <- lapply(object$partitionsX, function(x) as.integer(as.vector(x)))
  cols[vapply(cols, length, integer(1)) > 0L]
}

.wald_term_vars <- function(term){
  if(identical(term, "1")) character(0) else strsplit(term, ":", fixed=TRUE)[[1]]
}

.sparse_rank <- function(A){
  A <- as(as(A, "generalMatrix"), "CsparseMatrix")
  as.integer(Matrix::rankMatrix(A, method="qr"))
}

.fixed_rank <- function(object){
  .sparse_rank(object$W[, seq_len(nrow(object$b)), drop=FALSE])
}

# Term a contains term t structurally (factor sets) or implicitly (column space).
.wald_contains <- function(object, cols, a, t){
  if(a == t) return(FALSE)
  va <- .wald_term_vars(a); vt <- .wald_term_vars(t)
  if(all(vt %in% va)) return(TRUE)
  if(identical(t, "1")) return(FALSE)
  one <- Matrix::Matrix(1, nrow(object$W), 1, sparse=TRUE)
  base <- cbind(one, object$W[, cols[[a]], drop=FALSE])
  .sparse_rank(cbind(base, object$W[, cols[[t]], drop=FALSE])) == .sparse_rank(base)
}

# Rows of the Cholesky factor of the information matrix give the contrast L of
# an incremental Wald test in the requested order; L Phi L' = I by construction.
.wald_contrast <- function(M, cols, order, target){
  idx <- unlist(cols[order])
  U <- tryCatch(chol(M[idx, idx, drop=FALSE]), error=function(e) NULL)
  if(is.null(U)){
    stop("The fixed-effect information matrix is not positive definite (aliased fixed effects?).",
         call.=FALSE)
  }
  positions <- which(idx %in% cols[[target]])
  L <- matrix(0, length(positions), ncol(M))
  L[, idx] <- U[positions, , drop=FALSE]
  L
}

.mmes_rebuild_inputs <- function(object){
  if(is.null(object$inputArgs) || is.null(object$args$fixed)){
    stop("Denominator degrees of freedom require a fit from this sommer version; please refit the model.",
         call.=FALSE)
  }
  args <- c(list(fixed=object$args$fixed, rcov=object$args$rcov,
                 data=object$dataOriginal, returnParam=TRUE, verbose=FALSE,
                 dateWarning=FALSE),
            object$inputArgs[!vapply(object$inputArgs, is.null, logical(1))])
  if(!is.null(object$args$random)) args$random <- object$args$random
  inputs <- tryCatch(do.call(mmes, args), error=function(e){
    stop("Could not rebuild the model inputs: ", conditionMessage(e), call.=FALSE)
  })
  if(nrow(inputs$X) != nrow(object$W) || ncol(inputs$X) != nrow(object$b)){
    stop("The rebuilt model inputs do not match the fitted object; the data may have changed.",
         call.=FALSE)
  }
  inputs
}

.natural_to_working <- function(cs, theta){
  work <- theta
  work[1L] <- log(theta[1L])
  for(f in cs$factors){
    s <- as.integer(f$par_start); e <- as.integer(f$par_end)
    if(e < s) next
    report <- f$report
    if(is.null(report) || !identical(report$backend, "builtin")){
      stop("Denominator degrees of freedom are not available for covariance factors without a builtin report transform.",
           call.=FALSE)
    }
    for(k in seq_len(e - s + 1L)){
      idx <- s + k - 1L
      v <- theta[idx]
      work[idx] <- switch(report$transform[k],
        identity=v,
        exp=log(v),
        tanh=atanh(v),
        bounded_logit={
          lo <- report$lower[k]; hi <- report$upper[k]
          stats::qlogis((v - lo) / (hi - lo))
        },
        stop("Unknown report transform: ", report$transform[k], call.=FALSE))
    }
  }
  work
}

# Information of the fixed effects, X'V^{-1}X, at natural covariance parameters.
.mmes_fixed_information_at <- function(object, inputs, thetaList){
  cs <- inputs$covStruct
  for(i in seq_along(cs)){
    cs[[i]]$par <- .natural_to_working(cs[[i]], thetaList[[i]])
    cs[[i]]$free[] <- FALSE
  }
  p <- ncol(inputs$X)
  if(isTRUE(inputs$henderson)){
    solver <- if(identical(object$solver, "cholmod")) "cholmod" else "ldlt"
    r <- .Call("_sommer_ai_mme_sp2", PACKAGE="sommer",
               inputs$X, inputs$Z, inputs$Zind, inputs$Ai, inputs$yvar, inputs$W,
               inputs$useH, inputs$residualBlock, inputs$residualIndex,
               1L, 1e-4, 1e-4, inputs$tolParInv, cs, 0, 1, FALSE, 0L, solver,
               1e-8, 0L, 8L, 20L, isTRUE(inputs$REML),
               inputs$responsePrepared, inputs$preparedMean, inputs$preparedSd,
               inputs$preparedIntercept)
    D <- cbind(Matrix::Diagonal(p), Matrix::Matrix(0, p, nrow(r$C) - p, sparse=TRUE))
    D <- as(as(D, "generalMatrix"), "CsparseMatrix")
    Phi <- as.matrix(predict_mmes_vcov_cpp(r, D))
  }else{
    r <- .Call("_sommer_ai_reml_direct_sp2", PACKAGE="sommer",
               inputs$X, inputs$Z, inputs$Zind, inputs$Ai, inputs$yvar, inputs$W,
               inputs$useH, inputs$residualBlock, inputs$residualIndex,
               1L, 1e-4, 1e-4, inputs$tolParInv, cs, 0, 1, FALSE, 0L,
               isTRUE(inputs$REML), inputs$responsePrepared, inputs$preparedMean,
               inputs$preparedSd, inputs$preparedIntercept)
    Phi <- as.matrix(r$VarBeta)
  }
  solve(Phi)
}

# Numerical derivatives of M = X'V^{-1}X with respect to the free covariance
# parameters on their reported scale (the scale of theta_se).
.mmes_df_machinery <- function(object, second=FALSE, rel_step=1e-3){
  if(inherits(object, "mmes.glmm")){
    stop("Satterthwaite and Kenward-Roger degrees of freedom are not available for PQL fits.",
         call.=FALSE)
  }
  if(isFALSE(object$REML)){
    stop("Satterthwaite and Kenward-Roger degrees of freedom require a REML fit.", call.=FALSE)
  }
  inputs <- .mmes_rebuild_inputs(object)
  thetaList <- lapply(object$covPar, as.numeric)
  counts <- vapply(thetaList, length, integer(1))
  theta <- unlist(thetaList, use.names=FALSE)
  free <- unlist(lapply(object$covStruct, function(x) as.logical(x$free)), use.names=FALSE)
  if(length(free) != length(theta) || !all(dim(object$theta_se) == length(theta))){
    stop("Covariance parameters and theta_se do not conform.", call.=FALSE)
  }
  freeIdx <- which(free)
  Wm <- as.matrix(object$theta_se)[freeIdx, freeIdx, drop=FALSE]
  owner <- rep(seq_along(counts), counts)
  relist <- function(v) split(v, owner)

  evalM <- function(v){
    tryCatch(.mmes_fixed_information_at(object, inputs, relist(v)), error=function(e){
      stop("Evaluating X'V^{-1}X at perturbed covariance parameters failed: ",
           conditionMessage(e), call.=FALSE)
    })
  }
  M0 <- evalM(theta)
  h <- rel_step * pmax(abs(theta[freeIdx]), 1e-4)
  shift <- function(ks, signs){
    v <- theta
    v[freeIdx[ks]] <- v[freeIdx[ks]] + signs * h[ks]
    evalM(v)
  }
  nf <- length(freeIdx)
  Mplus <- lapply(seq_len(nf), function(k) shift(k, 1))
  Mminus <- lapply(seq_len(nf), function(k) shift(k, -1))
  P <- lapply(seq_len(nf), function(k) (Mplus[[k]] - Mminus[[k]]) / (2 * h[k]))
  H <- NULL
  if(second){
    H <- vector("list", nf * nf)
    dim(H) <- c(nf, nf)
    for(i in seq_len(nf)){
      H[[i, i]] <- (Mplus[[i]] - 2 * M0 + Mminus[[i]]) / h[i]^2
      for(j in seq_len(i - 1L)){
        H[[i, j]] <- (shift(c(i, j), c(1, 1)) - shift(c(i, j), c(1, -1)) -
                        shift(c(i, j), c(-1, 1)) + shift(c(i, j), c(-1, -1))) /
          (4 * h[i] * h[j])
        H[[j, i]] <- H[[i, j]]
      }
    }
  }
  Phi <- solve(M0)
  out <- list(Phi=Phi, P=P, H=H, W=Wm)
  if(second){
    # Kenward-Roger adjusted covariance; the V second-derivative term enters through H.
    S1 <- matrix(0, nrow(Phi), ncol(Phi))
    S2 <- S1
    for(i in seq_len(nf)) for(j in seq_len(nf)){
      S1 <- S1 + Wm[i, j] * H[[i, j]]
      S2 <- S2 + Wm[i, j] * (P[[i]] %*% Phi %*% P[[j]])
    }
    out$PhiA <- Phi + Phi %*% (S1 - 2 * S2) %*% Phi
    out$PhiA <- (out$PhiA + t(out$PhiA)) / 2
  }
  out
}

.satterthwaite_ddf <- function(L, mach){
  Phi <- mach$Phi
  VL <- L %*% Phi %*% t(L)
  e <- eigen((VL + t(VL)) / 2, symmetric=TRUE)
  q <- nrow(L)
  nu <- vapply(seq_len(q), function(k){
    l <- drop(t(e$vectors[, k]) %*% L)
    g <- vapply(mach$P, function(Pi) -drop(t(l) %*% Phi %*% Pi %*% Phi %*% l), numeric(1))
    2 * e$values[k]^2 / drop(t(g) %*% mach$W %*% g)
  }, numeric(1))
  if(q == 1L) return(nu)
  E <- sum(ifelse(nu > 2, nu / (nu - 2), 0))
  if(E > q) 2 * E / (E - q) else NA_real_
}

.kenward_roger <- function(L, b, mach){
  Phi <- mach$Phi
  q <- nrow(L)
  Theta <- t(L) %*% solve(L %*% Phi %*% t(L)) %*% L
  G <- lapply(mach$P, function(Pi) Theta %*% Phi %*% Pi %*% Phi)
  nf <- length(G)
  A1 <- 0; A2 <- 0
  for(i in seq_len(nf)) for(j in seq_len(nf)){
    A1 <- A1 + mach$W[i, j] * sum(diag(G[[i]])) * sum(diag(G[[j]]))
    A2 <- A2 + mach$W[i, j] * sum(G[[i]] * t(G[[j]]))
  }
  B <- (A1 + 6 * A2) / (2 * q)
  g <- ((q + 1) * A1 - (q + 4) * A2) / ((q + 2) * A2)
  den <- 3 * q + 2 * (1 - g)
  c1 <- g / den; c2 <- (q - g) / den; c3 <- (q + 2 - g) / den
  Estar <- 1 / (1 - A2 / q)
  Vstar <- (2 / q) * (1 + c1 * B) / ((1 - c2 * B)^2 * (1 - c3 * B))
  rho <- Vstar / (2 * Estar^2)
  m <- 4 + (q + 2) / (q * rho - 1)
  lambda <- m / (Estar * (m - 2))
  Lb <- L %*% b
  Fstat <- drop(t(Lb) %*% solve(L %*% mach$PhiA %*% t(L)) %*% Lb) / q
  list(F=lambda * Fstat, ddf=m, lambda=lambda)
}

wald_mmes <- function(object, ssType=c("incremental", "conditional"),
                      denDF=c("none", "residual", "satterthwaite", "kr"),
                      terms=NULL){
  if(!inherits(object, "mmes")) stop("object must inherit from class 'mmes'.", call.=FALSE)
  ssType <- match.arg(ssType)
  denDF <- match.arg(denDF)
  cols <- .wald_term_columns(object)
  termNames <- names(cols)
  if(is.null(terms)) terms <- termNames
  unknown <- setdiff(terms, termNames)
  if(length(unknown)){
    stop("Unknown fixed terms: ", paste(unknown, collapse=", "), ". Available: ",
         paste(termNames, collapse=", "), call.=FALSE)
  }
  b <- as.numeric(object$b)
  Phi <- .wald_fixed_vcov(object)
  M <- solve(Phi)

  contrasts <- lapply(terms, function(t){
    k <- match(t, termNames)
    order <- if(ssType == "incremental") seq_along(termNames) else {
      within <- vapply(termNames, function(a) .wald_contains(object, cols, a, t), logical(1))
      c(setdiff(which(!within), k), k, which(within))
    }
    .wald_contrast(M, cols, order, t)
  })

  mach <- NULL
  if(denDF %in% c("satterthwaite", "kr")){
    mach <- .mmes_df_machinery(object, second=(denDF == "kr"))
  }
  residualDF <- nrow(object$W) - .fixed_rank(object)

  rows <- lapply(seq_along(terms), function(i){
    L <- contrasts[[i]]
    q <- nrow(L)
    Lb <- L %*% b
    wald <- drop(t(Lb) %*% solve(L %*% Phi %*% t(L)) %*% Lb)
    Fval <- wald / q
    ddf <- switch(denDF, none=NA_real_, residual=residualDF,
                  satterthwaite=.satterthwaite_ddf(L, mach), kr=NA_real_)
    if(denDF == "kr"){
      kr <- .kenward_roger(L, b, mach)
      Fval <- kr$F
      ddf <- kr$ddf
    }
    p <- if(denDF == "none") stats::pchisq(wald, q, lower.tail=FALSE) else
      stats::pf(Fval, q, ddf, lower.tail=FALSE)
    data.frame(Df=q, denDF=ddf, Wald=wald, F.value=if(denDF == "none") NA_real_ else Fval,
               p.value=p)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- terms
  attr(out, "ssType") <- ssType
  attr(out, "denDF") <- denDF
  class(out) <- c("wald.mmes", "data.frame")
  out
}

print.wald.mmes <- function(x, digits=max(3, getOption("digits") - 3), ...){
  cat("Wald tests for fixed effects (", attr(x, "ssType"), ", denDF = ",
      attr(x, "denDF"), ")\n", sep="")
  y <- x
  class(y) <- "data.frame"
  if(identical(attr(x, "denDF"), "none")){
    y$denDF <- NULL
    y$F.value <- NULL
  }
  stats::printCoefmat(y, digits=digits, has.Pvalue=TRUE, P.values=TRUE,
                      cs.ind=NULL, zap.ind=integer(0), tst.ind=integer(0),
                      na.print="")
  invisible(x)
}
