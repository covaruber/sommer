relationship_marker_data <- function(){
  X <- matrix(c(-1, 0, 1, 0, -1, 1,
                -1, -1, 0, -1, -1, -1,
                0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1,
                .2, -.4, 1, -1, .6, 0), nrow=6)
  dimnames(X) <- list(paste0("id", 1:6), paste0("m", 1:5))
  X
}

relationship_marker_reference <- function(X, min.MAF=0){
  p <- colMeans(X + 1) / 2
  keep <- pmin(p, 1 - p) > min.MAF & apply(X, 2, var) > 0
  Y <- X[, keep, drop=FALSE]
  p <- p[keep]
  pq2 <- 2 * p * (1 - p)
  W <- sweep(Y, 2, colMeans(Y))
  A <- tcrossprod(W) / sum(pq2)
  endelman <- tcrossprod(W)
  endelman <- endelman / mean(diag(endelman))
  W <- sweep(Y, 2, 2 * p - 1, "*") - Y^2
  W <- sweep(W, 2, pq2, "+")
  D <- tcrossprod(W) / sum(pq2^2)
  W <- sweep(1 - abs(Y), 2, pq2)
  list(A=A, endelman=endelman, D=D,
       Su=tcrossprod(W) / sum(pq2 * (1 - pq2)))
}

relationship_hybrid_reference <- function(A, G, index, tau=1, omega=1, ridge=1e-6){
  invert <- function(M){
    if(inherits(try(chol(M), silent=TRUE), "try-error")){
      M <- M + diag(ridge, nrow(M))
    }
    solve(M)
  }
  i <- which(index == 1)
  j <- which(index == 0)
  A22 <- A[j, j, drop=FALSE]
  Ai <- invert(A22)
  H22 <- invert(tau * invert(G) + (1 - omega) * Ai)
  B <- A[i, j, drop=FALSE] %*% Ai
  H <- A
  H[i, i] <- H[i, i, drop=FALSE] + B %*% (H22 - A22) %*% t(B)
  H[i, j] <- H[i, j, drop=FALSE] + B %*% (H22 - A22)
  H[j, i] <- t(H[i, j, drop=FALSE])
  H[j, j] <- H22
  H
}

test_that("relationship kernels retain their normalization and fractional coding", {
  X <- relationship_marker_data()
  for(maf in c(0, .1, .3)){
    ref <- relationship_marker_reference(X, maf)
    expect_equal(A.mat(X, min.MAF=maf), ref$A, tolerance=1e-12)
    expect_equal(D.mat(X, min.MAF=maf), ref$D, tolerance=1e-12)
    expect_equal(D.mat(X, FALSE, min.MAF=maf), ref$Su, tolerance=1e-12)
    expect_equal(sommer:::amat(X, FALSE, maf), unname(ref$endelman), tolerance=1e-12)
  }
  # The second column has MAF exactly 1/12, so the strict cutoff removes it.
  expect_equal(A.mat(X, min.MAF=1/12), A.mat(X[, -2, drop=FALSE]), tolerance=1e-12)
  expect_equal(A.mat(X[, 1, drop=FALSE]),
               relationship_marker_reference(X[, 1, drop=FALSE])$A)
  tiny <- matrix(c(-1e-20, 0, 2e-20), 3)
  expect_equal(sommer:::amat(tiny, FALSE, 0),
               unname(relationship_marker_reference(tiny)$endelman), tolerance=1e-12)
})

test_that("marker wrappers mean-impute and preserve dimensions and names", {
  X <- relationship_marker_data()
  X[1:2, 1] <- NA
  expected <- X
  expected[1:2, 1] <- mean(X[, 1], na.rm=TRUE)
  capture.output(a <- A.mat(X, return.imputed=TRUE))
  capture.output(d <- D.mat(X, return.imputed=TRUE))
  expect_equal(a$X, expected)
  expect_equal(d$X, expected)
  expect_equal(a$A, A.mat(expected))
  expect_equal(d$D, D.mat(expected))
  capture.output(single <- A.mat(X[, 1, drop=FALSE], return.imputed=TRUE))
  expect_equal(single$X, expected[, 1, drop=FALSE])
  capture.output(tr <- sommer:::.amat_marker_transform(X, 0))
  expect_equal(tcrossprod(tr$Z) / tr$k, a$A, tolerance=1e-12)
})

test_that("epistatic kernels share the same imputed marker formulas", {
  X <- relationship_marker_data()
  X[1:2, 1] <- NA
  capture.output(A <- A.mat(X, min.MAF=.1))
  for(nishio in c(TRUE, FALSE)){
    capture.output(D <- D.mat(X, nishio, min.MAF=.1))
    expected <- list("A#A"=A*A, "A#D"=A*D, "D#D"=D*D)
    for(type in names(expected)){
      output <- capture.output(E <- E.mat(X, nishio, type, min.MAF=.1))
      expect_equal(E, expected[[type]], tolerance=1e-12)
      expect_length(grep("Imputing", output), 1L)
    }
  }
  expect_equal(sommer:::emat(A, A), unname(A*A))
  expect_error(sommer:::emat(diag(2), diag(3)), "identical dimensions")
  for(type in list("bad", NA_character_, character(), c("A#A", "D#D"))){
    expect_error(E.mat(X, type=type), "type must be")
  }
})

test_that("unusable marker inputs fail explicitly rather than returning NaNs", {
  for(fun in list(A.mat, D.mat, E.mat)){
    expect_error(fun(matrix(0, 4, 2)), "No variable markers")
    expect_error(fun(relationship_marker_data(), min.MAF=.5), "No variable markers")
    expect_error(fun(matrix(numeric(), 4, 0)), "at least two individuals")
    expect_error(fun(matrix(c(-1, 0, 1), 1)), "at least two individuals")
    for(maf in c(-.1, .6, NA_real_, Inf)){
      expect_error(fun(relationship_marker_data(), min.MAF=maf), "min.MAF")
    }
    for(bad in c(Inf, 2)){
      X <- relationship_marker_data()
      X[1, 1] <- bad
      expect_error(fun(X), "finite genotypes")
    }
    X <- relationship_marker_data()
    X[, 1] <- NA_real_
    expect_error(fun(X), "no observed genotypes")
  }
  expect_error(sommer:::amat(matrix(NA_real_, 3, 2), TRUE, 0), "finite genotypes")
})

test_that("hybrid kernels match block and precision formulas in any index order", {
  set.seed(3701)
  Z <- matrix(rnorm(64), 8)
  A <- tcrossprod(Z)/8 + diag(8)
  dimnames(A) <- list(paste0("id", 1:8), paste0("id", 1:8))
  for(j in list(c(2, 4, 7), 5L, 1:8)){
    g <- length(j)
    Z <- matrix(rnorm(g*g), g)
    G <- tcrossprod(Z)/g + diag(g)
    dimnames(G) <- list(rownames(A)[j], rownames(A)[j])
    index <- as.numeric(!seq_len(nrow(A)) %in% j)
    for(weights in list(c(1,1), c(1.4,1), c(1.4,.8), c(0,0))){
      tau <- weights[1]; omega <- weights[2]
      ref <- relationship_hybrid_reference(A, G, index, tau, omega)
      actual <- sommer:::hmat(A, G, index, 1e-6, tau, omega)
      expect_equal(actual, unname(ref), tolerance=1e-10)
      precision <- solve(A)
      precision[j,j] <- precision[j,j,drop=FALSE] +
        tau*solve(G) - omega*solve(A[j,j,drop=FALSE])
      expect_equal(actual, unname(solve(precision)), tolerance=1e-10)
      H <- H.mat(A, G, tau, omega)
      expect_equal(H[rownames(A),rownames(A)], ref, tolerance=1e-10)
      expect_equal(rownames(H), rev(unique(c(rownames(G), rownames(A)))))
      expect_equal(H, t(H), tolerance=1e-12)
    }
  }
})

test_that("hybrid Cholesky retries preserve the existing diagonal ridge", {
  for(singularA in c(FALSE, TRUE)){
    A <- diag(3)
    if(singularA){
      A[2:3,2:3] <- 1
      A[1,2:3] <- A[2:3,1] <- .1
    }
    for(singularG in c(FALSE, TRUE)){
      G <- if(singularG) matrix(1,2,2) else diag(2)
      for(weights in list(c(1,1), c(2,1), c(1.4,.8))){
        ref <- relationship_hybrid_reference(A,G,c(1,0,0),weights[1],weights[2])
        actual <- sommer:::hmat(A,G,c(1,0,0),1e-6,weights[1],weights[2])
        expect_equal(actual, ref, tolerance=1e-8)
      }
    }
  }
  # H22 precision is zero before its ridge retry.
  expect_equal(sommer:::hmat(diag(3),diag(2),c(1,0,0),1e-6,1,2),
               diag(c(1,1e6,1e6)), tolerance=1e-10)
  expect_error(sommer:::hmat(-diag(3),diag(2),c(1,0,0),1e-6,1,1),
               "A22 is not positive definite")
  expect_error(sommer:::hmat(diag(3),-diag(2),c(1,0,0),1e-6,1,1),
               "G22 is not positive definite")
  expect_error(sommer:::hmat(diag(3),diag(2),c(1,0,0),1e-6,1,3),
               "H22 precision is not positive definite")
})

test_that("hybrid inputs are checked before subsetting or factorization", {
  A <- diag(3)
  dimnames(A) <- list(letters[1:3], letters[1:3])
  G <- A[2:3,2:3]
  expect_error(H.mat(unname(A),G), "row and column names")
  bad <- G
  dimnames(bad) <- list(c("b","z"),c("b","z"))
  expect_error(H.mat(A,bad), "must also occur in A")
  dimnames(bad) <- list(c("b","b"),c("b","b"))
  expect_error(H.mat(A,bad), "unique row and column names")
  bad <- G
  colnames(bad) <- rev(colnames(bad))
  expect_error(H.mat(A,bad), "identical")
  bad <- G
  bad[1,2] <- .5
  expect_error(H.mat(A,bad), "symmetric")
  bad[1,2] <- Inf
  expect_error(H.mat(A,bad), "finite")
  expect_error(sommer:::hmat(A,G,c(0,0),1e-6,1,1), "index must")
  expect_error(sommer:::hmat(A,G,c(1,0,2),1e-6,1,1), "index must")
  expect_error(sommer:::hmat(A,G,c(0,0,0),1e-6,1,1), "dimensions must match")
  expect_error(H.mat(A,G,tolparinv=0), "tolparinv")
  expect_error(H.mat(A,G,tau=NA_real_), "tau and omega")
  expect_error(H.mat(A,G,tau=1e308,omega=-1e308), "finite")
  A <- A[1,1,drop=FALSE]
  expect_equal(H.mat(A,A), A)
})
