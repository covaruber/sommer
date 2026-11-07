shrink_marker_data <- function(){
  set.seed(2917)
  X <- matrix(sample(c(-1, 0, 1), 12 * 25, replace=TRUE), 12)
  X[1:3,] <- -1
  X[,24] <- 0
  X[,25] <- 1
  dimnames(X) <- list(paste0("id", 1:12), paste0("m", 1:25))
  X
}

shrink_reference <- function(X, min.MAF=0, ploidy=NULL){
  tr <- sommer:::.amat_marker_transform(X, min.MAF, ploidy)
  W <- tr$Z
  m <- ncol(W)
  mu <- rowMeans(W)
  Z <- sweep(W, 1, mu)
  S <- tcrossprod(Z) / m
  target <- diag(mean(diag(S)), nrow(S))
  b2 <- 0
  for(j in seq_len(m)){
    b2 <- b2 + sum((tcrossprod(Z[,j]) - S)^2) / m^2
  }
  d2 <- sum((S - target)^2)
  intensity <- if(d2 == 0) 0 else max(0, min(1, b2 / d2))
  A <- ((1 - intensity) * S + intensity * target + tcrossprod(mu)) * m / tr$k
  colnames(A) <- rownames(A) <- rownames(X)
  attr(A, "shrinkage") <- intensity
  list(A=A, intensity=intensity, X=tr$X)
}

test_that("native EJ shrinkage matches the covariance estimator and retains row means", {
  X <- shrink_marker_data()
  for(maf in c(0, .1, .3)){
    ref <- shrink_reference(X, min.MAF=maf)
    A <- A.mat(X, min.MAF=maf, shrink=TRUE)
    expect_equal(A, ref$A, tolerance=1e-12)
    expect_gt(attr(A, "shrinkage"), 0)
    expect_lt(attr(A, "shrinkage"), 1)
    expect_equal(mean(diag(A)), mean(diag(A.mat(X, min.MAF=maf))), tolerance=1e-12)
    expect_equal(A, t(A), ignore_attr=TRUE, tolerance=1e-12)
    expect_gte(min(eigen(A, symmetric=TRUE, only.values=TRUE)$values), -1e-12)
  }
  expect_identical(A.mat(X), A.mat(X, shrink=FALSE))
  expect_null(attr(A.mat(X), "shrinkage"))
  expect_identical(A.mat(X, 0, TRUE, NULL, FALSE), A.mat(X, return.imputed=TRUE))
  A <- A.mat(X, shrink=TRUE)
  perm <- rev(seq_len(nrow(X)))
  expect_equal(unname(A.mat(X[perm,,drop=FALSE], shrink=TRUE)),
               unname(A[perm,perm]), ignore_attr=TRUE, tolerance=1e-12)
  expect_equal(A.mat(-X, shrink=TRUE), A, tolerance=1e-12)
})

test_that("EJ shrinkage preserves mean-imputation and polyploid dosage scaling", {
  X <- shrink_marker_data()
  X[4:5,1] <- NA_real_
  for(ploidy in list(NULL, 2, 4, 6)){
    dosage <- if(is.null(ploidy)) X else (X + 1) * ploidy / 2
    capture.output(ref <- shrink_reference(dosage, min.MAF=.1, ploidy=ploidy))
    output <- capture.output(result <- A.mat(dosage, min.MAF=.1, ploidy=ploidy,
                                             shrink=TRUE, return.imputed=TRUE))
    expect_named(result, c("X", "A"))
    expect_equal(result$X, ref$X)
    expect_equal(result$A, ref$A, tolerance=1e-12)
    expect_identical(dimnames(result$X), dimnames(dosage))
    expect_length(grep("Imputing", output), 1L)
  }
  X <- shrink_marker_data()
  diploid <- A.mat(X, shrink=TRUE)
  for(ploidy in c(2, 3, 4, 8, 49)){
    actual <- A.mat((X + 1) * ploidy / 2, ploidy=ploidy, shrink=TRUE)
    expect_equal(attr(actual, "shrinkage"), attr(diploid, "shrinkage"), tolerance=1e-12)
    expect_equal(actual, diploid * ploidy / 2, tolerance=1e-12)
  }
})

test_that("EJ shrinkage handles zero covariance and clipped intensity", {
  X <- matrix(rep(c(-1, 0, 1), 4), nrow=3)
  expect_equal(attr(A.mat(X, shrink=TRUE), "shrinkage"), 0)
  expect_equal(A.mat(X, shrink=TRUE), A.mat(X), ignore_attr=TRUE, tolerance=1e-12)
  expect_equal(A.mat(X[,1,drop=FALSE], shrink=TRUE),
               A.mat(X[,1,drop=FALSE]), ignore_attr=TRUE, tolerance=1e-12)
  X <- matrix(rep(c(-1, -1, 0), 10), nrow=3)
  X[,10] <- c(-1, 1, -1)
  expect_equal(attr(A.mat(X, shrink=TRUE), "shrinkage"), 1)
  expect_equal(A.mat(X, shrink=TRUE), shrink_reference(X)$A, tolerance=1e-12)
  native <- sommer:::amat_shrink(matrix(c(1, 0, -1, 0, 0, 1, 0, -1), 2), 1)
  expect_equal(native$intensity, 0)
  expect_equal(native$A, 2 * diag(2))
})

test_that("EJ shrinkage matches rrBLUP for the same retained marker set", {
  skip_if_not_installed("rrBLUP")
  X <- shrink_marker_data()
  for(missing in c(FALSE, TRUE)){
    if(missing) X[4:5,1] <- NA_real_
    for(maf in c(0, .1, .3)){
      capture.output(tr <- sommer:::.amat_marker_transform(X, maf))
      capture.output(expected <- rrBLUP::A.mat(X[,tr$keep,drop=FALSE],
                                               min.MAF=0, max.missing=1, shrink=TRUE))
      capture.output(actual <- A.mat(X, min.MAF=maf, shrink=TRUE))
      expect_equal(actual, expected, ignore_attr=TRUE, tolerance=1e-12)
    }
  }
})

test_that("invalid shrinkage requests and unusable native inputs fail explicitly", {
  X <- shrink_marker_data()
  for(shrink in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "EJ",
                     list(method="EJ"), list(method="REG"))){
    expect_error(A.mat(X, shrink=shrink), "shrink must be TRUE or FALSE")
  }
  expect_error(A.mat(matrix(0, 4, 2), shrink=TRUE), "No variable markers")
  expect_error(A.mat(X, min.MAF=.5, shrink=TRUE), "No variable markers")
  X[1,1] <- Inf
  expect_error(A.mat(X, shrink=TRUE), "finite genotypes")
  expect_error(sommer:::amat_shrink(matrix(numeric(), 2, 0), 1), "finite matrix")
  expect_error(sommer:::amat_shrink(matrix(1, 1, 2), 1), "finite matrix")
  expect_error(sommer:::amat_shrink(matrix(NA_real_, 2, 2), 1), "finite matrix")
  for(k in c(0, -1, NA_real_, Inf)){
    expect_error(sommer:::amat_shrink(matrix(1, 2, 2), k), "normalization")
  }
})
