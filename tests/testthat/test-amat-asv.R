asv_marker_data <- function(){
  X <- matrix(c(-1, 0, 1, -1, 0, 1,
                -1, -1, -1, -1, -1, 0,
                0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1,
                .2, -.4, 1, -1, .6, 0), nrow=6)
  dimnames(X) <- list(paste0("id", 1:6), paste0("m", 1:5))
  X
}

asv_reference <- function(X, min.MAF=0, ploidy=NULL){
  for(j in seq_len(ncol(X))){
    X[is.na(X[,j]),j] <- mean(X[,j], na.rm=TRUE)
  }
  copies <- if(is.null(ploidy)) 2 else ploidy
  dosage <- if(is.null(ploidy)) X + 1 else X
  p <- colMeans(dosage) / copies
  keep <- pmin(p, 1 - p) > min.MAF & apply(X, 2, var) > 0
  markers <- X[,keep,drop=FALSE]
  W <- sweep(markers, 2, colMeans(markers))
  k <- sum(apply(markers, 2, var))
  A <- tcrossprod(W) / k
  colnames(A) <- rownames(A) <- rownames(X)
  list(A=A, X=X, W=W, k=k, keep=keep)
}

test_that("ASV implements equation (2) and unit average pairwise semivariance", {
  X <- asv_marker_data()
  for(maf in c(0, .1, .3)){
    ref <- asv_reference(X, maf)
    A <- A.mat(X, min.MAF=maf, method="ASV")
    n <- nrow(A)
    P <- diag(n) - matrix(1/n, n, n)
    semivariance <- (outer(diag(A), diag(A), "+") - 2*A) / 2
    expect_equal(A, ref$A, tolerance=1e-12)
    expect_equal(sum(diag(P %*% A)) / (n-1), 1, tolerance=1e-12)
    expect_equal(mean(semivariance[upper.tri(semivariance)]), 1, tolerance=1e-12)
    expect_equal(mean(diag(A)), (n-1)/n, tolerance=1e-12)
    expect_equal(rowSums(A), rep(0, n), ignore_attr=TRUE, tolerance=1e-12)
    expect_equal(A, t(A), tolerance=1e-12)
    expect_gte(min(eigen(A, symmetric=TRUE, only.values=TRUE)$values), -1e-12)
    vanraden <- A.mat(X, min.MAF=maf)
    expect_equal(A, vanraden * (n-1)/sum(diag(vanraden)), tolerance=1e-12)
    expect_null(attr(A, "shrinkage"))
    tr <- sommer:::.amat_marker_transform(X, maf, method="ASV")
    expect_equal(tr$k, ref$k, tolerance=1e-12)
    expect_equal(tcrossprod(tr$Z)/tr$k, A, tolerance=1e-12)
  }
  expect_identical(A.mat(X), A.mat(X, method="VanRaden"))
  expect_identical(A.mat(X, shrink=TRUE),
                   A.mat(X, shrink=TRUE, method="VanRaden"))
  expect_equal(A.mat(matrix(c(-1, 0, 1), 3), method="ASV"),
               tcrossprod(c(-1, 0, 1)))
  expect_equal(A.mat(matrix(c(-1, 1), 2), method="ASV"),
               matrix(c(.5, -.5, -.5, .5), 2))
})

test_that("ASV retains imputation and names and is invariant to equivalent dosage coding", {
  X <- asv_marker_data()
  X[1:2,1] <- NA_real_
  for(ploidy in list(NULL, 2, 3, 4, 8, 49)){
    dosage <- if(is.null(ploidy)) X else (X+1) * ploidy/2
    capture.output(ref <- asv_reference(dosage, .1, ploidy))
    output <- capture.output(result <- A.mat(dosage, min.MAF=.1, ploidy=ploidy,
                                             method="ASV", return.imputed=TRUE))
    expect_named(result, c("X", "A"))
    expect_equal(result$X, ref$X)
    expect_equal(result$A, ref$A, tolerance=1e-12)
    expect_identical(dimnames(result$X), dimnames(dosage))
    expect_identical(dimnames(result$A), list(rownames(X), rownames(X)))
    expect_length(grep("Imputing", output), 1L)
    capture.output(base <- A.mat(X, min.MAF=.1, method="ASV"))
    expect_equal(result$A, base, tolerance=1e-12)
  }
  X <- asv_marker_data()
  expect_equal(A.mat(-X, method="ASV"), A.mat(X, method="ASV"), tolerance=1e-12)
  order <- c(6, 2, 4, 1, 5, 3)
  expect_equal(A.mat(X[order,,drop=FALSE], method="ASV"),
               A.mat(X, method="ASV")[order,order], tolerance=1e-12)
  expect_equal(A.mat(X, min.MAF=1/12, method="ASV"),
               A.mat(X[,-2,drop=FALSE], method="ASV"), tolerance=1e-12)
})

test_that("real polyploid data meets ASV scaling and the empirical variance reference", {
  data("DT_polyploid", package="enhancer", envir=environment())
  capture.output(numo <- enhancer::atcg1234(GT_polyploid[,1:100], ploidy=4,
                                           silent=TRUE))
  common <- intersect(DT_polyploid$Name, rownames(numo$M))
  X <- numo$M[common,,drop=FALSE]
  A <- A.mat(X, ploidy=4, method="ASV")
  expect_equal(A, asv_reference(X, ploidy=4)$A, tolerance=1e-12)
  expect_equal(sum(diag(A)), nrow(A)-1, tolerance=1e-12)
})

test_that("ASV rejects unsupported combinations and unusable inputs explicitly", {
  X <- asv_marker_data()
  expect_error(A.mat(X, method="ASV", shrink=TRUE), "not supported")
  for(method in list("bad", NA_character_, character(), c("ASV", "VanRaden"))){
    expect_error(A.mat(X, method=method), "arg")
  }
  expect_error(A.mat(matrix(0, 4, 2), method="ASV"), "No variable markers")
  expect_error(A.mat(X, min.MAF=.5, method="ASV"), "No variable markers")
  expect_error(A.mat(X, min.MAF=-1, method="ASV"), "min.MAF")
  expect_error(A.mat(X, ploidy=1, method="ASV"), "ploidy")
  X[,1] <- NA_real_
  expect_error(A.mat(X, method="ASV"), "no observed genotypes")
  X <- asv_marker_data()
  X[1,1] <- Inf
  expect_error(A.mat(X, method="ASV"), "finite genotypes")
  expect_error(sommer:::amat_asv(matrix(0, 3, 2)), "ASV normalization")
  expect_error(sommer:::amat_asv(matrix(numeric(), 2, 0)), "finite centered matrix")
  expect_error(sommer:::amat_asv(matrix(1, 1, 2)), "finite centered matrix")
  expect_error(sommer:::amat_asv(matrix(NA_real_, 3, 2)), "finite centered matrix")
  W <- matrix(c(-1, 0, 1), 3)
  for(scale in c(1e-200, 1e200)){
    expect_equal(sommer:::amat_asv(W * scale), tcrossprod(W), tolerance=1e-12)
  }
})
