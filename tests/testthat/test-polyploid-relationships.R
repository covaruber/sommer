polyploid_relationship_data <- function(ploidy=4){
  X <- matrix(c(0, 1, 2, 3, 4, 2,
                0, 0, 0, 0, 0, 1,
                2, 2, 2, 2, 2, 2,
                4, 4, 4, 4, 4, 4,
                .4, 1.2, 2.6, 3.8, 2, 0), nrow=6) * ploidy/4
  dimnames(X) <- list(paste0("id", 1:6), paste0("m", 1:5))
  X
}

polyploid_dominance_reference <- function(X, ploidy, min.MAF=0){
  for(j in seq_len(ncol(X))){
    X[is.na(X[,j]),j] <- mean(X[,j], na.rm=TRUE)
  }
  p <- colMeans(X) / ploidy
  keep <- pmin(p, 1-p) > min.MAF & apply(X, 2, var) > 0
  Y <- X[,keep,drop=FALSE]
  p <- p[keep]
  # Start from heterozygous allele-pair counts, then remove the additive projection.
  W <- Y * (ploidy-Y) -
    sweep(Y, 2, (ploidy-1)*(1-2*p), "*")
  W <- sweep(W, 2, ploidy*(ploidy-1)*p^2, "-")
  D <- tcrossprod(W) / sum(2*ploidy*(ploidy-1)*(p*(1-p))^2)
  colnames(D) <- rownames(D) <- rownames(X)
  list(D=D, X=X, keep=keep)
}

test_that("polyploid dominance uses the orthogonal digenic kernel", {
  for(ploidy in c(2, 3, 4, 6, 8)){
    X <- polyploid_relationship_data(ploidy)
    for(maf in c(0, .1, .3)){
      ref <- polyploid_dominance_reference(X, ploidy, maf)
      D <- D.mat(X, min.MAF=maf, ploidy=ploidy)
      expect_equal(D, ref$D, tolerance=1e-12)
      expect_equal(D, t(D), tolerance=1e-12)
      expect_gte(min(eigen(D, symmetric=TRUE, only.values=TRUE)$values), -1e-12)
      expect_equal(D.mat(ploidy-X, min.MAF=maf, ploidy=ploidy), D, tolerance=1e-12)
      order <- c(6, 2, 4, 1, 5, 3)
      expect_equal(D.mat(X[order,,drop=FALSE], min.MAF=maf, ploidy=ploidy),
                   D[order,order], tolerance=1e-12)
    }
    for(p in c(.1, .3, .5, .8)){
      dosage <- 0:ploidy
      weights <- dbinom(dosage, ploidy, p)
      t <- dosage-ploidy*p
      w <- -t^2+(1-2*p)*t+ploidy*p*(1-p)
      expect_equal(sum(weights*w), 0, tolerance=1e-12)
      expect_equal(sum(weights*w*t), 0, tolerance=1e-12)
      expect_equal(sum(weights*w^2), 2*ploidy*(ploidy-1)*(p*(1-p))^2,
                   tolerance=1e-12)
    }
  }
  X <- polyploid_relationship_data()
  expect_equal(D.mat(X[,1,drop=FALSE], ploidy=4),
               polyploid_dominance_reference(X[,1,drop=FALSE], 4)$D, tolerance=1e-12)
  expect_equal(D.mat(X, min.MAF=1/24, ploidy=4),
               D.mat(X[,-2,drop=FALSE], ploidy=4), tolerance=1e-12)
})

test_that("explicit diploid dosages preserve both original dominance methods", {
  X <- polyploid_relationship_data(2)
  for(nishio in c(TRUE, FALSE)){
    expect_equal(D.mat(X, nishio, ploidy=2), D.mat(X-1, nishio), tolerance=1e-12)
    for(type in c("A#A", "A#D", "D#D")){
      expect_equal(E.mat(X, nishio, type, ploidy=2),
                   E.mat(X-1, nishio, type), tolerance=1e-12)
    }
  }
  expect_identical(D.mat(X-1), D.mat(X-1, ploidy=NULL))
  expect_identical(E.mat(X-1), E.mat(X-1, ploidy=NULL))
})

test_that("polyploid epistasis shares component filters and mean imputation", {
  for(ploidy in c(3, 4, 6, 8)){
    X <- polyploid_relationship_data(ploidy)
    X[1:2,1] <- NA_real_
    capture.output(result <- D.mat(X, min.MAF=.1, return.imputed=TRUE, ploidy=ploidy))
    ref <- polyploid_dominance_reference(X, ploidy, .1)
    expect_named(result, c("X", "D"))
    expect_equal(result$X, ref$X)
    expect_equal(result$D, ref$D, tolerance=1e-12)
    expect_identical(dimnames(result$X), dimnames(X))
    capture.output(A <- A.mat(X, min.MAF=.1, ploidy=ploidy))
    expected <- list("A#A"=A*A, "A#D"=A*result$D, "D#D"=result$D^2)
    for(type in names(expected)){
      output <- capture.output(E <- E.mat(X, type=type, min.MAF=.1, ploidy=ploidy))
      expect_equal(E, expected[[type]], tolerance=1e-12)
      expect_length(grep("Imputing", output), 1L)
      expect_identical(dimnames(E), dimnames(A))
    }
    capture.output(E <- E.mat(X, nishio=FALSE, type="A#A", min.MAF=.1, ploidy=ploidy))
    expect_equal(E, A*A, tolerance=1e-12)
  }
})

test_that("the supplied polyploid data supports dominance and all epistasis types", {
  data("DT_polyploid", package="enhancer", envir=environment())
  capture.output(numo <- enhancer::atcg1234(GT_polyploid[,1:100], ploidy=4,
                                           silent=TRUE))
  X <- numo$M[intersect(DT_polyploid$Name, rownames(numo$M)),,drop=FALSE]
  D <- D.mat(X, ploidy=4)
  expect_equal(D, polyploid_dominance_reference(X, 4)$D, tolerance=1e-12)
  A <- A.mat(X, ploidy=4, min.MAF=.02)
  D <- D.mat(X, ploidy=4, min.MAF=.02)
  for(type in c("A#A", "A#D", "D#D")){
    expected <- switch(type, "A#A"=A*A, "A#D"=A*D, "D#D"=D*D)
    expect_equal(E.mat(X, type=type, ploidy=4), expected, tolerance=1e-12)
  }
})

test_that("hybrid relationship algebra already accepts polyploid additive matrices", {
  X <- polyploid_relationship_data()
  A <- A.mat(X, ploidy=4) + diag(.2, nrow(X))
  for(j in list(c(2, 4, 5), 3L, seq_len(nrow(X)))){
    G <- if(length(j) == 1L) A[j,j,drop=FALSE] else
      A.mat(X[j,,drop=FALSE], ploidy=4) + diag(.1, length(j))
    for(weights in list(c(1, 1), c(1.4, .8))){
      precision <- solve(A)
      precision[j,j] <- precision[j,j,drop=FALSE] +
        weights[1]*solve(G) - weights[2]*solve(A[j,j,drop=FALSE])
      H <- H.mat(A, G, tau=weights[1], omega=weights[2])
      expect_equal(H[rownames(A),rownames(A)], solve(precision), tolerance=1e-10)
      expect_equal(H, t(H), tolerance=1e-12)
      expect_identical(rownames(H), rev(unique(c(rownames(G), rownames(A)))))
    }
  }
})

test_that("polyploid dominance rejects invalid dosages and unspecified models", {
  X <- polyploid_relationship_data()
  for(fun in list(D.mat, E.mat)){
    for(ploidy in list(1, 2.5, NA_real_, Inf, "4", numeric(), c(2,4))){
      expect_error(fun(X, ploidy=ploidy), "ploidy")
    }
    expect_error(fun(matrix(2, 4, 2), ploidy=4), "No variable markers")
    expect_error(fun(X, min.MAF=.5, ploidy=4), "No variable markers")
    for(maf in c(-.1, .6, NA_real_, Inf)){
      expect_error(fun(X, min.MAF=maf, ploidy=4), "min.MAF")
    }
    for(bad in c(-1, 5, Inf)){
      invalid <- X
      invalid[1,1] <- bad
      expect_error(fun(invalid, ploidy=4), "finite genotypes")
    }
    invalid <- X
    invalid[,1] <- NA_real_
    expect_error(fun(invalid, ploidy=4), "no observed genotypes")
  }
  for(nishio in list(NA, NULL, logical(), c(TRUE,FALSE), 1, "TRUE")){
    expect_error(D.mat(X, nishio=nishio, ploidy=4), "nishio must")
    expect_error(E.mat(X, nishio=nishio, type="A#D", ploidy=4), "nishio must")
  }
  expect_error(D.mat(X, nishio=FALSE, ploidy=4), "not supported")
  for(type in c("A#D", "D#D")){
    expect_error(E.mat(X, nishio=FALSE, type=type, ploidy=4), "not supported")
    expect_error(sommer:::emat_dosage(X, FALSE, type, 0, 4), "not supported")
  }
  expect_error(sommer:::dmat_dosage(X, FALSE, 0, 4), "not supported")
  expect_error(sommer:::dmat_dosage(X, TRUE, 0, 1), "ploidy")
  expect_error(sommer:::emat_dosage(X, TRUE, "bad", 0, 4), "type must")
  expect_error(sommer:::dmat_dosage(X-1, TRUE, 0, 4), "finite genotypes")
})
