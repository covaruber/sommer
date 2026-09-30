markerData <- function(){
  set.seed(7)
  n <- 60; m <- 120
  M <- matrix(sample(c(-1, 0, 1), n*m, TRUE, prob=c(.3, .4, .3)), n)
  rownames(M) <- paste0("g", 1:n); colnames(M) <- paste0("snp", 1:m)
  M[, 5] <- 1
  lambda <- 0.05
  G <- (1 - lambda) * A.mat(M) + lambda * diag(n)
  Gi <- solve(G); attr(Gi, "inverse") <- TRUE
  d <- data.frame(id=factor(rep(rownames(M), each=2), levels=rownames(M)),
                  env=factor(rep(c("A", "B"), n)))
  d$y <- 5 + ((M + 1) %*% rnorm(m, 0, 0.3))[as.character(d$id), 1] + rnorm(nrow(d))
  list(M=M, Gi=Gi, d=d, lambda=lambda)
}

test_that("marker effects back-transform the GBLUP exactly", {
  s <- markerData(); d <- s$d; Gi <- s$Gi
  fit <- mmes(y~1, random=~vsm(ism(id), Gu=Gi), data=d, verbose=FALSE)
  me <- meffects_mmes(fit, 1, s$M, blend=s$lambda)
  tr <- sommer:::.amat_marker_transform(s$M, 0)
  u <- fit$uList[[1]][rownames(s$M), 1]
  Ginv <- as.matrix(Gi)[rownames(s$M), rownames(s$M)]
  expect_equal(unname(drop(tr$Z %*% me$effect[tr$keep])),
               unname(u - s$lambda * drop(Ginv %*% u)), tolerance=1e-8)
  expect_true(is.na(me$effect[5]))
  expect_equal(length(attr(me, "markersKept")), 119L)
  expect_equal(me$effect, meffects_mmes(fit, 1, s$M, blend=s$lambda, Gu=Gi)$effect)
  custom <- meffects_mmes(fit, 1, method="custom", Z=tr$Z, scale=1/tr$k, blend=s$lambda)
  expect_equal(custom$effect, me$effect[tr$keep], tolerance=1e-12)
})

test_that("marker-effect standard errors match a dense calculation", {
  s <- markerData(); d <- s$d; Gi <- s$Gi
  fit <- mmes(y~1, random=~vsm(ism(id), Gu=Gi), data=d, verbose=FALSE, computeCi=2)
  me <- meffects_mmes(fit, 1, s$M, blend=s$lambda, se=TRUE, chunk=37)
  tr <- sommer:::.amat_marker_transform(s$M, 0)
  lev <- rownames(fit$uList[[1]])
  Ginv <- as.matrix(Gi)[lev, lev]
  H <- (1 - s$lambda) / tr$k * Ginv %*% tr$Z[lev, ]
  cols <- fit$partitions[[1]][1, 1]:fit$partitions[[1]][1, 2]
  Cuu <- as.matrix(fit$Ci)[cols, cols]
  V <- crossprod(H, (as.numeric(fit$theta[[1]]) * solve(Ginv) - Cuu) %*% H)
  expect_equal(me$se[tr$keep], unname(sqrt(pmax(diag(V), 0))), tolerance=1e-6)
  expect_equal(me$p.value[tr$keep], 2 * pnorm(-abs(me$effect[tr$keep] / me$se[tr$keep])))
})

test_that("marker effects are returned per coordinate for structured terms", {
  s <- markerData(); d <- s$d; Gi <- s$Gi
  fit <- mmes(y~env, random=~vsm(dsm(env), ism(id), Gu=Gi), data=d, verbose=FALSE,
              nIters=15)
  me <- meffects_mmes(fit, 1, s$M, blend=s$lambda)
  expect_equal(unique(me$coordinate), c("A", "B"))
  expect_equal(nrow(me), 2L * ncol(s$M))
  expect_error(meffects_mmes(fit, 1, s$M[-1, ]), "missing from M")
  expect_error(meffects_mmes(fit, "nope", s$M), "Unknown random term")
})
