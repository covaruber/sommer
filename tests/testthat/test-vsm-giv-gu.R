giv_fixture <- function(){
  set.seed(3)
  nInd <- 30
  M <- matrix(sample(c(-1, 0, 1), nInd * 120, replace=TRUE), nrow=nInd)
  G <- tcrossprod(scale(M)) / ncol(M) + diag(1e-2, nInd)
  ids <- paste0("id", 1:nInd)
  Ginv <- solve(G)
  Ginv <- (Ginv + t(Ginv)) / 2
  dimnames(Ginv) <- list(ids, ids)
  attr(Ginv, "inverse") <- TRUE
  DT <- data.frame(id=factor(rep(ids, 3), levels=ids))
  u <- as.numeric(t(chol(G)) %*% rnorm(nInd))
  DT$y <- 10 + 2 * u[as.integer(DT$id)] + rnorm(nrow(DT))
  DT$y2 <- 0.5 * DT$y + rnorm(nrow(DT))
  lo <- which(lower.tri(Ginv, diag=TRUE), arr.ind=TRUE)
  lo <- lo[order(lo[, 1], lo[, 2]), ]
  list(Ginv=Ginv, data=DT, ids=ids,
       lower=cbind(Row=lo[, 1], Column=lo[, 2], Ainverse=Ginv[lo]))
}

as_giv <- function(x, rowNames, df=FALSE){
  if(df) x <- as.data.frame(x)
  attr(x, "rowNames") <- rowNames
  x
}

fit_giv <- function(Gu, data, ...){
  suppressMessages(mmes(y~1, random=~vsm(ism(id), Gu=Gu), rcov=~units,
                        data=data, verbose=FALSE, dateWarning=FALSE, ...))
}

test_that("3-column Gu (lower, upper, both, data.frame, shuffled) matches full Gu", {
  fx <- giv_fixture()
  lo <- fx$lower
  up <- lo[, c(2, 1, 3)]
  both <- rbind(lo, up[up[, 1] != up[, 2], ])
  set.seed(9)
  cases <- list(
    lower=as_giv(lo, fx$ids),
    upper=as_giv(up, fx$ids),
    both=as_giv(both, fx$ids),
    df=as_giv(lo, fx$ids, df=TRUE),
    shuffled=as_giv(lo[sample(nrow(lo)), ], fx$ids)
  )
  ref <- fit_giv(fx$Ginv, fx$data, nIters=15)
  for(nm in names(cases)){
    fit <- fit_giv(cases[[nm]], fx$data, nIters=15)
    expect_equal(unlist(fit$theta), unlist(ref$theta), tolerance=1e-10, info=nm)
    expect_equal(fit$u, ref$u, tolerance=1e-10, info=nm)
  }
})

test_that("3-column Gu keeps extra attributes and unobserved levels", {
  fx <- giv_fixture()
  giv <- as_giv(fx$lower, fx$ids)
  attr(giv, "inbreeding") <- setNames(rep(0, length(fx$ids)), fx$ids)
  attr(giv, "logdet") <- -1
  attr(giv, "geneticGroups") <- c(0, 0)
  DTsub <- droplevels(fx$data[fx$data$id != "id30", ])
  ref <- fit_giv(fx$Ginv, DTsub, nIters=15)
  fit <- fit_giv(giv, DTsub, nIters=15)
  expect_equal(unlist(fit$theta), unlist(ref$theta), tolerance=1e-10)
  expect_true("id30" %in% rownames(fit$u))
})

test_that("3-column Gu works with multi-trait, direct engine and covm", {
  fx <- giv_fixture()
  giv <- as_giv(fx$lower, fx$ids)
  DL <- rbind(data.frame(id=fx$data$id, trait="y", v=fx$data$y),
              data.frame(id=fx$data$id, trait="y2", v=fx$data$y2))
  fm <- function(Gu) suppressMessages(mmes(
    v~trait, random=~vsm(usm(trait), ism(id), Gu=Gu),
    rcov=~vsm(dsm(trait), ism(units)), data=DL, nIters=10,
    verbose=FALSE, dateWarning=FALSE))
  expect_equal(unlist(fm(giv)$theta), unlist(fm(fx$Ginv)$theta), tolerance=1e-10)

  ref <- fit_giv(fx$Ginv, fx$data, nIters=10, henderson=FALSE)
  fit <- fit_giv(giv, fx$data, nIters=10, henderson=FALSE)
  expect_equal(unlist(fit$theta), unlist(ref$theta), tolerance=1e-10)

  DT <- fx$data
  DT$id2 <- DT$id
  fc <- function(Gu) suppressMessages(mmes(
    y~1, random=~covm(vsm(ism(id), Gu=Gu), vsm(ism(id2), Gu=Gu)),
    rcov=~units, data=DT, nIters=5, verbose=FALSE, dateWarning=FALSE))
  expect_equal(unlist(fc(giv)$theta), unlist(fc(fx$Ginv)$theta), tolerance=1e-10)
})

test_that("3-column Gu input is validated", {
  fx <- giv_fixture()
  lo <- fx$lower
  vs <- function(Gu) with(fx$data, vsm(ism(id), Gu=Gu))

  expect_error(vs(lo), "rowNames")
  noDiag <- lo[!(lo[, 1] == 5 & lo[, 2] == 5), ]
  expect_error(vs(as_giv(noDiag, fx$ids)), "missing diagonal.*id5")
  both <- rbind(lo, lo[lo[, 1] != lo[, 2], c(2, 1, 3)][1, , drop=FALSE])
  both[nrow(both), 3] <- both[nrow(both), 3] + 1
  expect_error(vs(as_giv(both, fx$ids)), "conflicting")
  bad <- lo; bad[1, 1] <- length(fx$ids) + 1
  expect_error(vs(as_giv(bad, fx$ids)), "between 1 and")
  bad <- lo; bad[1, 3] <- NA
  expect_error(vs(as_giv(bad, fx$ids)), "finite")
  expect_error(vs(as_giv(lo, rep("a", length(fx$ids)))), "unique")
  notInv <- as_giv(lo, fx$ids); attr(notInv, "INVERSE") <- FALSE
  expect_error(vs(notInv), "precision")
})

test_that("3x3 full Gu is not mistaken for 3-column input", {
  Gu <- diag(3) + 0.1
  dimnames(Gu) <- list(letters[1:3], letters[1:3])
  expect_false(sommer:::.is_giv(Gu))
  expect_true(sommer:::.is_giv(as_giv(Gu, letters[1:3])))
})

test_that("symmetric Gu is stored as one triangle and fits identically", {
  fx <- giv_fixture()
  v <- with(fx$data, vsm(ism(id), Gu=fx$Ginv))
  expect_s4_class(v$Gu, "dsCMatrix")
  expect_true(isTRUE(attr(v$Gu, "inverse")))
  gen <- methods::as(Matrix::Matrix(fx$Ginv, sparse=TRUE), "generalMatrix")
  attr(gen, "inverse") <- TRUE
  for(s in c("ldlt", "cholmod")){
    a <- fit_giv(fx$Ginv, fx$data, nIters=10, solver=s)
    b <- fit_giv(gen, fx$data, nIters=10, solver=s)
    expect_equal(unlist(a$theta), unlist(b$theta), tolerance=1e-12, info=s)
  }
})
