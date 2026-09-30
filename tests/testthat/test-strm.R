strmData <- function(){
  set.seed(3)
  nid <- 50
  ids <- paste0("i", seq_len(nid))
  A <- diag(nid)
  for(k in 1:30){
    a <- sample(nid, 2)
    A[a[1], a[2]] <- A[a[2], a[1]] <- 0.25
  }
  A <- A + diag(0.3, nid)
  dimnames(A) <- list(ids, ids)
  Ai <- solve(A)
  attr(Ai, "inverse") <- TRUE
  n <- 360
  d <- data.frame(id=factor(sample(ids, n, TRUE), levels=ids),
                  dam=factor(sample(ids[1:25], n, TRUE), levels=ids),
                  pe=factor(sample(ids[26:50], n, TRUE), levels=ids),
                  env=factor(sample(c("E1", "E2"), n, TRUE)))
  S <- matrix(c(4, 1.2, 0.5, 1.2, 2, 0.3, 0.5, 0.3, 1.5), 3)
  U <- t(chol(A)) %*% matrix(rnorm(nid*3), nid) %*% chol(S)
  rownames(U) <- ids
  d$y <- 10 + U[as.character(d$id), 1] + U[as.character(d$dam), 2] +
    U[as.character(d$pe), 3] + rnorm(n, sd=1.5)
  list(d=d, Ai=Ai)
}

strmFit <- function(random, d){
  mmes(y~1, random=random, data=d, verbose=FALSE, tolParConvLL=1e-10,
       tolParConvNorm=1e-10, nIters=100)
}

test_that("strm() with two terms reproduces covm()", {
  s <- strmData(); d <- s$d; Ai <- s$Ai
  mc <- strmFit(~covm(vsm(ism(id), Gu=Ai), vsm(ism(dam), Gu=Ai)), d)
  ms <- strmFit(~strm(vsm(ism(id)), vsm(ism(dam)), Gu=Ai), d)
  expect_equal(ms$theta[[1]], mc$theta[[1]], tolerance=1e-5, ignore_attr=TRUE)
  expect_equal(tail(ms$llik[1,], 1), tail(mc$llik[1,], 1), tolerance=1e-8)
})

test_that("strm() with a diagonal term covariance equals independent terms", {
  s <- strmData(); d <- s$d; Ai <- s$Ai
  md <- strmFit(~strm(vsm(ism(id)), vsm(ism(dam)), cov=dsm, Gu=Ai), d)
  mi <- strmFit(~vsm(ism(id), Gu=Ai) + vsm(ism(dam), Gu=Ai), d)
  expect_equal(unname(diag(md$theta[[1]])),
               unname(c(mi$theta[[1]], mi$theta[[2]])), tolerance=1e-5)
  expect_equal(md$theta[[1]][1, 2], 0)
})

test_that("strm() fits three correlated terms with labels", {
  s <- strmData(); d <- s$d; Ai <- s$Ai
  m <- strmFit(~strm(dir=vsm(ism(id)), mat=vsm(ism(dam)), pe=vsm(ism(pe)), Gu=Ai), d)
  expect_equal(dim(m$theta[[1]]), c(3L, 3L))
  expect_equal(colnames(m$uList[[1]]), c("dir", "mat", "pe"))
  expect_equal(nrow(m$uList[[1]]), 50L)
  native <- covparams_mmes_se(m, 1L)
  expect_true(all(c("variance[dir]", "covariance[mat,dir]") %in% native$parameter))
  expect_true(all(is.finite(native$StdError)))
})

test_that("strm() shares inner covariance factors across terms", {
  s <- strmData(); d <- s$d; Ai <- s$Ai
  m <- strmFit(~strm(vsm(dsm(env), ism(id)), vsm(dsm(env), ism(dam)), Gu=Ai), d)
  expect_equal(dim(m$theta[[1]]), c(4L, 4L))
  expect_equal(colnames(m$uList[[1]]), c("t1:E1", "t1:E2", "t2:E1", "t2:E2"))
  expect_error(strm(vsm(dsm(d$env), ism(d$id)), vsm(ism(d$dam))), "same covariance factors")
  expect_error(strm(vsm(ism(d$id))), "at least two")
})
