metricSites <- function(ns=40, reps=3){
  set.seed(1)
  co <- cbind(x=runif(ns, 0, 10), y=runif(ns, 0, 10))
  d <- data.frame(x=rep(co[,1], each=reps), y=rep(co[,2], each=reps), one=factor(1))
  K <- exp(-as.matrix(dist(co))/3)
  u <- drop(t(chol(K)) %*% rnorm(ns))
  d$yv <- 10 + 2*u[rep(seq_len(ns), each=reps)] + rnorm(nrow(d))
  list(coords=co, data=d)
}

test_that("metricm analytic derivatives match numerical derivatives", {
  co <- metricSites()$coords
  for(model in c("exponential", "power", "gaussian", "spherical", "circular")){
    for(anisotropy in c("none", "product", "geometric")){
      if(anisotropy == "product" && model %in% c("spherical", "circular")) next
      cf <- suppressWarnings(metricm(co, model=model, anisotropy=anisotropy))$covFactor
      p <- cf$par
      K <- cf$evaluator$fun(p)
      expect_equal(unname(diag(K)), rep(1, nrow(K)))
      for(k in seq_along(p)){
        h <- 1e-6
        up <- p; dn <- p
        up[k] <- p[k] + h; dn[k] <- p[k] - h
        numeric <- (cf$evaluator$fun(up) - cf$evaluator$fun(dn)) / (2*h)
        expect_lt(max(abs(numeric - cf$derivative$fun(p, k))), 1e-6)
      }
    }
  }
})

test_that("exponential, power and Matern(0.5) fits are equivalent", {
  s <- metricSites()
  ctl <- list(verbose=FALSE, tolParConvLL=1e-10, tolParConvNorm=1e-10, nIters=100)
  mExp <- do.call(mmes, c(list(yv~1, random=~vsm(metricm(cbind(x, y)), ism(one)),
                               data=s$data), ctl))
  mMat <- do.call(mmes, c(list(yv~1, random=~vsm(maternm(cbind(x, y), nu=0.5,
                               fixed=c(FALSE, TRUE)), ism(one)), data=s$data), ctl))
  mPow <- do.call(mmes, c(list(yv~1, random=~vsm(metricm(cbind(x, y), model="power"),
                               ism(one)), data=s$data), ctl))
  rangeExp <- mExp$covPar[[1]][2]
  expect_equal(rangeExp, mMat$covPar[[1]][2], tolerance=1e-5)
  expect_equal(rangeExp, -1/log(mPow$covPar[[1]][2]), tolerance=1e-5)
  expect_equal(tail(mExp$llik[1,], 1), tail(mPow$llik[1,], 1), tolerance=1e-6)
  native <- covparams_mmes_se(mPow, 1L)
  expect_equal(native$parameter, c("variance", "rho"))
  expect_true(all(is.finite(native$StdError)))
})

test_that("metricm validates dimensions and model combinations", {
  co3 <- cbind(runif(5), runif(5), runif(5))
  expect_error(metricm(co3, model="circular"), "two dimensions")
  expect_error(metricm(runif(5), anisotropy="product"), "two-dimensional")
  expect_error(metricm(cbind(runif(5), runif(5)), model="gaussian", metric="manhattan"),
               "manhattan")
  expect_error(metricm(cbind(runif(5), runif(5)), model="spherical", anisotropy="product"),
               "Product anisotropy")
  cf <- metricm(cbind(runif(6), runif(6)), anisotropy="geometric")$covFactor
  expect_equal(cf$par_names, c("range", "angle", "ratio"))
  cfm <- maternm(cbind(runif(6), runif(6)), anisotropy="geometric")$covFactor
  expect_equal(cfm$par_names, c("range", "nu", "angle", "ratio"))
})
