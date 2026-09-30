oatsVcc <- function(){
  data(DT_yatesoats, package="enhancer", envir=environment())
  DT_yatesoats
}

vccFit <- function(..., henderson=TRUE){
  mmes(..., verbose=FALSE, tolParConvLL=1e-11, tolParConvNorm=1e-11, nIters=100,
       henderson=henderson)
}

test_that("returnParam exposes the variance-parameter table", {
  oats <- oatsVcc()
  p <- mmes(Y ~ V*N, random=~B + B:MP, rcov=~units, data=oats, returnParam=TRUE)
  expect_equal(p$vcParams$term, c("vsm(ism(B))", "vsm(ism(B:MP))", "vsm(ism(units))"))
  expect_equal(p$vcParams$kind, rep("scale", 3))
  expect_equal(p$vcParams$index, 1:3)
})

test_that("equality and ratio constraints match exact reformulations in both engines", {
  oats <- oatsVcc()
  ZB <- model.matrix(~B-1, oats)
  ZM <- model.matrix(~B:MP-1, oats)
  ZM <- ZM[, colSums(ZM) > 0]
  Zeq <- cbind(ZB, ZM)
  Zr <- cbind(sqrt(2) * ZB, ZM)
  for(h in c(TRUE, FALSE)){
    eq <- vccFit(Y ~ V*N, random=~B + B:MP, rcov=~units, data=oats, henderson=h,
                 vcc=data.frame(parameter=c(1, 2), group=1))
    ref <- vccFit(Y ~ V*N, random=~vsm(ism(Zeq)), rcov=~units, data=oats, henderson=h)
    expect_equal(unname(unlist(eq$covPar)), unname(unlist(ref$covPar))[c(1, 1, 2)],
                 tolerance=1e-6)
    expect_equal(tail(c(eq$llik), 1), tail(c(ref$llik), 1), tolerance=1e-8)
    expect_equal(eq$theta_se[1, 1], ref$theta_se[1, 1], tolerance=1e-5)
    expect_equal(eq$theta_se[1, 1], eq$theta_se[2, 2], tolerance=1e-10)
    expect_equal(eq$vcParams$group, c(1L, 1L, NA))

    ratio <- vccFit(Y ~ V*N, random=~B + B:MP, rcov=~units, data=oats, henderson=h,
                    vcc=data.frame(parameter=c("vsm(ism(B:MP)):sigma2", "vsm(ism(B)):sigma2"),
                                   group=1, scale=c(1, 2)))
    refR <- vccFit(Y ~ V*N, random=~vsm(ism(Zr)), rcov=~units, data=oats, henderson=h)
    expect_equal(ratio$covPar[[1]], 2 * ratio$covPar[[2]], tolerance=1e-10,
                 ignore_attr=TRUE)
    expect_equal(unname(unlist(ratio$covPar))[2:3], unname(unlist(refR$covPar)),
                 tolerance=1e-6)
    expect_equal(tail(c(ratio$llik), 1), tail(c(refR$llik), 1), tolerance=1e-8)
  }
})

test_that("equated correlations reach the constrained REML optimum", {
  set.seed(11)
  d <- data.frame(g=factor(rep(1:30, each=8)), t=factor(rep(1:4, 60)),
                  h=factor(rep(1:40, each=6)))
  R <- 0.5^abs(outer(1:4, 1:4, "-"))
  eg <- matrix(rnorm(30*4), 30) %*% chol(R)
  eh <- matrix(rnorm(40*4), 40) %*% chol(R)
  d$y <- rnorm(nrow(d)) + eg[cbind(as.integer(d$g), as.integer(d$t))] +
    eh[cbind(as.integer(d$h), as.integer(d$t))]
  constrained <- vccFit(y~t, random=~vsm(ar1m(t), ism(g)) + vsm(ar1m(t), ism(h)),
                        data=d, vcc=data.frame(parameter=c(2, 4), group=1))
  rho <- constrained$covPar[[1]][2]
  expect_equal(rho, constrained$covPar[[2]][2], ignore_attr=TRUE, tolerance=1e-10)
  profile <- function(r){
    fit <- vccFit(y~t, random=~vsm(ar1m(t, rho=r, fixed=TRUE), ism(g)) +
                    vsm(ar1m(t, rho=r, fixed=TRUE), ism(h)), data=d)
    tail(c(fit$llik), 1)
  }
  best <- tail(c(constrained$llik), 1)
  step <- 0.3 * (1 - abs(rho))
  expect_equal(profile(rho), best, tolerance=1e-6)
  expect_lt(profile(rho + step), best)
  expect_lt(profile(rho - step), best)
})

test_that("vcc validates groups and propagates fixed parameters", {
  oats <- oatsVcc()
  expect_error(mmes(Y ~ V*N, random=~B + vsm(ar1m(N), ism(B)), rcov=~units, data=oats,
                    verbose=FALSE, vcc=data.frame(parameter=c(1, 3), group=1)),
               "different kinds")
  expect_error(mmes(Y ~ V*N, random=~B, rcov=~units, data=oats, verbose=FALSE,
                    vcc=data.frame(parameter=c(1, 9), group=1)), "Unknown parameters")
  expect_error(mmes(Y ~ V*N, random=~B, rcov=~units, data=oats, verbose=FALSE,
                    vcc=data.frame(parameter=1, group=1)), "single parameter")
  fx <- vccFit(Y ~ V*N, random=~vsm(ism(B), sigma2=50, fixedSigma2=TRUE) + B:MP,
               rcov=~units, data=oats, vcc=data.frame(parameter=c(1, 2), group=1))
  expect_equal(unname(unlist(fx$covPar))[1:2], c(50, 50), tolerance=1e-8)
})
