oatsWaldData <- function(drop=0L){
  data(DT_yatesoats, package="enhancer", envir=environment())
  d <- DT_yatesoats
  if(drop > 0L){
    set.seed(2)
    d <- d[-sample(nrow(d), drop), ]
  }
  d
}

tightFit <- function(...){
  mmes(..., verbose=FALSE, tolParConvLL=1e-12, tolParConvNorm=1e-12, nIters=100)
}

test_that("theta_se is the inverse REML information on the reported scale", {
  oats <- oatsWaldData()
  m0 <- tightFit(Y ~ V * N, rcov=~units, data=oats)
  s2 <- as.numeric(m0$covPar[[1]])
  expect_equal(as.numeric(m0$theta_se), 2 * s2^2 / (nrow(oats) - 12), tolerance=1e-6)

  mh <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=oats)
  md <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=oats, henderson=FALSE)
  expect_equal(mh$theta_se, md$theta_se, tolerance=1e-4)
})

test_that("incremental Wald without random effects reproduces lm anova", {
  oats <- oatsWaldData()
  m0 <- tightFit(Y ~ V * N, rcov=~units, data=oats)
  a <- anova(lm(Y ~ V * N, data=oats))
  w <- wald_mmes(m0, denDF="residual")
  expect_s3_class(w, "wald.mmes")
  expect_equal(w$F.value[-1], a$`F value`[1:3], tolerance=1e-8)
  expect_equal(w$denDF[-1], rep(a$Df[4], 3))
  expect_equal(w$p.value[-1], a$`Pr(>F)`[1:3], tolerance=1e-6)
  expect_equal(anova(m0, denDF="residual")$Wald, w$Wald)
})

test_that("Satterthwaite and Kenward-Roger df reproduce the split-plot strata", {
  skip_on_cran()
  oats <- oatsWaldData()
  m <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=oats)
  for(method in c("satterthwaite", "kr")){
    w <- wald_mmes(m, denDF=method)
    expect_equal(w$denDF, c(5, 10, 45, 45), tolerance=1e-4)
    expect_equal(w[c("V", "N", "V:N"), "F.value"], c(1.4853, 37.6857, 0.3028),
                 tolerance=1e-4)
  }
})

test_that("Kenward-Roger df and F on unbalanced data", {
  skip_on_cran()
  ob <- oatsWaldData(9L)
  m <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=ob)
  kr <- wald_mmes(m, denDF="kr")
  expect_equal(kr["V:N", "denDF"], 36.885, tolerance=0.01)
  expect_equal(kr["V:N", "F.value"], 0.2792, tolerance=0.01)
})

test_that("conditional Wald respects marginality", {
  ob <- oatsWaldData(9L)
  m <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=ob)
  inc <- wald_mmes(m)
  con <- wald_mmes(m, ssType="conditional")
  # The highest-order term and N (added after V) coincide; V is adjusted for N.
  expect_equal(con["V:N", "Wald"], inc["V:N", "Wald"], tolerance=1e-8)
  expect_equal(con["N", "Wald"], inc["N", "Wald"], tolerance=1e-8)
  expect_false(isTRUE(all.equal(con["V", "Wald"], inc["V", "Wald"])))
  m2 <- tightFit(Y ~ N * V, random=~B + B:MP, rcov=~units, data=ob)
  expect_equal(con["V", "Wald"], wald_mmes(m2)["V", "Wald"], tolerance=1e-6)
})

test_that("pairwise comparisons can use small-sample df for fixed contrasts", {
  skip_on_cran()
  oats <- oatsWaldData()
  m <- tightFit(Y ~ V * N, random=~B + B:MP, rcov=~units, data=oats)
  Dt <- m$Dtable
  Dt[c(1, 2), "average"] <- TRUE
  Dt[c(3, 4), "include"] <- TRUE
  Dt[4, "average"] <- TRUE
  p <- predict(m, Dtable=Dt, D="N", pairwise=TRUE, df="satterthwaite")
  expect_equal(p$pairwise$df, rep(45, 6), tolerance=1e-4)
  pr <- predict(m, D="B", pairwise=TRUE, df="kr")
  expect_true(all(is.infinite(pr$pairwise$df)))
})

test_that("Wald tests are available for PQL fits but small-sample df are not", {
  set.seed(4)
  d <- data.frame(x=rnorm(120), g=factor(rep(1:12, each=10)))
  d$y <- rpois(120, exp(0.3 + 0.4 * d$x))
  fit <- mmes(y~x, random=~g, data=d, family=poisson(), verbose=FALSE)
  expect_true(all(is.finite(wald_mmes(fit)$p.value)))
  expect_error(wald_mmes(fit, denDF="kr"), "PQL")
})
