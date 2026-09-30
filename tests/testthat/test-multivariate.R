multivariateFixture <- function(n=40L){
  set.seed(17)
  data <- data.frame(id=factor(rep(seq_len(n), each=3L)))
  genetic <- matrix(rnorm(n * 2L), n) %*% chol(matrix(c(1, 0.35, 0.35, 1.4), 2))
  data$trait_a <- genetic[data$id, 1L] + rnorm(nrow(data), sd=0.7)
  data$trait_b <- genetic[data$id, 2L] + rnorm(nrow(data), sd=1)
  data
}

test_that("stackTraits builds stable trait and record coordinates", {
  wide <- multivariateFixture(8L)
  long <- stackTraits(wide, c("trait_a", "trait_b"), keep="id")

  expect_equal(nrow(long), 2L * nrow(wide))
  expect_equal(levels(long$trait), c("trait_a", "trait_b"))
  expect_equal(long$record[seq_len(nrow(wide))],
               long$record[nrow(wide) + seq_len(nrow(wide))])
  expect_equal(long$value[seq_len(nrow(wide))], wide$trait_a)
  expect_equal(long$value[nrow(wide) + seq_len(nrow(wide))], wide$trait_b)
})

test_that("multi-trait residual pairing is explicit or warns when implicit", {
  wide <- multivariateFixture()
  long <- stackTraits(wide, c("trait_a", "trait_b"), keep="id")
  fitExplicit <- mmes(
    value~trait,
    random=~vsm(usm(trait), ism(id)),
    rcov=~vsm(usm(trait), ism(record)),
    data=long, verbose=FALSE, nIters=20
  )
  fitImplicit <- mmes(
    value~trait,
    random=~vsm(usm(trait), ism(id)),
    rcov=~vsm(usm(trait), ism(units)),
    data=long, verbose=FALSE, nIters=20
  )
  expect_equal(unname(unlist(fitExplicit$covPar)),
               unname(unlist(fitImplicit$covPar)), tolerance=1e-8)

  covariance <- covmatrix_mmes(fitExplicit, 1L)
  expect_true(isSymmetric(covariance$covariance))
  expect_equal(covariance$correlation[1L, 2L],
               covmatrix_mmes(fitExplicit, 1L, se=FALSE)$correlation[1L, 2L])
  nativeSE <- covparams_mmes_se(fitExplicit, 1L)$StdError
  expect_true(all(is.finite(covariance$covariance.se)))
  expect_gt(max(nativeSE), 0)
  expect_equal(unname(diag(covariance$correlation.se)), c(0, 0))
})

test_that("heterogeneous corgm and mixed-family PQL use trait dispersions", {
  set.seed(12)
  n <- 50L
  wide <- data.frame(id=factor(rep(seq_len(n), each=3L)))
  genetic <- matrix(rnorm(n * 2L), n) %*% chol(matrix(c(1, 0.4, 0.4, 0.8), 2))
  wide$disease <- rbinom(nrow(wide), 1, plogis(-0.5 + genetic[wide$id, 1L]))
  wide$yield <- 10 + genetic[wide$id, 2L] + rnorm(nrow(wide))
  long <- stackTraits(wide, c("disease", "yield"), keep="id")
  descriptor <- corgm(long$trait, variance="heterogeneous")$covFactor
  expect_equal(length(descriptor$par), length(descriptor$free))
  expect_equal(length(descriptor$par), length(descriptor$par_names))

  fit <- mmes(
    value~trait,
    random=~vsm(usm(trait), ism(id)),
    rcov=~vsm(corgm(trait, variance="heterogeneous"), ism(record)),
    family=familym("trait", disease=binomial(), yield=gaussian()),
    data=long, verbose=FALSE, nIters=20,
    pqlControl=list(maxit=6L)
  )
  expect_s3_class(fit, "mmes.glmm")
  expect_equal(unname(fit$dispersion["disease"]), 1)
  expect_true(is.finite(fit$dispersion["yield"]))
  expect_true(all(is.finite(covmatrix_mmes(fit, 1L, se=FALSE)$correlation)))
  expect_match(
    tryCatch({
      mmes(value~trait, random=~id,
           rcov=~vsm(usm(trait), ism(record)),
           family=familym("trait", disease=binomial(), yield=gaussian()),
           data=long, verbose=FALSE)
      ""
    }, error=conditionMessage),
    "trait-specific variances"
  )
  long$trait <- factor(long$trait, levels=c("yield", "disease"))
  expect_error(
    mmes(value~trait, random=~vsm(usm(trait), ism(id)),
         rcov=~vsm(corgm(trait, variance="heterogeneous"), ism(record)),
         family=familym("trait", disease=binomial(), yield=gaussian()),
         data=long, verbose=FALSE),
    "first level of 'trait'"
  )
})

test_that("weighted complete-block Kronecker residual agrees across engines", {
  wide <- multivariateFixture(30L)
  long <- stackTraits(wide, c("trait_a", "trait_b"), keep="id")
  weights <- runif(nrow(long), 1, 3)
  fit <- function(henderson, W=diag(weights)){
    args <- list(
      fixed=value~trait,
      random=~vsm(usm(trait), ism(id)),
      rcov=~vsm(usm(trait), ism(record)),
      data=long, verbose=FALSE, nIters=25,
      henderson=henderson, tolParConvLL=1e-10, tolParConvNorm=1e-9
    )
    if(!missing(W) && !is.null(W)) args$W <- W
    do.call(mmes, args)
  }
  henderson <- fit(TRUE)
  direct <- fit(FALSE)
  expect_equal(unlist(henderson$covPar), unlist(direct$covPar), tolerance=1e-5)
  weightedIdentity <- fit(TRUE, diag(length(weights)))
  unweighted <- fit(TRUE, NULL)
  expect_equal(unlist(weightedIdentity$covPar), unlist(unweighted$covPar),
               tolerance=1e-7)
})