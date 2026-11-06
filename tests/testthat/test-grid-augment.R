gridAugmentData <- function(seed=11, drop=25, na=20){
  set.seed(seed)
  dat <- expand.grid(trial=factor(c("A", "B")), range=factor(1:8), row=factor(1:10))
  dat$geno <- factor(sample(sprintf("g%02d", 1:20), nrow(dat), replace=TRUE))
  dat$y <- 10 + as.numeric(dat$trial) + rnorm(nrow(dat))
  dat$y[sample(nrow(dat), na)] <- NA
  dat[-sample(nrow(dat), drop), , drop=FALSE]
}

gridAugmentFit <- function(dat, mode, ..., rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial),
                           representation=list(), solve=list()){
  mmes(y ~ trial, random=~geno, rcov=rcov, data=dat,
       control=mmes_control(representation=utils::modifyList(list(gridAugment=mode), representation),
                            solve=solve, convergence=list(maxIters=25),
                            diagnostics=list(verbose=FALSE, dateWarning=FALSE)), ...)
}

expectSameFit <- function(a, b, tolerance=1e-6){
  expect_equal(tail(as.numeric(b$llik), 1), tail(as.numeric(a$llik), 1), tolerance=tolerance)
  expect_equal(unlist(b$theta), unlist(a$theta), tolerance=tolerance)
  expect_equal(as.numeric(b$b), as.numeric(a$b), tolerance=tolerance)
  expect_equal(rownames(b$b), rownames(a$b))
  expect_equal(as.numeric(b$u), as.numeric(a$u), tolerance=tolerance)
  expect_equal(as.numeric(fitted(b)), as.numeric(fitted(a)), tolerance=tolerance)
  expect_equal(as.numeric(residuals(b)), as.numeric(residuals(a)), tolerance=tolerance)
  expect_equal(b$AIC, a$AIC, tolerance=tolerance)
  expect_equal(b$BIC, a$BIC, tolerance=tolerance)
  expect_equal(dim(b$W), dim(a$W))
}

test_that("grid augmentation of dsumm AR1 x AR1 residuals is exact for REML", {
  dat <- gridAugmentData()
  off <- gridAugmentFit(dat, "off")
  auto <- gridAugmentFit(dat, "auto")

  info <- auto$engineDiagnostics$gridAugmentation
  expect_true(info$active)
  expect_equal(info$blocks, 2L)
  expect_equal(info$observed, sum(!is.na(dat$y)))
  expect_equal(info$cells, 2L * 8L * 10L - sum(!is.na(dat$y)))
  expect_null(off$engineDiagnostics$gridAugmentation$active)

  expectSameFit(off, auto)
  expect_equal(nrow(auto$Ci), nrow(off$Ci))

  pOff <- predict(off, D="trial")$pvals
  pAuto <- predict(auto, D="trial")$pvals
  expect_equal(pAuto$predicted.value, pOff$predicted.value, tolerance=1e-6)
  expect_equal(pAuto$std.error, pOff$std.error, tolerance=1e-6)

  vOff <- summary(off)$varcomp
  vAuto <- summary(auto)$varcomp
  expect_equal(vAuto$VarComp, vOff$VarComp, tolerance=1e-6)
  expect_equal(vAuto$VarCompSE, vOff$VarCompSE, tolerance=1e-5)
})

test_that("grid augmentation keeps a single symbolic analysis and sparse grid blocks", {
  dat <- gridAugmentData(seed=12)
  auto <- gridAugmentFit(dat, "auto")
  d <- auto$engineDiagnostics
  expect_true(d$gridAugmentation$active)
  expect_equal(d$CsymbolicAnalyses, 1)
  expect_equal(d$gridSparseBlocks, 2)
  expect_gt(d$gridSparseBlockUses, 0)
  expect_equal(d$gridBandFallbacks, 0)
})

test_that("grid augmentation matches with the ldlt solver and a single-trial residual", {
  dat <- droplevels(gridAugmentData(seed=13)[gridAugmentData(seed=13)$trial == "A", ])
  rcov <- ~vsm(ar1m(range), ar1m(row), ism(units))
  fit <- function(mode) mmes(y ~ 1, random=~geno, rcov=rcov, data=dat,
    control=mmes_control(representation=list(gridAugment=mode), solve=list(method="ldlt"),
                         convergence=list(maxIters=25), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  off <- fit("off")
  auto <- fit("auto")
  expect_true(auto$engineDiagnostics$gridAugmentation$active)
  expect_equal(auto$engineDiagnostics$gridAugmentation$blocks, 1L)
  expectSameFit(off, auto)
})

test_that("banded sparse grid blocks match the dense grid path", {
  dat <- expand.grid(trial=factor(c("A", "B")), range=factor(1:8), row=factor(1:10))
  set.seed(14)
  dat$geno <- factor(sample(sprintf("g%02d", 1:20), nrow(dat), replace=TRUE))
  dat$y <- rnorm(nrow(dat))
  sparse <- gridAugmentFit(dat, "off")
  dense <- gridAugmentFit(dat, "off", representation=list(sparseFactors="off"))
  expect_gt(sparse$engineDiagnostics$gridSparseBlocks, 0)
  expect_equal(dense$engineDiagnostics$gridSparseBlocks, 0)
  expect_equal(sparse$engineDiagnostics$CsymbolicAnalyses, 1)
  expectSameFit(dense, sparse)
})

test_that("grid augmentation is skipped quietly or errors when forced on ineligible models", {
  dat <- gridAugmentData(seed=15)
  ml <- mmes(y ~ trial, random=~geno, rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial),
             data=dat, control=mmes_control(task=list(objective="ml"), convergence=list(maxIters=10),
                                            diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(isTRUE(ml$engineDiagnostics$gridAugmentation$active))
  expect_error(
    mmes(y ~ trial, random=~geno, rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial),
         data=dat, control=mmes_control(task=list(objective="ml"), representation=list(gridAugment="force"),
                                        convergence=list(maxIters=10), diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "gridAugment='force'")

  diag <- gridAugmentFit(dat, "auto", rcov=~vsm(dsm(trial), ism(units)))
  expect_false(isTRUE(diag$engineDiagnostics$gridAugmentation$active))
  expect_error(gridAugmentFit(dat, "force", rcov=~vsm(dsm(trial), ism(units))), "gridAugment='force'")
})

test_that("the gridAugment option sets the default mode", {
  dat <- gridAugmentData(seed=16)
  old <- options(sommer.mme.gridAugment="off")
  on.exit(options(old), add=TRUE)
  fit <- mmes(y ~ trial, random=~geno, rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial),
              data=dat, control=mmes_control(convergence=list(maxIters=10), diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  expect_false(isTRUE(fit$engineDiagnostics$gridAugmentation$active))
})

test_that("reduced-rank AI uses direct residual traces on augmented grids", {
  set.seed(17)
  dat <- expand.grid(trial=factor(c("A", "B", "C")), range=factor(1:6), row=factor(1:8))
  dat$geno <- factor(sample(sprintf("g%02d", 1:15), nrow(dat), replace=TRUE))
  gEff <- rnorm(15)
  dat$y <- as.numeric(dat$trial) + gEff[as.integer(dat$geno)] * c(1, 0.8, 1.2)[as.integer(dat$trial)] + rnorm(nrow(dat))
  dat <- dat[-sample(nrow(dat), 15), ]
  fit <- function(mode) mmes(y ~ trial, random=~vsm(rrm(trial, 1), ism(geno)) + vsm(dsm(trial), ism(geno)),
    rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial), data=dat,
    control=mmes_control(representation=list(gridAugment=mode), convergence=list(maxIters=6),
                         diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  off <- suppressWarnings(fit("off"))
  auto <- suppressWarnings(fit("auto"))
  expect_true(auto$engineDiagnostics$gridAugmentation$active)
  expect_gt(auto$engineDiagnostics$directResidualTraces, 0)
  expect_equal(as.numeric(auto$llik), as.numeric(off$llik), tolerance=1e-6)
  expect_equal(unlist(auto$theta), unlist(off$theta), tolerance=1e-5)
  expect_equal(as.numeric(fitted(auto)), as.numeric(fitted(off)), tolerance=1e-5)
})
