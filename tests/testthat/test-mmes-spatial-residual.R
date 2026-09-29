test_that("complete separable spatial residuals stay descriptor-only", {
  set.seed(41)
  dat <- expand.grid(
    trial=factor(seq_len(2)),
    range=factor(seq_len(4)),
    row=factor(seq_len(5))
  )
  dat$y <- rnorm(nrow(dat))

  fit <- mmes(
    y ~ trial,
    rcov=~vsm(ar1m(range), ar1m(row), ism(units)),
    data=dat,
    nIters=3,
    verbose=FALSE
  )

  expect_true(fit$residualCovarianceCompact)
  expect_equal(dim(fit$theta[[1L]]), c(0L, 0L))
  expect_equal(fit$covParNative$parameter, c("variance", "rho", "rho"))
  expect_true(all(is.finite(fit$covParNative$estimate)))
})

test_that("ragged spatial residual blocks stay compact and match the dense engine", {
  set.seed(42)
  dat <- expand.grid(
    trial=factor(seq_len(2)),
    range=factor(seq_len(4)),
    row=factor(seq_len(5))
  )
  dat <- dat[-c(1L, 7L, 12L), , drop=FALSE]
  dat$plot <- paste0("p", seq_len(nrow(dat)))
  dat$y <- rnorm(nrow(dat))

  fitFor <- function(henderson){
    mmes(
      y ~ trial,
      rcov=~vsm(dsm(trial), ar1m(range), ar1m(row), ism(units)),
      data=dat,
      nIters=15,
      verbose=FALSE,
      henderson=henderson
    )
  }
  compact <- fitFor(TRUE)
  dense <- fitFor(FALSE)

  expect_true(compact$residualCovarianceCompact)
  expect_equal(dim(compact$theta[[1L]]), c(0L, 0L))
  expect_equal(unlist(compact$covPar), unlist(dense$covPar), tolerance=1e-4)

  layout <- mmes(
    y ~ trial,
    rcov=~vsm(dsm(trial), ar1m(range), ar1m(row), ism(units)),
    data=dat, returnParam=TRUE, verbose=FALSE
  )
  # A unique plot column must not split residual blocks.
  expect_equal(length(unique(layout$residualBlock)), 2L)
})

sectioned_trials <- function(seed){
  set.seed(seed)
  one <- function(trial, nr, nc, rhoRange, rhoRow, s2){
    g <- expand.grid(range=seq_len(nr), row=seq_len(nc))
    K <- s2 * kronecker(rhoRow^abs(outer(seq_len(nc), seq_len(nc), "-")),
                        rhoRange^abs(outer(seq_len(nr), seq_len(nr), "-")))
    g$y <- as.vector(t(chol(K)) %*% rnorm(nr * nc)) + 10
    g$trial <- trial
    g
  }
  dat <- rbind(one("A", 6, 7, 0.6, 0.1, 2), one("B", 5, 8, -0.2, 0.5, 4),
               one("C", 7, 6, 0.3, -0.3, 1))
  dat <- dat[-sample(nrow(dat), 12L), , drop=FALSE]
  dat$trial <- factor(dat$trial)
  dat$range <- factor(dat$range, levels=1:7)
  dat$row <- factor(dat$row, levels=1:8)
  dat
}

spatial_dsumm <- ~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial)

test_that("dsumm() estimates section-specific parameters exactly", {
  dat <- sectioned_trials(7)
  fitFor <- function(henderson){
    mmes(y ~ trial, rcov=spatial_dsumm, data=dat, nIters=40,
         verbose=FALSE, henderson=henderson)
  }
  grid <- fitFor(TRUE)
  dense <- fitFor(FALSE)

  expect_true(grid$residualCovarianceCompact)
  expect_true(dense$residualCovarianceCompact)
  expect_equal(dim(dense$theta[[1L]]), c(0L, 0L))
  expect_equal(unlist(grid$covPar), unlist(dense$covPar), tolerance=1e-5)

  native <- covparams_mmes(grid)
  expect_equal(native$factor[native$section == "B"],
               c("trial", "ar1m(range)", "ar1m(row)"))
  expect_equal(native$parameter[native$section == "B"], c("variance", "rho", "rho"))
  expect_setequal(unique(native$section), c("A", "B", "C"))
})

test_that("dsumm() separates into independent per-section fits", {
  dat <- sectioned_trials(8)
  joint <- mmes(y ~ trial, rcov=spatial_dsumm, data=dat, nIters=60,
                tolParConvLL=1e-9, verbose=FALSE)
  native <- covparams_mmes(joint)

  for(section in levels(dat$trial)){
    sub <- droplevels(dat[dat$trial == section, ])
    alone <- mmes(y ~ 1, rcov=~vsm(ar1m(range), ar1m(row), ism(units)),
                  data=sub, nIters=60, tolParConvLL=1e-9, verbose=FALSE)
    expect_equal(native$estimate[native$section %in% section],
                 covparams_mmes(alone)$estimate, tolerance=1e-3)
  }
})

test_that("dsumm(levels=) gives unlisted sections an iid residual", {
  dat <- sectioned_trials(10)
  # Section C has no spatial coordinates; it must still be kept as iid.
  dat$range[dat$trial == "C"] <- NA
  dat$row[dat$trial == "C"] <- NA

  fitFor <- function(henderson){
    mmes(y ~ trial,
         rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial, levels=c("A", "B")),
         data=dat, nIters=80, tolParConvLL=1e-10, tolParConvNorm=1e-10,
         verbose=FALSE, henderson=henderson)
  }
  grid <- fitFor(TRUE)
  dense <- fitFor(FALSE)
  expect_equal(nrow(grid$y), nrow(dat))
  expect_equal(unlist(grid$covPar), unlist(dense$covPar), tolerance=1e-6)

  native <- covparams_mmes(grid)
  expect_equal(native$parameter[native$section == "C"], "variance")
  expect_equal(native$estimate[native$section == "C"],
               stats::var(dat$y[dat$trial == "C"]), tolerance=1e-4)
})

test_that("dsumm() supports multi-parameter inner structures", {
  set.seed(12)
  dat <- expand.grid(site=factor(c("s1", "s2")), env=factor(c("e1", "e2", "e3")),
                     id=factor(seq_len(15)))
  dat$y <- rnorm(nrow(dat))
  fitFor <- function(henderson){
    mmes(y ~ site, rcov=~dsumm(vsm(usm(env), ism(units)), by=site),
         data=dat, nIters=80, tolParConvLL=1e-10, tolParConvNorm=1e-10,
         verbose=FALSE, henderson=henderson)
  }
  grid <- fitFor(TRUE)
  dense <- fitFor(FALSE)
  expect_equal(unlist(grid$covPar), unlist(dense$covPar), tolerance=1e-6)
  expect_equal(sum(covparams_mmes(grid)$section == "s2"), 7L)
})

test_that("dsumm() validates its input", {
  dat <- sectioned_trials(9)
  expect_error(
    mmes(y ~ trial, rcov=~dsumm(vsm(ar1m(range), ar1m(row), ism(units)), by=trial,
                                levels="Z"), data=dat, verbose=FALSE),
    "Unknown levels"
  )
  expect_error(
    mmes(y ~ trial, random=~dsumm(vsm(ar1m(range), ism(row)), by=trial),
         data=dat, verbose=FALSE),
    "residual vsm"
  )
  expect_error(
    mmes(y ~ trial, rcov=~dsumm(vsm(ar1m(range), ism(units)), by=trial[-1]),
         data=dat, verbose=FALSE),
    "one value per observation"
  )
})