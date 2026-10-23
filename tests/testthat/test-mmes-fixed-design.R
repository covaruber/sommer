test_that("zero columns after response filtering have actionable diagnostics", {
  maturityRecords <- data.frame(
    maturity=c(10, 12, 14, 16, NA, NA),
    locale=factor(rep(c("north", "south", "unused"), each=2L)),
    covariate=c(0, 0, 0, 0, 1, 2)
  )
  quietControl <- mmes_control(diagnostics=list(verbose=FALSE))
  expect_error(
    mmes(maturity ~ locale + covariate, data=maturityRecords, control=quietControl),
    "1 of 3 columns are all zero \\(4 retained observations\\)"
  )
  expect_error(
    mmes(maturity ~ locale + covariate, data=maturityRecords, control=quietControl),
    "Affected fixed terms: covariate.*"
  )
  expect_error(
    mmes(maturity ~ locale + covariate, data=maturityRecords, control=quietControl),
    "droplevels\\(\\)"
  )
})

test_that("unobserved interaction cells are explained without dense rank analysis", {
  maturityRecords <- data.frame(
    maturity=c(10, 12, 20, 22),
    locale=factor(c("north", "north", "south", "south")),
    cultivar=factor(c("early", "early", "late", "late"))
  )
  expect_error(
    mmes(maturity ~ locale + locale:cultivar, data=maturityRecords,
         control=mmes_control(diagnostics=list(verbose=FALSE))),
    "Example zero columns:.*"
  )
  expect_error(
    mmes(maturity ~ locale + locale:cultivar, data=maturityRecords,
         control=mmes_control(diagnostics=list(verbose=FALSE))),
    "interaction\\(\\.\\.\\., drop=TRUE\\)"
  )
})

test_that("parent and cell effects report aliasing even without zero columns", {
  maturityRecords <- expand.grid(
    locale=factor(c("north", "south")),
    cultivar=factor(c("early", "late")),
    replicate=seq_len(3L)
  )
  maturityRecords$maturity <- seq_len(nrow(maturityRecords))
  maturityRecords$cell <- interaction(maturityRecords$locale,
                                      maturityRecords$cultivar, drop=TRUE)
  quietControl <- mmes_control(diagnostics=list(verbose=FALSE))
  expect_error(
    mmes(maturity ~ locale + cell, data=maturityRecords, control=quietControl),
    "rank 4 < 5 columns"
  )
  expect_error(
    mmes(maturity ~ locale + cell, data=maturityRecords, control=quietControl),
    "control=mmes_control\\(representation=list\\(fixedRank='drop'\\)\\)"
  )
  maturityRecords$covariate <- seq_len(nrow(maturityRecords))
  maturityRecords$duplicate <- maturityRecords$covariate
  expect_error(
    mmes(maturity ~ covariate + duplicate, data=maturityRecords, control=quietControl),
    "rank 2 < 3 columns"
  )
})

test_that("rank-reduction advice precedes long diagnostic labels and truncation", {
  for(hasZeroColumns in c(TRUE, FALSE)){
    diagnosticDesign <- Matrix::Matrix(cbind(diag(2L),
      if(hasZeroColumns) c(0, 0) else c(1, 1)), sparse=TRUE)
    colnames(diagnosticDesign) <- paste0(
      c("first", "second", "third"), strrep("long_coefficient_label", 100L))
    attr(diagnosticDesign, "assign") <- rep(1L, 3L)
    diagnosticError <- tryCatch(
      sommer:::.mmes_check_fixed_design(diagnosticDesign,
        strrep("long_term_label", 100L)),
      error=identity)
    expect_s3_class(diagnosticError, "error")
    visibleMessage <- substr(conditionMessage(diagnosticError), 1L, 1000L)
    expect_match(visibleMessage,
      "control=mmes_control(representation=list(fixedRank='drop'))", fixed=TRUE)
    expect_match(visibleMessage, "Omitted coefficients are constrained to zero", fixed=TRUE)
    expect_match(visibleMessage, "supported Gaussian models", fixed=TRUE)
    expect_lt(regexpr("control=mmes_control", visibleMessage, fixed=TRUE)[[1L]],
      regexpr("long_term_label", conditionMessage(diagnosticError), fixed=TRUE)[[1L]])
  }
})

test_that("wide nonzero designs report rank deficiency without QR warnings", {
  wideDesign <- Matrix::Matrix(cbind(diag(2L), c(1, 1)), sparse=TRUE)
  colnames(wideDesign) <- c("first", "second", "third")
  attr(wideDesign, "assign") <- rep(1L, 3L)
  expect_error(
    sommer:::.mmes_check_fixed_design(wideDesign, "cell"),
    "rank 2 < 3 columns"
  )
})

test_that("full-rank cell means still fit and setup preserves aliased designs", {
  maturityRecords <- data.frame(
    maturity=c(10, 12, 14, 16, 20, 22, 24, 26),
    locale=factor(rep(c("north", "south"), each=4L)),
    cultivar=factor(rep(c("early", "late"), each=2L, times=2L))
  )
  maturityRecords$cell <- interaction(maturityRecords$locale,
                                      maturityRecords$cultivar, drop=TRUE)
  cellFit <- mmes(maturity ~ 0 + cell, data=maturityRecords,
    control=mmes_control(diagnostics=list(verbose=FALSE)))
  expect_true(cellFit$convergence)
  expect_equal(as.numeric(fitted(cellFit)),
    as.numeric(ave(maturityRecords$maturity, maturityRecords$cell)),
    tolerance=1e-8)
  expect_equal(as.numeric(cellFit$theta[[1L]]), 2, tolerance=1e-4)

  aliasedSetup <- mmes(maturity ~ locale + cell, data=maturityRecords,
    control=mmes_control(output=list(setup=TRUE),
                         diagnostics=list(verbose=FALSE)))
  expect_equal(ncol(aliasedSetup$X), 5L)
  expect_equal(sommer:::.sparse_rank(aliasedSetup$X), 4L)
})
