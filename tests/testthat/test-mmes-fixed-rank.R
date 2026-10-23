.rank_control <- function(formulation="henderson", objective="reml", inverse=0L){
  mmes_control(formulation=formulation, representation=list(fixedRank="drop"),
    task=list(objective=objective), output=list(inverse=inverse),
    convergence=list(maxIters=60L, likelihoodTol=1e-7, parameterTol=1e-7),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
}

.rank_records <- function(){
  set.seed(805)
  records <- expand.grid(locale=factor(c("north", "south")),
    cultivar=factor(c("early", "late")), lineage=factor(seq_len(8L)),
    replicate=seq_len(2L))
  records$cell <- interaction(records$locale, records$cultivar, drop=TRUE)
  lineageEffects <- rnorm(8L, sd=2)
  records$response <- 10 + as.integer(records$cell) +
    lineageEffects[as.integer(records$lineage)] + rnorm(nrow(records))
  records
}

test_that("fixed-rank control is opt-in and strictly validated", {
  expect_identical(mmes_control()$representation$fixedRank, "error")
  expect_identical(.rank_control("direct")$representation$fixedRank, "drop")
  expect_error(mmes_control(representation=list(fixedRank=TRUE)),
    "representation\\$fixedRank")
  expect_error(mmes_control(representation=list(fixedRank="automatic")),
    "representation\\$fixedRank")
  expect_error(mmes(response ~ locale + cell, data=.rank_records(),
    control=mmes_control(diagnostics=list(verbose=FALSE))),
    "rank deficient")
})

test_that("native sparse basis reconstructs zero, aliased and wide designs", {
  designExamples <- list(
    cbind(1, seq_len(8L), 2 * seq_len(8L), 0),
    cbind(diag(3L), c(1, 1, 1), 0),
    cbind(1, seq_len(8L) * 1e-10, seq_len(8L) * 1e10),
    cbind(1, seq_len(8L), sin(seq_len(8L)),
          seq_len(8L) + 2 * sin(seq_len(8L)))
  )
  for(designValues in designExamples){
    originalDesign <- as(as(as(Matrix::Matrix(designValues, sparse=TRUE),
      "dMatrix"), "generalMatrix"), "CsparseMatrix")
    colnames(originalDesign) <- paste0("column", seq_len(ncol(originalDesign)))
    basis <- sommer:::fixed_design_basis_cpp(originalDesign)
    expect_equal(as.matrix(originalDesign[, basis$retained, drop=FALSE] %*%
      basis$relation), as.matrix(originalDesign), tolerance=1e-12,
      ignore_attr=TRUE)
    expect_equal(length(basis$retained), basis$rank)
    expect_equal(sort(c(basis$retained, basis$omitted)), seq_len(ncol(originalDesign)))
    expect_true(all(diff(basis$retained) > 0L))
    expect_s4_class(basis$relation, "dgCMatrix")
  }
  zeroDesign <- as(as(Matrix::Matrix(0, 5, 2, sparse=TRUE), "generalMatrix"),
                   "CsparseMatrix")
  expect_error(sommer:::fixed_design_basis_cpp(zeroDesign), "No nonzero fixed-effect")
  nonfiniteDesign <- as(Matrix::Matrix(matrix(c(1, Inf, 2, 3), 2), sparse=TRUE),
                        "CsparseMatrix")
  expect_error(sommer:::fixed_design_basis_cpp(nonfiniteDesign), "non-finite")
})

test_that("both engines preserve fits and random partitions under ML and REML", {
  records <- .rank_records()
  for(formulationName in c("henderson", "direct")){
    for(objectiveName in c("reml", "ml")){
      expect_warning(
        reducedFit <- mmes(response ~ locale + cell, random=~lineage, data=records,
          control=.rank_control(formulationName, objectiveName)),
        "rank reduction retained 4 of 5"
      )
      referenceFit <- mmes(response ~ 0 + cell, random=~lineage, data=records,
        control=mmes_control(formulation=formulationName, task=list(objective=objectiveName),
          convergence=list(maxIters=60L, likelihoodTol=1e-7, parameterTol=1e-7),
          diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      expect_true(reducedFit$convergence)
      expect_true(referenceFit$convergence)
      expect_equal(unlist(reducedFit$theta), unlist(referenceFit$theta), tolerance=1e-5)
      referenceDesign <- cbind(model.matrix(~0 + cell, records),
                               model.matrix(~0 + lineage, records))
      referenceFitted <- referenceDesign %*% referenceFit$bu
      expect_equal(as.numeric(fitted(reducedFit)), as.numeric(referenceFitted),
                   tolerance=1e-5)
      expect_equal(as.numeric(reducedFit$u), as.numeric(referenceFit$u), tolerance=1e-5)
      expect_equal(nrow(reducedFit$b), 5L)
      expect_equal(as.numeric(reducedFit$b[reducedFit$fixedRank$omitted, ]), 0)
      expect_equal(reducedFit$partitions[[1L]][1L, 1L], 6)
      expect_equal(ncol(reducedFit$W), nrow(reducedFit$bu))
      expect_equal(nrow(reducedFit$coefficientMap), nrow(reducedFit$bu))
      expect_equal(reducedFit$AIC + 2 * tail(as.numeric(reducedFit$llik), 1L), 8,
                   tolerance=1e-8)
      fixedSummary <- summary(reducedFit)$betas
      expect_true(all(is.na(fixedSummary$Std.Error[reducedFit$fixedRank$omitted])))
      expect_true(all(is.na(fixedSummary$t.value[reducedFit$fixedRank$omitted])))
      expect_match(fixedSummary$Status[reducedFit$fixedRank$omitted], "aliased")
      if(formulationName == "direct") expect_equal(dim(reducedFit$VarBeta), c(4L, 4L))
    }
  }
})

test_that("zero interaction columns and setup use the same full-rank subset", {
  records <- .rank_records()
  records <- records[records$locale == "north" | records$cultivar == "late", ]
  expect_warning(
    reducedFit <- mmes(response ~ locale + locale:cultivar, data=records,
      control=.rank_control()), "rank reduction"
  )
  setupControl <- .rank_control()
  setupControl$output$setup <- TRUE
  setup <- mmes(response ~ locale + locale:cultivar, data=records, control=setupControl)
  expect_equal(setup$fixedRank$retained, reducedFit$fixedRank$retained)
  expect_equal(ncol(setup$X), 3L)
  expect_equal(nrow(reducedFit$b), 4L)
  expect_equal(as.numeric(fitted(reducedFit)),
    as.numeric(ave(records$response, records$cell)), tolerance=1e-7)
})

test_that("rank selection uses the rows retained after missing-response filtering", {
  records <- data.frame(response=c(10, 12, 20, 22, NA, NA),
    locale=factor(rep(c("north", "south", "unused"), each=2L)),
    covariate=c(0, 0, 0, 0, 1, 2))
  expect_warning(
    reducedFit <- mmes(response ~ locale + covariate, data=records,
      control=.rank_control()), "rank reduction retained 2 of 3")
  expect_equal(nrow(reducedFit$W), 4L)
  expect_equal(reducedFit$fixedRank$columnNames[reducedFit$fixedRank$omitted],
               "covariate")
  expect_equal(as.numeric(fitted(reducedFit)), c(11, 11, 21, 21), tolerance=1e-8)
  expect_equal(as.numeric(reducedFit$theta[[1L]]), 2, tolerance=1e-5)
})

test_that("Henderson prediction distinguishes estimable means from constrained coefficients", {
  records <- .rank_records()
  expect_warning(
    reducedFit <- mmes(response ~ locale + cell, random=~lineage, data=records,
      control=.rank_control(inverse=2L)), "rank reduction"
  )
  fixedDesign <- model.matrix(~locale + cell, records)
  observedContrasts <- cbind(fixedDesign[seq_len(4L), , drop=FALSE],
                             matrix(0, 4L, nrow(reducedFit$u)))
  prediction <- predict(reducedFit, D=observedContrasts, VarU=FALSE,
                         df="residual", pairwise=TRUE)
  expect_equal(prediction$pvals$predicted.value,
    as.numeric(observedContrasts %*% reducedFit$bu))
  expect_true(all(is.finite(prediction$pvals$std.error)))
  defaultPrediction <- predict(reducedFit, D=observedContrasts)
  expect_equal(defaultPrediction$vcov, prediction$vcov, tolerance=1e-10)
  expect_true(all(is.finite(defaultPrediction$pvals$var.u)))
  expect_equal(sommer:::.fixed_rank(reducedFit), 4L)
  constrainedContrast <- matrix(0, 1L, nrow(reducedFit$bu))
  constrainedContrast[1L, reducedFit$fixedRank$omitted[1L]] <- 1
  expect_error(predict(reducedFit, D=constrainedContrast, VarU=FALSE), "not estimable")
  expect_error(wald_mmes(reducedFit), "Wald term tests.*not available")
  expect_error(anova(reducedFit, reducedFit), "Likelihood-ratio comparisons.*not available")
  expect_error(sommer:::.mmes_df_machinery(reducedFit), "not available")
  postFit <- postPEV(reducedFit, mode=2L)
  expect_equal(as.matrix(postFit$Ci), as.matrix(reducedFit$Ci), tolerance=1e-8)
  diagonalFit <- postPEV(reducedFit, mode=1L)
  expect_equal(Matrix::diag(diagonalFit$Ci), Matrix::diag(reducedFit$Ci),
               tolerance=1e-8)
  expect_true(all(is.na(summary(postFit)$betas$Std.Error[reducedFit$fixedRank$omitted])))
  restoredFit <- unserialize(serialize(reducedFit, NULL))
  expect_equal(predict(restoredFit, D=observedContrasts, VarU=FALSE)$vcov,
               prediction$vcov, tolerance=1e-10)
  if(requireNamespace("emmeans", quietly=TRUE)){
    expect_error(emmeans::emmeans(reducedFit, ~locale), "emmeans.*not available")
  }
})

test_that("full-rank opt-in leaves outputs unchanged without warnings", {
  records <- .rank_records()
  controlDefault <- mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  controlReduced <- mmes_control(representation=list(fixedRank="drop"),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  referenceFit <- mmes(response ~ locale, random=~lineage, data=records, control=controlDefault)
  expect_warning(
    optedFit <- mmes(response ~ locale, random=~lineage, data=records, control=controlReduced),
    NA
  )
  expect_null(optedFit$fixedRank)
  expect_equal(optedFit$b, referenceFit$b, tolerance=1e-10)
  expect_equal(optedFit$theta, referenceFit$theta, tolerance=1e-10)
})

test_that("rank reduction composes with REML omissions and latent RR mappings", {
  records <- .rank_records()
  capturedWarnings <- character()
  reducedFit <- withCallingHandlers(
    mmes(response ~ locale + cell,
      random=~vsm(rrm(locale, 1L, loadings=matrix(c(1, 0.7), 2L, 1L),
                       fixed=rep(TRUE, 2L)), ism(lineage)) + vsm(ism(locale)),
      data=records, control=.rank_control(inverse=2L)),
    warning=function(fitWarning){
      capturedWarnings <<- c(capturedWarnings, conditionMessage(fitWarning))
      invokeRestart("muffleWarning")
    })
  expect_length(capturedWarnings, 2L)
  expect_true(any(grepl("rank reduction", capturedWarnings)))
  expect_true(any(grepl("Omitting 1 structurally unestimable", capturedWarnings)))
  expect_true(reducedFit$engineDiagnostics$latentReducedRank)
  expect_equal(ncol(reducedFit$coefficientMap), nrow(reducedFit$C))
  expect_equal(nrow(reducedFit$coefficientMap), nrow(reducedFit$bu))
  expect_true(all(reducedFit$uList[[2L]] == 0))
  expect_equal(as.numeric(fitted(reducedFit)),
               as.numeric(reducedFit$W %*% reducedFit$bu))
  observedContrasts <- as.matrix(reducedFit$W[seq_len(4L), , drop=FALSE])
  predictions <- predict(reducedFit, D=observedContrasts, VarU=FALSE)
  expect_true(all(is.finite(predictions$pvals$std.error)))
  postFit <- postPEV(reducedFit, mode=2L)
  expect_equal(as.matrix(postFit$Ci), as.matrix(reducedFit$Ci), tolerance=1e-8)
})

test_that("alternate Henderson solvers use the same selected basis", {
  records <- .rank_records()
  referenceFit <- NULL
  for(solverName in c("ldlt", "cholmod", "pcg")){
    solverControl <- mmes_control(representation=list(fixedRank="drop"),
      solve=list(method=solverName), convergence=list(maxIters=1L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))
    expect_warning(
      solverFit <- mmes(response ~ locale + cell, random=~lineage,
        data=records, control=solverControl), "rank reduction")
    if(is.null(referenceFit)) referenceFit <- solverFit
    expect_equal(solverFit$fixedRank$retained, referenceFit$fixedRank$retained)
    expect_equal(as.numeric(fitted(solverFit)), as.numeric(fitted(referenceFit)),
                 tolerance=1e-6)
    observedContrasts <- as.matrix(solverFit$W[seq_len(4L), , drop=FALSE])
    expect_true(all(is.finite(predict(solverFit, D=observedContrasts,
                                     VarU=FALSE)$pvals$std.error)))
  }
  saddleControl <- mmes_control(
    representation=list(fixedRank="drop", equations="saddle-point"),
    update=list(blendWeight=0), convergence=list(maxIters=1L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  expect_warning(
    saddleFit <- mmes(response ~ locale + cell, random=~lineage,
      data=records, control=saddleControl), "rank reduction")
  expect_equal(saddleFit$fixedRank$retained, referenceFit$fixedRank$retained)
  expect_equal(as.numeric(fitted(saddleFit)), as.numeric(fitted(referenceFit)),
               tolerance=1e-6)
})

test_that("unsupported workflows fail explicitly rather than ignore rank policy", {
  records <- .rank_records()
  expect_error(mmes(response ~ locale + cell, data=records, family=poisson(),
    control=.rank_control()), "requires a Gaussian identity-link")
  expect_error(mmes(response ~ locale + cell, data=records,
    control=mmes_control(representation=list(fixedRank="drop", compression="groups"))),
    "without PQL.*compression")
  expect_error(mmes_control(task=list(mode="solve"),
    representation=list(fixedRank="drop")), "representation\\$fixedRank")
  expect_error(mmes(response ~ locale + cell, random=~lineage, data=records,
    control=mmes_control(representation=list(fixedRank="drop", factorScores="fixed-shape"))),
    "without PQL.*factor scores")
  lineagePrecision <- diag(nlevels(records$lineage))
  dimnames(lineagePrecision) <- rep(list(levels(records$lineage)), 2L)
  attr(lineagePrecision, "inverse") <- TRUE
  expect_error(mmes(response ~ locale + cell,
    random=~vsm(ism(lineage), Gu=lineagePrecision, rotation=TRUE),
    data=records, control=.rank_control()),
    "does not currently support rotation=TRUE")
})
