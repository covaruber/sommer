.observed_prediction_fit <- function(){
  records <- data.frame(
    response=c(10, 12, 20, 22, 24, 30, 32, NA),
    locale=factor(c("north", "north", "north", "north", "north", "south", "south", "south")),
    cultivar=factor(c("early", "early", "late", "late", "late", "late", "late", "early"))
  )
  expect_warning(
    fittedModel <- mmes(response ~ locale + locale:cultivar, data=records,
      control=mmes_control(representation=list(fixedRank="drop"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "rank reduction"
  )
  fittedModel
}

test_that("rank-reduced classify prediction stops before constructing any grid", {
  fittedModel <- .observed_prediction_fit()
  local_mocked_bindings(
    sparse.model.matrix=function(...) stop("Cartesian allocation reached"),
    .mmes_vcov_solver=function(...) stop("Covariance allocation reached"),
    .package="sommer")
  expect_error(predict(fittedModel, D="locale:cultivar"), "before constructing.*Cartesian")
  expect_error(predict(fittedModel, D="locale:cultivar"), "observed=TRUE")
  expect_error(predict(fittedModel, D="locale:cultivar", levels=list("locale:cultivar"="north:early")),
    "before constructing.*Cartesian")
})

test_that("observed predictions match estimable full contrasts without grid construction", {
  fittedModel <- .observed_prediction_fit()
  local_mocked_bindings(
    sparse.model.matrix=function(...) stop("Cartesian allocation reached"),
    .package="sommer")
  predictedCells <- predict(fittedModel, D="locale:cultivar", observed=TRUE)
  expect_equal(predictedCells$pvals[[1L]], c("north:early", "north:late", "south:late"))
  expect_equal(predictedCells$pvals$predicted.value, c(11, 22, 31), tolerance=1e-8)
  expect_equal(predictedCells$observed$observations, c(2L, 3L, 2L))
  expect_match(predictedCells$observed$weighting, "observation-weighted")
  expect_equal(Matrix::rowSums(abs(predictedCells$D)), c(1, 2, 3),
               ignore_attr=TRUE)
  explicitPrediction <- predict(fittedModel, D=predictedCells$D)
  expect_equal(predictedCells$vcov, explicitPrediction$vcov, tolerance=1e-10)
  expect_equal(predictedCells$pvals$std.error, explicitPrediction$pvals$std.error)
  expect_equal(diag(predictedCells$vcov),
    as.numeric(fittedModel$theta[[1L]]) / c(2, 3, 2), tolerance=1e-7)
  expect_equal(predictedCells$VarU, matrix(0, 3L, 3L))
  expect_equal(predictedCells$sampling.vcov, predictedCells$vcov)
  expect_equal(dim(predictedCells$D), c(3L, nrow(fittedModel$bu)))
  expect_s4_class(predictedCells$D, "sparseMatrix")
})

test_that("observed margins use retained observation weights, not balanced-grid weights", {
  fittedModel <- .observed_prediction_fit()
  predictedLocales <- predict(fittedModel, D="locale", observed=TRUE, VarU=FALSE)
  expect_equal(predictedLocales$pvals$predicted.value, c(17.6, 31), tolerance=1e-8)
  expect_equal(predictedLocales$observed$observations, c(5L, 2L))
  reversedClassify <- predict(fittedModel, D="cultivar:locale", observed=TRUE)
  expect_equal(reversedClassify$pvals[[1L]], c("early:north", "late:north", "late:south"))
  expect_equal(reversedClassify$pvals$predicted.value, c(11, 22, 31), tolerance=1e-8)
})

test_that("observed levels restrict and order existing combinations before allocations", {
  fittedModel <- .observed_prediction_fit()
  selectedPredictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
    levels=list("locale:cultivar"=c("south:late", "north:early")),
    sed=TRUE, pairwise=TRUE, df="residual")
  expect_equal(selectedPredictions$pvals[[1L]], c("south:late", "north:early"))
  expect_equal(selectedPredictions$pvals$predicted.value, c(31, 11), tolerance=1e-8)
  expect_equal(selectedPredictions$observed$observations, c(2L, 2L))
  expect_equal(nrow(selectedPredictions$pairwise), 1L)
  expect_true(all(is.finite(selectedPredictions$sed)))
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    levels=list("locale:cultivar"="south:early")), "not observed.*south:early")
})

test_that("memory guards run before grouped design and covariance allocation", {
  fittedModel <- .observed_prediction_fit()
  local_mocked_bindings(
    .mmes_vcov_solver=function(...) stop("Covariance allocation reached"),
    .package="sommer")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE, maxMemoryMB=1e-8),
    "exceeding maxMemoryMB")
  fittedModel$control$execution$denseMemoryMB <- 1e-8
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE),
    "exceeding maxMemoryMB")
  expect_error(predict(fittedModel, D=fittedModel$W, PEV=FALSE),
    "PEV=FALSE does not remove")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
                       maxMemoryMB=0), "maxMemoryMB must be")
})

test_that("observed prediction does not silently reinterpret custom rules", {
  fittedModel <- .observed_prediction_fit()
  expect_error(predict(fittedModel, D="locale:cultivar", observed=NA), "observed must be")
  expect_error(predict(fittedModel, D=fittedModel$W, observed=TRUE), "explicit D matrices")
  expect_equal(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=fittedModel$Dtable)$pvals$predicted.value, c(11, 22, 31), tolerance=1e-8)
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    levels=list(locale="north")), "levels must be a named list")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    levels=list("locale:cultivar"=c("north:early", "north:early"))),
    "unique observed")
  expect_error(predict(fittedModel, D="missing", observed=TRUE), "fitted fixed formula")
  fittedModel$fixedRank <- NULL
  expect_error(predict(fittedModel, D="locale", observed=TRUE), "requires a rank-reduced fit")
})

.observed_random_fit <- function(structured=FALSE){
  set.seed(317)
  records <- expand.grid(locale=factor(c("north", "south")),
    cultivar=factor(c("early", "late")), block=factor(c("b1", "b2", "b3")),
    replicate=seq_len(3L))
  records <- records[records$locale == "north" | records$cultivar == "late", ]
  records <- records[!(records$locale == "north" & records$cultivar == "early" &
    records$block == "b3") & !(records$block == "b2" & records$replicate > 1L), ]
  records$response <- 10 + as.integer(records$locale) +
    c(-2, 1, 3)[as.integer(records$block)] + rnorm(nrow(records), sd=0.4)
  expect_warning(
    fittedModel <- mmes(response ~ locale + locale:cultivar,
      random=if(structured) ~vsm(dsm(locale), ism(block)) else ~block,
      data=records, control=mmes_control(representation=list(fixedRank="drop"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "rank reduction")
  fittedModel
}

test_that("observed Dtable includes empirical random contributions and their covariance", {
  fittedModel <- .observed_random_fit()
  randomRules <- fittedModel$Dtable
  randomRules$include[randomRules$type == "random"] <- TRUE
  local_mocked_bindings(sparse.model.matrix=function(...) stop("Cartesian allocation reached"),
    .package="sommer")
  predictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
                         Dtable=randomRules, pairwise=TRUE, df="residual")
  rowLabels <- paste(fittedModel$data$locale, fittedModel$data$cultivar, sep=":")
  expectedValues <- tapply(as.numeric(fitted(fittedModel)), rowLabels, mean)
  expect_equal(predictions$pvals$predicted.value,
    as.numeric(expectedValues[predictions$pvals[[1L]]]), tolerance=1e-8)
  explicit <- predict(fittedModel, D=predictions$D, pairwise=TRUE, df="residual")
  for(field in c("vcov", "VarU", "sampling.vcov", "pairwise")){
    expect_equal(predictions[[field]], explicit[[field]], tolerance=1e-10)
  }
  expect_equal(predictions$observed$randomTerms,
               randomRules$term[randomRules$type == "random"])
  expect_true(any(predictions$pvals$var.u > 0))
  expect_true(all(predictions$Dtable$include[predictions$Dtable$type == "random"]))
  reordered <- randomRules[rev(seq_len(nrow(randomRules))), ]
  expect_equal(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=reordered)$D, predictions$D)
  expect_equal(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=predictions$Dtable)$D, predictions$D)
})

test_that("observed random include-and-average uses equally weighted co-occurring levels", {
  fittedModel <- .observed_random_fit()
  randomRules <- fittedModel$Dtable
  randomRules$include[randomRules$type == "random"] <- TRUE
  randomRules$average[randomRules$type == "random"] <- TRUE
  predictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE, Dtable=randomRules)
  randomColumns <- nrow(fittedModel$b) + seq_len(nrow(fittedModel$u))
  expect_equal(as.numeric(predictions$D[1L, randomColumns]), c(0.5, 0.5, 0))
  expect_equal(as.numeric(predictions$D[2L, randomColumns]), rep(1/3, 3L))
  empiricalRules <- randomRules
  empiricalRules$average[empiricalRules$type == "random"] <- FALSE
  empirical <- predict(fittedModel, D="locale:cultivar", observed=TRUE, Dtable=empiricalRules)
  expect_equal(as.numeric(empirical$D[1L, randomColumns]), c(0.75, 0.25, 0))
  expect_false(isTRUE(all.equal(predictions$D, empirical$D)))
})

test_that("observed random averages support selected fitted levels and structured terms", {
  for(structured in c(FALSE, TRUE)){
    fittedModel <- .observed_random_fit(structured)
    randomRules <- fittedModel$Dtable
    randomRow <- which(randomRules$type == "random")
    randomRules$average[randomRow] <- TRUE
    randomTerm <- randomRules$term[randomRow]
    selected <- setNames(list("b3"), randomTerm)
    predictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
      Dtable=randomRules, levels=selected)
    randomColumns <- nrow(fittedModel$b) + seq_len(nrow(fittedModel$u))
    expectedWeights <- rep(0, length(randomColumns))
    expectedWeights[rep(rownames(fittedModel$uList[[1L]]),
      times=ncol(fittedModel$uList[[1L]])) == "b3"] <- if(structured) 0.5 else 1
    expect_equal(as.numeric(predictions$D[1L, randomColumns]), expectedWeights)
    expect_equal(predictions$levels[[randomTerm]], "b3")
    expect_false("levels" %in% names(predictions$Dtable))
    explicit <- predict(fittedModel, D=predictions$D)
    expect_equal(predictions$vcov, explicit$vcov, tolerance=1e-10)
    if(structured){
      cellSelection <- setNames(list("north:b3"), randomTerm)
      cellPredictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
        Dtable=randomRules, levels=cellSelection)
      expect_equal(sum(cellPredictions$D[1L, randomColumns]), 1)
      expect_equal(sum(cellPredictions$D[1L, randomColumns] != 0), 1L)
      randomRules$include[randomRow] <- TRUE
      randomRules$average[randomRow] <- FALSE
      empirical <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
        Dtable=randomRules)
      groupLabels <- paste(fittedModel$data$locale, fittedModel$data$cultivar, sep=":")
      firstGroup <- groupLabels == empirical$pvals[[1L]][1L]
      expect_equal(as.numeric(empirical$D[1L, randomColumns]),
        as.numeric(Matrix::colMeans(fittedModel$W[firstGroup, randomColumns, drop=FALSE])))
      randomRules$average[randomRow] <- TRUE
      equalWeights <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
        Dtable=randomRules)
      expect_equal(sum(equalWeights$D[1L, randomColumns] != 0), 2L)
      expect_equal(sum(equalWeights$D[1L, randomColumns]), 1)
    }
  }
})

test_that("observed classify subsets and random restrictions are independently reusable", {
  fittedModel <- .observed_random_fit()
  randomRules <- fittedModel$Dtable
  randomRow <- which(randomRules$type == "random")
  randomRules$average[randomRow] <- TRUE
  restrictions <- c(list("locale:cultivar"=c("south:late", "north:early")),
    setNames(list(c("b1", "b3")), randomRules$term[randomRow]))
  predictions <- predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=randomRules, levels=restrictions)
  expect_equal(predictions$pvals[[1L]], c("south:late", "north:early"))
  randomColumns <- nrow(fittedModel$b) + seq_len(nrow(fittedModel$u))
  expect_equal(as.numeric(predictions$D[1L, randomColumns]), c(0.5, 0, 0.5))
  expect_equal(predictions$levels, restrictions)
  expect_equal(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=predictions$Dtable, levels=predictions$levels)$D, predictions$D)
})

test_that("observed Dtable rejects invalid and non-estimable rules before solves", {
  fittedModel <- .observed_random_fit()
  randomRules <- fittedModel$Dtable
  randomRow <- which(randomRules$type == "random")
  randomRules$include[randomRow] <- TRUE
  randomRules$average[randomRow] <- TRUE
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=randomRules, levels=setNames(list("b3"), randomRules$term[randomRow])),
    "No selected random levels co-occur")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=randomRules, levels=setNames(list("absent"), randomRules$term[randomRow])),
    "Unknown levels")
  malformedRules <- randomRules[-1L, ]
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=malformedRules), "each fitted term exactly once")
  malformedRules <- randomRules
  malformedRules$include[1L] <- TRUE
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=malformedRules), "Only random-term flags")
  malformedRules <- randomRules
  malformedRules$average[randomRow] <- NA
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=malformedRules), "nonmissing logical")
  local_mocked_bindings(.mmes_vcov_solver=function(...) stop("Covariance allocation reached"),
    .package="sommer")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE,
    Dtable=randomRules, maxMemoryMB=1e-8), "exceeding maxMemoryMB")
})

test_that("observed fixed predictions exclude random effects and retain numeric covariates", {
  set.seed(274)
  records <- expand.grid(locale=factor(c("north", "south")),
    cultivar=factor(c("early", "late")), lineage=factor(seq_len(6L)),
    replicate=seq_len(3L))
  records <- records[records$locale == "north" | records$cultivar == "late", ]
  records$temperature <- rnorm(nrow(records))
  records$response <- 10 + as.integer(records$locale) +
    2 * records$temperature + rep(rnorm(6L), length.out=nrow(records)) +
    rnorm(nrow(records))
  expect_warning(
    fittedModel <- mmes(response ~ locale + locale:cultivar + temperature,
      random=~lineage, data=records,
      control=mmes_control(representation=list(fixedRank="drop"),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE))),
    "rank reduction"
  )
  predictedCells <- predict(fittedModel, D="locale:cultivar", observed=TRUE)
  fixedDesign <- fittedModel$W[, seq_len(nrow(fittedModel$b)), drop=FALSE]
  fixedPredictor <- as.numeric(fixedDesign %*% fittedModel$b)
  cellLabels <- paste(records$locale, records$cultivar, sep=":")
  expectedMeans <- tapply(fixedPredictor, cellLabels, mean)
  expect_equal(predictedCells$pvals$predicted.value,
    as.numeric(expectedMeans[predictedCells$pvals[[1L]]]), tolerance=1e-8)
  expect_true(all(predictedCells$D[, nrow(fittedModel$b) +
    seq_len(nrow(fittedModel$u)), drop=FALSE] == 0))
  expect_equal(predictedCells$VarU, matrix(0, 3L, 3L))
  expect_error(predict(fittedModel, D="temperature", observed=TRUE),
               "nonmissing categorical")
})

test_that("ambiguous observed labels cannot merge distinct combinations silently", {
  fittedModel <- .observed_prediction_fit()
  fittedModel$data$locale <- c("a:b", "a:b", "a", "a", "a", "other", "other")
  fittedModel$data$cultivar <- c("c", "c", "b:c", "b:c", "b:c", "last", "last")
  expect_error(predict(fittedModel, D="locale:cultivar", observed=TRUE),
               "labels are ambiguous")
})
