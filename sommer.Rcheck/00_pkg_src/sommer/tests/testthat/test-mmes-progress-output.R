test_that("EM progress reports direction diagnostics and matches quiet fitting", {
  data(DT_example)
  fit <- function(verbose){
    mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=3), update=list(blendWeight=1), proposal=list(stepWeight=1), solve=list(method="ldlt"), diagnostics=list(verbose=verbose, dateWarning=FALSE)))
  }
  output <- capture.output(verbose_fit <- fit(TRUE))
  quiet_fit <- fit(FALSE)

  header <- grep("^ *iter +LogLik", output, value=TRUE)
  expect_length(header, 1L)
  expect_match(header, "deltaLL +update +policy +trust +line +notPD +time\\(s\\)")
  expect_lte(nchar(header), 100L)

  rows <- grep("^ +[0-9]+ +", output, value=TRUE)
  expect_length(rows, ncol(verbose_fit$llik))
  tokens <- strsplit(trimws(c(header, rows)), " +")
  expect_true(all(lengths(tokens) == length(tokens[[1L]])))
  column_ends <- function(line){
    positions <- gregexpr("\\S+", line)[[1L]]
    as.integer(positions + attr(positions, "match.length") - 1L)
  }
  for(row in rows){
    expect_identical(column_ends(row), column_ends(header))
  }

  expect_equal(verbose_fit$monitor, quiet_fit$monitor)
  expect_equal(verbose_fit$llik, quiet_fit$llik)
})

test_that("pure-AI progress retains numerical safeguard diagnostics", {
  data(DT_example)
  for(solver in c("ldlt", "cholmod")){
    fit <- function(verbose){
      mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example, control=mmes_control(convergence=list(maxIters=3), update=list(blendWeight=0), solve=list(method=solver), diagnostics=list(verbose=verbose, dateWarning=FALSE)))
    }
    output <- capture.output(verbose_fit <- fit(TRUE))
    quiet_fit <- fit(FALSE)
    header <- grep("^ *iter +LogLik", output, value=TRUE)
    expect_length(header, 1L)
    expect_match(header, "deltaLL +update +policy +trust +line +notPD +time\\(s\\)")
    expect_false(grepl("pivot|wall|cpu\\(sec\\)", header))
    expect_false(grepl("Qgain|stepNorm", header))

    rows <- grep("^ +[0-9]+ +", output, value=TRUE)
    fields <- strsplit(trimws(rows), " +")
    expect_length(rows, ncol(verbose_fit$llik))
    expect_true(all(lengths(fields) == 9L))
    expect_true(all(vapply(fields, `[[`, character(1), 4L) == "AI"))
    expect_true(all(vapply(fields, `[[`, character(1), 5L) == "off"))
    expect_equal(verbose_fit$monitor, quiet_fit$monitor)
    expect_equal(verbose_fit$llik, quiet_fit$llik)
  }
})

test_that("pure-AI progress reports failed structures and trust scaling", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(6L))),
    genotype=factor(paste0("G", seq_len(24L)))
  )
  trial_effect <- c(45, -30, -12, 52, 18, -20)
  genotype_effect <- sin(seq_len(24L) / 3) * 8
  data$BLUEs <- 100 + trial_effect[data$trial] +
    genotype_effect[data$genotype] + cos(seq_len(nrow(data)) / 5) * 3

  output <- capture.output(fit <- mmes(BLUEs~trial, random=~vsm(fam(trial, 2), ism(genotype)), rcov=~units, data=data, control=mmes_control(convergence=list(maxIters=12), solve=list(method="ldlt"), diagnostics=list(verbose=TRUE, dateWarning=FALSE), update=list(blendWeight=0))))
  rows <- grep("^ +[0-9]+ +", output, value=TRUE)
  fields <- strsplit(trimws(rows), " +")
  expect_length(rows, ncol(fit$llik))
  expect_true(all(lengths(fields) == 9L))
  header <- grep("^ *iter +LogLik", output, value=TRUE)
  column_ends <- function(line){
    positions <- gregexpr("\\S+", line)[[1L]]
    as.integer(positions + attr(positions, "match.length") - 1L)
  }
  for(row in rows){
    expect_identical(column_ends(row), column_ends(header))
  }
  not_pd <- vapply(fields, `[[`, character(1), 8L)
  alpha <- as.numeric(vapply(fields, `[[`, character(1), 6L))
  expect_true(all(grepl("^(-|[0-9]+(,[0-9]+)*)$", not_pd)))
  expect_true(any(grepl("^[0-9]+$", not_pd)))
  expect_true(all(is.finite(alpha) & alpha >= 0 & alpha <= 1))
  expect_true(any(alpha < 1))
  line_alpha <- as.numeric(vapply(fields, `[[`, character(1), 7L))
  expect_true(all(is.finite(line_alpha) & line_alpha >= 0 & line_alpha <= 1))
  expect_equal(line_alpha[1L], 1)
  expect_identical(any(line_alpha < 1), fit$engineDiagnostics$lineSearchHalvings > 0L)
  expect_false(any(grepl(paste0("local fallback|trust scaling applied|update rejected|",
    "REML likelihood decreased|Global likelihood line search accepted"), output)))
})

test_that("all mmes progress paths share compact columns without changing fits", {
  data(DT_example)
  dataset <- DT_example
  progressControl <- function(formulation="henderson", objective="reml",
                              solver="auto", recovery="em", blendWeight=NULL,
                              curvature="exact-ai", equations="normal"){
    mmes_control(formulation=formulation, task=list(objective=objective),
      solve=list(method=solver), update=list(recovery=recovery,
        blendWeight=blendWeight), curvature=list(method=curvature),
      representation=list(equations=equations),
      convergence=list(maxIters=3L),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  }
  cases <- list(
    progressControl(),
    progressControl(objective="ml"),
    progressControl(solver="cholmod"),
    progressControl(solver="pcg", blendWeight=0),
    progressControl(blendWeight=0, curvature="hybrid-ai"),
    progressControl(blendWeight=c(1, 0.4, 0)),
    progressControl(recovery="pxem", blendWeight=c(1, 0.4, 0)),
    progressControl(blendWeight=0, equations="saddle-point"),
    progressControl(formulation="direct", blendWeight=0),
    progressControl(formulation="direct", blendWeight=1),
    progressControl(formulation="direct", recovery="pxem",
      blendWeight=c(1, 0.4, 0)),
    progressControl(formulation="direct"))
  expected <- NULL
  for(control in cases){
    fit <- function(verbose){
      fitControl <- do.call(mmes_control, utils::modifyList(
        unclass(control), list(diagnostics=list(verbose=verbose)),
        keep.null=TRUE))
      mmes(fixed=Yield~Env, random=~Name, data=dataset, control=fitControl)
    }
    output <- capture.output(loud <- fit(TRUE))
    quiet <- fit(FALSE)
    header <- grep("^ *iter +LogLik", output, value=TRUE)
    expect_length(header, 1L)
    if(is.null(expected)) expected <- header
    expect_identical(header, expected)
    rows <- strsplit(trimws(grep("^ +[0-9]+ +", output, value=TRUE)), " +")
    expect_length(rows, length(loud$llik))
    expect_true(all(lengths(rows) == 9L))
    expect_equal(as.numeric(vapply(rows, `[[`, character(1), 2L)),
      as.numeric(loud$llik), tolerance=1e-6)
    expect_identical(rows[[1L]][3L], "-")
    elapsed <- as.numeric(vapply(rows, `[[`, character(1), 9L))
    expect_true(all(diff(elapsed) >= 0))
    for(field in c("llik", "monitor", "theta", "bu", "InfMat", "dLu")){
      expect_equal(loud[[field]], quiet[[field]], info=field)
    }
  }
})

test_that("recovery progress names the selected method and trigger", {
  data(DT_example)
  for(reml in c(TRUE, FALSE)) for(method in c("em", "pxem")){
    output <- capture.output(fit <- mmes(Yield~Env, random=~Name, data=DT_example, control=mmes_control(task=list(objective=if(reml) "reml" else "ml"), update=list(recovery=method), convergence=list(maxIters=6L), proposal=list(stepWeight=1e-7), diagnostics=list(verbose=TRUE, dateWarning=FALSE))))
    rows <- strsplit(trimws(grep("^ +[0-9]+ +", output, value=TRUE)), " +")
    expect_identical(rows[[1L]][4L], "AI")
    expect_identical(rows[[1L]][5L], paste0("auto:", method))
    recovery <- fit$engineDiagnostics$recovery$iterations
    expect_identical(rows[[recovery]][4L], toupper(method))
    expect_identical(rows[[recovery]][5L], "low-gain")
    expect_identical(rows[[recovery + 1L]][4L], "AI")
    expect_identical(rows[[recovery + 1L]][5L], paste0("auto:", method))
  }
  dataset <- DT_example
  dataset$Name2 <- dataset$Name
  output <- capture.output(fit <- mmes(Yield~Env, random=~Name+Name2, data=dataset, control=mmes_control(convergence=list(maxIters=3L), diagnostics=list(verbose=TRUE, dateWarning=FALSE))))
  expect_true(any(grepl("AI-reg +info-solve", output)))
  expect_equal(fit$engineDiagnostics$recovery$steps, 0L)
})
