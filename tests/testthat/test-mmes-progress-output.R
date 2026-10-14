test_that("Henderson progress rows include PD and trust diagnostics", {
  data(DT_example)
  for(solver in c("ldlt", "cholmod")){
    fit <- function(verbose){
      mmes(Yield~Env, random=~Name, rcov=~units, data=DT_example,
        nIters=3, solver=solver, verbose=verbose, dateWarning=FALSE)
    }
    output <- capture.output(verbose_fit <- fit(TRUE))
    quiet_fit <- fit(FALSE)
    header <- grep("^ *iteration", output, value=TRUE)
    expect_length(header, 1L)
    expect_match(header, "EMweight +notPD +trustAlpha +lineAlpha")
    expect_identical(grepl("pivot", header), solver == "ldlt")

    rows <- grep("^ +[0-9]+ +", output, value=TRUE)
    fields <- strsplit(trimws(rows), " +")
    expect_length(rows, ncol(verbose_fit$llik))
    expect_true(all(lengths(fields) == if(solver == "ldlt") 10L else 9L))
    column_ends <- function(line){
      positions <- gregexpr("\\S+", line)[[1L]]
      as.integer(positions + attr(positions, "match.length") - 1L)
    }
    for(row in rows){
      expect_identical(column_ends(row), column_ends(header))
    }
    for(values in fields){
      expect_match(values[3L], "^[0-9]{2}:[0-9]{2}:[0-9]{2}$")
      numeric_columns <- c(2L, 4L, 6L, 8L, 9L, if(solver == "ldlt") 10L)
      expect_true(all(grepl("^-?[0-9]+\\.[0-9]{3}(e[+-][0-9]+)?$",
        values[numeric_columns])))
    }
    expect_identical(vapply(fields, `[[`, character(1), 7L), rep("-", length(rows)))
    expect_equal(as.numeric(vapply(fields, `[[`, character(1), 8L)), rep(1, length(rows)))
    expect_equal(as.numeric(vapply(fields, `[[`, character(1), 9L)), rep(1, length(rows)))
    expect_false(any(grepl(paste0("local fallback|trust scaling applied|update rejected|",
      "REML likelihood decreased|Global likelihood line search accepted"), output)))
    expect_equal(verbose_fit$monitor, quiet_fit$monitor)
    expect_equal(verbose_fit$llik, quiet_fit$llik)
  }
})

test_that("Henderson progress rows collect failed structures and applied trust scaling", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(6L))),
    genotype=factor(paste0("G", seq_len(24L)))
  )
  trial_effect <- c(45, -30, -12, 52, 18, -20)
  genotype_effect <- sin(seq_len(24L) / 3) * 8
  data$BLUEs <- 100 + trial_effect[data$trial] +
    genotype_effect[data$genotype] + cos(seq_len(nrow(data)) / 5) * 3

  output <- capture.output(fit <- mmes(
    BLUEs~trial, random=~vsm(fam(trial, 2), ism(genotype)), rcov=~units,
    data=data, nIters=12, solver="ldlt", verbose=TRUE, dateWarning=FALSE,
    emWeight=c(exp(seq(log(1), log(0.03), length.out=8L)), rep(0.03, 4L))
  ))
  rows <- grep("^ +[0-9]+ +", output, value=TRUE)
  fields <- strsplit(trimws(rows), " +")
  expect_length(rows, ncol(fit$llik))
  expect_true(all(lengths(fields) == 10L))
  expect_gte(length(rows), 10L)
  header <- grep("^ *iteration", output, value=TRUE)
  column_ends <- function(line){
    positions <- gregexpr("\\S+", line)[[1L]]
    as.integer(positions + attr(positions, "match.length") - 1L)
  }
  for(row in rows){
    expect_identical(column_ends(row), column_ends(header))
  }
  not_pd <- vapply(fields, `[[`, character(1), 7L)
  alpha <- as.numeric(vapply(fields, `[[`, character(1), 8L))
  expect_true(all(grepl("^(-|[0-9]+(,[0-9]+)*)$", not_pd)))
  expect_true(any(not_pd == "-"))
  expect_true(any(grepl("^[0-9]+$", not_pd)))
  expect_true(any(not_pd == "1,2"))
  expect_true(all(is.finite(alpha) & alpha >= 0 & alpha <= 1))
  expect_true(any(alpha < 1))
  line_alpha <- as.numeric(vapply(fields, `[[`, character(1), 9L))
  expect_true(all(is.finite(line_alpha) & line_alpha >= 0 & line_alpha <= 1))
  expect_equal(line_alpha[1L], 1)
  expect_identical(any(line_alpha < 1), fit$engineDiagnostics$lineSearchHalvings > 0L)
  expect_false(any(grepl(paste0("local fallback|trust scaling applied|update rejected|",
    "REML likelihood decreased|Global likelihood line search accepted"), output)))
})