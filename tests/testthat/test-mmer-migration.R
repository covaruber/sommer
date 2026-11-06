legacy_migration <- function(fixed=y~1, random=NULL, rcov=~units,
                             modelCall=quote(mmer(y~1, data=dat))){
  sommer:::.mmer_mmes_message(fixed, random, rcov, modelCall)
}

migration_call <- function(text){
  start <- regexpr("mmes\\(", text)[1L]
  end <- regexpr("\nplease switch!", text)[1L]
  parse(text=substr(text, start, end - 1L))[[1L]]
}

test_that("migration prints fixed, random and default residual formulas", {
  text <- legacy_migration(y~env, ~id + block)
  expect_match(text, "Equivalent mmes model", fixed=TRUE)
  expect_match(text, "please switch!", fixed=TRUE)
  model <- migration_call(text)
  expect_identical(model$fixed, quote(y~env))
  expect_identical(model$random, quote(~vsm(ism(id)) + vsm(ism(block))))
  expect_identical(model$rcov, quote(~vsm(ism(units))))
  expect_identical(model$data, quote(dat))
  expect_null(migration_call(legacy_migration())$random)
})

test_that("migration translates covariance calls without evaluating matrices", {
  text <- legacy_migration(random=
    ~sommer::vsr(usr(env), id, Gu=unavailable_matrix) +
      vs(dsr(env), id) + vsr(atr(env, "a"), id),
    rcov=~vsr(dsr(env), units))
  model <- migration_call(text)
  expect_identical(model$random, quote(
    ~vsm(usm(env), ism(id),
         Gu=structure(solve(unavailable_matrix), inverse=TRUE)) +
      vsm(dsm(env), ism(id)) + vsm(atm(env, "a"), ism(id))))
  expect_identical(model$rcov, quote(~vsm(dsm(env), ism(units))))
  expect_identical(migration_call(legacy_migration(random=~vsr(id, Gu=NULL)))$random,
                   quote(~vsm(ism(id))))
})

test_that("formula objects and data expressions are preserved", {
  fixed <- y~I(x + 1)
  random <- ~id
  text <- legacy_migration(fixed, random, modelCall=
    quote(mmer(fixed_model, random=random_model, data=subset(dat, keep))))
  model <- migration_call(text)
  expect_identical(model$fixed, quote(y~I(x + 1)))
  expect_identical(model$data, quote(subset(dat, keep)))
  expect_identical(fixed, y~I(x + 1))
  expect_identical(random, ~id)
})

test_that("non-equivalent legacy options are explicitly flagged", {
  text <- legacy_migration(random=~vsr(csr(env, custom), id, Gtc=constraint),
    modelCall=quote(mmer(y~1, data=dat, weights=w, init=initial)))
  expect_match(text, "manual adjustments required", fixed=TRUE)
  expect_match(text, "Custom covariance", fixed=TRUE)
  expect_match(text, "Gti/Gtc", fixed=TRUE)
  expect_match(text, "init/constraints", fixed=TRUE)
  expect_match(text, "weights/W", fixed=TRUE)
  expect_false(grepl("Equivalent mmes model", text, fixed=TRUE))
})

test_that("multivariate guidance uses stacked traits and explicit record pairing", {
  text <- legacy_migration(cbind(y, z)~env, ~id)
  expect_match(text, 'stackTraits(dat, traits = c("y", "z"))', fixed=TRUE)
  expect_match(text, "fixed = value ~ 0 + trait + trait:(env)", fixed=TRUE)
  expect_match(text, "vsm(usm(trait), ism(id))", fixed=TRUE)
  expect_match(text, "vsm(usm(trait), ism(record))", fixed=TRUE)
  expect_match(text, "missing-response filtering", fixed=TRUE)
  intercept <- legacy_migration(cbind(y, z)~1)
  expect_match(intercept, "fixed = value ~ 0 + trait", fixed=TRUE)
  expect_false(grepl("trait:()", intercept, fixed=TRUE))
  spaced <- legacy_migration(cbind(`trait one`, `trait two`)~1)
  expect_match(spaced, 'traits = c("trait one", "trait two")', fixed=TRUE)
})

test_that("printed equivalent can prepare a matching covariance model", {
  dat <- expand.grid(id=factor(letters[1:4]), env=factor(c("a", "b")),
                     rep=seq_len(3))
  dat$y <- sin(seq_len(nrow(dat)))
  K <- diag(4) + 0.1
  dimnames(K) <- list(levels(dat$id), levels(dat$id))
  model <- migration_call(legacy_migration(y~env,
    ~vsr(dsr(env), id, Gu=K), ~vsr(dsr(env), units)))
  model$control <- quote(mmes_control(output=list(setup=TRUE),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  setup <- eval(model)
  legacy <- suppressMessages(mmer(y~env,
    random=~vsr(dsr(env), id, Gu=K), rcov=~vsr(dsr(env), units),
    data=dat, returnParam=TRUE, verbose=FALSE))
  expect_equal(as.matrix(setup$X), do.call(cbind, legacy$X), ignore_attr=TRUE)
  expect_equal(as.numeric(setup$yvar), as.numeric(legacy$yvar))
  modernKernel <- Reduce(`+`, Map(function(Z, Ai){
    as.matrix(Z) %*% solve(as.matrix(Ai), t(as.matrix(Z)))
  }, setup$Z, setup$Ai))
  legacyKernel <- Reduce(`+`, Map(function(Z, K) Z %*% K %*% t(Z),
                                 legacy$Z, legacy$K))
  expect_equal(modernKernel, as.matrix(legacyKernel), ignore_attr=TRUE,
               tolerance=1e-10)
})

test_that("mmer emits guidance after fitting even with verbose disabled", {
  dat <- data.frame(y=sin(seq_len(12)), id=factor(rep(letters[1:4], 3)))
  fitted <- FALSE
  reshape <- sommer:::reshape_mmer
  local_mocked_bindings(reshape_mmer=function(...){
    result <- reshape(...)
    fitted <<- TRUE
    result
  }, .package="sommer")
  messages <- character()
  fit <- withCallingHandlers(
    mmer(y~1, random=~id, data=dat, nIters=2, verbose=FALSE),
    message=function(m){
      expect_true(fitted)
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(messages[1L], "deprecated")
  expect_match(messages[2L], "Equivalent mmes model", fixed=TRUE)
  expect_match(messages[2L], "please switch!", fixed=TRUE)
  expect_s3_class(fit, "mmer")
})
