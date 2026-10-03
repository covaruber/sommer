test_that("rank-two factor-analytic updates remain finite and positive definite", {
  data <- expand.grid(
    trial=factor(paste0("C", seq_len(6L))),
    genotype=factor(paste0("G", seq_len(24L)))
  )

  trial_effect <- c(45, -30, -12, 52, 18, -20)
  genotype_effect <- sin(seq_len(24L) / 3) * 8
  interaction <- cos(seq_len(nrow(data)) / 5) * 3
  data$BLUEs <- 100 + trial_effect[data$trial] +
    genotype_effect[data$genotype] + interaction

  fit <- mmes(
    BLUEs ~ trial,
    random=~vsm(fam(trial, 2), ism(genotype)),
    rcov=~units,
    data=data,
    nIters=20,
    verbose=FALSE
  )

  expect_s3_class(fit, "mmes")
  expect_true(all(is.finite(fit$monitor)))
  expect_true(all(vapply(fit$covPar, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(vapply(
    fit$theta,
    function(x) min(eigen(x, symmetric=TRUE, only.values=TRUE)$values) > 0,
    logical(1)
  )))
})

test_that("FA/RR Woodbury and near-boundary fallback match generic covariance factors", {
  set.seed(314)
  data <- expand.grid(trial=factor(seq_len(8L)), genotype=factor(seq_len(12L)))
  data$BLUEs <- rnorm(nrow(data))
  loadings <- matrix(runif(16L, -0.3, 0.3), 8L, 2L)
  loadings[1L, 2L] <- 0
  loadings[1L, 1L] <- loadings[2L, 2L] <- 0.6
  precision <- Matrix::bandSparse(12L, k=c(-1L, 0L, 1L),
    diagonals=list(rep(-0.1, 11L), rep(1, 12L), rep(-0.1, 11L)))
  dimnames(precision) <- list(levels(data$genotype), levels(data$genotype))
  attr(precision, "inverse") <- TRUE
  for(model in c("fa", "rr", "boundary")){
    specific <- if(model == "rr") rep(1, 8L) else rep(0.7, 8L)
    if(model == "boundary") specific[8L] <- 1e-10
    factor <- if(model == "rr") rrm(data$trial, 2L, loadings=loadings) else
      fam(data$trial, 2L, loadings=loadings, specific=specific)
    expect_identical(factor$covFactor$precision$backend, "woodbury")
    expect_identical(factor$covFactor$precision$kind, if(model == "rr") "rr" else "fa")
    factor$covFactor$free[] <- FALSE
    generic <- ownm(data$trial, K=tcrossprod(loadings) + diag(specific))
    generic$covFactor$model <- "fa"
    for(reml in c(TRUE, FALSE)){
      fits <- lapply(list(factor, generic), function(shape){
        mmes(BLUEs~trial, random=~vsm(shape, ism(genotype), Gu=precision),
          rcov=~vsm(dsm(trial), ism(units)), data=data, nIters=4,
          REML=reml, solver="ldlt", verbose=FALSE, dateWarning=FALSE)
      })
      expect_equal(as.numeric(fits[[1L]]$bu), as.numeric(fits[[2L]]$bu), tolerance=1e-7)
      expect_equal(as.numeric(fits[[1L]]$llik), as.numeric(fits[[2L]]$llik), tolerance=1e-7)
      expect_equal(unname(fits[[1L]]$theta[[1L]]), unname(fits[[2L]]$theta[[1L]]), tolerance=1e-7)
    }
  }
})

test_that("analytic FA/RR derivatives match centered finite differences", {
  set.seed(271)
  data <- expand.grid(trial=factor(seq_len(5L)), genotype=factor(seq_len(10L)))
  data$y <- rnorm(nrow(data))
  compare_model <- function(model){
    structured <- if(model == "fa") fam(data$trial, 2L) else rrm(data$trial, 2L)
    factor <- structured$covFactor
    rows <- if(model == "fa") factor$fa_row else factor$rr_row
    cols <- if(model == "fa") factor$fa_col else factor$rr_col
    diagonal <- if(model == "fa") factor$fa_diag else factor$rr_diag
    nload <- length(rows)
    order <- factor$order
    covariance <- function(par){
      loading <- matrix(0, nlevels(data$trial), order)
      for(index in seq_len(nload)){
        loading[rows[index], cols[index]] <- if(diagonal[index]) exp(par[index]) else par[index]
      }
      specific <- if(model == "fa") c(1, exp(par[nload + seq_len(nlevels(data$trial)-1L)])) else
        rep(1, nlevels(data$trial))
      value <- tcrossprod(loading) + diag(specific)
      value / value[1L, 1L]
    }
    fit <- function(shape){
      mmes(y~1, random=~vsm(shape, ism(genotype)), rcov=~units,
        data=data, nIters=12, tolParConvLL=1e-8,
        verbose=FALSE, dateWarning=FALSE)
    }
    numeric <- structured
    numeric$covFactor$derivative <- list(backend="numeric", rel_step=1e-6)
    list(native=fit(structured), numeric=fit(numeric))
  }
  for(model in c("fa", "rr")){
    fit <- compare_model(model)
    expect_equal(unlist(fit$native$theta), unlist(fit$numeric$theta), tolerance=2e-3)
    expect_equal(as.numeric(fit$native$llik), as.numeric(fit$numeric$llik), tolerance=2e-3)
    expect_equal(as.numeric(fit$native$bu), as.numeric(fit$numeric$bu), tolerance=2e-3)
  }
})
