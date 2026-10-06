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

test_that("compound-symmetry Woodbury matches its generic covariance fallback", {
  set.seed(132)
  data <- expand.grid(trial=factor(seq_len(8L)), genotype=factor(seq_len(12L)))
  data$BLUEs <- rnorm(nrow(data))
  precision <- Matrix::Diagonal(12L)
  dimnames(precision) <- list(levels(data$genotype), levels(data$genotype))
  attr(precision, "inverse") <- TRUE

  for(rho in c(0.3, -0.1)){
    factor <- csm(data$trial, rho=rho)
    factor$covFactor$free[] <- FALSE
    covariance <- matrix(rho, 8L, 8L)
    diag(covariance) <- 1
    generic <- ownm(data$trial, K=covariance)
    fits <- lapply(list(factor, generic), function(shape){
      mmes(BLUEs~trial, random=~vsm(shape, ism(genotype), Gu=precision),
        rcov=~units, data=data, nIters=4, solver="ldlt", verbose=FALSE,
        dateWarning=FALSE)
    })
    expect_equal(as.numeric(fits[[1L]]$bu), as.numeric(fits[[2L]]$bu),
                 tolerance=1e-7)
    expect_equal(as.numeric(fits[[1L]]$llik), as.numeric(fits[[2L]]$llik),
                 tolerance=1e-7)
    expect_equal(unname(fits[[1L]]$theta[[1L]]),
                 unname(fits[[2L]]$theta[[1L]]), tolerance=1e-7)
  }
})

test_that("stationary AR precision matches dense covariance factors", {
  set.seed(93)
  q <- 7L
  data <- expand.grid(trial=factor(seq_len(q)), genotype=factor(seq_len(8L)))
  data$BLUEs <- rnorm(nrow(data))
  precision <- Matrix::Diagonal(8L)
  dimnames(precision) <- list(levels(data$genotype), levels(data$genotype))
  attr(precision, "inverse") <- TRUE
  cases <- list(c(order=1L, heterogeneous=TRUE),
                c(order=2L, heterogeneous=FALSE),
                c(order=2L, heterogeneous=TRUE),
                c(order=3L, heterogeneous=FALSE),
                c(order=3L, heterogeneous=TRUE))

  for(case in cases){
    order <- as.integer(case[["order"]])
    heterogeneous <- as.logical(case[["heterogeneous"]])
    fixed <- rep(TRUE, order + if(heterogeneous) q - 1L else 0L)
    shape <- if(order == 1L){
      ar1m(data$trial, variance="heterogeneous", fixed=fixed)
    }else if(order == 2L){
      ar2m(data$trial, variance=if(heterogeneous) "heterogeneous" else "homogeneous",
           fixed=fixed)
    }else{
      ar3m(data$trial, variance=if(heterogeneous) "heterogeneous" else "homogeneous",
           fixed=fixed)
    }
    factor <- shape$covFactor
    correlation <- matrix(
      sommer:::.ar_correlation_from_pacf(tanh(factor$par[seq_len(order)]), q),
      q, q
    )
    if(heterogeneous){
      standardDeviation <- c(1, exp(factor$par[order + seq_len(q - 1L)] / 2))
      covariance <- diag(standardDeviation) %*% correlation %*% diag(standardDeviation)
    }else{
      covariance <- correlation
    }
    generic <- ownm(data$trial, K=covariance)
    fits <- lapply(list(shape, generic), function(term){
      mmes(BLUEs~trial, random=~vsm(term, ism(genotype), Gu=precision),
        rcov=~units, data=data, nIters=4, solver="ldlt", verbose=FALSE,
        dateWarning=FALSE)
    })
    expect_equal(as.numeric(fits[[1L]]$bu), as.numeric(fits[[2L]]$bu),
                 tolerance=1e-7, info=paste("AR order", order, "heterogeneous", heterogeneous))
    expect_equal(as.numeric(fits[[1L]]$llik), as.numeric(fits[[2L]]$llik),
                 tolerance=1e-7, info=paste("AR order", order, "heterogeneous", heterogeneous))
  }
})

test_that("compound-symmetry and AR precision are shared with residual Kronecker blocks", {
  set.seed(28)
  q <- 6L
  data <- expand.grid(trial=factor(seq_len(q)), id=factor(seq_len(8L)))
  data$y <- rnorm(nrow(data))
  shapes <- list(list(shape=csm(data$trial, rho=0.25), rho=0.25),
                 list(shape=csm(data$trial, rho=-0.1), rho=-0.1),
                 list(shape=ar2m(data$trial, fixed=rep(TRUE, 2L)), rho=NULL))

  for(item in shapes){
    shape <- item$shape
    factor <- shape$covFactor
    covariance <- if(identical(factor$representation$kind, "compound_symmetry")){
      matrix(item$rho, q, q)
    }else{
      matrix(sommer:::.ar_correlation_from_pacf(tanh(factor$par), q), q, q)
    }
    diag(covariance) <- 1
    factor$free[] <- FALSE
    shape$covFactor <- factor
    generic <- ownm(data$trial, K=covariance)
    fits <- lapply(list(shape, generic), function(term){
      mmes(y~trial, random=~id, rcov=~vsm(term, ism(id)), data=data,
        nIters=4, solver="ldlt", verbose=FALSE, dateWarning=FALSE)
    })
    expect_equal(as.numeric(fits[[1L]]$llik), as.numeric(fits[[2L]]$llik),
                 tolerance=1e-7, info=factor$model)
    expect_equal(as.numeric(fits[[1L]]$bu), as.numeric(fits[[2L]]$bu),
                 tolerance=1e-7, info=factor$model)
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
