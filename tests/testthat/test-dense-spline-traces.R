.dense_spline_data <- function(){
  dataset <- expand.grid(col=seq_len(7L), row=seq_len(8L))
  dataset$geno <- factor(rep(seq_len(8L), length.out=nrow(dataset)))
  dataset$section <- factor(rep(c("A", "B"), length.out=nrow(dataset)))
  dataset$y <- 3 + sin(dataset$col) + cos(dataset$row / 2) +
    as.integer(dataset$geno) / 5 + sin(seq_len(nrow(dataset)) * 3)
  spl2Dmats("col", "row", dataset, nsegments=c(5L, 6L))
}

test_that("dense spline residual traces agree with observation-space fitting", {
  M <- .dense_spline_data()
  random <- ~vsm(ism(M$fC)) + vsm(ism(M$fR)) +
    vsm(ism(M$fC.R)) + vsm(ism(M$C.fR)) + vsm(ism(M$fC.fR))
  for(objective in c("reml", "ml")){
    for(sectioned in c(FALSE, TRUE)){
      residual <- if(sectioned) ~vsm(dsm(section), ism(units)) else ~units
      weights <- if(sectioned) Matrix::Diagonal(nrow(M$data),
        x=seq(0.7, 1.8, length.out=nrow(M$data))) else NULL
      fit <- function(formulation, solver="ldlt"){
        args <- list(fixed=y~geno, random=random, rcov=residual, data=M$data,
          control=mmes_control(formulation=formulation,
            task=list(objective=objective),
            solve=if(formulation == "direct") list() else list(method=solver),
            score=if(formulation == "direct") list() else
              list(denseCrossproduct="force"),
            update=list(blendWeight=0),
            convergence=list(maxIters=2L),
            diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
        if(!is.null(weights)) args$W <- weights
        do.call(mmes, args)
      }
      reference <- fit("direct")
      for(solver in c("ldlt", "cholmod")){
        candidate <- fit("henderson", solver)
        expect_gt(candidate$engineDiagnostics$denseScoreCrossproducts, 0L)
        expect_equal(unlist(candidate$covPar), unlist(reference$covPar),
          tolerance=1e-6)
        expect_equal(as.numeric(candidate$bu), as.numeric(reference$bu),
          tolerance=1e-6)
        expect_equal(diff(as.numeric(candidate$llik)),
          diff(as.numeric(reference$llik)), tolerance=1e-6)
        if(!sectioned){
          kernels <- c(lapply(M[c("fC", "fR", "fC.R", "C.fR", "fC.fR")],
            tcrossprod), list(diag(nrow(M$data))))
          derivatives <- Map(`*`, kernels, as.numeric(unlist(candidate$covPar)))
          Vi <- chol2inv(chol(Reduce(`+`, derivatives)))
          X <- model.matrix(~geno, M$data)
          VX <- Vi %*% X
          P <- Vi - VX %*% solve(crossprod(X, VX), t(VX))
          Py <- P %*% M$data$y
          traceOperator <- if(objective == "reml") P else Vi
          expected <- vapply(derivatives, function(derivative)
            sum(traceOperator * derivative) - sum(Py * (derivative %*% Py)),
            numeric(1))
          expect_equal(as.numeric(candidate$dLu), unname(expected),
            tolerance=1e-6)
        }
      }
    }
  }
})

test_that("dense score policy and weighted cross-products obey their boundaries", {
  .native_run("test-mme-dense-score", c("mme_dense_score.h", "mme_control.h"), eigen=TRUE)
})

test_that("auto dense score dispatch is structural and off preserves the legacy path", {
  dataset <- data.frame(y=sin(seq_len(200L)), id=factor(rep(1:20, 10)))
  design <- sin(outer(seq_len(200L), seq_len(80L), `*`))
  fit <- function(mode, dense=TRUE, budget=1000, policy="speed"){
    mmes(y~1, random=if(dense) ~vsm(ism(design)) else ~id, data=dataset,
      control=mmes_control(score=list(denseCrossproduct=mode),
        execution=list(denseMemoryMB=budget, memoryPolicy=policy),
        update=list(blendWeight=0), convergence=list(maxIters=1L),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  automatic <- fit("auto")
  forced <- fit("force")
  legacy <- fit("off")
  expect_gt(automatic$engineDiagnostics$denseScoreCrossproducts, 0L)
  expect_identical(legacy$engineDiagnostics$denseScoreCrossproducts, 0L)
  expect_gt(legacy$engineDiagnostics$legacyScoreCrossproducts, 0L)
  expect_equal(automatic[c("dLu", "InfMat", "bu", "covPar")],
    legacy[c("dLu", "InfMat", "bu", "covPar")], tolerance=1e-8)
  expect_equal(automatic[c("dLu", "InfMat", "bu", "covPar")],
    forced[c("dLu", "InfMat", "bu", "covPar")], tolerance=1e-8)
  sparse <- fit("auto", FALSE)
  expect_identical(sparse$engineDiagnostics$denseScoreCrossproducts, 0L)
  expect_gt(sparse$engineDiagnostics$legacyScoreCrossproducts, 0L)
  limited <- fit("auto", budget=0.1, policy="limit")
  expect_identical(limited$engineDiagnostics$denseScoreCrossproducts, 0L)
  expect_gt(limited$engineDiagnostics$denseScoreMemoryFallbacks, 0L)
  expect_equal(limited[c("dLu", "InfMat", "bu", "covPar")],
    legacy[c("dLu", "InfMat", "bu", "covPar")], tolerance=1e-8)
  forcedLimited <- fit("force", budget=0.1, policy="limit")
  expect_identical(forcedLimited$engineDiagnostics$denseScoreCrossproducts, 0L)
  expect_gt(forcedLimited$engineDiagnostics$denseScoreMemoryFallbacks, 0L)
  expect_equal(forcedLimited[c("dLu", "InfMat", "bu", "covPar")],
    legacy[c("dLu", "InfMat", "bu", "covPar")], tolerance=1e-8)
  expect_error(fit("force", budget=0.1), "Dense-design score cross-products.*denseMemoryMB")
})

test_that("forced dense score requests reject unsupported equation domains", {
  expect_error(mmes_control(formulation="direct",
    score=list(denseCrossproduct="force")), "direct")
  expect_error(mmes_control(solve=list(method="pcg"),
    score=list(denseCrossproduct="force")), "deterministic normal Henderson")
  expect_error(mmes_control(representation=list(equations="saddle-point"),
    score=list(denseCrossproduct="force")), "deterministic normal Henderson")
  expect_error(mmes_control(representation=list(compression="groups"),
    score=list(denseCrossproduct="force")), "without compression")
  expect_error(mmes_control(score=list(denseCrossproduct="typo")), "Experimental toggle")
  expect_identical(mmes_control()$score$denseCrossproduct, "auto")
  dataset <- data.frame(y=rep(c(0, 1, 2), 4L), id=factor(rep(1:4, each=3L)))
  control <- mmes_control(score=list(denseCrossproduct="force"),
    convergence=list(maxIters=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE))
  expect_error(mmes(y~1, random=~id, family=poisson(), data=dataset,
    control=control), "Gaussian identity-link")
  weights <- diag(12L)
  weights[1L, 2L] <- weights[2L, 1L] <- 0.1
  expect_error(mmes(y~1, random=~id, data=dataset, W=weights, control=control),
    "Forced dense score cross-products require diagonal")
})

test_that("established covariance paths retain their original score assembly", {
  for(structure in c("ar2", "us", "fa")){
    workload <- generate_mme_workload(structure=structure, n_ind=6L,
      n_levels=4L, seed=814L)
    for(solver in c("ldlt", "cholmod")){
      fit <- function(mode){
        mmes(workload$fixed, random=workload$random, data=workload$data,
          control=mmes_control(solve=list(method=solver),
            score=list(denseCrossproduct=mode),
            convergence=list(maxIters=2L),
            diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
      }
      automatic <- fit("auto")
      legacy <- fit("off")
      expect_equal(automatic$engineDiagnostics$denseScoreCrossproducts, 0)
      expect_equal(automatic[c("dLu", "InfMat", "bu", "covPar", "llik")],
        legacy[c("dLu", "InfMat", "bu", "covPar", "llik")], tolerance=1e-10)
      for(name in c("blockSchurActive", "blockChainActive", "blockBandedActive",
          "factorSchurActive", "residualEngine")){
        expect_true(name %in% names(automatic$engineDiagnostics))
        expect_identical(automatic$engineDiagnostics[[name]],
          legacy$engineDiagnostics[[name]])
      }
    }
  }
})

test_that("dense score selection is visible through AI recovery and native progress", {
  dataset <- data.frame(y=sin(seq_len(200L)))
  design <- sin(outer(seq_len(200L), seq_len(80L), `*`))
  for(recovery in c(FALSE, TRUE)){
    text <- capture.output(fit <- mmes(y~1, random=~vsm(ism(design)),
      data=dataset, control=mmes_control(
        update=if(recovery) list() else list(blendWeight=0),
        convergence=list(maxIters=2L),
        diagnostics=list(verbose=TRUE, dateWarning=FALSE))))
    expect_match(paste(text, collapse="\n"), "score-crossproduct=dense-design",
      fixed=TRUE)
    expect_gte(fit$engineDiagnostics$denseScoreCrossproducts, length(fit$llik))
    expect_gte(fit$engineDiagnostics$denseScoreRowProductsAvoided, 2e6)
  }
})
