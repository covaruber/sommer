.augmented_data <- function(){
  set.seed(1931)
  dataset <- expand.grid(id=factor(seq_len(8L)), level=factor(seq_len(3L)),
    replicate=factor(seq_len(3L)))
  dataset$unit <- interaction(dataset$id, dataset$replicate)
  dataset$h <- factor(rep(seq_len(5L), length.out=nrow(dataset)))
  dataset$y <- as.integer(dataset$level) / 3 +
    rnorm(8L)[dataset$id] + rnorm(nrow(dataset), sd=0.7)
  dataset
}

.augmented_control <- function(assembly="augmented-working", ...){
  values <- utils::modifyList(list(
    curvature=list(assembly=assembly), convergence=list(maxIters=1L),
    update=list(blendWeight=0),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE, profile="force")),
    list(...), keep.null=TRUE)
  do.call(mmes_control, values)
}

.augmented_native <- function(setup, assembly="augmented-working", tangent=NULL){
  covariance <- setup$covStruct
  attr(covariance, "aiAssembly") <- assembly
  attr(covariance, "directionBlend") <- FALSE
  attr(covariance, "emInformationOnly") <- TRUE
  attr(covariance, "covarianceInformation") <- "none"
  if(!is.null(tangent)) attr(covariance, "aiTangent") <- tangent
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    0, 1, FALSE, solver=setup$solver, reml=setup$REML,
    responsePrepared=setup$responsePrepared, preparedMean=setup$preparedMean,
    preparedSd=setup$preparedSd, preparedIntercept=setup$preparedIntercept)
}

test_that("augmented AI algebra passes independent native checks", {
  .native_run("test-mme-ai-augmented", "mme_ai_augmented.h", eigen=TRUE)
})

test_that("AI assembly is strictly validated and transported", {
  expect_identical(mmes_control()$curvature$assembly, "sensitivity")
  control <- mmes_control(curvature=list(assembly="AUGMENTED-WORKING"))
  expect_identical(control$curvature$assembly, "augmented-working")
  expect_identical(attr(control, "native")$aiAssembly, "augmented-working")
  for(value in list("", NA_character_, NULL, TRUE, c("sensitivity", "augmented-working"))){
    expect_error(mmes_control(curvature=list(assembly=value)), "curvature\\$assembly")
  }
  for(settings in list(list(formulation="direct"),
      list(solve=list(method="pcg")), list(task=list(mode="solve")),
      list(representation=list(equations="conditional"), output=list(covariance="none",
        coefficientMatrix=FALSE)),
      list(representation=list(factorScores="fixed-shape")))){
    expect_error(do.call(mmes_control, c(settings,
      list(curvature=list(assembly="augmented-working")))), "curvature\\$assembly")
  }
  dataset <- .augmented_data()
  control$curvature$assembly <- "typo"
  expect_error(mmes(y~level, random=~id, data=dataset, control=control),
    "curvature\\$assembly")
  control <- .augmented_control()
  expect_error(mmes(y~level, random=~id, data=dataset, family=poisson(),
    control=control), "curvature\\$assembly")
  expect_error(mmes(y~level, random=~id, data=dataset,
    .internal=list(pqlInner=TRUE), control=control), "curvature\\$assembly")
  expect_error(mmes(y~level, random=~vsm(rrm(level, 1L), ism(id)),
    data=dataset, control=control), "curvature\\$assembly")
  for(control in list(mmes_control(formulation="direct"),
      mmes_control(solve=list(method="pcg")))){
    control$curvature$assembly <- "augmented-working"
    expect_error(mmes(y~level, random=~id, data=dataset, control=control),
      "curvature\\$assembly")
  }
  setup <- mmes(y~level, random=~id, data=dataset,
    control=.augmented_control(output=list(setup=TRUE)))
  expect_identical(attr(setup$covStruct, "aiAssembly"), "augmented-working")
  expect_error(.augmented_native(setup, "typo"), "aiAssembly|assembly")
  fit <- .augmented_native(setup)
  expect_identical(fit$engineDiagnostics$aiAssembly$selected, "augmented-working")
  expect_equal(fit$engineDiagnostics$aiAssembly$randomObservationColumns, 0)
  setup$solver <- "pcg"
  expect_error(.augmented_native(setup), "deterministic|PCG|pcg")
  latent <- mmes(y~level, random=~vsm(rrm(level, 1L), ism(id)),
    data=dataset, control=.augmented_control("sensitivity", output=list(setup=TRUE)))
  expect_error(.augmented_native(latent), "latent|reduced.rank")
})

test_that("same-point augmented AI matches sensitivity across covariance structures", {
  dataset <- .augmented_data()
  relationship <- 0.3^abs(outer(seq_len(8L), seq_len(8L), "-"))
  precision <- solve(relationship)
  dimnames(precision) <- rep(list(levels(dataset$id)), 2L)
  attr(precision, "inverse") <- TRUE
  models <- list(
    list(random=~id+h),
    list(random=~vsm(ar1m(level, rho=0.25), ism(id), Gu=precision)+h),
    list(random=~vsm(usm(level), ism(id), Gu=precision)+h),
    list(random=~vsm(fam(level, 1L), ism(id), Gu=precision)+h),
    list(random=~vsm(csm(level, rho=0.2), ism(id))+h),
    list(random=~vsm(ar1m(level), ism(id))+h,
      rcov=~vsm(ar1m(level, rho=-0.2), ism(unit))),
    list(random=~id+h, rcov=~vsm(dsm(level), ism(units))),
    list(rcov=~vsm(dsm(level), ism(units))),
    list(random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE)+h))
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(index in seq_along(models)){
      setup <- do.call(mmes, c(list(fixed=y~level, data=dataset,
        control=.augmented_control(solve=list(method=solver),
          task=list(objective=if(reml) "reml" else "ml"),
          output=list(setup=TRUE))), models[[index]]))
      reference <- .augmented_native(setup, "sensitivity")
      candidate <- .augmented_native(setup)
      for(field in c("InfMat", "dLu", "bu", "llik", "covPar")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8,
          info=paste(index, solver, reml, field))
      }
      expect_equal(candidate$engineDiagnostics$CsymbolicAnalyses,
        reference$engineDiagnostics$CsymbolicAnalyses)
      expect_equal(candidate$engineDiagnostics$profile$counters$numericFactorCalls,
        reference$engineDiagnostics$profile$counters$numericFactorCalls)
      expect_equal(candidate$engineDiagnostics$aiAssembly$randomObservationColumns, 0)
    }
  }
})

test_that("budgeted augmented panels preserve cross-panel information", {
  dataset <- .augmented_data()
  for(solver in c("ldlt", "cholmod")){
    setup <- mmes(y~level, random=~vsm(ar1m(level), ism(id))+h,
      data=dataset, control=.augmented_control(solve=list(method=solver),
        output=list(setup=TRUE)))
    full <- .augmented_native(setup)
    diagnostic <- full$engineDiagnostics$aiAssembly
    directions <- nrow(full$InfMat)
    expect_false(diagnostic$tiled)
    expect_equal(diagnostic$panelColumns, directions)
    panelBytes <- 8 * (8 * length(full$bu) + directions)
    fixedBytes <- diagnostic$workspaceBytes - directions * panelBytes
    budget <- fixedBytes + 2.25 * panelBytes
    attr(setup$covStruct, "denseMemoryMB") <- budget / 1e6
    tiled <- .augmented_native(setup)
    diagnostic <- tiled$engineDiagnostics$aiAssembly
    expect_true(diagnostic$tiled)
    expect_equal(diagnostic$panelColumns, 2)
    expect_equal(diagnostic$solveCalls, ceiling(directions / 2))
    expect_lte(diagnostic$workspaceBytes, budget)
    expect_equal(diagnostic$randomObservationColumns, 0)
    for(field in c("InfMat", "dLu", "bu", "llik", "covPar")){
      expect_equal(tiled[[field]], full[[field]], tolerance=1e-8)
    }
    attr(setup$covStruct, "denseMemoryMB") <-
      (fixedBytes + panelBytes - 16) / 1e6
    expect_error(.augmented_native(setup), "workspace|budget|memory")
  }
})

test_that("augmented AI matches an independent dense derivative-product formula", {
  dataset <- .augmented_data()
  for(reml in c(TRUE, FALSE)){
    setup <- mmes(y~level, random=~id+h, data=dataset,
      control=.augmented_control(output=list(setup=TRUE),
        task=list(objective=if(reml) "reml" else "ml")))
    fit <- .augmented_native(setup)
    designs <- lapply(setup$Z, as.matrix)
    variances <- vapply(fit$covPar, function(value) as.numeric(value)[1L], numeric(1)) /
      fit$Cscale
    derivatives <- c(lapply(seq_along(designs), function(index)
      variances[index] * tcrossprod(designs[[index]])),
      list(tail(variances, 1L) * diag(nrow(dataset))))
    V <- Reduce(`+`, derivatives)
    Vi <- solve(V)
    X <- as.matrix(setup$X)
    P <- Vi - Vi %*% X %*% solve(crossprod(X, Vi %*% X), t(X) %*% Vi)
    a <- P %*% as.numeric(fit$yWorking)
    working <- do.call(cbind, lapply(derivatives, function(value) value %*% a))
    expected <- crossprod(working, P %*% working)
    expect_equal(unname(fit$InfMat), unname(expected), tolerance=1e-9)
    score <- vapply(derivatives, function(value)
      drop(crossprod(a, value %*% a)) - sum((if(reml) P else Vi) * value),
      numeric(1))
    expect_equal(as.numeric(fit$dLu), -score, tolerance=1e-8)
  }
})

test_that("weights, incomplete residual blocks and compressed statistics retain AI", {
  dataset <- .augmented_data()
  diagonal <- diag(seq(0.7, 1.4, length.out=nrow(dataset)))
  correlated <- Matrix::bandSparse(nrow(dataset), k=c(-1L, 0L, 1L),
    diagonals=list(rep(0.03, nrow(dataset)-1L), rep(1.2, nrow(dataset)),
      rep(0.03, nrow(dataset)-1L)))
  for(solver in c("ldlt", "cholmod")) for(weight in list(NULL, diagonal, correlated)){
    for(incomplete in if(is.null(weight)) c(FALSE, TRUE) else FALSE){
      selected <- if(incomplete) setdiff(seq_len(nrow(dataset)), c(2L, 21L)) else
        seq_len(nrow(dataset))
      arguments <- list(fixed=y~level, random=~vsm(ar1m(level), ism(id))+h,
        rcov=~vsm(ar1m(level, rho=0.25), ism(unit)),
        data=dataset[selected, ],
        control=.augmented_control(solve=list(method=solver), output=list(setup=TRUE)))
      if(!is.null(weight)) arguments$W <- weight[selected, selected]
      setup <- do.call(mmes, arguments)
      reference <- .augmented_native(setup, "sensitivity")
      candidate <- .augmented_native(setup)
      for(field in c("InfMat", "dLu", "bu", "llik")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
      }
    }
  }
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    model <- function(compression, assembly){
      mmes(y~level, random=~id+h, rcov=~vsm(dsm(level), ism(units)),
        data=dataset, W=diagonal,
        control=.augmented_control(assembly, solve=list(method=solver),
          representation=list(compression=compression),
          task=list(objective=if(reml) "reml" else "ml")))
    }
    reference <- model("off", "sensitivity")
    for(compression in c("off", "groups", "classes")){
      candidate <- model(compression, "augmented-working")
      for(field in c("InfMat", "dLu", "bu", "llik", "theta_se")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
      }
    }
  }
})

test_that("mixed and empty direction maps preserve full information pullbacks", {
  dataset <- .augmented_data()
  for(solver in c("ldlt", "cholmod")){
    setup <- mmes(y~level, random=~vsm(ar1m(level), ism(id))+h,
      data=dataset, control=.augmented_control(solve=list(method=solver),
        output=list(setup=TRUE)))
    full <- .augmented_native(setup)
    count <- nrow(full$InfMat)
    tangent <- cbind(seq_len(count) / count, rep(c(0, 1), length.out=count))
    for(map in list(tangent, matrix(0, count, 0L), diag(count)[, 1L, drop=FALSE])){
      candidate <- .augmented_native(setup, tangent=map)
      expect_equal(candidate$InfMat, crossprod(map, full$InfMat %*% map),
        tolerance=1e-8)
      expect_equal(candidate$dLu, full$dLu)
    }
  }
})

test_that("grid, sectioned and generic sparse residuals retain augmented information", {
  set.seed(1932)
  dataset <- expand.grid(site=factor(letters[1:2]), row=factor(seq_len(3L)),
    column=factor(seq_len(4L)))
  dataset$id <- factor(rep(seq_len(5L), length.out=nrow(dataset)))
  dataset$y <- rnorm(nrow(dataset))
  for(incomplete in c(FALSE, TRUE)) for(solver in c("ldlt", "cholmod")){
    selected <- if(incomplete) dataset[-c(2L, 17L), ] else dataset
    for(residual in list(
        ~vsm(ar1m(row, rho=0.2), ar1m(column, rho=-0.15), ism(units)),
        ~dsumm(vsm(ar1m(row, rho=0.2), ar1m(column, rho=-0.15), ism(units)),
          by=site),
        ~dsumm(vsm(usm(row), ism(units)), by=site))){
      setup <- mmes(y~site, random=~id, rcov=residual, data=selected,
        control=.augmented_control(output=list(setup=TRUE),
          solve=list(method=solver)))
      candidate <- .augmented_native(setup)
      reference <- .augmented_native(setup, "sensitivity")
      for(field in c("InfMat", "dLu", "llik", "bu")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
      }
    }
  }
})

test_that("saddle and rotated coefficient representations use the same AI contract", {
  dataset <- .augmented_data()
  for(reml in c(TRUE, FALSE)){
    setup <- mmes(y~level, random=~id+h,
      rcov=~vsm(ar1m(level, rho=0.2), ism(unit)), data=dataset,
      control=.augmented_control(output=list(setup=TRUE),
        representation=list(equations="saddle-point"),
        task=list(objective=if(reml) "reml" else "ml")))
    candidate <- .augmented_native(setup)
    reference <- .augmented_native(setup, "sensitivity")
    for(field in c("InfMat", "dLu", "llik", "bu")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
    }
  }
  selected <- droplevels(dataset[dataset$replicate == levels(dataset$replicate)[1L], ])
  precision <- solve(0.4^abs(outer(seq_len(8L), seq_len(8L), "-")))
  dimnames(precision) <- rep(list(levels(selected$id)), 2L)
  attr(precision, "inverse") <- TRUE
  for(solver in c("ldlt", "cholmod")){
    model <- function(assembly){
      mmes(y~level, random=~vsm(ism(id), Gu=precision, rotation=TRUE),
        data=selected, control=.augmented_control(assembly,
          solve=list(method=solver), convergence=list(maxIters=3L)))
    }
    candidate <- model("augmented-working")
    reference <- model("sensitivity")
    for(field in c("InfMat", "dLu", "llik", "bu")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
    }
  }
})

test_that("recursive recovery, profiling and reduced fits preserve assembly selection", {
  dataset <- .augmented_data()
  configurations <- list(
    list(update=list(blendWeight=NULL)),
    list(update=list(blendWeight=0.3)),
    list(update=list(blendWeight=NULL, recovery="pxem")),
    list(update=list(blendWeight=NULL, commonScale="force")),
    list(curvature=list(method="hybrid-ai")),
    list(proposal=list(extrapolation="aitken")),
    list(curvature=list(reduced="force"), update=list(blendWeight=NULL),
      vcc=data.frame(parameter=1:2, group=1L)))
  for(solver in c("ldlt", "cholmod")) for(index in seq_along(configurations)){
    settings <- configurations[[index]]
    constraints <- settings$vcc
    settings$vcc <- NULL
    model <- function(assembly){
      args <- utils::modifyList(list(solve=list(method=solver),
        convergence=list(maxIters=5L)), settings, keep.null=TRUE)
      control <- do.call(.augmented_control, c(list(assembly=assembly), args))
      mmes(y~level, random=~id+h, data=dataset, vcc=constraints, control=control)
    }
    reference <- model("sensitivity")
    candidate <- model("augmented-working")
    for(field in c("bu", "theta", "InfMat", "theta_se", "dLu")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-6,
        info=paste(index, solver, field))
    }
    expect_identical(candidate$engineDiagnostics$aiAssembly$selected,
      "augmented-working")
    expect_gt(candidate$engineDiagnostics$aiAssembly$refreshes, 0)
    expect_gte(candidate$engineDiagnostics$aiAssembly$solveCalls,
      candidate$engineDiagnostics$aiAssembly$refreshes)
    expect_true(candidate$executionPlan$admission$aiAssembly$active)
  }
})

test_that("complete fits and post-fit coefficient uncertainty retain their contracts", {
  dataset <- .augmented_data()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)) for(inverse in 0:2){
    model <- function(assembly, information="full"){
      mmes(y~level, random=~id, data=dataset,
        control=.augmented_control(assembly, solve=list(method=solver),
          task=list(objective=if(reml) "reml" else "ml"),
          output=list(inverse=inverse, covariance=information),
          update=list(blendWeight=NULL),
          convergence=list(maxIters=60L, likelihoodTol=1e-9, parameterTol=1e-8)))
    }
    reference <- model("sensitivity")
    candidate <- model("augmented-working")
    expect_true(candidate$convergence)
    expect_equal(tail(c(candidate$llik), 1L), tail(c(reference$llik), 1L),
      tolerance=1e-7)
    for(field in c("bu", "theta", "theta_se", "uPevList")){
      expect_equal(candidate[[field]], reference[[field]], tolerance=1e-6)
    }
    expect_equal(dim(candidate$C), rep(length(candidate$bu), 2L))
    if(inverse == 2L){
      expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-6)
    }
    if(inverse == 0L){
      omitted <- model("augmented-working", "none")
      expect_null(omitted$theta_se)
      expect_null(omitted$InfMat)
      expect_equal(postPEV(candidate, mode=2L)$Ci, postPEV(reference, mode=2L)$Ci,
        tolerance=1e-6)
      expect_equal(predict(candidate, D="id")$pvals,
        predict(reference, D="id")$pvals, tolerance=1e-6)
    }
  }
})

test_that("new evaluations cannot reuse stale numerical borders", {
  dataset <- .augmented_data()
  model <- function(data, assembly, W=NULL){
    arguments <- list(fixed=y~level, random=~id+h, data=data,
      control=.augmented_control(assembly, solve=list(method="cholmod")))
    if(!is.null(W)) arguments$W <- W
    do.call(mmes, arguments)
  }
  for(change in c("response", "design", "weights")){
    changed <- dataset
    weights <- NULL
    if(change == "response") changed$y <- changed$y + sin(seq_len(nrow(changed)))
    if(change == "design") changed$h <- rev(changed$h)
    if(change == "weights") weights <- diag(seq(0.6, 1.4, length.out=nrow(changed)))
    model(dataset, "augmented-working")
    candidate <- model(changed, "augmented-working", weights)
    reference <- model(changed, "sensitivity", weights)
    expect_equal(candidate$InfMat, reference$InfMat, tolerance=1e-8)
    expect_equal(candidate$dLu, reference$dLu, tolerance=1e-8)
  }
})

test_that("fixed-rank reduction and fixed covariance parameters preserve reporting", {
  dataset <- .augmented_data()
  dataset$cell <- interaction(dataset$level, dataset$replicate)
  model <- function(assembly){
    mmes(y~level+cell, random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE),
      rcov=~vsm(ism(units), sigma2=0.8, fixedSigma2=TRUE), data=dataset,
      control=.augmented_control(assembly, representation=list(fixedRank="drop")))
  }
  expect_warning(reference <- model("sensitivity"), "rank reduction")
  expect_warning(candidate <- model("augmented-working"), "rank reduction")
  expect_identical(candidate$convergence, reference$convergence)
  for(field in c("bu", "theta", "InfMat", "theta_se", "fixedRank")){
    expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
  }
  expect_equal(length(candidate$bu), length(reference$bu))
})
