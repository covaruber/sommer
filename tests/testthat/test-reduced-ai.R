.reduced_ai_data <- function(){
  set.seed(331)
  dataset <- expand.grid(id=factor(seq_len(18L)), level=factor(seq_len(3L)),
                         replicate=seq_len(3L))
  dataset$h <- factor(rep(seq_len(9L), length.out=nrow(dataset)))
  dataset$unit <- interaction(dataset$id, dataset$replicate)
  dataset$y <- rnorm(nrow(dataset)) + rep(rnorm(18L), 9L) +
    as.integer(dataset$level) / 3
  dataset
}

.reduced_ai_native <- function(setup, tangent=NULL){
  covariance <- setup$covStruct
  attr(covariance, "directionBlend") <- FALSE
  attr(covariance, "emInformationOnly") <- TRUE
  attr(covariance, "covarianceInformation") <- "none"
  if(!is.null(tangent)) attr(covariance, "aiTangent") <- tangent
  sommer:::ai_mme_sp2(setup$X, setup$Z, setup$Zind, setup$Ai, setup$yvar,
    setup$W, setup$useH, setup$residualBlock, setup$residualIndex,
    1L, setup$tolParConvLL, setup$tolParConvNorm, setup$tolParInv, covariance,
    0, 1, FALSE, solver=setup$solver, reml=setup$REML)
}

test_that("empty and random-only AI tangents preserve likelihoods and information shapes", {
  dataset <- .reduced_ai_data()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(compression in c("off", "groups")){
      setup <- mmes(y~level, random=~id+h, data=dataset,
        control=mmes_control(solve=list(method=solver),
          task=list(objective=if(reml) "reml" else "ml"),
          representation=list(compression=compression),
          output=list(setup=TRUE), diagnostics=list(verbose=FALSE)))
      full <- .reduced_ai_native(setup)
      for(tangent in list(matrix(numeric(), 3L, 0L),
                          matrix(c(1, 0, 0), 3L, 1L))){
        reduced <- .reduced_ai_native(setup, tangent)
        expect_equal(dim(reduced$InfMat), rep(ncol(tangent), 2L))
        if(ncol(tangent) > 0L){
          expect_equal(reduced$InfMat, crossprod(tangent, full$InfMat %*% tangent),
            tolerance=1e-8)
        }
        for(field in c("llik", "bu", "dLu")){
          expect_equal(reduced[[field]], full[[field]], tolerance=1e-10)
        }
      }
    }
  }
})

test_that("reduced AI uses curvature module and revalidates mutated requests", {
  expect_identical(mmes_control()$curvature$reduced, "off")
  expect_identical(mmes_control(curvature=list(reduced="force"))$curvature$reduced,
    "force")
  expect_identical(mmes_control(curvature=list(reduced="auto"))$curvature$reduced,
    "auto")
  expect_error(mmes_control(curvature=list(reduced="invalid")),
    "Experimental toggle")
  control <- mmes_control()
  control$curvature$reduced <- "invalid"
  expect_error(mmes(y~level, random=~id, data=.reduced_ai_data(), control=control),
    "Experimental toggle")
})

test_that("reduced information preserves boundary tangents and regularized trust units", {
  helper <- .native_ai_recovery()
  information <- diag(c(1e-30, 2, 3, 4))
  score <- c(-1e-8, 3, 2, 1)
  basis <- cbind(c(1, 0, 0, 0), c(0, 2, 1, 0))
  for(feasible in list(basis, basis[, 2L, drop=FALSE])){
    probe <- helper$reduced_ai_recovery_probe(information, score, basis,
      feasible, c(1, 0.5, 2, 1))
    expect_equal(probe$projected, crossprod(feasible, information %*% feasible),
      tolerance=1e-12)
    expect_identical(probe$packedReliable, probe$fullReliable)
    expect_equal(probe$packed, probe$full, tolerance=1e-12)
  }
  expect_error(helper$reduced_ai_recovery_probe(information, score, basis,
    diag(4), rep(1, 4)), "outside the reduced information tangent")
})

test_that("packed native AI equals independent contraction of full information", {
  dataset <- .reduced_ai_data()
  models <- list(
    list(random=~id+h, vcc=data.frame(parameter=c(1L, 2L), group=1L, scale=c(2, 1))),
    list(random=~vsm(ar1m(level), ism(id))+vsm(ar1m(level), ism(h)),
      vcc=data.frame(parameter=1:4, group=c(1, 2, 1, 2))),
    list(random=~vsm(ar1m(level), ism(id))+vsm(ar1m(level), ism(h)),
      vcc=data.frame(parameter=c(1L, 5L, 2L, 4L), group=c(1, 1, 2, 2))),
    list(random=~vsm(dsm(level), ism(id))+h,
      vcc=data.frame(parameter=c(1L, 5L), group=1L)),
    list(random=~vsm(usm(level), ism(id))+h,
      vcc=data.frame(parameter=c(1L, 8L), group=1L)),
    list(random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE)+h,
      vcc=data.frame(parameter=1:2, group=1L)),
    list(random=~vsm(ar1m(level), ism(id))+h,
      rcov=~vsm(ar1m(level), ism(unit)),
      vcc=data.frame(parameter=c(1L, 4L, 2L, 5L), group=c(1, 1, 2, 2))),
    list(random=~id+h, vcc=data.frame(parameter=1:2, group=1L),
      W=diag(rep(c(0.75, 1, 1.25), length.out=nrow(dataset))))
  )
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(index in seq_along(models)){
      setup <- do.call(mmes, c(list(fixed=y~level, data=dataset,
        control=mmes_control(solve=list(method=solver),
          task=list(objective=if(reml) "reml" else "ml"),
          output=list(setup=TRUE),
          diagnostics=list(verbose=FALSE))), models[[index]]))
      tangent <- attr(setup$covStruct, "vcmap")$T
      full <- .reduced_ai_native(setup)
      reduced <- .reduced_ai_native(setup, tangent)
      expect_equal(reduced$InfMat, crossprod(tangent, full$InfMat %*% tangent),
        tolerance=1e-8, info=paste(solver, reml, index))
      for(field in c("dLu", "bu", "llik", "covPar")){
        expect_identical(reduced[[field]], full[[field]], info=field)
      }
      diagnostics <- reduced$engineDiagnostics
      expect_equal(diagnostics$nativeAISensitivityColumns, ncol(tangent))
      expect_equal(diagnostics$nativeAISensitivityEntries,
        length(full$bu) * ncol(tangent))
      expect_equal(diagnostics$nativeAIAssemblyEntries, ncol(tangent)^2)
      expect_null(reduced$theta_se)
    }
  }
})

test_that("reduced tied fits preserve optimizer, reporting and coefficient uncertainty", {
  dataset <- .reduced_ai_data()
  for(solver in c("ldlt", "cholmod")) for(reml in c(TRUE, FALSE)){
    for(information in c("full", "none")){
      model <- function(mode){
        mmes(y~level, random=~id+h, data=dataset,
          vcc=data.frame(parameter=1:2, group=1L, scale=c(2, 1)),
          control=mmes_control(output=list(covariance=information, inverse=2L),
            diagnostics=list(profile="force", verbose=FALSE, dateWarning=FALSE),
            curvature=list(reduced=mode),
            convergence=list(maxIters=12L, likelihoodTol=1e-10,
              parameterTol=1e-8),
            task=list(objective=if(reml) "reml" else "ml")))
      }
      full <- model("off")
      reduced <- model("force")
      for(field in c("llik", "theta", "covPar", "bu", "InfMat", "theta_se",
                    "Ci", "uPevList", "convergence")){
        expect_equal(reduced[[field]], full[[field]], tolerance=1e-7, info=field)
      }
      diagnostics <- reduced$engineDiagnostics$reducedAI
      expect_true(diagnostics$active)
      expect_gt(diagnostics$evaluations, 0L)
      expect_equal(diagnostics$sensitivityColumns, 2 * diagnostics$evaluations)
      expect_equal(diagnostics$sensitivityColumnsSkipped, diagnostics$evaluations)
      expect_equal(diagnostics$maxAssemblyEntries, 4)
      expect_false(".aiTangent" %in% names(reduced))
      expect_lte(reduced$engineDiagnostics$recovery$finalProjectedScoreNorm,
        reduced$engineDiagnostics$recovery$scoreTolerance)
      if(information == "full"){
        expect_equal(diagnostics$fullInformationEvaluations, 1L)
        expect_equal(covparams_mmes_se(reduced), covparams_mmes_se(full), tolerance=1e-7)
      }else{
        expect_equal(diagnostics$fullInformationEvaluations, 0L)
        expect_null(reduced$InfMat)
        expect_null(reduced$theta_se)
        expect_error(vpredict(reduced, ~V1/(V1+V3)),
          "refit with output=list\\(covariance='full'\\)")
      }
    }
  }
})

test_that("reduced certification handles iteration limits, fixed residuals and empty tangents", {
  dataset <- .reduced_ai_data()
  for(iterations in c(1L, 2L)){
    for(fixedResidual in c(FALSE, TRUE)){
      for(information in c("full", "none")){
        model <- function(mode){
          mmes(y~level, random=~id+h, data=dataset, rcov=if(fixedResidual) ~vsm(ism(units), sigma2=1, fixedSigma2=TRUE)
              else ~units, vcc=data.frame(parameter=1:2, group=1L),
            control=mmes_control(output=list(covariance=information),
              curvature=list(reduced=mode),
              convergence=list(maxIters=iterations),
              diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
        }
        full <- model("off")
        reduced <- model("force")
        for(field in c("llik", "theta", "bu", "InfMat", "theta_se", "convergence")){
          expect_equal(reduced[[field]], full[[field]], tolerance=1e-8, info=field)
        }
        expect_equal(reduced$engineDiagnostics$reducedAI$fullInformationEvaluations,
          as.integer(information == "full"))
      }
    }
  }
  for(information in c("full", "none")){
    fit <- mmes(y~level, random=~vsm(ism(id), sigma2=1, fixedSigma2=TRUE)+
        vsm(ism(h), sigma2=1, fixedSigma2=TRUE), rcov=~vsm(ism(units), sigma2=1, fixedSigma2=TRUE), data=dataset, vcc=data.frame(parameter=1:2, group=1L), control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(curvature=list(reduced="force"))), list(convergence=list(maxIters=1L), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))
    expect_true(fit$convergence)
    expect_true(all(is.finite(unlist(fit$covPar))))
  }
})

test_that("reduced assembly rejects unsupported force requests and reports auto fallback", {
  dataset <- .reduced_ai_data()
  base <- list(fixed=y~level, random=~id+h, data=dataset,
    vcc=data.frame(parameter=1:2, group=1L))
  control <- mmes_control(convergence=list(maxIters=2L),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE),
    curvature=list(reduced="force"))
  for(settings in list(
      list(control=list(update=list(blendWeight=0))),
      list(control=list(solve=list(method="pcg"))),
      list(control=list(formulation="direct")),
      list(control=list(task=list(mode="solve"))),
      list(family=poisson()),
      list(control=list(curvature=list(method="hybrid-ai"))))){
    override <- settings$control
    settings$control <- NULL
    if(is.null(override)) override <- list()
    merged <- utils::modifyList(unclass(control), override, keep.null=TRUE)
    expect_error(do.call(mmes, c(base, settings, list(control=do.call(mmes_control, merged)))),
      "curvature\\$reduced")
  }
  expect_error(mmes(y~level, random=~id, data=dataset, control=mmes_control(curvature=list(reduced="force"), output=list(setup=TRUE))), "requires vcc")
  fallbackControl <- do.call(mmes_control, utils::modifyList(
    unclass(control), list(curvature=list(reduced="auto"),
      update=list(blendWeight=0)), keep.null=TRUE))
  fallback <- do.call(mmes, c(base, list(control=fallbackControl)))
  expect_false(fallback$engineDiagnostics$reducedAI$active)
  expect_match(fallback$engineDiagnostics$reducedAI$fallbackReason, "legacy optimizer")
})
