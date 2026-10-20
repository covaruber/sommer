test_that("compression controls validate in representation module", {
  expect_identical(mmes_control()$representation$compression, "off")
  for(value in list(NA_character_, "", TRUE, character(), c("off", "groups"), 1)){
    expect_error(mmes_control(representation=list(compression=value)),
      "compression")
  }
  for(solver in c("auto", "ldlt", "cholmod", "pcg")){
    control <- mmes_control(solve=list(method=solver),
      representation=list(compression="CLASSES"))
    expect_identical(control$representation$compression, "classes")
  }
  expect_error(mmes_control(representation=list(unknownCompression="groups")),
    "Unknown representation settings")
})

.mme_compression_data <- function(){
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(seq_len(3L)),
    replicate=seq_len(4L))
  dataset$response <- sin(seq_len(nrow(dataset)) * 0.7) +
    as.integer(dataset$id) / 4 + as.integer(dataset$environment) / 3
  dataset$x <- as.integer(dataset$id) / 8
  dataset
}

.mme_compression_fit <- function(solver, method, reml=TRUE, inverse=0L,
                                 weights=NULL, iterations=1L, relationship=NULL,
                                 acceleration="none", dFirst="auto"){
  dataset <- .mme_compression_data()
  arguments <- list(fixed=response~environment+x,
    random=if(is.null(relationship)) ~id else ~vsm(ism(id), Gu=relationship),
    rcov=~vsm(dsm(environment), ism(units)), data=dataset,
    control=mmes_control(solve=c(list(method=solver, dFirstML=dFirst),
        if(solver == "pcg") list(tol=1e-11)),
      task=list(objective=if(reml) "reml" else "ml"),
      convergence=list(maxIters=iterations),
      output=list(inverse=inverse),
      update=list(blendWeight=0),
      representation=list(compression=method),
      proposal=list(extrapolation=if(acceleration == "aitken") "aitken" else "none"),
      curvature=list(method=if(acceleration == "hybrid-ai") "hybrid-ai" else "exact-ai"),
      score=if(solver == "pcg") list(traceProbes=12L) else list(),
      logdet=if(solver == "pcg") list(lanczosSteps=40L) else list(),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  if(!is.null(weights)) arguments$W <- weights
  suppressMessages(do.call(mmes, arguments))
}

test_that("all deterministic solvers preserve compressed likelihood, exact traces, AI and PEV", {
  weights <- Matrix::Diagonal(x=seq(0.5, 2, length.out=nrow(.mme_compression_data())))
  for(solver in c("ldlt", "cholmod", "auto")){
    for(reml in c(TRUE, FALSE)){
      for(inverse in 0:2){
        for(weight in list(NULL, weights)){
          reference <- .mme_compression_fit(solver, "off", reml, inverse, weight)
          for(method in c("groups", "classes")){
            candidate <- .mme_compression_fit(solver, method, reml, inverse, weight)
            for(field in c("bu", "C", "Ci", "llik", "dLu", "InfMat", "theta",
                          "theta_se", "monitor", "uPevList", "yWorking")){
              actual <- candidate[[field]]
              expected <- reference[[field]]
              if(field == "Ci" && inverse == 2L){
                actual <- as.matrix(actual)
                expected <- as.matrix(expected)
              }
              expect_equal(actual, expected, tolerance=1e-8,
                info=paste(solver, reml, inverse, method, field))
            }
            expect_equal(candidate$W, reference$W)
            expect_equal(candidate$y, reference$y)
            expect_equal(fitted(candidate), fitted(reference), tolerance=1e-8)
            expect_equal(residuals(candidate), residuals(reference), tolerance=1e-8)
            diagnostics <- candidate$engineDiagnostics$compression
            expect_identical(diagnostics$method, method)
            expect_equal(diagnostics$groups, 18)
            expect_equal(diagnostics$observations, 72)
            expect_null(candidate$engineDiagnostics$pcg)
          }
        }
      }
    }
  }
})

test_that("shared compression rejects unsupported models on deterministic solvers", {
  dataset <- .mme_compression_data()
  for(solver in c("ldlt", "cholmod", "auto")){
    control <- mmes_control(solve=list(method=solver),
      representation=list(compression="groups"))
    expect_error(suppressMessages(mmes(response~environment+x, random=~id, rcov=~vsm(usm(environment), ism(units)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(control), list(diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))), "independent residuals")
    expect_error(suppressMessages(mmes(response~environment+x, random=~id, data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(control), list(formulation=if(FALSE) "henderson" else "direct", diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))), "formulation='direct'")
    expect_error(suppressMessages(mmes(response~environment+x, random=~id, data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(control), list(task=list(mode=if(TRUE) "solve" else "fit"), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))), "task\\$mode='solve'")
    expect_error(suppressMessages(mmes(response~environment+x, random=~id, data=dataset, family=poisson(), control=do.call(mmes_control, utils::modifyList(unclass(control), list(diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))), "Gaussian identity")
  }
})

test_that("compression preserves fixed-only ML and nonreplicated continuous designs", {
  dataset <- .mme_compression_data()
  dataset$x <- seq_len(nrow(dataset)) / nrow(dataset)
  for(solver in c("ldlt", "cholmod")){
    fit <- function(method, reml) suppressMessages(mmes(response~environment+x, rcov=~vsm(dsm(environment), ism(units)), data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(representation=list(compression=method))), list(task=list(objective=if(reml) "reml" else "ml"), convergence=list(maxIters=2L), update=list(blendWeight=0), diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE))))
    for(reml in c(TRUE, FALSE)){
      reference <- fit("off", reml)
      for(method in c("groups", "classes")){
        candidate <- fit(method, reml)
        expect_equal(candidate$engineDiagnostics$compression$groups, nrow(dataset))
        for(field in c("bu", "llik", "dLu", "InfMat")){
          expect_equal(candidate[[field]], reference[[field]], tolerance=1e-8)
        }
      }
    }
    expect_error(suppressMessages(mmes(response~environment+x, random=~id, data=dataset, control=do.call(mmes_control, utils::modifyList(unclass(mmes_control(representation=list(compression="classes"), execution=list(denseMemoryMB=0.0001))), list(diagnostics=list(verbose=FALSE, dateWarning=FALSE)), keep.null=TRUE)))), "denseMemoryMB")
  }
})
