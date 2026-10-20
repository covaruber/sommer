test_that("shared sparse factorization and AI recovery policies preserve their thresholds", {
  .native_run("test-mme-sparse",
    c("mme_sparse.h", "mme_recovery.h", "mme_posterior.h"), eigen=TRUE)
})

test_that("joint and random-only solvers compose with irregular residual solves", {
  dataset <- expand.grid(id=factor(seq_len(6L)), environment=factor(seq_len(3L)),
    replicate=factor(seq_len(2L)))
  dataset$record <- interaction(dataset$id, dataset$replicate)
  dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / 5
  dataset <- dataset[-c(1L, 7L), ]
  weights <- diag(nrow(dataset)) + matrix(0.02, nrow(dataset), nrow(dataset))
  fitModel <- function(solver, reml, dFirst, inverseMode){
    mmes(response~environment, random=~id, rcov=~vsm(usm(environment), ism(record)), W=weights,
      data=dataset, control=mmes_control(solve=list(method=solver,
          blockSchur="off", dFirstML=dFirst),
        task=list(objective=if(reml) "reml" else "ml"),
        convergence=list(maxIters=3L),
        output=list(inverse=inverseMode),
        update=list(blendWeight=0),
        diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
  }
  for(reml in c(TRUE, FALSE)){
    for(mode in 0:2){
      reference <- fitModel("ldlt", reml, "off", mode)
      candidate <- fitModel("cholmod", reml, "off", mode)
      for(field in c("llik", "theta", "covPar", "bu", "InfMat", "dLu", "uPevList")){
        expect_equal(candidate[[field]], reference[[field]], tolerance=1e-7)
      }
      if(mode > 0L){
        expect_equal(as.matrix(candidate$Ci), as.matrix(reference$Ci), tolerance=1e-7)
      }
      if(!reml && mode != 1L){
        for(solver in c("ldlt", "cholmod")){
          dFirst <- fitModel(solver, reml, "force", mode)
          expect_true(dFirst$engineDiagnostics$dFirstMLActive)
          for(field in c("llik", "theta", "bu", "InfMat", "dLu")){
            expect_equal(dFirst[[field]], reference[[field]], tolerance=1e-7)
          }
        }
      }
    }
  }
})
