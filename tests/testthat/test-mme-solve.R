library(sommer)

.solve_reference_bu <- function(fit, args){
  inputs <- do.call(mmes, c(args, list(returnParam=TRUE, verbose=FALSE, dateWarning=FALSE)))
  cs <- inputs$covStruct
  for(i in seq_along(cs)){
    cs[[i]]$par <- fit$covParWorking[[i]]
    cs[[i]]$free[] <- FALSE
  }
  r <- .Call("_sommer_ai_mme_sp2", PACKAGE="sommer",
             inputs$X, inputs$Z, inputs$Zind, inputs$Ai, inputs$yvar, inputs$W,
             inputs$useH, inputs$residualBlock, inputs$residualIndex,
             1L, 1e-4, 1e-4, inputs$tolParInv, cs, 0, 1, FALSE, 0L, "ldlt",
             1e-8, 0L, 8L, 20L, TRUE, inputs$responsePrepared,
             inputs$preparedMean, inputs$preparedSd, inputs$preparedIntercept)
  as.numeric(r$bu)
}

.solve_check <- function(args, nIters=15L){
  fit <- suppressWarnings(suppressMessages(do.call(mmes,
    c(args, list(nIters=nIters, verbose=FALSE, dateWarning=FALSE)))))
  sol <- suppressMessages(do.call(mmes,
    c(args, list(solveOnly=TRUE, covPar=fit, pcgTol=1e-12, verbose=FALSE, dateWarning=FALSE))))
  ref <- .solve_reference_bu(fit, args)
  expect_s3_class(sol, "mmesSolve")
  expect_true(sol$convergence)
  expect_equal(as.numeric(sol$bu), ref, tolerance=1e-7)
  expect_equal(rownames(sol$bu), rownames(fit$bu))
  invisible(list(fit=fit, sol=sol))
}

test_that("solveOnly matches the engine for diagonal residuals and pedigree-like Gu", {
  data(DT_example)
  .solve_check(list(fixed=Yield~Env, random=~Name, rcov=~units, data=DT_example))

  set.seed(1)
  ids <- levels(DT_example$Name)
  A <- diag(length(ids)) * 1.2
  A[cbind(2:length(ids), 1:(length(ids) - 1L))] <- -0.3
  A[cbind(1:(length(ids) - 1L), 2:length(ids))] <- -0.3
  dimnames(A) <- list(ids, ids)
  Ai <- Matrix::Matrix(A, sparse=TRUE)
  attr(Ai, "inverse") <- TRUE
  .solve_check(list(fixed=Yield~Env, random=~vsm(ism(Name), Gu=Ai), rcov=~units,
                    data=DT_example))
})

test_that("solveOnly matches the engine for structured random and residual covariances", {
  data(DT_example)
  .solve_check(list(fixed=Yield~Env, random=~vsm(usm(Env), ism(Name)),
                    rcov=~vsm(dsm(Env), ism(units)), data=DT_example))

  long <- stackTraits(DT_example, traits=c("Yield", "Weight"))
  long$value <- as.numeric(scale(long$value))
  .solve_check(list(fixed=value~trait+trait:Env, random=~vsm(usm(trait), ism(Name)),
                    rcov=~vsm(usm(trait), ism(record)), data=long))
})

test_that("solveOnly matches the engine with spatial residuals and diagonal weights", {
  data(DT_cpdata)
  DT <- DT_cpdata
  DT$rowf <- factor(DT$Row)
  DT$colf <- factor(DT$Col)
  .solve_check(list(fixed=Yield~1, random=~id,
                    rcov=~vsm(ar1m(rowf), ar1m(colf)), data=DT), nIters=8L)

  w <- seq(0.5, 1.5, length.out=nrow(DT))
  .solve_check(list(fixed=Yield~1, random=~id, rcov=~units, data=DT,
                    W=Matrix::Diagonal(x=w)))
})

test_that("covPar accepts lists and parameter tables and validates them", {
  data(DT_example)
  args <- list(fixed=Yield~Env, random=~vsm(usm(Env), ism(Name)), rcov=~units,
               data=DT_example, verbose=FALSE, dateWarning=FALSE)
  fit <- suppressMessages(do.call(mmes, c(args, list(nIters=15L))))
  first <- suppressMessages(do.call(mmes, c(args, list(solveOnly=TRUE, covPar=fit, pcgTol=1e-12))))

  fromList <- suppressMessages(do.call(mmes,
    c(args, list(solveOnly=TRUE, covPar=fit$covPar, pcgTol=1e-12))))
  expect_equal(fromList$bu, first$bu, tolerance=1e-6)

  fromTable <- suppressMessages(do.call(mmes,
    c(args, list(solveOnly=TRUE, covPar=first$vcParams, pcgTol=1e-12))))
  expect_equal(fromTable$bu, first$bu, tolerance=1e-10)

  bad <- first$vcParams[-1L, ]
  expect_error(suppressMessages(do.call(mmes, c(args, list(solveOnly=TRUE, covPar=bad)))),
               "missing parameters")
  expect_error(suppressMessages(do.call(mmes, c(args, list(solveOnly=TRUE)))),
               "needs every covariance parameter")
  expect_error(suppressMessages(do.call(mmes, c(args, list(covPar=fit)))),
               "only used with solveOnly")

  fixedArgs <- list(fixed=Yield~Env,
                    random=~vsm(ism(Name), sigma2=4, fixedSigma2=TRUE),
                    rcov=~vsm(ism(units), sigma2=8, fixedSigma2=TRUE),
                    data=DT_example, verbose=FALSE, dateWarning=FALSE)
  inFormula <- suppressMessages(do.call(mmes, c(fixedArgs, list(solveOnly=TRUE, pcgTol=1e-12))))
  table <- inFormula$vcParams
  expect_equal(table$value, c(4, 8))
  # BLUP depends only on the variance ratio.
  table$value <- table$value * 3
  scaled <- suppressMessages(do.call(mmes, c(fixedArgs, list(solveOnly=TRUE, covPar=table, pcgTol=1e-12))))
  expect_equal(scaled$bu, inFormula$bu, tolerance=1e-8)

  G0 <- fit$theta[[1]]
  fromMatrix <- suppressMessages(do.call(mmes,
    c(args, list(solveOnly=TRUE, covPar=list(G0, fit$theta[[2]]), pcgTol=1e-12))))
  expect_equal(fromMatrix$bu, first$bu, tolerance=1e-6)
  expect_equal(unname(fromMatrix$theta[[1]]), unname(G0), tolerance=1e-10)
})
