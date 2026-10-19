## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>")

## ----identify-runtime---------------------------------------------------------
R.version.string
R.home()
.libPaths()
extSoftVersion()["BLAS"]
sessionInfo()[c("BLAS", "LAPACK")]

## ----install-flexiblas, eval=FALSE--------------------------------------------
# install.packages(c("flexiblas", "parallelly"),
#                  repos="https://cloud.r-project.org")

## ----verify-flexiblas, eval=FALSE---------------------------------------------
# library(flexiblas)
# stopifnot(flexiblas_avail())
# flexiblas_current_backend()
# flexiblas_get_num_threads()
# parallelly::availableCores()

## ----select-flexiblas-backend, eval=FALSE-------------------------------------
# flexiblas::flexiblas_list()
# flexiblas::flexiblas_list_loaded()
# 
# # Only use OPENBLAS if that name appears in the available backend list.
# backend <- flexiblas::flexiblas_load_backend("OPENBLAS")
# flexiblas::flexiblas_switch(backend)
# flexiblas::flexiblas_current_backend()
# flexiblas::flexiblas_get_num_threads()

## ----install-blasctl, eval=FALSE----------------------------------------------
# install.packages(c("RhpcBLASctl", "parallelly"),
#                  repos="https://cloud.r-project.org")

## ----verify-blasctl, eval=FALSE-----------------------------------------------
# library(RhpcBLASctl)
# extSoftVersion()["BLAS"]
# blas_get_num_procs()
# parallelly::availableCores()

## ----verify-sommer, eval=FALSE------------------------------------------------
# library(sommer)
# packageVersion("sommer")
# find.package("sommer")
# stopifnot("control" %in% names(formals(mmes)))

## ----verification-fit---------------------------------------------------------
blas <- unname(extSoftVersion()["BLAS"])
hasFlexiBLAS <- requireNamespace("flexiblas", quietly=TRUE) &&
  isTRUE(tryCatch(flexiblas::flexiblas_avail(), error=function(error) FALSE))
hasDirectControl <- !is.na(blas) &&
  !grepl("flexiblas", blas, ignore.case=TRUE) &&
  grepl("openblas|mkl|acml|goto", blas, ignore.case=TRUE) &&
  requireNamespace("RhpcBLASctl", quietly=TRUE)

if(requireNamespace("sommer", quietly=TRUE) &&
   "control" %in% names(formals(sommer::mmes)) &&
   requireNamespace("parallelly", quietly=TRUE) &&
   (hasFlexiBLAS || hasDirectControl)){
  local({
    suppressPackageStartupMessages(library(sommer))
    getThreads <- if(hasFlexiBLAS) flexiblas::flexiblas_get_num_threads else
      RhpcBLASctl::blas_get_num_procs
    before <- getThreads()
    exampleData <- expand.grid(Name=factor(seq_len(12L)),
      Env=factor(letters[1:3]), replicate=seq_len(2L))
    exampleData$Yield <- as.integer(exampleData$Env) +
      sin(as.integer(exampleData$Name)) + cos(seq_len(nrow(exampleData)))
    for(request in list(1L, "auto")){
      fit <- sommer::mmes(Yield~Env,
  random=~Name,
  rcov=~units,
  data=exampleData,
  control = mmes_control(solve = list(method = "cholmod"),
  convergence = list(maxIters = 2),
  diagnostics = list(dateWarning = FALSE, verbose = FALSE),
  execution = list(blasThreads = request)))
      print(fit$threadingDiagnostics[c("blasBackend", "controller",
        "blasThreadsRequested", "blasThreadsConfigured", "availableCPUs",
        "openmpMaxThreads", "numericalPath", "averageCPUs")])
      stopifnot(getThreads() == before)
    }
  })
}else{
  message("Verification skipped: install the matching controller and parallelly, and load a sommer build with blasThreads.")
}

## ----user-model, eval=FALSE---------------------------------------------------
# library(sommer)
# mixm2 <- mmes(value~envf,
#   random=~vsm(fam(envf, 1), ism(id), Gu=MMTi) +
#            vsm(dsm(envf), ism(repf)),
#   rcov=~vsm(dsm(envf), ism(units)),
#   data=DT,
#   control = mmes_control(convergence = list(maxIters = 5),
#   execution = list(blasThreads = "auto")))
# 
# mixm2$threadingDiagnostics

