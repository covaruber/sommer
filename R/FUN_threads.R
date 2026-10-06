.mmes_blas_controller <- function(){
  if(requireNamespace("flexiblas", quietly=TRUE) &&
     isTRUE(tryCatch(flexiblas::flexiblas_avail(), error=function(error) FALSE))){
    current <- tryCatch(flexiblas::flexiblas_get_num_threads(),
      error=function(error) NA_integer_)
    if(length(current) == 1L && is.finite(current) && current > 0){
      return(list(provider="flexiblas",
        backend=flexiblas::flexiblas_current_backend(),
        get=flexiblas::flexiblas_get_num_threads,
        set=flexiblas::flexiblas_set_num_threads))
    }
  }
  blas <- unname(extSoftVersion()["BLAS"])
  if(!is.na(blas) && grepl("flexiblas", blas, ignore.case=TRUE)) return(NULL)
    if(!is.na(blas) && grepl("openblas|mkl|acml|goto", blas, ignore.case=TRUE) &&
      requireNamespace("RhpcBLASctl", quietly=TRUE)){
    current <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
      error=function(error) NA_integer_)
    if(length(current) == 1L && is.finite(current) && current > 0){
      return(list(provider="RhpcBLASctl",
        backend=blas,
        get=RhpcBLASctl::blas_get_num_procs,
        set=RhpcBLASctl::blas_set_num_threads))
    }
  }
  NULL
}

.mmes_available_cpus <- function(){
  if(requireNamespace("parallelly", quietly=TRUE)){
    return(as.integer(parallelly::availableCores()))
  }
  NA_integer_
}

.mmes_blas_scope <- function(requested=NULL, controller=.mmes_blas_controller()){
  if(is.character(requested) && identical(requested, "auto")){
    requested <- .mmes_available_cpus()
    if(is.na(requested)){
      stop("blasThreads='auto' requires the optional parallelly package; alternatively supply a positive integer.",
        call.=FALSE)
    }
  }
  if(!is.null(requested) && (!is.numeric(requested) || length(requested) != 1L ||
     is.na(requested) || !is.finite(requested) || requested < 1 ||
     requested > .Machine$integer.max || requested != floor(requested))){
    stop("blasThreads must be NULL, 'auto', or one positive integer.", call.=FALSE)
  }
  previous <- if(is.null(controller)) NA_integer_ else as.integer(controller$get())
  restore <- function() invisible(NULL)
  if(!is.null(requested)){
    if(is.null(controller)){
      stop("Runtime BLAS control is unavailable. Install flexiblas for a FlexiBLAS R installation, or RhpcBLASctl for a supported BLAS; alternatively configure BLAS before launching R and leave blasThreads=NULL.",
        call.=FALSE)
    }
    if(length(previous) != 1L || is.na(previous) || previous < 1L){
      stop("Cannot safely restore the current BLAS thread setting.", call.=FALSE)
    }
    restore <- function() controller$set(previous)
    applied <- FALSE
    on.exit(if(!applied) restore(), add=TRUE)
    controller$set(as.integer(requested))
  }
  configured <- if(is.null(controller)) NA_integer_ else as.integer(controller$get())
  if(!is.null(requested) && !identical(configured, as.integer(requested))){
    stop("The BLAS backend did not accept the requested blasThreads setting.", call.=FALSE)
  }
  scope <- list(restore=restore, info=list(
    blas=unname(extSoftVersion()["BLAS"]),
    blasBackend=if(is.null(controller)) NA_character_ else controller$backend,
    controller=if(is.null(controller)) "unavailable" else controller$provider,
    blasThreadsRequested=requested,
    blasThreadsConfigured=configured,
    availableCPUs=.mmes_available_cpus(),
    environment=Sys.getenv(c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
      "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "SLURM_CPUS_PER_TASK"))))
  scope$started <- proc.time()
  if(!is.null(requested)) applied <- TRUE
  scope
}

.mmes_thread_report <- function(scope){
  formatCount <- function(value) if(length(value) != 1L || is.na(value)) "unknown" else value
  cat("BLAS backend: ", if(length(scope$info$blasBackend) == 1L &&
    !is.na(scope$info$blasBackend)) scope$info$blasBackend else scope$info$blas, ".\n", sep="")
  cat("BLAS threads configured: ", formatCount(scope$info$blasThreadsConfigured),
    " (", scope$info$controller, "); allocation-aware CPUs: ",
    formatCount(scope$info$availableCPUs), ".\n", sep="")
}

.mmes_thread_finish <- function(result, scope, verbose=FALSE){
  elapsed <- proc.time() - scope$started
  info <- scope$info
  info$elapsedSeconds <- unname(elapsed["elapsed"])
  info$cpuSeconds <- unname(elapsed["user.self"] + elapsed["sys.self"])
  info$averageCPUs <- if(info$elapsedSeconds > 0) info$cpuSeconds / info$elapsedSeconds else NA_real_
  native <- result$factorizationDiagnostics
  if(!is.null(native)){
    info$openmpMaxThreads <- native$openmpMaxThreads
    info$factorizations <- native$timings
  }
  engine <- result$engineDiagnostics
  info$numericalPath <- if(isTRUE(engine$factorSchurActive)) "latent-factor Schur" else
    if(isTRUE(engine$blockChainActive)) "block-chain Schur" else
    if(isTRUE(engine$blockSchurActive)) "block Schur" else
    if(length(result$solver) == 1L) result$solver else "unknown"
  result$threadingDiagnostics <- info
  if(verbose){
    cat(sprintf("Fit CPU usage: %.3f CPU-s / %.3f elapsed-s = %.2f average CPUs (not an active-thread count).\n",
      info$cpuSeconds, info$elapsedSeconds, info$averageCPUs))
    for(path in names(info$factorizations)){
      timing <- info$factorizations[[path]]
      if(timing$calls > 0L){
        cat(sprintf("%s usage (%d calls, %s): %.3f CPU-s / %.3f elapsed-s = %.2f average CPUs.\n",
          path, timing$calls, timing$measurement, timing$cpuSeconds, timing$elapsedSeconds, timing$averageCPUs))
      }
    }
  }
  result
}