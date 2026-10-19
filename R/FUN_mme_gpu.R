.mmes_gpu_hardware <- function(){
  visible <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA_character_)
  if(!is.na(visible) && visible %in% c("", "-1")) return(FALSE)
  if(.Platform$OS.type == "unix" && dir.exists("/sys/bus/pci/devices")){
    devices <- list.files("/sys/bus/pci/devices", full.names = TRUE)
    readValue <- function(path){
      tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(error) "")
    }
    for(device in devices){
      if(identical(readValue(file.path(device, "vendor")), "0x10de") &&
         any(grepl("^0x03", readValue(file.path(device, "class"))))){
        return(TRUE)
      }
    }
  }
  command <- unname(Sys.which("nvidia-smi"))
  if(!nzchar(command)) return(FALSE)
  output <- tryCatch(suppressWarnings(system2(command,
    c("--query-gpu=name", "--format=csv,noheader"), stdout = TRUE, stderr = TRUE,
    timeout = 2L)), error = function(error) character())
  length(output) > 0L && is.null(attr(output, "status")) &&
    !any(grepl("no devices|failed|error", output, ignore.case = TRUE))
}

.mmes_gpu_status <- function(hardware = .mmes_gpu_hardware,
                             libraries = function()
                               .Call("_sommer_mme_gpu_runtime_status", PACKAGE = "sommer")){
  status <- list(hardware = isTRUE(hardware()), libraries = NA,
    compiled = NA, available = FALSE, reason = "No accessible NVIDIA GPU detected.")
  if(!status$hardware) return(status)
  runtime <- tryCatch(libraries(), error = function(error)
    list(libraries = FALSE, compiled = NA, available = FALSE,
      reason = conditionMessage(error)))
  utils::modifyList(status, runtime)
}

.mmes_gpu_message <- function(status){
  if(!isTRUE(status$hardware)){
    return("sommer: No accessible NVIDIA GPU detected; only CPU options are available.")
  }
  if(isTRUE(status$available)){
    return(paste0("sommer: NVIDIA GPU and CUDA are ready. Enable eligible GPU kernels with ",
      "options(sommer.mme.gpu=TRUE) or mmes_control(execution=list(gpu=TRUE)). CPU remains the default."))
  }
  if(!isTRUE(status$libraries)){
    return(paste0("sommer: NVIDIA GPU detected, but CUDA runtime, cuBLAS or cuSOLVER cannot be loaded. ",
      "Install compatible CUDA libraries and set SOMMER_CUDA_HOME or the loader path. Only CPU options are available. ",
      status$reason))
  }
  if(!isTRUE(status$compiled)){
    return(paste0("sommer: NVIDIA GPU and CUDA libraries detected, but this build has no GPU kernels. ",
      "Reinstall with SOMMER_CUDA_HOME=/path/to/cuda R CMD INSTALL --preclean . before enabling sommer.mme.gpu."))
  }
  paste0("sommer: NVIDIA GPU and CUDA libraries detected, but GPU initialization failed. ",
    "Only CPU options are available. ", status$reason)
}

.mmes_gpu_require <- function(){
  status <- .mmes_gpu_status()
  if(!isTRUE(status$available)){
    stop("gpu=TRUE cannot run: ", .mmes_gpu_message(status), call. = FALSE)
  }
  invisible(status)
}