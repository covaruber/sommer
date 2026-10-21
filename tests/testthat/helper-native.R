.native_require <- function(condition, reason){
  if(isTRUE(condition)) return(invisible(NULL))
  if(identical(Sys.getenv("SOMMER_REQUIRE_NATIVE"), "true")){
    stop(reason, call.=FALSE)
  }
  testthat::skip(reason)
}

.native_ai_recovery <- function(){
  source <- testthat::test_path("..", "native", "test-mme-ai-recovery.cpp")
  header <- testthat::test_path("..", "..", "src", "mme_ai_recovery.h")
  .native_require(all(file.exists(c(source, header))),
    "Native checks require the package source tree")
  .native_require(nzchar(Sys.which("c++")), "C++ compiler unavailable")
  .native_require(requireNamespace("RcppArmadillo", quietly=TRUE),
    "RcppArmadillo unavailable")
  helper <- new.env()
  Rcpp::sourceCpp(source, env=helper)
  helper
}

.native_run <- function(name, headers, eigen=FALSE, strict=FALSE){
  root <- testthat::test_path("..", "..")
  source <- file.path(root, "tests", "native", paste0(name, ".cpp"))
  .native_require(all(file.exists(c(source, file.path(root, "src", headers)))),
    "Native checks require the package source tree")
  compiler <- Sys.which("c++")
  .native_require(nzchar(compiler), "C++ compiler unavailable")
  include <- if(eigen) system.file("include", package="RcppEigen") else character()
  if(eigen) .native_require(nzchar(include), "RcppEigen headers unavailable")
  binary <- tempfile(paste0("sommer-", name, "-"),
    fileext=if(.Platform$OS.type == "windows") ".exe" else "")
  on.exit(unlink(binary), add=TRUE)
  flags <- if(strict) c("-Wall", "-Wextra", "-pedantic") else "-w"
  args <- c("-std=c++11", "-O0", flags,
    if(eigen) paste0("-I", shQuote(include)), shQuote(source), "-o", shQuote(binary))
  output <- suppressWarnings(system2(compiler, args, stdout=TRUE, stderr=TRUE))
  if(!is.null(attr(output, "status"))){
    stop("Native compilation failed for ", name, ":\n",
         paste(output, collapse="\n"), call.=FALSE)
  }
  output <- suppressWarnings(system2(binary, stdout=TRUE, stderr=TRUE))
  if(!is.null(attr(output, "status"))){
    stop("Native test failed for ", name, ":\n",
         paste(output, collapse="\n"), call.=FALSE)
  }
  testthat::succeed()
}
