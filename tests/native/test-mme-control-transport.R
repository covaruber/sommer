# Run from the repository root. Set TMPDIR to a directory inside the repository
# when the environment restricts compiler scratch files.
(function(){
  cache <- tempfile("sommer-control-cache-")
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  Rcpp::sourceCpp("tests/native/test-mme-control-transport.cpp",
                  cacheDir = cache, showOutput = FALSE)
  old <- options("sommer.mme.testMode", "sommer.mme.memoryPolicy")
  on.exit(options(old), add = TRUE)
  descriptor <- list()
  fails <- function(expr, message) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), identical(conditionMessage(error), message))
  }
  options(sommer.mme.testMode = NULL)
  stopifnot(identical(control_mode_probe(descriptor), "auto"),
            identical(control_mode_probe(descriptor, "off"), "off"))
  options(sommer.mme.testMode = TRUE)
  stopifnot(identical(control_mode_probe(descriptor), "force"))
  attr(descriptor, "testMode") <- "FaLsE"
  stopifnot(identical(control_mode_probe(descriptor), "off"))
  attr(descriptor, "testMode") <- NA
  stopifnot(identical(control_mode_probe(descriptor), "auto"),
            identical(control_mode_probe(descriptor, "off"), "off"))
  attr(descriptor, "testMode") <- 1L
  stopifnot(identical(control_mode_probe(descriptor), "force"))
  attr(descriptor, "testMode") <- c(TRUE, FALSE)
  stopifnot(identical(control_mode_probe(descriptor), "force"))
  options(sommer.mme.testMode = FALSE)
  stopifnot(identical(control_mode_probe(descriptor), "off"))
  options(sommer.mme.testMode = "force")
  fails(control_mode_probe(descriptor), "sommer.mme.testMode must be TRUE or FALSE.")
  attr(descriptor, "testMode") <- NA
  stopifnot(identical(control_mode_probe(descriptor), "auto"))
  attr(descriptor, "testMode") <- "auto"
  stopifnot(identical(control_mode_probe(descriptor), "auto"))
  attr(descriptor, "testMode") <- "invalid"
  fails(control_mode_probe(descriptor),
        "testMode must be 'auto', 'off', 'force', TRUE, or FALSE.")
  attr(descriptor, "testMode") <- NULL
  options(sommer.mme.testMode = NA)
  fails(control_mode_probe(descriptor), "sommer.mme.testMode must be TRUE or FALSE.")
  options(sommer.mme.testMode = c(TRUE, FALSE))
  fails(control_mode_probe(descriptor), "sommer.mme.testMode must be TRUE or FALSE.")
  options(sommer.mme.testMode = 1L)
  fails(control_mode_probe(descriptor), "sommer.mme.testMode must be TRUE or FALSE.")
  options(sommer.mme.testMode = NULL)
  stopifnot(identical(control_proposal_probe(descriptor),
    list(hybrid = FALSE, aitken = FALSE, accelerated = FALSE, blend = "em", adaptive = FALSE)))
  attr(descriptor, "acceleration") <- "hybrid-ai"
  attr(descriptor, "blend") <- "pxem"
  attr(descriptor, "adaptiveRecovery") <- TRUE
  stopifnot(identical(control_proposal_probe(descriptor),
    list(hybrid = TRUE, aitken = FALSE, accelerated = TRUE, blend = "pxem", adaptive = TRUE)))
  # Recursive wrappers replace attributes; never cache their earlier policy.
  attr(descriptor, "acceleration") <- "none"
  attr(descriptor, "adaptiveRecovery") <- FALSE
  stopifnot(!control_proposal_probe(descriptor)$hybrid,
            !control_proposal_probe(descriptor)$accelerated,
            !control_proposal_probe(descriptor)$adaptive)
  attr(descriptor, "blend") <- "invalid"
  fails(control_proposal_probe(descriptor), "Unknown EM blending method.")
  options(sommer.mme.memoryPolicy = NULL)
  stopifnot(identical(control_memory_probe(descriptor), "speed"))
  options(sommer.mme.memoryPolicy = "limit")
  stopifnot(identical(control_memory_probe(descriptor), "limit"))
  attr(descriptor, "memoryPolicy") <- "speed"
  stopifnot(identical(control_memory_probe(descriptor), "speed"))
  attr(descriptor, "memoryPolicy") <- "Speed"
  fails(control_memory_probe(descriptor), "memoryPolicy must be 'speed' or 'limit'.")
  attr(descriptor, "memoryPolicy") <- NULL
  options(sommer.mme.memoryPolicy = TRUE)
  fails(control_memory_probe(descriptor), "memoryPolicy must be 'speed' or 'limit'.")
  cat("Native control transport tests passed.\n")
})()
