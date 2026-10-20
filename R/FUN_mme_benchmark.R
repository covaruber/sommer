# Experimental Benchmarking and Crossover Analysis Framework for sommer MME Solvers
# Evaluates experimental methods against standard baselines across problem dimensions,
# levels of factors, and covariance structures to calibrate crossover decision paths.

#' Generate Parameterized MME Workloads for Solver Evaluation
#'
#' Instantiates reproducible synthetic mixed model datasets and formulas spanning
#' various factor dimensions, covariance structures, relationship densities, and tasks.
#'
#' @param name Descriptive identifier for the workload.
#' @param n_ind Number of unique individuals / genotypes. Default is 60.
#' @param n_levels Number of levels in the primary factor (e.g. environments, timepoints, or sites). Default is 5.
#' @param structure Covariance structure type. Supported: \code{"ar1"}, \code{"ar2"},
#'   \code{"us"}, \code{"diag"}, \code{"fa"}, \code{"rr"}, \code{"spatial"}, \code{"sar"}, \code{"car"}.
#' @param relationship Relationship matrix type for individuals: \code{"none"} (identity),
#'   \code{"sparse"} (tridiagonal pedigree-like inverse), or \code{"dense"} (genomic relationship inverse).
#' @param missing_rate Proportion of missing responses (0.0 to 0.5). Default is 0.0.
#' @param reml Logical; whether to fit by REML (TRUE) or ML (FALSE). Default is TRUE.
#' @param computeCi Integer; final inverse demand: 0 (none), 1 (Takahashi selected inverse), 2 (full inverse). Default is 0.
#' @param seed Random seed for reproducibility. Default is 42.
#'
#' @return A list of class \code{"mme_workload"} containing \code{data}, \code{fixed}, \code{random},
#'   \code{rcov}, \code{REML}, \code{computeCi}, and metadata.
#' @export
generate_mme_workload <- function(name = NULL,
                                 n_ind = 60L,
                                 n_levels = 5L,
                                 structure = c("ar1", "ar2", "us", "diag", "fa", "rr", "spatial", "sar", "car"),
                                 relationship = c("none", "dense", "sparse"),
                                 missing_rate = 0.0,
                                 reml = TRUE,
                                 computeCi = 0L,
                                 seed = 42L){
  structure <- match.arg(structure)
  relationship <- match.arg(relationship)
  n_ind <- as.integer(n_ind)
  n_levels <- as.integer(n_levels)
  if(n_ind < 4L) stop("n_ind must be >= 4.", call.=FALSE)
  if(n_levels < 2L) stop("n_levels must be >= 2.", call.=FALSE)

  set.seed(seed)
  ind_ids <- paste0("ID_", sprintf("%04d", seq_len(n_ind)))
  level_ids <- paste0("L_", sprintf("%03d", seq_len(n_levels)))

  # Build relationship matrix Gu if requested
  Gu <- NULL
  if(relationship == "dense"){
    M <- matrix(stats::rnorm(n_ind * 40L), n_ind, 40L)
    G <- stats::cor(t(M)) + diag(0.05, n_ind)
    dimnames(G) <- list(ind_ids, ind_ids)
    Gu <- solve(G)
    attr(Gu, "inverse") <- TRUE
  } else if(relationship == "sparse"){
    # Tridiagonal precision matrix
    Gu <- Matrix::bandSparse(n_ind, k=c(0, 1), diagonals=list(rep(2.0, n_ind), rep(-0.4, n_ind - 1L)), symmetric=TRUE)
    dimnames(Gu) <- list(ind_ids, ind_ids)
    attr(Gu, "inverse") <- TRUE
  }

  # Build dataset
  df <- expand.grid(id = factor(ind_ids, levels=ind_ids),
                    level = factor(level_ids, levels=level_ids),
                    rep = seq_len(2L))
  n_obs <- nrow(df)

  # Underlying signal simulation
  u_ind <- stats::rnorm(n_ind, sd=1.5)
  names(u_ind) <- ind_ids
  beta_lvl <- stats::rnorm(n_levels, sd=2.0)
  names(beta_lvl) <- level_ids

  mu <- 10.0 + beta_lvl[as.character(df$level)] + u_ind[as.character(df$id)]
  noise <- stats::rnorm(n_obs, sd=1.0)
  df$response <- mu + noise

  if(missing_rate > 0.0){
    n_miss <- min(n_obs - 10L, floor(missing_rate * n_obs))
    if(n_miss > 0L){
      miss_idx <- sample(n_obs, n_miss)
      df$response[miss_idx] <- NA_real_
    }
  }

  fixed_form <- response ~ level
  rcov_form <- ~ units

  # Construct random effect formula by covariance structure
  random_form <- switch(
    structure,
    "ar1" = if(is.null(Gu)) ~ vsm(ar1m(level), ism(id)) else ~ vsm(ar1m(level), ism(id), Gu=Gu),
    "ar2" = if(is.null(Gu)) ~ vsm(ar2m(level), ism(id)) else ~ vsm(ar2m(level), ism(id), Gu=Gu),
    "us"  = if(is.null(Gu)) ~ vsm(usm(level), ism(id)) else ~ vsm(usm(level), ism(id), Gu=Gu),
    "diag" = if(is.null(Gu)) ~ vsm(dsm(level), ism(id)) else ~ vsm(dsm(level), ism(id), Gu=Gu),
    "fa"  = if(is.null(Gu)) ~ vsm(fam(level, 1L), ism(id)) else ~ vsm(fam(level, 1L), ism(id), Gu=Gu),
    "rr"  = if(is.null(Gu)) ~ vsm(rrm(level, 1L), ism(id)) else ~ vsm(rrm(level, 1L), ism(id), Gu=Gu),
    "spatial" = {
      # 2D spatial grid
      df$row <- factor(rep(seq_len(n_levels), length.out=n_obs))
      df$col <- factor(rep(seq_len(ceiling(n_obs / n_levels)), each=n_levels, length.out=n_obs))
      rcov_form <- ~ vsm(ar1m(row), ar1m(col))
      if(is.null(Gu)) ~ ism(id) else ~ vsm(ism(id), Gu=Gu)
    },
    "sar" = {
      # Synthetic spatial proximity matrix
      W_prox <- Matrix::bandSparse(n_levels, k=1, diagonals=list(rep(1.0, n_levels - 1L)), symmetric=TRUE)
      dimnames(W_prox) <- list(level_ids, level_ids)
      if(is.null(Gu)) ~ vsm(sar(level, W_prox), ism(id)) else ~ vsm(sar(level, W_prox), ism(id), Gu=Gu)
    },
    "car" = {
      W_prox <- Matrix::bandSparse(n_levels, k=1, diagonals=list(rep(1.0, n_levels - 1L)), symmetric=TRUE)
      dimnames(W_prox) <- list(level_ids, level_ids)
      if(is.null(Gu)) ~ vsm(car(level, W_prox), ism(id)) else ~ vsm(car(level, W_prox), ism(id), Gu=Gu)
    }
  )

  if(is.null(name)){
    name <- sprintf("%s_%dL_%dID_%s_%s", structure, n_levels, n_ind, relationship, if(reml) "REML" else "ML")
  }

  workload <- list(
    name = name,
    data = df,
    fixed = fixed_form,
    random = random_form,
    rcov = rcov_form,
    REML = reml,
    computeCi = as.integer(computeCi),
    dimensions = list(
      n_obs = n_obs,
      n_ind = n_ind,
      n_levels = n_levels,
      structure = structure,
      relationship = relationship,
      missing_rate = missing_rate
    )
  )
  class(workload) <- "mme_workload"
  workload
}

#' Generate a Factorial Grid of MME Workloads
#'
#' @param structures Vector of covariance structures to test.
#' @param n_ind Vector of individual counts.
#' @param n_levels Vector of factor level counts.
#' @param relationships Vector of relationship types.
#' @param reml Logical; whether to fit by REML.
#' @param computeCi Vector of computeCi modes.
#'
#' @return A named list of \code{mme_workload} objects.
#' @export
mme_workload_grid <- function(structures = c("ar1", "us", "diag"),
                              n_ind = c(40L, 100L),
                              n_levels = c(4L, 10L),
                              relationships = c("none", "dense"),
                              reml = TRUE,
                              computeCi = 0L){
  grid <- expand.grid(
    structure = structures,
    n_ind = n_ind,
    n_levels = n_levels,
    relationship = relationships,
    computeCi = computeCi,
    stringsAsFactors = FALSE
  )
  workloads <- list()
  for(i in seq_len(nrow(grid))){
    row <- grid[i, ]
    wl <- generate_mme_workload(
      structure = row$structure,
      n_ind = row$n_ind,
      n_levels = row$n_levels,
      relationship = row$relationship,
      reml = reml,
      computeCi = row$computeCi
    )
    if(length(unique(computeCi)) > 1L){
      wl$name <- paste0(wl$name, "_Ci", row$computeCi)
    }
    workloads[[wl$name]] <- wl
  }
  workloads
}

.mme_benchmark_row <- function(workload, name, engine, fit=NULL,
                                elapsed=NA_real_, parity=NA, difference=NA_real_,
                                error=NA_character_){
  diagnostics <- fit$engineDiagnostics
  data.frame(workload=name, structure=workload$dimensions$structure,
    n_levels=workload$dimensions$n_levels, n_ind=workload$dimensions$n_ind,
    relationship=workload$dimensions$relationship, engine=engine,
    elapsed_sec=elapsed,
    iterations=if(is.null(fit$monitor)) NA_integer_ else ncol(fit$monitor),
    converged=isTRUE(fit$convergence),
    blockSchurActive=isTRUE(diagnostics$blockSchurActive),
    blockBandedActive=isTRUE(diagnostics$blockBandedActive),
    factorSchurActive=isTRUE(diagnostics$factorSchurActive),
    parity_ok=parity, max_rel_diff=difference, error=error,
    .mme_profile_columns(diagnostics), stringsAsFactors=FALSE)
}

.mme_benchmark_difference <- function(fit, reference){
  relative <- function(candidate, baseline){
    if(!length(candidate) || length(candidate) != length(baseline) ||
       any(!is.finite(candidate)) || any(!is.finite(baseline))) return(Inf)
    max(abs(candidate - baseline) / (abs(baseline) + 1))
  }
  max(relative(tail(as.numeric(fit$llik), 1L), tail(as.numeric(reference$llik), 1L)),
      relative(unlist(fit$theta, use.names=FALSE), unlist(reference$theta, use.names=FALSE)),
      relative(as.numeric(fit$bu), as.numeric(reference$bu)))
}

#' Run Experimental MME Benchmark Suite
#'
#' Executes workloads across solver engines, recording elapsed time, convergence,
#' active numerical paths and numerical parity.
#' @param workloads A list of \code{mme_workload} objects (or a single workload).
#' @param engines A named list of \code{\link{mmes_control}} objects.
#' @param replications Positive number of timed repetitions, excluding warmup.
#' @param verifyParity Whether to compare candidates against the first engine.
#' @param tolerance Positive relative tolerance for parity verification.
#' @param nIters Maximum REML/ML iterations for benchmark fits.
#' @param verbose Whether to print progress.
#' @return A \code{sommer_mme_benchmark} object. Failed runs retain the same
#'   result columns with an error message. Parity is \code{NA} when not checked
#'   or when the reference failed; comparisons use final likelihoods and estimates.
#' @export
run_mme_benchmark <- function(workloads,
                              engines = NULL,
                              replications = 3L,
                              verifyParity = TRUE,
                              tolerance = 1.0e-5,
                              nIters = 5L,
                              verbose = TRUE){
  replications <- .mme_positive(replications, "replications", integer=TRUE)
  tolerance <- .mme_positive(tolerance, "tolerance")
  if(inherits(workloads, "mme_workload")){
    workloads <- list(workloads)
    names(workloads) <- workloads[[1L]]$name
  }
  mmeControl <- get("mmes_control", mode="function")
  mmesFn <- get("mmes", mode="function")
  if(is.null(engines)){
    engines <- list(
      baseline = mmeControl(),
      banded_force = mmeControl(solve=list(blockBanded="force")),
      latent_force = mmeControl(solve=list(factorSchurMultiFactor="force"))
    )
  }

  results <- list()

  for(wl_name in names(workloads)){
    wl <- workloads[[wl_name]]
    if(verbose) message(crayon::bold(paste0("\n--- Workload: ", wl_name, " ---")))

    ref_fit <- NULL
    ref_engine_name <- names(engines)[1L]

    for(eng_name in names(engines)){
      ctrl <- engines[[eng_name]]
      ctrl$task$objective <- if(wl$REML) "reml" else "ml"
      ctrl$output$inverse <- wl$computeCi
      ctrl$convergence$maxIters <- nIters
      ctrl$diagnostics$verbose <- FALSE
      ctrl$diagnostics$dateWarning <- FALSE
      if(verbose) message(sprintf("  Running engine: %-18s ... ", eng_name), appendLF=FALSE)

      # Warmup run
      warm_fit <- tryCatch({
        mmesFn(fixed = wl$fixed, random = wl$random, rcov = wl$rcov, data = wl$data,
               control = ctrl)
      }, error = function(e){
        list(error = e$message)
      })

      if(!is.null(warm_fit$error)){
        if(verbose) message(crayon::red(paste("FAILED:", warm_fit$error)))
        results[[length(results) + 1L]] <-
          .mme_benchmark_row(wl, wl_name, eng_name, error=warm_fit$error)
        next
      }

      # Timed replications
      reps <- replications
      times <- numeric(reps)
      last_fit <- warm_fit

      for(r in seq_len(reps)){
        t0 <- proc.time()
        last_fit <- tryCatch(
          mmesFn(fixed=wl$fixed, random=wl$random, rcov=wl$rcov, data=wl$data,
                 control=ctrl),
          error=function(e) list(error=conditionMessage(e)))
        if(!is.null(last_fit$error)) break
        t1 <- proc.time()
        times[r] <- (t1 - t0)[3L]
      }
      if(!is.null(last_fit$error)){
        if(verbose) message(crayon::red(paste("FAILED:", last_fit$error)))
        results[[length(results) + 1L]] <-
          .mme_benchmark_row(wl, wl_name, eng_name, error=last_fit$error)
        next
      }

      med_time <- stats::median(times)
      diag <- last_fit$engineDiagnostics

      if(verbose){
        active_paths <- c()
        if(isTRUE(diag$blockSchurActive)) active_paths <- c(active_paths, "Schur")
        if(isTRUE(diag$blockChainActive)) active_paths <- c(active_paths, "Chain")
        if(isTRUE(diag$blockBandedActive)) active_paths <- c(active_paths, paste0("Banded(", diag$blockBandedWidth, ")"))
        if(isTRUE(diag$factorSchurActive)) active_paths <- c(active_paths, "Latent")
        if(isTRUE(diag$dFirstMLActive)) active_paths <- c(active_paths, "DFirst")
        paths_str <- if(length(active_paths)) paste(active_paths, collapse=", ") else "standard"
        message(crayon::green(sprintf("%.4fs (med of %d) [%s]", med_time, reps, paths_str)))
      }

      # Parity check against reference engine
      parity_ok <- NA
      max_rel_diff <- NA_real_
      if(eng_name == ref_engine_name){
        ref_fit <- last_fit
        max_rel_diff <- .mme_benchmark_difference(last_fit, ref_fit)
        parity_ok <- max_rel_diff <= tolerance
      } else if(verifyParity && !is.null(ref_fit)){
        max_rel_diff <- .mme_benchmark_difference(last_fit, ref_fit)
        parity_ok <- (max_rel_diff <= tolerance)
        if(!parity_ok && verbose){
          message(crayon::yellow(sprintf("    Warning: Parity deviation = %.2e (tolerance = %.2e)",
                                         max_rel_diff, tolerance)))
        }
      }

      results[[length(results) + 1L]] <- .mme_benchmark_row(
        wl, wl_name, eng_name, last_fit, med_time, parity_ok, max_rel_diff)
    }
  }

  res_df <- do.call(rbind, results)
  out <- list(
    results = res_df,
    workloads = workloads,
    engines = engines,
    replications = replications,
    tolerance = tolerance
  )
  class(out) <- "sommer_mme_benchmark"
  out
}

.mme_phase_names <- c("setup", "factorEval", "assembly", "factorization", "acceptance",
                      "tracePrep", "sensitivities", "aiAssembly", "scoreTraces",
                      "update", "constraints", "proposal", "final")

# One-row data frame of phase seconds and workload counters; NA when profiling is off.
.mme_profile_columns <- function(diag){
  prof <- if(is.list(diag)) diag$profile else NULL
  cols <- as.list(rep(NA_real_, length(.mme_phase_names)))
  names(cols) <- paste0("phase_", .mme_phase_names)
  counters <- list(nEffects = NA_real_, nnzC = NA_real_, freeParams = NA_real_,
                   numericFactorSeconds = NA_real_, sensitivityRHSColumns = NA_real_,
                   exactAIRefreshes = NA_real_)
  dominant <- NA_character_
  if(is.list(prof)){
    secs <- prof$seconds
    cols[paste0("phase_", names(secs))] <- as.list(unname(secs))
    sourceNames <- c(nEffects="nEffects", nnzC="nnzC",
      freeParams="freeCovarianceParameters", numericFactorSeconds="numericFactorSeconds",
      sensitivityRHSColumns="sensitivityRHSColumns", exactAIRefreshes="exactAIRefreshes")
    for(name in names(sourceNames)){
      value <- prof$counters[[sourceNames[[name]]]]
      if(!is.null(value)) counters[[name]] <- value
    }
    if(length(secs)) dominant <- names(which.max(secs))
  }
  as.data.frame(c(cols, counters, list(dominantPhase = dominant)), stringsAsFactors = FALSE)
}

#' Profile an MME Workload by Phase
#'
#' Fits one workload with phase profiling enabled and tabulates where the time goes.
#'
#' @param workload An \code{mme_workload} from \code{\link{generate_mme_workload}}.
#' @param engine A \code{\link{mmes_control}} object; profiling is forced on.
#' @param nIters Maximum REML/ML iterations. Default is 5.
#' @param dominantShare Share of fit time above which a phase is flagged dominant. Default is 0.3.
#'
#' @return A list of class \code{"sommer_mme_profile"} with a per-phase \code{table},
#'   workload \code{counters}, per-iteration seconds, and the total elapsed time.
#' @export
profile_mme_workload <- function(workload, engine = mmes_control(), nIters = 5L,
                                 dominantShare = 0.3){
  if(!inherits(workload, "mme_workload")){
    stop("workload must be of class 'mme_workload'.", call.=FALSE)
  }
  engine$diagnostics$profile <- "force"
  engine$task$objective <- if(workload$REML) "reml" else "ml"
  engine$output$inverse <- workload$computeCi
  engine$convergence$maxIters <- nIters
  engine$diagnostics$verbose <- FALSE
  engine$diagnostics$dateWarning <- FALSE
  t0 <- proc.time()
  fit <- mmes(fixed = workload$fixed, random = workload$random, rcov = workload$rcov,
              data = workload$data, control = engine)
  elapsed <- unname((proc.time() - t0)[3L])
  prof <- fit$engineDiagnostics$profile
  if(!is.list(prof)) stop("The fit did not return a profile.", call.=FALSE)
  secs <- prof$seconds
  table <- data.frame(phase = names(secs), seconds = unname(secs),
                      share = unname(secs) / sum(secs), stringsAsFactors = FALSE)
  table$dominant <- table$share >= dominantShare
  out <- list(workload = workload$name, table = table, counters = prof$counters,
              secondsByIteration = prof$secondsByIteration,
              elapsed_sec = elapsed, engine = engine)
  class(out) <- "sommer_mme_profile"
  out
}

#' Print method for sommer_mme_profile
#' @param x An object of class \code{sommer_mme_profile}.
#' @param ... Ignored.
#' @method print sommer_mme_profile
#' @export
print.sommer_mme_profile <- function(x, ...){
  cat("-- sommer MME Phase Profile:", x$workload, "--\n")
  tab <- x$table
  tab$seconds <- round(tab$seconds, 4)
  tab$share <- round(tab$share, 3)
  print(tab, row.names = FALSE)
  cat(sprintf("Elapsed %.3fs | nEffects %d | nnz(C) %.0f | free covariance parameters %d | solver %s\n",
              x$elapsed_sec, x$counters$nEffects, x$counters$nnzC,
              x$counters$freeCovarianceParameters, x$counters$solver))
  invisible(x)
}

#' Analyze Crossover Frontiers from MME Benchmarks
#'
#' Identifies problem configurations (factors, levels, dimensions) where experimental solver
#' paths provide speedups over baseline methods, fitting empirical crossover conditions.
#'
#' @param benchmark An object of class \code{sommer_mme_benchmark}.
#' @param baseline_engine Name of the reference baseline engine. Default is the first engine tested.
#' @param candidate_engine Name of the experimental engine to evaluate.
#'
#' @return An object of class \code{"sommer_mme_crossover"}.
#' @export
analyze_mme_crossovers <- function(benchmark,
                                   baseline_engine = NULL,
                                   candidate_engine = NULL){
  if(!inherits(benchmark, "sommer_mme_benchmark")){
    stop("benchmark must be of class 'sommer_mme_benchmark'.", call.=FALSE)
  }
  df <- benchmark$results
  engines_tested <- unique(df$engine)
  if(length(engines_tested) < 2L){
    stop("Benchmark must contain at least 2 engines to analyze crossovers.", call.=FALSE)
  }

  if(is.null(baseline_engine)) baseline_engine <- engines_tested[1L]
  if(is.null(candidate_engine)) candidate_engine <- engines_tested[2L]

  base_df <- df[df$engine == baseline_engine, ]
  cand_df <- df[df$engine == candidate_engine, ]

  merged <- merge(base_df, cand_df, by=c("workload", "structure", "n_levels", "n_ind", "relationship"),
                  suffixes=c(".base", ".cand"))

  merged$speedup <- merged$elapsed_sec.base / merged$elapsed_sec.cand
  merged$advantageous <- is.finite(merged$speedup) & merged$speedup > 1.05 &
    !is.na(merged$parity_ok.cand) & merged$parity_ok.cand

  # Categorize by structure and levels
  crossover_summary <- list()
  for(struct in unique(merged$structure)){
    sub <- merged[merged$structure == struct, ]
    win_cases <- sub[sub$advantageous, ]
    loss_cases <- sub[!sub$advantageous & !is.na(sub$speedup), ]

    crossover_summary[[struct]] <- list(
      structure = struct,
      total_cases = nrow(sub),
      wins = nrow(win_cases),
      losses = nrow(loss_cases),
      min_speedup = if(nrow(sub)) min(sub$speedup, na.rm=TRUE) else NA_real_,
      max_speedup = if(nrow(sub)) max(sub$speedup, na.rm=TRUE) else NA_real_,
      median_speedup = if(nrow(sub)) stats::median(sub$speedup, na.rm=TRUE) else NA_real_,
      recommendation = if(nrow(win_cases) == 0){
        "Keep off by default; no winning regime detected."
      } else if(nrow(loss_cases) == 0){
        "Safe to promote to default auto-dispatch across all tested dimensions."
      } else {
        paste0("Conditionally advantageous: win when n_levels >= ", min(win_cases$n_levels),
               " and n_ind >= ", min(win_cases$n_ind), ".")
      }
    )
  }

  out <- list(
    comparison = paste(candidate_engine, "vs", baseline_engine),
    baseline_engine = baseline_engine,
    candidate_engine = candidate_engine,
    data = merged,
    summaries = crossover_summary
  )
  class(out) <- "sommer_mme_crossover"
  out
}

#' Print method for sommer_mme_benchmark
#' @param x An object of class \code{sommer_mme_benchmark}.
#' @param ... Ignored.
#' @method print sommer_mme_benchmark
#' @export
print.sommer_mme_benchmark <- function(x, ...){
  cat("-- sommer MME Benchmark Results --\n")
  cat("Workloads tested :", length(x$workloads), "\n")
  cat("Engines evaluated:", paste(names(x$engines), collapse=", "), "\n\n")

  print_cols <- c("workload", "engine", "elapsed_sec", "iterations", "parity_ok", "max_rel_diff")
  print(x$results[, print_cols], row.names=FALSE)
  invisible(x)
}

#' Print method for sommer_mme_crossover
#' @param x An object of class \code{sommer_mme_crossover}.
#' @param ... Ignored.
#' @method print sommer_mme_crossover
#' @export
print.sommer_mme_crossover <- function(x, ...){
  cat("-- sommer MME Crossover Analysis --\n")
  cat("Candidate:", x$candidate_engine, "| Baseline:", x$baseline_engine, "\n\n")

  for(struct in names(x$summaries)){
    s <- x$summaries[[struct]]
    cat(sprintf("Structure [%s]:\n", struct))
    cat(sprintf("  Cases: %d | Wins: %d | Median Speedup: %.2fx (range: %.2f - %.2fx)\n",
                s$total_cases, s$wins, s$median_speedup, s$min_speedup, s$max_speedup))
    cat("  Recommendation:", s$recommendation, "\n\n")
  }
  invisible(x)
}
