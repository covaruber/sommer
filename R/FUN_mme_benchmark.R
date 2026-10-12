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
    workloads[[wl$name]] <- wl
  }
  workloads
}

#' Run Experimental MME Benchmark Suite
#'
#' Executes a battery of workloads across candidate solver engines, measuring wall-clock time,
#' iterations, peak memory (if available), active numerical paths, and asserting numerical parity.
#'
#' @param workloads A list of \code{mme_workload} objects (or a single workload).
#' @param engines A named list of \code{\link{mmes_control}} objects. Defaults to comparing
#'   \code{baseline} (standard auto solver) against candidate experimental toggles.
#' @param replications Number of timed repetitions (warmup run excluded). Default is 3.
#' @param verifyParity Logical; whether to assert that experimental candidates match the reference
#'   engine in likelihood, variance components, and BLUPs. Default is TRUE.
#' @param tolerance Numerical tolerance for parity verification. Default is 1e-5.
#' @param nIters Maximum REML/ML iterations for benchmark fits. Default is 5.
#' @param verbose Logical; whether to print progress. Default is TRUE.
#'
#' @return An object of class \code{"sommer_mme_benchmark"}.
#' @export
run_mme_benchmark <- function(workloads,
                              engines = NULL,
                              replications = 3L,
                              verifyParity = TRUE,
                              tolerance = 1.0e-5,
                              nIters = 5L,
                              verbose = TRUE){
  if(inherits(workloads, "mme_workload")){
    workloads <- list(workloads)
    names(workloads) <- workloads[[1L]]$name
  }
  mmeControl <- get("mmes_control", mode="function")
  mmesFn <- get("mmes", mode="function")
  if(is.null(engines)){
    engines <- list(
      baseline = mmeControl(solver="auto"),
      banded_force = mmeControl(experimental=list(blockBanded="force")),
      latent_force = mmeControl(experimental=list(factorSchurMultiFactor="force"))
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
      if(verbose) message(sprintf("  Running engine: %-18s ... ", eng_name), appendLF=FALSE)

      # Warmup run
      warm_fit <- tryCatch({
        mmesFn(fixed = wl$fixed, random = wl$random, rcov = wl$rcov, data = wl$data,
               REML = wl$REML, computeCi = wl$computeCi, nIters = nIters,
               engine = ctrl, verbose = FALSE, dateWarning = FALSE)
      }, error = function(e){
        list(error = e$message)
      })

      if(!is.null(warm_fit$error)){
        if(verbose) message(crayon::red(paste("FAILED:", warm_fit$error)))
        results[[length(results) + 1L]] <- data.frame(
          workload = wl_name,
          structure = wl$dimensions$structure,
          n_levels = wl$dimensions$n_levels,
          n_ind = wl$dimensions$n_ind,
          relationship = wl$dimensions$relationship,
          engine = eng_name,
          elapsed_sec = NA_real_,
          iterations = NA_integer_,
          converged = FALSE,
          error = warm_fit$error,
          stringsAsFactors = FALSE
        )
        next
      }

      # Timed replications
      reps <- max(1L, as.integer(replications))
      times <- numeric(reps)
      last_fit <- warm_fit

      for(r in seq_len(reps)){
        t0 <- proc.time()
        last_fit <- mmesFn(fixed = wl$fixed, random = wl$random, rcov = wl$rcov, data = wl$data,
                           REML = wl$REML, computeCi = wl$computeCi, nIters = nIters,
                           engine = ctrl, verbose = FALSE, dateWarning = FALSE)
        t1 <- proc.time()
        times[r] <- (t1 - t0)[3L]
      }

      med_time <- stats::median(times)
      n_conv <- if(!is.null(last_fit$monitor)) ncol(last_fit$monitor) else NA_integer_
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
      parity_ok <- TRUE
      max_rel_diff <- 0.0
      if(eng_name == ref_engine_name){
        ref_fit <- last_fit
      } else if(verifyParity && !is.null(ref_fit)){
        # Compare likelihood
        rel_llik <- abs(last_fit$llik - ref_fit$llik) / (abs(ref_fit$llik) + 1.0)
        # Compare variance components
        v_last <- unlist(last_fit$theta)
        v_ref <- unlist(ref_fit$theta)
        rel_theta <- if(length(v_last) == length(v_ref)){
          max(abs(v_last - v_ref) / (abs(v_ref) + 1.0))
        } else 1.0
        # Compare BLUPs
        u_last <- as.numeric(last_fit$bu)
        u_ref <- as.numeric(ref_fit$bu)
        rel_bu <- if(length(u_last) == length(u_ref)){
          max(abs(u_last - u_ref) / (abs(u_ref) + 1.0))
        } else 1.0

        max_rel_diff <- max(rel_llik, rel_theta, rel_bu, na.rm=TRUE)
        parity_ok <- (max_rel_diff <= tolerance)
        if(!parity_ok && verbose){
          message(crayon::yellow(sprintf("    Warning: Parity deviation = %.2e (tolerance = %.2e)",
                                         max_rel_diff, tolerance)))
        }
      }

      results[[length(results) + 1L]] <- data.frame(
        workload = wl_name,
        structure = wl$dimensions$structure,
        n_levels = wl$dimensions$n_levels,
        n_ind = wl$dimensions$n_ind,
        relationship = wl$dimensions$relationship,
        engine = eng_name,
        elapsed_sec = med_time,
        iterations = n_conv,
        converged = isTRUE(last_fit$converged),
        blockSchurActive = isTRUE(diag$blockSchurActive),
        blockBandedActive = isTRUE(diag$blockBandedActive),
        factorSchurActive = isTRUE(diag$factorSchurActive),
        parity_ok = parity_ok,
        max_rel_diff = max_rel_diff,
        error = NA_character_,
        stringsAsFactors = FALSE
      )
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
  merged$advantageous <- merged$speedup > 1.05 & merged$parity_ok.cand

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
