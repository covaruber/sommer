# Sweep Axis A: Block-Banded Cholesky vs Baseline Generic CHOLMOD
.libPaths(c("/tmp/sommer-review-implementation-lib", .libPaths()))
library(sommer)

cat("=== Axis A Sweep: Block-Banded vs Generic CHOLMOD ===\n")

levels_vec <- c(8L, 20L, 50L, 100L)
inds_vec   <- c(20L, 60L)
structures <- c("ar1", "ar2")

results <- list()

for(struct in structures){
  for(n_l in levels_vec){
    for(n_i in inds_vec){
      wl_name <- sprintf("%s_%dL_%dID", struct, n_l, n_i)
      cat(sprintf("Running %-20s ... ", wl_name))
      
      wl <- generate_mme_workload(
        name = wl_name,
        n_ind = n_i,
        n_levels = n_l,
        structure = struct,
        relationship = "none",
        reml = TRUE,
        computeCi = 0L,
        seed = 42L
      )
      
      engines <- list(
        cholmod_base = mmes_control(
          solver = "cholmod",
          experimental = list(blockBanded = "off")
        ),
        banded_force = mmes_control(
          solver = "cholmod",
          experimental = list(blockBanded = "force")
        )
      )
      
      bm <- tryCatch({
        run_mme_benchmark(
          workloads = wl,
          engines = engines,
          replications = 3L,
          verifyParity = TRUE,
          tolerance = 1e-5,
          nIters = 4L,
          verbose = FALSE
        )
      }, error = function(e){
        cat("ERROR:", e$message, "\n")
        NULL
      })
      
      if(!is.null(bm)){
        res <- bm$results
        t_base <- res$elapsed_sec[res$engine == "cholmod_base"]
        t_band <- res$elapsed_sec[res$engine == "banded_force"]
        active <- res$blockBandedActive[res$engine == "banded_force"]
        parity <- res$parity_ok[res$engine == "banded_force"]
        speedup <- t_base / t_band
        cat(sprintf("Base: %.4fs | Banded: %.4fs | Speedup: %.2fx | Active: %s | Parity: %s\n",
                    t_base, t_band, speedup, active, parity))
        results[[length(results) + 1L]] <- data.frame(
          structure = struct,
          levels = n_l,
          n_ind = n_i,
          t_base = t_base,
          t_band = t_band,
          speedup = speedup,
          active = active,
          parity = parity
        )
      }
    }
  }
}

df_res <- do.call(rbind, results)
print(df_res)
saveRDS(df_res, "/tmp/sweep_axis_a.rds")
