# Sweep Axis B: Multi-Factor Latent Woodbury Schur vs Baseline Generic CHOLMOD
.libPaths(c("/tmp/sommer-review-implementation-lib", .libPaths()))
library(sommer)

cat("=== Axis B Sweep: Multi-Factor Latent Woodbury Schur vs Generic CHOLMOD ===\n")

env_vec  <- c(4L, 8L)
rep_vec  <- c(2L, 3L)
id_vec   <- c(30L, 80L, 160L)

results <- list()

for(n_env in env_vec){
  for(n_rep in rep_vec){
    for(n_id in id_vec){
      fixture_name <- sprintf("FA1_MF_%dEnv_%dRep_%dID", n_env, n_rep, n_id)
      cat(sprintf("Running %-28s ... ", fixture_name))
      
      set.seed(42L)
      dataset <- expand.grid(
        id = factor(seq_len(n_id)),
        environment = factor(seq_len(n_env)),
        replicate = factor(seq_len(n_rep))
      )
      dataset$response <- sin(seq_len(nrow(dataset))) + as.integer(dataset$id) / n_id + stats::rnorm(nrow(dataset))
      
      M <- matrix(stats::rnorm(n_id * 30L), n_id, 30L)
      G <- stats::cor(t(M)) + diag(0.1, n_id)
      dimnames(G) <- list(levels(dataset$id), levels(dataset$id))
      precision <- solve(G)
      attr(precision, "inverse") <- TRUE
      
      fit_base <- tryCatch({
        t0 <- proc.time()
        m_base <- mmes(
          response ~ environment,
          random = ~ vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu = precision),
          data = dataset,
          solver = "cholmod",
          nIters = 4L,
          computeCi = 0L,
          getPEV = FALSE,
          verbose = FALSE,
          dateWarning = FALSE,
          engine = mmes_control(solver = "cholmod", experimental = list(factorSchurMultiFactor = "off"))
        )
        t1 <- proc.time()
        list(time = (t1 - t0)[3L], fit = m_base)
      }, error = function(e) list(error = e$message))
      
      fit_latent <- tryCatch({
        t0 <- proc.time()
        m_lat <- mmes(
          response ~ environment,
          random = ~ vsm(fam(environment, 1L), dsm(replicate), ism(id), Gu = precision),
          data = dataset,
          solver = "cholmod",
          nIters = 4L,
          computeCi = 0L,
          getPEV = FALSE,
          verbose = FALSE,
          dateWarning = FALSE,
          engine = mmes_control(solver = "cholmod", experimental = list(factorSchurMultiFactor = "force"))
        )
        t1 <- proc.time()
        list(time = (t1 - t0)[3L], fit = m_lat)
      }, error = function(e) list(error = e$message))
      
      if(is.null(fit_base$error) && is.null(fit_latent$error)){
        t_base <- fit_base$time
        t_lat  <- fit_latent$time
        speedup <- t_base / t_lat
        active <- fit_latent$fit$engineDiagnostics$factorSchurActive
        
        # Parity check
        rel_llik <- abs(tail(fit_latent$fit$llik, 1L) - tail(fit_base$fit$llik, 1L)) / (abs(tail(fit_base$fit$llik, 1L)) + 1.0)
        rel_bu <- max(abs(fit_latent$fit$bu - fit_base$fit$bu) / (abs(fit_base$fit$bu) + 1.0))
        parity <- (rel_llik < 1e-4 && rel_bu < 1e-4)
        
        cat(sprintf("Base: %.4fs | Latent: %.4fs | Speedup: %.2fx | Active: %s | Parity: %s\n",
                    t_base, t_lat, speedup, active, parity))
        
        results[[length(results) + 1L]] <- data.frame(
          env = n_env,
          rep = n_rep,
          n_id = n_id,
          total_effects = n_env * n_rep * n_id,
          t_base = t_base,
          t_lat = t_lat,
          speedup = speedup,
          active = active,
          parity = parity
        )
      } else {
        cat("ERROR:", fit_base$error, fit_latent$error, "\n")
      }
    }
  }
}

df_res_b <- do.call(rbind, results)
print(df_res_b)
saveRDS(df_res_b, "/tmp/sweep_axis_b.rds")
