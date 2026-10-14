# Fit the Schur-versus-sparse decision on FILL quantities (not on the generator's input knob).
# Usage: Rscript benchmarks/fit_fill_rule.R /tmp/calibrate_medium.csv [more.csv ...]
# Features, all from the symbolic Cholesky of the coefficient matrix (AMD order):
#   fill_fraction = nnz(L) / (n (n + 1) / 2)   how dense the sparse factor becomes
#   sparse_flops  = sum of squared column counts of L
#   schur_cost    = dense block plan cost (BlockSchurPolicy::cost); recomputed here from the
#                   recorded group count and border when the run could not record it (chain engine)
args <- commandArgs(trailingOnly=TRUE)
if(!length(args)) stop("usage: fit_fill_rule.R results.csv [...]", call.=FALSE)
data <- do.call(rbind, lapply(args, function(path){
  d <- read.csv(path, stringsAsFactors=FALSE)
  if("density" %in% names(d)) names(d)[names(d) == "density"] <- "input_nonzero"
  if(!"pattern" %in% names(d)) d$pattern <- "random"
  d[, intersect(c("pattern", "seed", "engine_subtype", "contrast_valid", "parity_relative", "group_size", "n_groups", "input_nonzero", "border", "time_on", "time_off",
                  "sparseNnzL", "sparseFlops", "schurCost", "blockSchurGroups", "blockSchurBorder", "nEffects"), names(d))]
}))

if("contrast_valid" %in% names(data)){
  data <- data[data$contrast_valid & data$parity_relative <= 1e-6, ]
}
data$fill_fraction <- data$sparseNnzL / (data$nEffects * (data$nEffects + 1) / 2)
data$log_ratio <- log(data$time_off / data$time_on)
data$log_sparse <- log(data$sparseFlops)
data$log_schur <- log(data$schurCost)
data$log_effects <- log(data$nEffects)
if("seed" %in% names(data)){
  usable <- is.finite(data$log_sparse) & is.finite(data$log_ratio) & data$fill_fraction > 0
  for(subtype in unique(data$engine_subtype[usable])){
    subset <- data[usable & data$engine_subtype == subtype, ]
    training <- subset[subset$seed != 303L, ]
    validation <- subset[subset$seed == 303L, ]
    if(nrow(training) < 8L || !nrow(validation)) next
    model <- lm(log_ratio ~ log_sparse + log_effects + log(fill_fraction) +
                  log(blockSchurBorder + 1) + log(blockSchurGroups + 1), data=training)
    predicted <- predict(model, validation)
    cat("\nSubtype:", subtype, "training cells:", nrow(training),
        "held-out seed-303 cells:", nrow(validation), "\n")
    print(summary(model)$coefficients)
    for(margin in c(0, log(1.2), log(1.5))){
      admitted <- is.finite(predicted) & predicted > margin
      cat("margin", round(margin, 3), "held-out admits", sum(admitted),
          "worst admitted ratio", if(any(admitted)) exp(min(validation$log_ratio[admitted])) else NA,
          "missed wins over 1.2x", sum(!admitted & validation$log_ratio > log(1.2)), "\n")
    }
  }
}
ok <- is.finite(data$log_sparse) & is.finite(data$log_schur) & is.finite(data$log_ratio)
d <- data[ok, ]
cat(sprintf("Dense-plan cost comparisons available: %d of %d cells; no costs imputed.\n", nrow(d), nrow(data)))
if(nrow(d) < 8L) quit(status=0)

cat("\nFill versus the generator's input knob (median fill fraction by input nonzero fraction and pattern):\n")
print(aggregate(fill_fraction ~ pattern + input_nonzero, d, median), row.names=FALSE)

models <- list(
  ratio_only = log_ratio ~ I(log_sparse - log_schur),
  two_term = log_ratio ~ log_sparse + log_schur,
  two_term_fill = log_ratio ~ log_sparse + log_schur + fill_fraction,
  two_term_size = log_ratio ~ log_sparse + log_schur + log_effects)
for(name in names(models)){
  fit <- lm(models[[name]], data=d)
  cat(sprintf("\n[%s] R2=%.2f adjR2=%.2f residual sd=%.2f (log scale)\n", name,
              summary(fit)$r.squared, summary(fit)$adj.r.squared, summary(fit)$sigma))
  print(round(summary(fit)$coefficients, 3))
}

# Decision quality of the best simple rule: predict Schur wins when predicted log ratio > margin.
best <- lm(log_ratio ~ log_sparse + log_schur, data=d)
d$pred <- predict(best, d)
for(margin in c(0, 0.1, 0.2, 0.3)){
  admit <- d$pred > margin
  cat(sprintf("margin %.1f: admits %d of %d cells; worst admitted measured ratio %.2f; worst rejected missed win %.2f\n",
              margin, sum(admit), nrow(d), if(any(admit)) exp(min(d$log_ratio[admit])) else NA,
              if(any(!admit)) exp(max(d$log_ratio[!admit])) else NA))
}
cat("\nCells where the two-term fit is most wrong (measured vs predicted ratio):\n")
d$measured <- exp(d$log_ratio)
d$predicted <- exp(d$pred)
worst <- head(d[order(-abs(d$log_ratio - d$pred)), c("pattern", "group_size", "n_groups", "input_nonzero", "border",
                                                     "fill_fraction", "measured", "predicted")], 8)
print(worst, row.names=FALSE, digits=3)
