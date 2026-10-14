args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)) args[1L] else "benchmarks/results/fill-current.csv"
stem <- if(length(args) > 1L) args[2L] else "benchmarks/results/fill-factorial"
data <- read.csv(input, stringsAsFactors=FALSE)
data <- data[is.finite(data$time_on) & is.finite(data$time_off) &
               data$time_on > 0 & data$time_off > 0 & data$parity_relative <= 1e-6, ]
names(data)[names(data) == "density"] <- "input_nonzero"
design_factors <- c("group_size", "n_groups", "input_nonzero", "border", "pattern")
for(variable in c(design_factors, "seed")) data[[variable]] <- factor(data[[variable]])
data$log_speedup <- log(data$time_off / data$time_on)
data$log_fill <- log(data$sparseNnzL)
data$geometry <- do.call(paste, c(data[design_factors], sep="_"))
dir.create(dirname(stem), recursive=TRUE, showWarnings=FALSE)
factor_terms <- paste(design_factors, collapse=" + ")
main_formula <- as.formula(paste("log_speedup ~ seed +", factor_terms))
interaction_formula <- as.formula(paste("log_speedup ~ seed + (", factor_terms, ")^2"))
main_model <- lm(main_formula, data=data)
interaction_model <- lm(interaction_formula, data=data)
cat("Available cells:", nrow(data), "of 360 planned; valid active contrasts:", sum(data$contrast_valid), "\n")
cat("Geometry cells:", length(unique(data$geometry)), "of 120 planned.\n")
cat("Primary analysis includes inactive force attempts as no-op controls; engine subtype is an outcome, not a randomized factor.\n")
cat("Incomplete randomized harvest: use nested partial F tests, not sequential sums of squares.\n")
cat("Timing inference is conditional on this machine and worker contention; graph seeds are blocks, timing repetitions are not independent graphs.\n\n")
cat("Main effects versus two-factor interactions:\n")
print(anova(main_model, interaction_model))
cat(sprintf("Main-effects R2=%.3f; interaction R2=%.3f; interaction adjusted R2=%.3f\n",
            summary(main_model)$r.squared, summary(interaction_model)$r.squared,
            summary(interaction_model)$adj.r.squared))
cat("Interaction-model rank:", interaction_model$rank, "of", ncol(model.matrix(interaction_model)), "\n")
if(interaction_model$rank != ncol(model.matrix(interaction_model))){
  stop("Factorial model is rank-deficient; do not derive rules from it.")
}

terms <- attr(terms(interaction_model), "term.labels")
factor_tests <- lapply(design_factors, function(variable){
  keep <- terms[!vapply(strsplit(terms, ":", fixed=TRUE), function(parts) variable %in% parts, logical(1))]
  reduced <- lm(reformulate(keep, response="log_speedup"), data=data)
  comparison <- anova(reduced, interaction_model)
  data.frame(factor=variable, df=comparison$Df[2L], partial_ss=comparison$`Sum of Sq`[2L],
             F=comparison$F[2L], p=comparison$`Pr(>F)`[2L])
})
factor_tests <- do.call(rbind, factor_tests)
cat("\nOmnibus contribution of each factor including its interactions:\n")
print(factor_tests[order(-factor_tests$partial_ss), ], row.names=FALSE)
interaction_tests <- drop1(interaction_model, test="F")
interaction_tests <- interaction_tests[grepl(":", rownames(interaction_tests)), , drop=FALSE]
cat("\nConditional two-factor interaction tests:\n")
print(interaction_tests[order(interaction_tests$`Sum of Sq`, decreasing=TRUE), , drop=FALSE])

grid <- expand.grid(lapply(data[c(design_factors, "seed")], levels), stringsAsFactors=FALSE)
for(variable in names(grid)) grid[[variable]] <- factor(grid[[variable]], levels=levels(data[[variable]]))
grid$predicted <- predict(interaction_model, grid)
marginals <- do.call(rbind, lapply(design_factors, function(variable){
  values <- aggregate(grid$predicted, list(level=grid[[variable]]), mean)
  data.frame(factor=variable, level=values$level, geometric_speedup=exp(values$x))
}))
cat("\nEqual-weight estimated marginal speedups over the full crossed grid:\n")
print(marginals, row.names=FALSE)
write.csv(marginals, paste0(stem, "-marginals.csv"), row.names=FALSE)

training <- droplevels(data[data$seed != "303", ])
validation <- data[data$seed == "303", ]
predict_seed_average <- function(model, newdata){
  estimates <- lapply(levels(training$seed), function(graph_seed){
    changed <- newdata
    changed$seed <- factor(graph_seed, levels=levels(training$seed))
    predict(model, changed)
  })
  rowMeans(do.call(cbind, estimates))
}
fit_fill <- lm(log_speedup ~ seed + log_fill, data=training)
fit_main <- lm(main_formula, data=training)
fit_interactions <- lm(interaction_formula, data=training)
if(fit_interactions$rank != ncol(model.matrix(fit_interactions))) stop("Training model is rank-deficient")
models <- list(fill_only=fit_fill, main_effects=fit_main, two_factor=fit_interactions)
metrics <- do.call(rbind, lapply(names(models), function(name){
  prediction <- predict_seed_average(models[[name]], validation)
  validation[[paste0("prediction_", name)]] <<- exp(prediction)
  admission <- prediction >= log(1.2)
  data.frame(model=name, heldout_cells=nrow(validation),
             log_RMSE=sqrt(mean((prediction - validation$log_speedup)^2)),
             admitted=sum(admission),
             admitted_losses_over_10pct=sum(admission & validation$log_speedup < log(0.9)),
             worst_admitted=if(any(admission)) exp(min(validation$log_speedup[admission])) else NA_real_,
             missed_wins_over_20pct=sum(!admission & validation$log_speedup > log(1.2)))
}))
cat("\nSeed-303 validation; decision threshold is predicted speedup >= 1.2:\n")
print(metrics, row.names=FALSE)
write.csv(metrics, paste0(stem, "-validation.csv"), row.names=FALSE)
write.csv(validation, paste0(stem, "-heldout.csv"), row.names=FALSE)
saveRDS(list(main=main_model, interactions=interaction_model, factor_tests=factor_tests,
             interaction_tests=interaction_tests, validation=metrics), paste0(stem, ".rds"))

png(paste0(stem, "-interactions.png"), width=1800, height=1000, res=130)
par(mfrow=c(2, 2), mar=c(5, 5, 4, 1), las=1)
for(condition in levels(grid$border)) for(graph_pattern in levels(grid$pattern)){
  panel <- grid[grid$border == condition & grid$pattern == graph_pattern, ]
  means <- aggregate(predicted ~ group_size + n_groups, panel, mean)
  colors <- c("#00796b", "#b64527")
  plot(NA, xlim=range(as.numeric(levels(panel$group_size))),
       ylim=range(exp(means$predicted), 1), log="xy", xlab="Block size",
       ylab="Adjusted speedup: off / on", main=paste(graph_pattern, "pattern; border", condition))
  abline(h=1, lty=3, col="#777777")
  for(index in seq_along(levels(panel$n_groups))){
    values <- means[means$n_groups == levels(panel$n_groups)[index], ]
    lines(as.numeric(as.character(values$group_size)), exp(values$predicted),
          type="b", pch=16, col=colors[index], lwd=2)
  }
  legend("topleft", paste(levels(panel$n_groups), "groups"), col=colors, lty=1, pch=16, bty="n")
}
invisible(dev.off())