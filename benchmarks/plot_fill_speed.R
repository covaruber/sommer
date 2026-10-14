args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)) args[1L] else "benchmarks/results/fill-current.csv"
stem <- if(length(args) > 1L) args[2L] else "benchmarks/results/fill-speed"
data <- read.csv(input, stringsAsFactors=FALSE)
data <- data[data$contrast_valid & data$parity_relative <= 1e-6 &
               is.finite(data$sparseNnzL) & data$sparseNnzL > 0 &
               data$time_on > 0 & data$time_off > 0, ]
colors <- c(on="#00796b", off="#b64527")
subtypes <- intersect(c("dense", "chain"), unique(data$engine_subtype))
dir.create(dirname(stem), recursive=TRUE, showWarnings=FALSE)

draw_speed <- function(){
  par(mfrow=c(1, length(subtypes)), mar=c(5, 5, 4, 1), las=1)
  for(subtype in subtypes){
    subset <- data[data$engine_subtype == subtype, ]
    fill_grid <- exp(seq(log(min(subset$sparseNnzL)), log(max(subset$sparseNnzL)), length.out=200L))
    plot(NA, xlim=range(subset$sparseNnzL), ylim=range(1 / subset$time_on, 1 / subset$time_off),
         log="xy", xlab="Predicted factor entries: nnz(L) [log scale]",
         ylab="Fits per second [log scale; higher is faster]",
         main=paste(if(subtype == "dense") "Dense Schur" else "Block chain", "vs generic sparse"))
    grid(col="#e5e5e5")
    for(mode in c("on", "off")){
      speed <- 1 / subset[[paste0("time_", mode)]]
      points(subset$sparseNnzL, speed, pch=if(mode == "on") 16 else 1,
             col=adjustcolor(colors[mode], alpha.f=0.45), cex=0.8)
      fit <- lm(log(speed) ~ log(subset$sparseNnzL))
      lines(fill_grid, exp(coef(fit)[1L] + coef(fit)[2L] * log(fill_grid)),
            col=colors[mode], lwd=2, lty=if(mode == "on") 1 else 2)
    }
    legend("topright", c("Block engine enabled", "Block engine disabled"),
           col=colors, pch=c(16, 1), lty=c(1, 2), lwd=2, bty="n", cex=0.85)
    mtext(sprintf("%d paired cells; descriptive log-log fits", nrow(subset)), side=3, line=0.3, cex=0.8)
  }
}

draw_ratio <- function(){
  par(mfrow=c(1, length(subtypes)), mar=c(5, 5, 4, 1), las=1)
  for(subtype in subtypes){
    subset <- data[data$engine_subtype == subtype, ]
    subset$ratio <- subset$time_off / subset$time_on
    plot(subset$sparseNnzL, subset$ratio, log="xy", pch=16, ylim=range(subset$ratio, 1, 1.2),
         col=adjustcolor(colors["on"], alpha.f=0.5),
         xlab="Predicted factor entries: nnz(L) [log scale]",
         ylab="Paired speedup: off / on [higher favors block engine]",
         main=if(subtype == "dense") "Dense Schur: paired speedup" else "Block chain: paired speedup")
    abline(h=1, lty=2, col="#333333")
    abline(h=1.2, lty=3, col="#777777")
    model <- lm(log(ratio) ~ log(sparseNnzL), data=subset)
    fill_grid <- exp(seq(log(min(subset$sparseNnzL)), log(max(subset$sparseNnzL)), length.out=200L))
    predictions <- predict(model, newdata=data.frame(sparseNnzL=fill_grid), interval="confidence")
    lines(fill_grid, exp(predictions[, "fit"]), col=colors["on"], lwd=2)
    lines(fill_grid, exp(predictions[, "lwr"]), col=colors["on"], lty=3)
    lines(fill_grid, exp(predictions[, "upr"]), col=colors["on"], lty=3)
    crossing <- exp(-coef(model)[1L] / coef(model)[2L])
    inside <- is.finite(crossing) && crossing >= min(fill_grid) && crossing <= max(fill_grid)
    cat(sprintf("%s: n=%d, paired log-log slope=%.4f, R2=%.3f, fitted equality at nnz(L)=%.0f (%s observed range); wins >1.2x=%d\n",
                subtype, nrow(subset), coef(model)[2L], summary(model)$r.squared,
                crossing, if(inside) "inside" else "outside", sum(subset$ratio > 1.2)))
    mtext("Equality line = 1; dotted horizontal line = 1.2x", side=3, line=0.3, cex=0.8)
  }
}

png(paste0(stem, ".png"), width=1800, height=850, res=130)
draw_speed()
invisible(dev.off())
pdf(paste0(stem, ".pdf"), width=14, height=6.5)
draw_speed()
invisible(dev.off())
png(paste0(stem, "-paired.png"), width=1800, height=850, res=130)
draw_ratio()
invisible(dev.off())