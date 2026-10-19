
pmonitor <- function(object,...){
  x <- object$monitor
  mmin <- min(x)
  mmax <- max(x)
  plot(x[1,], type="l", ylim=c(mmin,mmax))
  if(nrow(x) > 1){
    for(i in 2:nrow(x)){
      par(new=TRUE)
      plot(x[i,], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="l",...)
    }
  }
  legend("topright",legend = rownames(x), col=1:(nrow(x)), bty="n", lty=1)
}




