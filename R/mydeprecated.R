mmer2 <- function(fixed, random, rcov, data, weights, 
                  iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                  init=NULL, constraints=NULL, method="NR", 
                  getPEV=TRUE,
                  na.method.X="exclude",
                  na.method.Y="exclude",
                  return.param=FALSE, 
                  date.warning=TRUE,
                  verbose=TRUE,reshape.output=TRUE){
  .Deprecated("mmer",
              msg="This function is deprecated. Use 'mmer' instead.\nNow the mmer function can run both types of models; \nformula-based and matrix-based models. Type ?mmer")
  stop("DEPRECATED FUNCTION",call. = FALSE)
}

GWAS2 <- function(fixed, random, rcov, data, weights, 
                 iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                 init=NULL, constraints=NULL, method="NR", 
                 getPEV=TRUE,
                 na.method.X="exclude",
                 na.method.Y="exclude",
                 return.param=FALSE, 
                 date.warning=TRUE,
                 verbose=TRUE,
                 M=NULL, gTerm=NULL, n.PC = 0, min.MAF = 0.05, 
                 n.core=1, P3D = TRUE){
  
  .Deprecated("GWAS",
              msg="This function is deprecated. Use 'GWAS' instead.\nNow the GWAS function can run both types of GWAS models; \nformula-based and matrix-based models. Type ?GWAS")
  stop("DEPRECATED FUNCTION",call. = FALSE)
  
}

