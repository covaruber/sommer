#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction

"predict.mmec" <- function(object, Dtable=NULL, D=NULL,...){
  
  
  if(is.null(D)){ ## if user don't provide a D but a Dtable
    ## get the names of terms used for each effect from the Dtable
    termNames <- apply(data.frame(Dtable[,"term"]),1,function(x){
      y<-paste(all.vars(as.formula(paste("~",x))),collapse = ":")
      z<-intersect(strsplit(y, split=":")[[1]], colnames(object$data))
      return(z)
    })
    ## get from the Dtable which term will be used for the hypertable classification
    forD <- Dtable[which(Dtable$D),"term"]
    forD <- paste(all.vars(as.formula(paste("~",forD))),collapse = ":")
    classify <- intersect(strsplit(forD, split=":")[[1]], colnames(object$data))
    classifyL <- as.list(classify)
    object$data[,"Intercept"] <- "1"
    object$data[,"1"] <- "Intercept"
    levelsClassifyL <- lapply(classifyL, function(x){sort(unique(object$data[,x]))})
    ## build a dataset for building the D matrix
    dataPredict <- do.call(expand.grid, levelsClassifyL )
    colnames(dataPredict) <- classify
    ## D matrix
    Dformed <- Matrix(0, nrow(dataPredict), length(object$b)+length(object$u))
    colnames(Dformed) <- c(rownames(object$b),rownames(object$u))
    
    ## for each partition of X and Z we move one by one
    for(k in 1:length(object$partitionsAll)){ 
      # which fixed or random effects belong to this partition 
      # will be connected to the Dtable
      buToUse <- object$partitionsAll[[k]] 
      # now for each row of the dataPredict or row in the D matrix we...
      for(i in 1:nrow(dataPredict)){
        tk <- termNames[[k]] # obtain the term name
        isPresent <- intersect(tk,names(dataPredict)) # check which are present in the dataset
        if(length(isPresent) > 0){ # if the term is present
          prov <- as.matrix(dataPredict[i,isPresent]) # get what values are we looking for
          myMatch <- list() # move column by column or value by value to see where is the match in this paritition
          for(j in 1:ncol(prov)){
            myMatch[[j]] <- buToUse[grep(prov[,j], colnames(Dformed)[buToUse] , fixed = TRUE)]
          }
          ## the number of 1's in each row needs to be equal to the number of classify levels
          if(length(isPresent) > 1){ # if there was more than one term
            fill <- do.call(intersect,myMatch)
          }else{ # if it was a single term
            fill <- unlist(myMatch)
          }
        }else{ # if the term for partition k is not part of the dataPredict
          # very likely is going to be ignored so assume all bu's can be used
          fill <- as.vector(buToUse) 
        }
        # if there is fixed/random effects to be used
        if(length(fill) > 0){
          # check if should be included or averaged
          if(Dtable[k,"include"]){
            Dformed[i,fill] = 1  # include
          }# end for including
          if(Dtable[k,"average"]){
            Dformed[i,fill] = 1 / (length(buToUse[1,]))
            # Dformed[i,fill] = 1 / length(unlist(myMatch)) # averaged
          }# end for averaging
          # if we're averaging and we identify an intercept it should be 1 for everyone
        }# end for if there's something to fill
        int <- grep("Intercept",colnames(Dformed))
        int2 <- intersect(int,buToUse)
        if(length(int2) > 0){
          Dformed[i,int] = 1
        }# end for adding intercept
      } # end of for each row
    } # end of for each partition
    
  }else{ Dformed <- D}
  
  Ci <- object$Ci
  
  vcov <- Dformed%*%Ci%*%t(Dformed)
  se <- as.vector(sqrt(diag(vcov))) 
  
  predicted.value <- as.vector(Dformed %*% object$bu)
  
  dataPredict$predicted.value <- predicted.value
  dataPredict$se <- se
  
  return(list(pvals=dataPredict, vcov=vcov, D=Dformed, Ci=Ci))
  
}

"print.predict.mmec"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
                 The predictions are obtained by averaging/aggregating across
                 the hypertable calculated from model terms constructed solely
                 from factors in the include sets. You can customize the model
                 terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}