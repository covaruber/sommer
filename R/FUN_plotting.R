bathy.colors <- function(n, alpha=1){
  return(rgb(seq(0.9,0,len=n), seq(0.9,0,len=n), 1, alpha))
}
big.peaks.col <-function(x, tre){
  r <- rle(x)
  v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
  pos <- v[which(x[v] > tre)] #position of the real peaks
  hei <- x[pos]
  out <- list(pos=pos,hei=hei)
  return(out)
}
hits <- function(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE, plotting=TRUE){
  
  bagi <- function(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE, plotting=TRUE){
    #####
    
    if(is.null(gwasm$W)){
      #cat()
      stop("A GWAS model is needed to create the design matrix for bagging",call.=FALSE)
    }else{ ######## IF MODELS WAS A GWAS MODEL ###########
      
      if(pick){ ############# 'PICK=TRUE' ARGUMENT #################
        cat("Please move the cursor to the GWAS plot, \nclick over the marker dots you want in the design matrix \nand press 'Esc' key when done")
        plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20)
        
        picked <- locator()$x
        marker <- round(picked)
        res <- big.peaks.col(gwasm$W.scores$additive, tre=threshold)
        marker2 <- numeric(length(marker))
        for(t in 1:length(marker)){
          ## in the neighbour of the peak selected
          ff <- which(res$pos %in% c(c(marker[t]-50):c(marker[t]+50)))
          ##
          marker2[t] <- res$pos[ff[which(res$hei[ff] == max(res$hei[ff]))[1]]]
        }
        
        abline(v=marker2, lty=3, col="red")
        legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
        
        XX <- as.matrix(gwasm$W[,marker])
        X1 <- as.data.frame(apply(XX,2, as.factor))
        xvar <- colnames(X1)
        X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
        
      }else{ ############# 'PICK=FALSE' ARGUMENT #################
        big.peaks.col <- function(x, tre){
          r <- rle(x)
          v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
          pos <- v[which(x[v] > tre)] #position of the real peaks
          hei <- x[pos]
          out <- list(pos=pos,hei=hei)
          return(out)
        }
        make.full <- function(X) {
          svd.X <- svd(X)
          r <- max(which(svd.X$d > 1e-08))
          return(as.matrix(svd.X$u[, 1:r]))
        }
        ########################
        # plot(transfft(gwasm$W.scores$additive,.02), type="l")
        # abline(v=res$pos, lty=3,col="red")
        #res <- big.peaks.col(transfft(gwasm$W.scores$additive,.02),0) # smoothed
        #hh <- which(res$hei %in% sort(res$hei, decreasing=TRUE)[1:nmar])
        #res <- list(pos=res$pos[hh], hei=res$hei[hh])
        
        condicion <- gwasm$map
        if(!is.null(condicion)){ # if map is present
          res <- big.peaks.col(as.vector(gwasm$map$p.val), tre=threshold)
        }else{
          res <- big.peaks.col(as.vector(gwasm$W.scores$additive), tre=threshold)
        }
        
        
        if(length(res$pos) >= nmar){ # if enough markers significant
          ## build the number of cluster plus 5 to make sure you have repeatability
          if(length(res$pos) <= nmar+15){
            if(length(res$pos) < nmar){
              res$clus <- kmeans(res$pos,length(res$pos))$cluster 
            }else{
              res$clus <- kmeans(res$pos,nmar-1)$cluster 
            }
            
          }else{
            res$clus <- kmeans(res$pos,nmar+15)$cluster 
          }
          
          heights <- numeric()
          for(i in 1:max(res$clus)){
            vv <- which(res$clus == i)
            heights[i] <- (res$hei[vv][which(res$hei[vv] == max(res$hei[vv]))])[1]
          }
          ## the top clusters (tallest peaks)
          maxo <- which(heights %in% sort(heights, decreasing = TRUE)[1:nmar])
          ## subset the 'res' object
          good <- which(res$clus %in% maxo)
          res2 <- list(pos=res$pos[good], hei=res$hei[good], clus=res$clus[good])
          
          marker <- numeric()
          for(i in unique(res2$clus)){
            vv <- which(res2$clus == i)
            marker[i] <- (res2$pos[vv][which(res2$hei[vv] == max(res2$hei[vv]))])[1]
          }
          marker <- as.numeric(na.omit(marker))
          #marker <- res$pos
          #########
          #########
          
          if(!is.null(condicion)){ #if map was present
            if(plotting){
              col.scheme <- rep((transp(c("cadetblue","red"))),30)
              plot(gwasm$map$p.val, bty="n", col=col.scheme[factor(gwasm$map$Chrom, levels = unique(gwasm$map$Chrom, na.rm=TRUE))], xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), pch=20, cex=2, las=2)
              init.mrks <- apply(data.frame(unique(gwasm$map$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=gwasm$map$Chrom)
              fin.mrks <- apply(data.frame(unique(gwasm$map$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=gwasm$map$Chrom)
              inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
              axis(side=1, at=inter.mrks, labels=paste("Chr",unique(gwasm$map$Chrom),sep=""), cex.axis=.5)
              abline(v=marker, lty=3, col="red")
              legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
            }
            marker <- as.character(gwasm$map$Locus[marker])
          }else{ # if map was not present
            if(plotting){
              plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20, cex=1.3)
              abline(v=marker, lty=3, col="red")
              legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
            }
          }
          
          
          
          if(method=="cluster"){
            XX <- as.matrix(gwasm$W[,marker])
          }
          #########################
          if(method == "maximum"){
            n.mark=nmar
            pval <- gwasm$W.scores$additive
            names(pval) <- colnames(gwasm$W)
            marker <- order(pval)[1:n.mark]  
            XX <- as.matrix(gwasm$W[,marker])
          }
          
          X1 <- as.data.frame(apply(XX,2, as.factor))
          
          dada <- data.frame(y=(gwasm$fitted.y[,1]), X1)
          fit <- lm(as.formula(paste("y~ ", paste(c(colnames(dada)[-1]), collapse="+"))), data=dada)
          step <- MASS::stepAIC(fit,direction="both",trace=FALSE)
          # good markers after stepwise
          xvar <-(as.character(attr(summary(step)$terms, "variables")))[-c(1:2)]
          #head(X1)
          #xvar <- colnames(X1)
          #cat(("Markers selected"))
          #cat(paste(xvar))
          X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
          #X2 <- make.full(X2)
          #X2 <- as.data.frame(as.matrix(X2))
          #fit <- lm(gwasm$fitted.y~ -1 + as.matrix(X2))
          #fit <- lm(make.formula(colnames(X1)),data=data.frame(X1,y=gwasm$fitted.y)) #lm using 100 marks
          
          #step <- stepAIC(fit,direction="both",trace=FALSE)  #forward-backward stepwise regression
          
        }else{ # if not enough markers
          #cat()
          stop("Not enough significant markers in your model to create a design matrix \nwith the number of markers specified by you. Please lower the 'threshold' \nargument or the number of markers in the 'nmar' argument",call. = FALSE)
        } #### end of if enough markers
        
      } ############# END OF 'PICK' ARGUMENT #################
    }
    if(only.mark){
      X2 <- colnames(X1)
    }
    return(X2) #X2
  }
  
  #univariate model
  if(names(gwasm)[1] == "var.comp"){ 
    XXX <- bagi(gwasm=gwasm, nmar=nmar, threshold=threshold, pick=pick, method=method, only.mark=only.mark, plotting=plotting)
  }else{ # univariate in parallel
    XXX <- lapply(gwasm, bagi, nmar=nmar, threshold=threshold, pick=pick, method=method, only.mark=only.mark, plotting=plotting)
  }
  
  return(XXX)
}
jet.colors <- function(n, alpha=1) {
  if(n > 0) {
    if(length(alpha) != 1 & length(alpha) != n) {
      print('Warning: using only first alpha value')
      alpha <- alpha[1]
    }
    if(length(alpha) == 1) {
      alpha <- rep(alpha, n)
    }
    ## TODO Include alpha values
    return(colorRampPalette(c('#000066', 'blue', 'cyan', 'yellow',
                              'red', '#660000'))(n))
  } else {
    ## Return an empty character string if they requested nothing.
    character()
  }
}
map.plot <- function(data, trait=NULL, trait.scale="same", col.chr=NULL, col.trait=NULL, type="hist", cex=0.4, lwd=1, cex.axis=0.4, cex.trait=0.8, jump=5 ){
  sasa <- which(colnames(data) == "Locus")
  if(length(sasa) > 0){
    data$Locus <- as.character(data$Locus) 
  }
  ## transparent function
  ## data needs to have 2 columns; LG and Position
  ## trait needs to indicate the name to plot in the chromosome
  ## the trait can be expressed as "dot", "line" or "hist"
  ## the trait scale can be "same" or "ind", which is same for all or individual
  ## cex is only the cex for the ruler of cM
  ## cex.axis is for the titles
  transp <- function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
  }
  #####
  draw.arc <- function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180, 
                        angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05, col = NA, 
                        lwd = NA, ...) {
    getYmult<-function () {
      if (dev.cur() == 1) {
        warning("No graphics device open.")
        ymult <- 1
      }
      else {
        xyasp <- par("pin")
        xycr <- diff(par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
      }
      return(ymult)
    }
    if (all(is.na(col))) 
      col <- par("col")
    if (all(is.na(lwd))) 
      lwd <- par("lwd")
    xylim <- par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, 
                           lwd, ...) {
      delta.angle <- (angle2 - angle1)
      if (n != as.integer(n)) 
        n <- as.integer(1 + delta.angle/n)
      delta.angle <- delta.angle/n
      angleS <- angle1 + seq(0, length = n) * delta.angle
      angleE <- c(angleS[-1], angle2)
      if (n > 1) {
        half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
        adj.angle = delta.angle * half.lwd.user/(2 * (radius + 
                                                        half.lwd.user))
        angleS[2:n] = angleS[2:n] - adj.angle
        angleE[1:(n - 1)] = angleE[1:(n - 1)] + adj.angle
      }
      p1x <- x + radius * cos(angleS)
      p1y <- y + radius * sin(angleS) * ymult
      p2x <- x + radius * cos(angleE)
      p2y <- y + radius * sin(angleE) * ymult
      segments(p1x, p1y, p2x, p2y, col = col, lwd = lwd, ...)
    }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, radius, angle1, angle2, n, col, 
                       lwd, stringsAsFactors = FALSE)
    for (i in 1:nrow(args)) do.call("draw.arc.0", c(args[i, ], 
                                                    ...))
    invisible(args)
  }
  ## this function takes a dataframe with 2 basic columns:
  ## LG containint the linkage group
  ## Position, the position in cM
  len <- numeric()
  for(i in 1:max(unique(data$LG))){
    len[i] <- max(data[which(data$LG == i),"Position"])
  }
  coree <- which(len == max(len))
  coree1 <- data[which(data$LG == coree),]
  linesss <- coree1$Position / max(coree1$Position)
  # if one trait wants to be plotted or not
  if(!is.null(trait)){
    cols <- 1#(dim(data)[2] - 2)
  }else{cols <- 0}
  extra <- length(len) + cols *length(len)
  fact <- 1/ extra
  fact2 <- fact + (fact*cols) # real separation between chromosomes
  ## colors for traits
  if(!is.null(col.trait)){
    col.trait <- col.trait
  }else{col.trait <- c(1:6,1:6)}
  # colors for chromosomes
  
  
  plot.new()
  for(j in 1:max(unique(data$LG))){
    
    prov <- data[which(data$LG == j),] # extract the jth LG
    ## add zeros to the beggining and end of the LG so the curves look good
    dddd <- prov[1,]
    dddd2 <- prov[dim(prov)[1],]
    dddd[1,which(names(dddd) != "LG")] <- 0
    dddd2[1,which(names(dddd) != "LG" & names(dddd) != "Position")] <- 0
    prov <- rbind(dddd,prov,dddd2)
    ##-----------------------------------------------------------
    ## CHROMOSOME CHUNK
    chr <- prov$Position / max(coree1$Position)
    ### -------------------------------------------------
    ### depending if is a genetic map or physical map we decide how the ruler will be
    #if(max(prov$Position) < 1000){mark <- 10}
    #if(max(prov$Position) > 1000 & max(prov$Position) < 10000){mark <- 100}
    #if(max(prov$Position) > 10000 & max(prov$Position) < 100000){mark <- 1000}
    #if(max(prov$Position) > 100000 & max(prov$Position) < 1000000){mark <- 10000}
    #if(max(prov$Position) < 1000000){mark <- 100000}
    
    ruler <- 1 - (c(seq(0,max(prov$Position), by=jump), round(max(prov$Position),0 )) / max(coree1$Position) ); ruler2 <- c(seq(0,max(prov$Position), by=jump),round(max(prov$Position),0 ))
    if(!is.null(trait)){
      sss <- (fact2*j)-fact# + j*cols 
    }else{sss <- (fact2*j)}# + j*cols}
    
    
    ## heatmap fr density
    dd2 <- density(chr, n=length(chr))$y # regular density
    dd <- sort(density(chr, n=length(chr))$y, decreasing=T)
    ##
    if(!is.null(col.chr)){
      hc <- colorRampPalette(c(col.chr[1], col.chr[2]))( length(dd) )
    }else{hc <- gray.colors(n=length(dd), start = 0, end = 0.6, gamma = 2.2, alpha = NULL)}
    ###### for each chromosome
    for(k in 1:length(dd2)){
      # for each putative point which is the closest position
      ooo <- which(dd == dd2[k])
      lines(y=c(1-chr[k],1-chr[k]), x=c(sss,sss-(fact/3)), lwd=lwd, col=hc[ooo])
    }
    lines(y=c(1,1-max(chr)), x=c(sss,sss), lwd=3)
    lines(y=c(1,1-max(chr)), x=c(sss-(fact/3),sss-(fact/3)), lwd=3)
    text(x=sss-(fact/1.6), y=ruler, labels=ruler2, cex=cex) # cex of the cM ruler
    axis(3,at=(sss-(fact/3)) , labels=paste("LG",j, sep=""), cex.axis=cex.axis, font=2) # cex of the name of LGs
    #axis(2,at=0.275, labels="Density")
    ## --------------------------------------------------------
    draw.arc((sss + sss-(fact/3))/2, 1- max(chr), (sss - (sss-(fact/3)))/2, deg1=180, deg2=360, col="black", lwd=2, lend=1)
    draw.arc((sss + sss-(fact/3))/2, 1, (sss - (sss-(fact/3)))/2, deg1=0, deg2=180, col="black", lwd=2, lend=1)
    ## ----------------
    ## ----------------
    if(!is.null(trait)){
      
      chuy <- which(colnames(prov) %in% trait)
      if(length(chuy)==0){stop("The column indicated as trait is not present in the data provided\nPlease double check your data\n",call. = FALSE)}
      ## ---------------------------
      ## IF TRAIT IS A NUMERIC TRAIT
      ## ---------------------------
      if(is.numeric(prov[,trait])){
        w1 <- which(names(data) == trait) # column of trait to plot dots
        
        if(trait.scale == "same"){
          bobo <- max(data[,trait], na.rm=TRUE)
        }else{bobo <- max(prov[,trait], na.rm=TRUE)}
        
        
        dotss <- fact * ( prov[,trait]/ bobo )
        dotss2 <- sss + (fact/8) + dotss
        ####### for ablines in the trait selected
        sections <- (bobo - min(prov[,trait], na.rm=TRUE))/5
        sections2 <- seq(min(prov[,trait], na.rm=TRUE), bobo, by=sections)
        sections3 <- fact * ( sections2/ bobo )
        sections4 <- sss + (fact/8) + sections3
        ####### if trait is true
        # ablines for different values
        for(d in 1:length(sections4)){
          lines( x=c(sections4[d], sections4[d]), y=c(1,1-max(chr)), col="black", lty=3,lwd=0.5)
          text(x=sections4[d], y=1, labels=round(sections2[d],1), cex=0.4, srt=270)
        }
        
        # plot the trait values
        if(type=="dot"){
          points(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
        }
        if(type == "line"){
          polygon(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.4))
          lines(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
          # density()
          
        }
        if(type == "hist"){
          for(l in 1:length(dotss2)){
            # y is the position in the chromosome
            # x is how long is the line
            lines( x=c(sss + (fact/8),dotss2[l]), y=c(1-chr[l],1-chr[l]),lwd=cex.trait, col=transp(col.trait[j],0.8))
          }
        }
        axis(3,at=sss + (fact/2), labels=trait, cex.axis=cex.axis) # cex of the scale of the trait
      }
      ## ---------------------------
      ## IF TRAIT IS A FACTOR TRAIT
      ## ---------------------------
      if(is.factor(prov[,trait])){
        riel <- sss + (fact/2)
        lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel+(fact/2), y=yy, labels=prov[ww2,trait], cex=0.3)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex) # cex of the scale of the trait (letters)
      }
      ## ---------------------------
      ## ---------------------------
      
      ## ---------------------------
      ## IF TRAIT IS A character TRAIT
      ## ---------------------------
      if(is.character(prov[,trait])){
        riel <- sss + (fact/1.8)
        #lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          #points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel, y=yy, labels=prov[ww2,trait], cex=0.17)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex.axis)
      }
      ## ---------------------------
      ## ---------------------------
    }
    ################
    
  }
}
maxi.qtl <- function(sma, distan=10, no.qtl=5, q95=2.5, LOD.int=FALSE, LODdrop=2, allCI=TRUE){
  param=NULL
  sma <- data.frame(sma)
  rnn <- rownames(sma)
  param <- vector("list", length(unique(sma$chr)))
  param <- lapply(param, function(x){x <- c(q95,no.qtl)})
  ## sma is the dataframe
  ## param is a list with 2 values, lod threshold and #of qtls
  sma <- apply(sma, 2, FUN=as.numeric)
  sma <- as.data.frame(sma)
  rownames(sma) <- rnn
  if(is.null(param)){
    param <- vector("list", length(unique(sma$chr)))
    param <- lapply(param, function(x){x <- c(1,1)})
  }
  
  big.peaks.col <- function(x, tre){
    r <- rle(x)
    v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
    pos <- v[which(x[v] > tre)] #position of the real peaks
    hei <- x[pos]
    out <- list(pos=pos,hei=hei)
    return(out)
  }
  
  result <- list(NA)
  for(j in 1:max(sma$chr,na.rm=TRUE)){
    prov <- sma[which(sma$chr == j),]
    para <- param[[j]]
    roxy <- big.peaks.col(prov$lod, tre=para[1]) # first element is the lod threshold, 2nd is the #of qtls
    if(length(roxy$pos) > 0){
      prov2.1 <- prov[roxy$pos,] # to keep
      prov2 <- prov[roxy$pos,] # to discard
      
      w <- numeric()
      res <- data.frame(matrix(NA,nrow=para[2], ncol=3)) # 3 columns because always rqtl gives that
      names(res) <- names(prov2)
      
      GOOD <- TRUE
      k=1
      while(GOOD){
        #for(k in 1:para[2]){# for the peaks required
        mm <- max(prov2$lod, na.rm=TRUE)
        w[k] <- which(prov2$lod == mm)[1]
        provi <- prov2$pos[w[k]]
        
        res[k,c(1:3)] <- prov2[w[k],c(1:3)]
        rownames(res)[k] <- rownames(prov2)[w[k]] #$$$$$$$$$$$$$$$$$$$$
        to.elim <- distan#abs(prov2$pos[w[k]] - 10)
        zz <- setdiff(which(prov2$pos > (provi - to.elim) & prov2$pos < (provi + to.elim) ), w[k])# [-w[k]]
        if(length(zz) > 0){
          prov2 <- prov2[-zz,]
        }else{prov2 <- prov2}
        selected <- which(prov2$lod == mm)[1]
        prov2 <- prov2[-selected,]
        #rownames(prov2.1)[-selected]
        k=k+1
        if(length(which(is.na(prov2))) > 0 | dim(prov2)[1] == 0 | k > no.qtl){GOOD<-FALSE}
        #}
      }
      
      #maxqtls <- sort(prov2$lod, decreasing = TRUE)
      result[[j]] <- res#prov2[which(prov2$lod %in% maxqtls[1:para[2]]),]
    }
    
  }# end for each chromosome
  #############
  result <- lapply(result, function(x){unique(x)})
  good <- which(unlist(lapply(result, is.null)) == FALSE)
  result <- result[good]
  
  good <- which(unlist(lapply(result, function(m){is.null(dim(m))})) == FALSE)
  result <- result[good]
  
  res2 <- do.call("rbind", result)
  
  res2<- res2[which(!is.na(res2$chr)),]
  
  res3 <- res2
  if(LOD.int){
    
    #print(res3)
    
    lod2 <-list() 
    for(i in 1:dim(res3)[1]){
      #apply(res3,1,function(x,sma){
      
      x <- res3[i,]
      #print(x)
      mpp <- sma[which(sma$chr == as.numeric(x[1])),]
      
      babo <- which(rownames(mpp) == rownames(x))#which(mpp$chr == as.numeric(x[1]) & mpp$pos == as.numeric(x[2]))
      toch <- mpp[babo,]
      baba <- babo-1
      babu <- babo+1
      rt=0
      #print(toch)
      while((rt < LODdrop) & (baba > 1)){ # 2-lod interval
        rt <- abs(as.numeric(mpp[baba,] - toch)[3])
        baba <- baba - 1
      }
      ## baba bow tell us what is the lower bound
      mpp[baba,]
      
      rt=0
      while((rt < LODdrop) & (babu < dim(mpp)[1])){
        rt <- abs(as.numeric(mpp[babu,] - toch)[3])
        babu <- babu + 1
      }
      mpp[babu,]
      if(allCI){
        resx <- rbind( mpp[baba:(babo-1),], toch,  mpp[(babo+1):babu,])
      }else{
        resx <- rbind( mpp[baba,], toch,  mpp[babu,])
      }
      
      
      lod2[[i]] <- resx
    }
    
    res2 <- do.call("rbind",lod2)
  }
  return(res2)
  
}
manhattan <- function(map, col=NULL, fdr.level=0.05, show.fdr=TRUE, PVCN=NULL, ylim=NULL,...){
  
  if(!is.null(PVCN)){
    colnames(map)[which(colnames(map)==PVCN)] <- "p.val"
  }
  
  required.names <- c("Chrom","Position","p.val")
  
  che <- which(names(map)%in%required.names)
  if(length(che) < 3){
    stop("Column names; 'Chrom','Position' and 'p.val' need 
         to be present in the data frame provided.",call. = FALSE)
  }
  
  map <- map[with(map, order(Chrom, Position)), ]
  
  
  yylim <- ceiling(max(map$p.val,na.rm=TRUE))
  if(is.null(col)){
    col.scheme <- rep((transp(c("cadetblue","red"))),30)
  }else{
    col.scheme <- rep(col,30)
  }
  ffr <- fdr(map$p.val, fdr.level=fdr.level)$fdr.10
  
  if(!is.null(ylim)){
    yyylim <- ylim
  }else{
    yyylim <- c(0,yylim)
  }
  
  plot(map$p.val, bty="n", 
       col=col.scheme[factor(map$Chrom, levels = unique(map$Chrom, na.rm=TRUE))], 
       xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), 
       las=2, 
       ylim=yyylim,...)
  init.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=map$Chrom)
  fin.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=map$Chrom)
  inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
  axis(side=1, at=inter.mrks, labels=paste("Chr",unique(map$Chrom),sep=""), cex.axis=.5)
  if(show.fdr){
    abline(h=ffr, col="slateblue4", lty=3, lwd=2)
    legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr,2), sep=""), 
           bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
  }
  #abline(v=marker, lty=3, col="red")
  #legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
  #marker <- as.character(map$Locus[marker])
}

spatPlots <- function(object, by=NULL, colfunc=NULL,row="ROW",range="RANGE", wire=FALSE, terms=NULL){
  
  nameslist <- lapply(as.list(names(object$U)),function(x){strsplit(x,":")[[1]]})
  if(is.null(terms)){
    tolook <- 1:length(object$U)
  }else{
    tolook <- which(unlist(lapply(nameslist,function(x){y <- which(terms%in%x); if(length(y) > 0){return(TRUE)}else{return(FALSE)}})))
    # tolook <- grep(terms,names(object$U))
  }
  pp <- predict(object, RtermsToForce = tolook)
  fits <- pp$fitted
  resp <- colnames(fits)[grep("predicted.value",colnames(fits))]
  if(is.null(by)){
    form <- as.formula(paste(resp,"~",row,"*",range))
  }else{
    form <- as.formula(paste(resp,"~",row,"*",range,"|",by))
  }
  
  if(is.null(colfunc)){
    colfunc <- colorRampPalette(c("gold","springgreen","steelblue4"))
  }
  
  maint <- paste(terms,collapse = ",")
  
  if(wire){ # wireframe
    print(wireframe(form, data=fits,  
                    aspect=c(61/87,0.4), drape=TRUE,
                    light.source=c(10,0,10), #=colfunc,
                    main=maint)
    )
  }else{ # levelplot
    print(levelplot(form, data=fits, main=maint, col.regions = colfunc))
  }
  return(fits)
}
transp <- function (col, alpha = 0.5) {
  res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                c[3]/255, alpha))
  return(res)
}


