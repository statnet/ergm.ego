degreedist.egodata <- function(egodata, freq = FALSE, prob = !freq, 
                               by = NULL, brgmod = FALSE){
  if (class(egodata) != "egodata"){
    stop("The egodata object passed to degreedist.egodata must be of class egodata.")
  }
  color <- "#83B6E1"
  beside <- FALSE
  maxdeg <- max(table(egodata$alters["vertex.names"]))
  deg.ego <- summary(egodata ~ degree(0:maxdeg, by = by))
  names(deg.ego) <- 0:maxdeg
  if(!is.null(by)){
    deg.ego <- matrix(deg.ego, ncol = maxdeg + 1, byrow = TRUE)
    colnames(deg.ego) <- 0:maxdeg
    ncolors <- dim(deg.ego)[1]
    if(ncolors == 2){
      color <- c("#eff3ff", "#377FBC")
    } else if(ncolors < 10){
      color <- RColorBrewer::brewer.pal(ncolors,"Blues")
    } else if(ncolors >= 10){
      color <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(ncolors)
    }
    
    ltext <- sort(unique(egodata$egos[by]))
    lfill <- c(color, 0)
    lborder <- c(rep("black", times = ncolors), 0)
    ltitle <- by
  }
  if(prob){
    if(is.null(by)){
      deg.ego <- deg.ego/sum(deg.ego)
    } else {
      deg.ego <- t(deg.ego)/colSums(deg.ego)
      deg.ego <- t(deg.ego)
      beside <- TRUE
    }
  }
  maxfreq <- ifelse(is.null(by), yes = max(deg.ego, na.rm = TRUE), 
                    no = max(colSums(deg.ego), na.rm = TRUE))


  if (brgmod) {
    brgdraws <- simulate.ergm.ego(ergm.ego(egodata ~ edges), nsim = 50)
    deg.brg <- summary(brgdraws ~ degree(0:maxdeg))
    brgmeans <- apply(deg.brg, MARGIN = 2, FUN = mean)
    brgsd <- apply(deg.brg, MARGIN = 2, FUN = sd)
    upper <- brgmeans + 2 * brgsd
    lower <- brgmeans - 2 * brgsd
    maxfreq <- max(maxfreq, upper, na.rm = TRUE)
  }
  
  baraxis <- barplot(deg.ego, xlab = "Degree", ylab = NULL,
                     col = color, beside = beside, plot = TRUE,
                     ylim = c(0, maxfreq))
  
  if(brgmod){
    points(x = baraxis-.15, y = brgmeans, col = "firebrick",
           lwd = 1, pch = 18, cex = 1.25)
    suppressWarnings(arrows(x0 = baraxis-.15, y0 = upper,
                            x1 = baraxis-.15, y1 = lower,
                            code = 3, length = 0.1, 
                            angle = 90, col = "firebrick"))
  } 
  if(!is.null(by)){
    legend(x="topright", legend = ltext, title = ltitle, fill = lfill, 
           border = lborder, bty = "n")
  }
}


mixingmatrix.egodata <- function(egodata, attrname, rowprob = FALSE){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties <- merge(egos[c(egoIDcol,attrname)], alters[c(egoIDcol,attrname)], 
                by = egoIDcol, suffixes = c(".ego",".alter"))
  mxmat <- table(ties[,paste0(attrname, c(".ego", ".alter"))])
  dimnames(mxmat) <- list(ego = levs,  
                          alter = levs)
  if(rowprob){
    mxmat <- mxmat/rowSums(mxmat)
    mxmat <- round(mxmat, digits = 3)
  }
  mxmat
}
  
  