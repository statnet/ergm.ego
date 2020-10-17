#  File R/degreedist.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################


#' Plotting the degree distribution of an egocentric dataset
#' 
#' A [degreedist()] method for [`egodata`] objects: plot a histogram
#' of the degree distribution of actors in the egocentric dataset,
#' optionally broken down by group and/or compared with a Bernoulli
#' graph.
#' 
#' @aliases degreedist
#' @param object A \code{\link{egodata}} object.
#' @param freq,prob Whether to plot the raw frequencies or the conditional
#' proportions of the degree values. Defaults to the latter.
#' @param by A character vector giving the name of a vertex attribute; if
#' given, plots the frequences broken down by that attribute.
#' @param brgmod Plot the range of predicted frequencies/probabilities
#' according to a Bernoulli graph having the same expected density as the
#' observed.
#' @param main Main title of the plot.
#' @param plot Whether to plot the histogram; if `FALSE`, graphical
#'   parameters and `bgrmod` have no effect.
#' @param weight Whether sampling weights should be incorporated into
#'   the calculation (`TRUE`, the default) or ignored (`FALSE`).
#' @param ... Additional arguments to [simulate.ergm.ego()].
#'
#' @return Returns either a vector of degree frequencies/proportions
#'   if `by=NULL` or a matrix with a row for each category if not. If
#'   \code{plot==TRUE} returns invisibly.
#' 
#' @seealso \code{\link{degreedist}},
#' \code{\link[ergm:summary_formula]{summary}}
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egodata(faux.mesa.high)
#' 
#' degreedist(fmh.ego,by="Grade",brgmod=TRUE)
#' # Compare:
#' degreedist(faux.mesa.high)
#' 
#' @importFrom graphics arrows barplot legend points
#' @export
degreedist.egodata <- function(object, freq = FALSE, prob = !freq, 
                               by = NULL, brgmod = FALSE, main = NULL, plot = TRUE, weight = TRUE, ...){
  egodata <- object
  if(!weight) egodata$egoWt[] <- 1

  color <- "#83B6E1"
  beside <- TRUE

  ylabel <- if(prob) "Proportion" else "Frequency"
  if(!is.null(by)) ylabel <- paste(ylabel, "(within attr level)")

  egoIDcol <- egodata$egoIDcol
  
  egodata$egos[[egoIDcol]] <- factor(egodata$egos[[egoIDcol]])
  egodata$alters[[egoIDcol]] <- factor(egodata$alters[[egoIDcol]], levels=levels(egodata$egos[[egoIDcol]])) 
  egodata$egos[[egoIDcol]] <- as.integer(egodata$egos[[egoIDcol]])
  egodata$alters[[egoIDcol]] <- as.integer(egodata$alters[[egoIDcol]])
  
  degtable <- rep(0, nrow(egodata$egos))
  degtable[as.numeric(names(table(egodata$alters[egoIDcol])))] <- table(egodata$alters[egoIDcol])
  
  if(is.null(by)){
    deg.ego <- xtabs(egodata$egoWt~degtable)
    names(dimnames(deg.ego)) <- "degree"
    degrees <- as.integer(names(deg.ego))
  }else{
    deg.ego <- xtabs(egodata$egoWt~egodata$egos[[by]]+degtable)
    names(dimnames(deg.ego)) <- c(by, "degree")
    levs <- rownames(deg.ego)
    degrees <- as.integer(colnames(deg.ego))
    ncolors <- dim(deg.ego)[1]
    if(ncolors == 2){
      color <- c("#eff3ff", "#377FBC")
    } else if(ncolors < 10){
      color <- RColorBrewer::brewer.pal(ncolors,"Blues")
    } else if(ncolors >= 10){
      color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(ncolors)
    }
    
    ltext <- levs
    lfill <- c(color, 0)
    ltitle <- by
    maxfreq <- max(colSums(deg.ego))
  }
  if(prob){
    if(is.null(by)){
      scaledeg <- sum(deg.ego)
      deg.ego <- deg.ego/scaledeg
      maxfreq <- max(deg.ego, na.rm = TRUE)
    } else {
      scaledeg <- rowSums(deg.ego)
      deg.ego <- deg.ego/scaledeg
      deg.ego <- deg.ego
      maxfreq <- max(max(deg.ego, na.rm = TRUE))
      beside <- TRUE
    }
  }

  if(plot){
    if(brgmod) {
      ppopsize.mul <- max(egodata$egoWt)/min(egodata$egoWt)
      brgdraws <- simulate(suppressMessages(ergm.ego(egodata ~ edges, control=control.ergm.ego(ppopsize=nrow(egodata$egos)*ppopsize.mul))), nsim = 50, ...)
      deg.brg <- summary(brgdraws ~ degree(degrees))/ppopsize.mul
      brgmeans <- apply(deg.brg, MARGIN = 2, FUN = mean)
      brgsd <- apply(deg.brg, MARGIN = 2, FUN = sd)
      upper <- brgmeans + 2 * brgsd
      lower <- brgmeans - 2 * brgsd
      
      if(prob){
        if(is.null(by)){
          brgmeans <- brgmeans/scaledeg
          upper <- upper/scaledeg
          lower <- lower/scaledeg
        } else {
          upper <- upper/sum(brgmeans)
          lower <- lower/sum(brgmeans)
          brgmeans <- brgmeans/sum(brgmeans)
        }
        
      }
      maxfreq <- max(maxfreq, upper, na.rm = TRUE)
    }
    
    baraxis <- barplot(deg.ego, xlab = "Degree", ylab = ylabel,
                       col = color, beside = beside, plot = TRUE,
                       ylim = c(0, maxfreq), main = main)
    
    if(brgmod){
      baraxis <- if(is.null(by)){
                   baraxis - 0.15
                 } else {
                   colMeans(baraxis)
                 }
      points(x = baraxis, y = brgmeans, col = "firebrick",
             lwd = 1, pch = 18, cex = 1.25)
      suppressWarnings(arrows(x0 = baraxis, y0 = upper,
                              x1 = baraxis, y1 = lower,
                              code = 3, length = 0.1, 
                              angle = 90, col = "firebrick"))
    } 
    if(!is.null(by)){
      legend(x="top", legend = ltext, title = ltitle, fill = lfill, bg="white")
    }
  }
  
  if(plot) invisible(deg.ego) else deg.ego
}




#' Summarizing the mixing among groups in an egocentric dataset
#' 
#' A \code{\link[network]{mixingmatrix}} method for
#' \code{\link{egodata}} objects, to return counts of how often a ego
#' of each group nominates an alter of each group.
#' 
#' 
#' @aliases mixingmatrix
#' @param object A \code{\link{egodata}} object.
#' @param attrname A character vector containing the name of the network
#' attribute whose mixing matrix is wanted.
#' @param rowprob Whether the counts should be normalized by row sums. That is,
#' whether they should be proportions conditional on the ego's group.
#' @param weight Whether sampling weights should be incorporated into
#'   the calculation (`TRUE`, the default) or ignored (`FALSE`).
#' @param ... Additional arguments, currently unused.
#' @return A matrix with a row and a column for each level of \code{attrname}.
#' 
#' Note that, unlike \code{\link[network]{mixingmatrix}}, what is counted are
#' \emph{nominations}, not ties. This means that under an egocentric census,
#' the diagonal of \code{mixingmatrix.egodata} will be twice that returned by
#' \code{\link[network]{mixingmatrix}} for the original undirected network.
#' @seealso \code{\link[network]{mixingmatrix}}, \code{\link[ergm]{nodemix}},
#' \code{\link[ergm.ego]{summary}} method for egocentric data
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egodata(faux.mesa.high)
#' 
#' (mm <- mixingmatrix(faux.mesa.high,"Grade"))
#' (mm.ego <- mixingmatrix(fmh.ego,"Grade"))
#' 
#' stopifnot(isTRUE(all.equal({tmp<-unclass(mm$matrix); diag(tmp) <- diag(tmp)*2;
#' tmp}, mm.ego, check.attributes=FALSE)))
#' 
#' @export
mixingmatrix.egodata <- function(object, attrname, rowprob = FALSE, weight = TRUE, ...){
  egodata <- object
  if(!weight) egodata$egoWt[] <- 1

  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]], alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties <- merge(egos[c(egoIDcol,attrname)], alters[c(egoIDcol,attrname)], 
                by = egoIDcol, suffixes = c(".ego",".alter"))
  ties$wt <- object$egoWt[match(ties[[egoIDcol]],egos[[egoIDcol]])]
  mxmat <- matrix(0, nrow = length(levs), ncol = length(levs))
  
  for(i in 1:length(levs)){
    for(j in 1:length(levs)){
      mxmat[i,j] <- sum(ties$wt[ties[,2]==i & ties[,3]==j])
    }
  }
  dimnames(mxmat) <- list(ego = levs,  
                          alter = levs)
  if(rowprob){
    mxmat <- mxmat/rowSums(mxmat)
  }
  mxmat
}
  
  
