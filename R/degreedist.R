#  File R/degreedist.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################


#' Plotting the degree distribution of an egocentric dataset
#' 
#' A [degreedist()] method for [`egodata`] objects: plot a histogram
#' of the degree distribution of actors in the egocentric dataset,
#' optionally broken down by group and/or compared with a Bernoulli
#' graph.
#' 
#' @aliases degreedist
#' @param object A [`egor`] object.
#' @param freq,prob Whether to plot the raw frequencies or the conditional
#' proportions of the degree values. Defaults to the latter.
#' @param by A character vector giving the name of a vertex attribute; if
#' given, plots the frequences broken down by that attribute.
#' @param brgmod Plot the range of predicted frequencies/probabilities
#' according to a Bernoulli graph having the same expected density as the
#' observed.
#' @param main Main title of the plot.
#' @param plot Whether to plot the histogram; defaults to the same
#'   value as`brgmod`, i.e., `FALSE`.
#' @param weight Whether sampling weights should be incorporated into
#'   the calculation (`TRUE`, the default) or ignored (`FALSE`).
#' @param ... Additional arguments to [simulate.ergm.ego()].
#'
#' @return Returns either a vector of degree frequencies/proportions
#'   if `by=NULL` or a matrix with a row for each category if not. If
#'   \code{plot==TRUE} returns invisibly.
#' 
#' @seealso [degreedist()], \code{\link[ergm:summary_formula]{summary}} for formulas.
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' 
#' degreedist(fmh.ego,by="Grade",brgmod=TRUE)
#' # Compare:
#' degreedist(faux.mesa.high)
#' 
#' @importFrom graphics arrows barplot legend points
#' @importFrom methods is
#' @export
degreedist.egor <- function(object, freq = FALSE, prob = !freq, 
                            by = NULL, brgmod = FALSE, main = NULL, plot = brgmod, weight = TRUE, ...){
  color <- "#83B6E1"
  beside <- TRUE

  ylabel <- if(prob) "Proportion" else "Frequency"
  if(!is.null(by)) ylabel <- paste(ylabel, "(within attr level)")

  degtable <- .degreeseq(object)

  w <- if(weight) weights(object) else rep(1, nrow(object$ego))
  
  if(is.null(by)){
    deg.ego <- xtabs(w~degtable)
    names(dimnames(deg.ego)) <- "degree"
    degrees <- as.integer(names(deg.ego))
  }else{
    deg.ego <- xtabs(w~as_tibble(object$ego)[[by]]+degtable)
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
    } else {
      scaledeg <- rowSums(deg.ego)
      beside <- TRUE
    }
    deg.ego <- deg.ego/scaledeg
  }
  maxfreq <- max(deg.ego, na.rm=TRUE)
  
  if(plot){
    if(brgmod) {
      ppopsize.mul <- max(w)/min(w)
      brgdraws <- simulate(suppressMessages(ergm.ego(object ~ edges, control=control.ergm.ego(ppopsize=nrow(object$ego)*ppopsize.mul))), nsim = 50, ...)
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
#' [`egor`] objects, to return counts of how often a ego
#' of each group nominates an alter of each group.
#' 
#' 
#' @aliases mixingmatrix
#' @param object A [`egor`] object.
#' @param attrname A character vector containing the name of the network
#' attribute whose mixing matrix is wanted.
#' @param rowprob Whether the counts should be normalized by row sums. That is,
#' whether they should be proportions conditional on the ego's group.
#' @param weight Whether sampling weights should be incorporated into
#'   the calculation (`TRUE`, the default) or ignored (`FALSE`).
#' @param ... Additional arguments, currently unused.
#' @return A matrix with a row and a column for each level of \code{attrname}.
#' 
#' Note that, unlike [network::mixingmatrix()], what is counted are
#' \emph{nominations}, not ties. This means that under an egocentric census,
#' the diagonal of \code{mixingmatrix.egor} will be twice that returned by
#' [network::mixingmatrix()] for the original undirected network.
#' @seealso [network::mixingmatrix()], \code{\link[ergm:nodemix-ergmTerm]{nodemix}} ERGM term,
#' \code{\link[ergm.ego]{summary}} method for egocentric data
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' 
#' (mm <- mixingmatrix(faux.mesa.high,"Grade"))
#' (mm.ego <- mixingmatrix(fmh.ego,"Grade"))
#' 
#' @export
mixingmatrix.egor <- function(object, attrname, rowprob = FALSE, weight = TRUE, ...){
  ds <- .degreeseq(object)
  egos <- rep(.unfactor(as_tibble(object$ego)[[attrname]]), ds)
  alters <- .unfactor(object$alter[[attrname]])
  levs <- sort(unique(c(egos,alters)))

  w <- if(weight) rep(weights(object),ds) else rep(1, nrow(object$alter))

  mxmat <- outer(levs, levs, Vectorize(function(l1, l2) sum(w[egos==l1&alters==l2])))

  dimnames(mxmat) <- list(ego = levs,  
                          alter = levs)
  if(rowprob){
    mxmat <- mxmat/rowSums(mxmat)
  }
  structure(mxmat, class = c("mixingmatrix", "table"), directed = FALSE,
            bipartite = FALSE)
}
