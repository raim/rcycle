## PLOT FUNCTIONS



#' Add arrows for segments in the phases object.
#' @export
arrows.phases <- function(x, types='shoulder', phase='phi',
                          y0, dy, col, labels, ticks=FALSE, verb=0, ...) {

    if ( missing(y0) ) y0 <- mean(par('usr')[3:4])
    if ( missing(dy) ) dy <- diff(par('usr')[3:4])/10 #length(types)

    miss <- which(!types%in%names(x) & !types%in%names(x$x))
    if ( length(miss) ) {
        warning('omitting types not found in data: ',
                paste(types[miss], sep=';'))
        types <- types[-miss]
    }

    for ( i in seq_along(types) ) {

        type <- types[i]

        if ( type %in% names(x) ) {
    
            segs <- x[[type]]
        
            n <- nrow(segs)
            starts <- ends <- segs[,phase]
            starts <- c(starts[n]-2*pi, starts)
            ends <- c(ends, ends[1]+2*pi)

            ## arrow colors
            ids <- segs$ID
            scol <- setNames(1:nrow(segs), ids)
            if ( !missing(col) )
                if ( col%in%colnames(segs) )
                    scol <- setNames(segs[[col]], ids)

            ## TODO: shadow.arrows
            arrows(x0=starts, x1=ends, y0=y0, code=3, length=.05, col=scol, ...)
            ## add names
            if ( !missing(labels) )
                shadowtext(x=(starts+ends)/2, y=rep(y0, length(starts)),
                           labels=segs[[labels]],
                           col=scol)

        } else if ( type %in% colnames(x$x) ) {

            ## arrow colors
            ids <- x$x$ID
            scol <- setNames(1:nrow(x$x), ids)
            if ( !missing(col) )
                if ( col%in%colnames(x$x) )
                    scol <- setNames(x$x[,col], ids)

            shadowtext(x=x$x[,type], y=rep(y0, nrow(x$x)),
                       labels=ids, col=scol, ...)
        } 
        if ( ticks ) axis(4, at=y0, labels=type, las=2)
        
        y0 <- y0 + dy
    
    }
}


## plot state time series along phase
## TODO: why is the moving average in some cases, Urea, below the curve?

#' Plot genes of interest time series
#' @export
plotGOI <- function(phases, counts, goi, names, col,
                    win=.05,
                    xlab=expression(pseudophase*phi),
                    ylab="mov. avg. of counts/total") {

    ## phase
    phs <- phases[,"phase"]
    pord <- order(phs)

    if ( missing(col) )
        col <- 1:length(goi)

    if ( missing(names) )
        names <- rownames(counts)[goi]
    
    ## moving average
    Ncells <- ncol(counts)*win 

    pho <- apply(counts[goi, pord], 1, ma,  n=Ncells, circular=TRUE)

    
    matplot(x=phs[pord], pho, col=col, type="l", lty=1, axes=FALSE,
            xlab=xlab, ylab=ylab)
    
    legend("topright", names, 
           col=col, lty=1, lwd=2, bg="#ffffff99",
           box.col=NA, inset=c(-.1,-.1), xpd=TRUE, y.intersp=.7,
           seg.len=.7, cex=.8)
    
    circ.axis(1)
    axis(2)

}

## TODO:
## * fuse plotStates and plotGOI, and align with segmenTools::plotClusters
##   and clusterAverage
## * polygon or 'each' style,
## * polygon: moving average +/- quantiles, sd, etc.

#' Plot cohort state time series.
#' @export
plotStates <- function(phase, states, cls.srt, cls.col,
                       win=.01, norm=FALSE, lines=TRUE, ma=TRUE,
                       sid="", legend=TRUE, leg.nrow=1,
                       xtype='phase', xlab=expression(phase~phi),
                       ylab, ...) {

    ## TODO: allow phases object

    ## order by x-value
    ord <- order(phase)

    if ( norm )    
        states <- log2(states/apply(states, 1, mean))

    
    if ( missing(cls.srt) )
        cls.srt <- rownames(states)
    if ( missing(cls.col) )
        cls.col <- setNames(1:length(cls.srt), nm=cls.srt)
    
    Ncells <- ncol(states)*win 

    tlim <- c(states[cls.srt,])
    tlim <- range(tlim[is.finite(tlim)],na.rm=TRUE)

    if ( missing(ylab) )
        ylab <- ifelse(norm, expression(log[2](sample/mean)),
                       expression(mean~counts))

    xlim <- range(phase)
    
    plot(1, xlim=xlim, ylim=tlim, col=NA,
         main=NA, axes=FALSE, xlab=xlab, ylab=ylab, ...)
    if ( xtype%in%c("angle","phase") ) 
        circ.axis(1)
    else axis(1) 
    axis(2)
    if ( lines ) 
        for ( k in seq_along(cls.srt) ) 
            lines(phase[ord], states[cls.srt[k], ord],
                  col=paste0(cls.col[cls.srt[k]],77))
    if ( ma )
        for ( k in seq_along(cls.srt) ) 
            lines(phase[ord], ma(states[cls.srt[k], ord],
                                 n=Ncells, circular=TRUE),
                  col=paste0(cls.col[cls.srt[k]]), lwd=2)

    figlabel(paste0(sid), pos=ifelse(legend,"topleft","top"), cex=1.2, font=2)
    if ( legend )
        legend("topright", sub(".*_", "", cls.srt),
               col=cls.col[cls.srt],
               pch=15, pt.cex=1, cex=.8, inset=c(-0.1,-0.1),
               bg="#ffffff00", ncol=ceiling(length(cls.srt)/leg.nrow),
               x.intersp=.6,
               seg.len=0, xpd=TRUE, box.col=NA)

}

## simple plot of two eigenvectors from PCA
## TODO: * this is biplot-like but shows arrows for
## for pca$x instead of pca$rotation; find proper description for this?

#' Plot PCA-based circle and state vectors.
#' @export
plotPC <- function(phases, x=1, y=2, col, 
                   expand=TRUE,
                   eigen.axis=FALSE,
                   cohorts=TRUE, ccol, txt.cex=1, ...) {

    if ( !inherits(phases, "phases") )
        warning("phases must be an object of class 'phases', ",
                "as returned by get_pseudophase")

    xs <- paste0('PC', x)
    ys <- paste0('PC', y)
    xv <- paste0('EV', x)
    yv <- paste0('EV', y)

    ##xlab <- xs
    ##ylab <- ys
    
    ## proportion of variance
    if ( !'summary'%in%names(phases) )
        phases$summary <- summary(phases)$importance
    varp <- round(phases$summary['Proportion of Variance',]*100,1)

    xlab <- paste0(xs, " (", varp[xs], "%)")
    ylab <- paste0(ys, " (", varp[ys], "%)")
    xvlab <- paste0(xv, " (", varp[xs], "%)")
    yvlab <- paste0(yv, " (", varp[ys], "%)")


        
    xlim <- range(phases$rotation[,xs])
    ylim <- range(phases$rotation[,ys])

    if ( expand ) {
        mx <- max(abs(c(ylim, xlim)))
        xlim <- ylim <- c(-mx, mx)
    }
    
    if ( missing(col) ) 
        dense2d(phases$rotation[,xs],
                phases$rotation[,ys],
                xlim=xlim, ylim=ylim,
                xlab=NA, ylab=NA, axes=FALSE, ...)
    else 
        plot(phases$rotation[,xs], phases$rotation[,ys],
             xlim=xlim, ylim=ylim,
             xlab=NA, ylab=NA, col=col, axes=FALSE, ...)
    if ( eigen.axis ) {
        axis(3);axis(4)
        mtext(xvlab, 3, par('mgp')[1])
        mtext(yvlab, 4, par('mgp')[1])
    }
    
    if ( cohorts ) {

        cohorts <- phases$x

        if ( missing(ccol) ) {
            ccol <- rep('#000000', nrow(cohorts))
            names(ccol) <- rownames(cohorts)
        }
            
        xlim <- range(cohorts[,xs])
        ylim <- range(cohorts[,ys])
        
        if ( expand ) {
            mx <- max(abs(c(ylim, xlim)))
            xlim <- ylim <- c(-mx, mx)
        } else stop("expand=TRUE required for aligned cohort arrows")

        par(new=TRUE)
        plot(cohorts[,xs], cohorts[,ys],
             xlim=xlim, ylim=ylim,
             axes=FALSE, col=NA, xlab=xlab, ylab=ylab)
        arrows(x0=0,y0=0, x1=cohorts[,xs], y1=cohorts[,ys],
               col="white", lwd=4, length=.05)
        arrows(x0=0,y0=0, x1=cohorts[,xs], y1=cohorts[,ys],
               col=ccol[rownames(cohorts)], lwd=2, length=.05)
        shadowtext(cohorts[,xs], cohorts[,ys],
                   labels=sub(".*_","",rownames(cohorts)),
                   col=ccol[rownames(cohorts)],
                   cex=txt.cex, font=2, xpd=TRUE, r=.1)
        axis(1);axis(2)
    }
    
}

## overlaid phase histograms of cell classes 
#' Plot a phase histogram by cell classes.
#' @export
phaseHist <- function(phase, cls, cls.srt, cls.col, sid="", leg.nrow=1) {

    if ( missing(cls.srt) ) cls.srt <- unique(cls)

    ## TODO: use -pi,pi for centered data
    brks <- seq(min(phase)-pi/18, max(phase)+pi/18, 2*pi/36)
    
    hist(phs, breaks=brks, col = NA, border="#00000000", axes=FALSE,
         main=paste0(sid), xlab="pseudophase")
    for ( pid in unique(cls) ) {
        clcol <- "#CCCCCC"
        if ( pid %in% names(cls.col) )
            clcol <- cls.col[pid]
        hist(phase[cls==pid], breaks=brks, add=TRUE,
             border=NA, col=paste0(clcol,"77"))
    }
    circ.axis(1)
    axis(2)
    legend("topright", sub(".*_", "", cls.srt),
           col=cls.col[cls.srt],
           pch=15, pt.cex=2, cex=.8,
           bg="#ffffff77", ncol=ceiling(length(cls.srt)/leg.nrow), x.intersp=.6,
           seg.len=0, xpd=TRUE, box.col=NA)
}

## small legend for phase angle coloring
phcol.legend <- function(leg.pos="topright",
                         legend= c(expression(-pi),
                                   expression(0),
                                   expression(pi)),
                         title=expression(phase~phi), ...) {
    pcol <- num2col(c(-pi,0,pi), limits=c(-pi,pi))
    legend(leg.pos, legend=legend,
           title=title, col=pcol, pch=19, ...)
}

#' generate circular plot axis in radian
#' @export
circ.axis <- function(x,
                      at=c(-2*pi,  -pi,
                           -pi/2,    0,
                           pi/2,    pi,
                           1.5*pi,2*pi,
                           2.5*pi,3*pi),
                      labels=expression(-2*pi,-pi,
                                        -pi/2,  0,
                                        pi/2,  pi,
                                        '',  2*pi,
                                        '',  3*pi),...)
    for ( ax in x )
        axis(ax, at=at,
             labels=labels,
             ...)
