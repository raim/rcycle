## PLOT FUNCTIONS



#' Add arrows for segments in the phases object.
#' @export
arrows.phases <- function(x, types='shoulder', phase='phi',
                          y0, dy, col, labels, labels.top, lxpd=par('xpd'),
                          pos=NULL,
                          ticks=FALSE, verb=0, ...) {

    if ( missing(y0) ) y0 <- mean(par('usr')[3:4])
    if ( missing(dy) ) dy <- diff(par('usr')[3:4])/10 #length(types)

    miss <- which(!types%in%names(x) & !types%in%names(x$x.phases))
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

            ## TODO: fuse adjacent equal class labels!
            if ( !missing(labels) ) {

                labs <- segs[[labels]]

                ## AMPLITUDE CUTOFF
                ## only show labels for top X amplitudes 
                if ( !missing(labels.top) ) 
                    labs[rank(-segs$amp) > labels.top] <- ''
                
            }
            
            ## arrow colors
            ids <- segs$ID
            scol <- setNames(1:nrow(segs), ids)
            if ( !missing(col) )
                if ( col%in%colnames(segs) )
                    scol <- setNames(segs[[col]], ids)

            ## TODO: shadow.arrows
            arrows(x0=starts, x1=ends, y0=y0, code=3, length=.05, col=scol, ...)
            ## add names
            if ( !missing(labels) ) {
                if ( lxpd ) {
                    oxpd <- par('xpd')
                    par(xpd=TRUE)
                }
                shadowtext(x=(starts+ends)/2, y=rep(y0, length(starts)),
                           labels=labs, col=scol, pos=pos, xpd=lxpd)
                if ( lxpd ) par(xpd=oxpd)
            }

        } else if ( type %in% colnames(x$x.phases) ) {

            ## labels
            ids <- x$x.phases$ID

            ## AMPLITUDE CUTOFF
            ## only show labels for top X amplitudes
            labs <- ids
            if ( !missing(labels.top) ) 
                labs[rank(-x$amp) <= labels.top] <- ''

            ## arrow colors
            scol <- setNames(1:nrow(x$x.phases), ids)
            if ( !missing(col) )
                if ( col%in%colnames(x$x.phases) )
                    scol <- setNames(x$x.phases[,col], ids)

            if ( lxpd ) {
                oxpd <- par('xpd')
                par(xpd=TRUE)
            }
            shadowtext(x=x$x.phases[,type], y=rep(y0, nrow(x$x.phases)),
                       labels=labs, col=scol, ...)
            
            if ( lxpd ) par(xpd=oxpd)
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
stop('this needs to be udpated')
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
                       lwd=2, axes=TRUE, ylab, ylim, ...) {

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

    if ( missing(ylim) ) {
        ylim <- c(states[cls.srt,])
        ylim <- range(ylim[is.finite(ylim)],na.rm=TRUE)
    } 
    
    if ( missing(ylab) )
        ylab <- ifelse(norm, expression(log[2](sample/mean)),
                       expression(mean~counts))
    
    xlim <- range(phase)
    
    plot(1, xlim=xlim, ylim=ylim, col=NA,
         main=NA, axes=FALSE, xlab=xlab, ylab=ylab, ...)
    if ( axes ) {
        if ( xtype%in%c("angle","phase") ) 
            circ.axis(1)
        else axis(1) 
        axis(2)
    }
    if ( lines ) 
        for ( k in seq_along(cls.srt) ) 
            lines(phase[ord], states[cls.srt[k], ord],
                  col=paste0(cls.col[cls.srt[k]],77))
    if ( ma )
        for ( k in seq_along(cls.srt) ) 
            lines(phase[ord], ma(states[cls.srt[k], ord],
                                 n=Ncells, circular=TRUE),
                  col=paste0(cls.col[cls.srt[k]]), lwd=lwd)

    figlabel(paste0(sid), pos=ifelse(legend,"topleft","top"), cex=1.2, font=2)
    if ( legend )
        legend("topright", sub(".*_", "", cls.srt),
               col=cls.col[cls.srt],
               pch=15, pt.cex=1, cex=.8, inset=c(-0.1,-0.1),
               bg="#ffffff00", ncol=ceiling(length(cls.srt)/leg.nrow),
               x.intersp=.6,
               seg.len=0, xpd=TRUE, box.col=NA)

}

## plots both rotation and x of a PCA object, used from plotPC
## TODO: revise color selection via named vectors
monoplot <- function(x, type='rotation',
                     xs, ys, lines=FALSE, arrows=FALSE,
                     labels=FALSE, labels.top,
                     axis=FALSE,
                     colf = NULL, col = NULL, lwd=1, pch=1, cex=1, txt.cex=1,
                     xlim, ylim, ax=c(1,2), xlab, ylab, ...) {

    ## get matrix from PCA object to plot
    xy <- x[[type]]

    ## PCA phase object from rcycle
    phases <- NULL
    typep <- paste0(type,".phase")
    if ( typep %in% names(x) )
        phases <- x[[typep]]
    
    ## expand single color
    if ( length(col)==1 )
        col <- setNames(rep(col, nrow(xy)), rownames(xy))
    

    ## take color from PCA decoration
    if ( is.null(col) & !is.null(phases) )  
        if ( 'col' %in% colnames(phases) )
            col <- phases$col

    ## TODO: avoid this?
    if ( !is.null(col) ) 
        if ( !is.null(names(col)) )
            if ( all(names(col)%in%rownames(xy)) )
                col <- col[rownames(xy)]

    ## colored points or density plot?
    if ( is.null(col) )  { # density plot!
        if ( is.null(colf) )
            colf <- function(n) grey.colors(n, start=0, end=1)
        dense2d(xy[,xs],
                xy[,ys],
                xlim=xlim, ylim=ylim,
                colf = colf, pch=pch, cex=cex,
                xlab=NA, ylab=NA, axes=FALSE, ...)
        ## NOTE: black arrows with density function
        ## TODO: provide explicit density argument, not just default!
        col <- setNames(rep(1, nrow(xy)), rownames(xy))
    } else 
        plot(xy[,xs], xy[,ys], # scatter plot
             xlim=xlim, ylim=ylim,
             col=col, pch=pch, cex=cex, 
             xlab=NA, ylab=NA, axes=FALSE, ...)
    
    if ( labels ) {

        ## TODO: solve this cleaner and independent of _ convention
        labs <- sub(".*_","",rownames(xy))

        ## AMPLITUDE CUTOFF
        ## only show labels for top X amplitudes 
        if ( !missing(labels.top) & !is.null(phases) ) 
            labs[rank(-phases$amp) > labels.top] <- ''
                
        shadowtext(xy[,xs], xy[,ys],
                   labels=labs,
                   col=col,
                   cex=txt.cex, font=2, xpd=TRUE, r=.1)
    }
    if ( arrows ) {

        ## background color: black or white? 
        bgn <- rep(NA, length(col))
        bg.thresh <- 0.75
        for (i in 1:length(col)) {
            crgb <- col2rgb(col[i])/255
            L <- 0.2126 * crgb[1, 1] + 0.7152 * crgb[2, 1] + 
                0.0722 * crgb[3, 1]
            bgn[i] <- ifelse(L > bg.thresh, "#000000", "#FFFFFF")
        }
        arrows(x0=0,y0=0, x1=xy[,xs], y1=xy[,ys],
               col=bgn, lwd=lwd+2, length=.05)
        arrows(x0=0,y0=0, x1=xy[,xs], y1=xy[,ys],
               col=col, lwd=lwd, length=.05)
    }
    if ( lines )
        lines(xy[,xs], xy[,ys], col=col[1], type='b', pch=NA, lwd=lwd)
    
    
    if ( axis ) {
        axis(ax[1])
        axis(ax[2])
        mtext(xlab, ax[1], par('mgp')[1])
        mtext(ylab, ax[2], par('mgp')[1])
    }
}
    

#' Plot PCA-based circle and state vectors.
#' 
#' NOTE: scale = 1, pc.biplot = FALSE is equivalent to the defaults
#' of R's biplot.prcomp.
#' 
#' @param z use the zth PC component for coloring
#' between quantiles given in z.q
#' @param scale scale parameter, if provided it converts the plot
#' to a true biplot 
#' @export
plotPC <- function(phases, x=1, y=2,
                   z, z.q = c(.05,.95), z.legend = FALSE, # color by PCz

                   vectors = TRUE, # eigenvector colors, cells
                   vlines=FALSE, varrows = FALSE,
                   vlabels = FALSE, vlabels.top=Inf,
                   colf = NULL, col = NULL, lwd = 1,
                   pch=19, cex=.5, vaxis = vectors,
                   
                   scores = TRUE, # scores: rotated data/PCs
                   slines = FALSE, sarrows = FALSE,
                   slabels = FALSE, slabels.top, 
                   scolf = NULL, scol = NULL, slwd=1,
                   spch=1, scex=1, saxis = scores, 

                   txt.cex = 1,
                   
                   scale = 1, # true biplot scaling as in biplot(scale=1)
                   pc.biplot = FALSE, # as in biplot
                   xlim, ylim, expand = 1, 

                   arcsinh = FALSE, # useful to emphasize crowded data, TODO: avoid
                   zero.axis = FALSE, zero.axis.label = zero.axis,
                   pc.arrows = FALSE, pc.lwd = 1,
                   show.var=TRUE,
                   ...) {


    ## PC names to interface data and for plot labels
    xs <- paste0('PC', x) # Rotated data: cohorts
    ys <- paste0('PC', y)
    xv <- paste0('EV', x) # Eigenvectors: cells
    yv <- paste0('EV', y)

    alls <- grep("^PC",colnames(phases$x))
    allv <- grep("^PC",colnames(phases$rotation))

    ## proportion of variance
    if ( !'summary'%in%names(phases) )
        phases$summary <- summary(phases)$importance
    varp <- round(phases$summary['Proportion of Variance',]*100,1)

    phases$var_prct <- varp

    ## TRUE biplot:
    ## scale data by eigenvalues as in biplot: matrices used for biplot
    ## ... should, when multiplied together, approximate X.
    ## https://stats.stackexchange.com/questions/66926/what-are-the-four-axes-on-pca-biplot
    ## https://stats.stackexchange.com/questions/141085/positioning-the-arrows-on-a-pca-biplot
    ## see biplot.prcomp code: src/library/stats/R/biplot.R

    ## test if this object was already scaled and refuse plotting
    ## until solved: put scaling into separate function and record!
    is.scaled <- FALSE # were matrices already scaled by previous call?
    if ( "SCALED" %in% names(phases) )
        is.scaled <- phases$SCALED 
    if ( is.scaled ) {
        ## until this is cleanly solved
        stop("the object returned by plotPC can not be used for replotting")
    }
    
    lam <- phases$sdev
    n <- NROW(phases$x)
    lam <- lam * sqrt(n) # singular values! TODO: n+1 ?
    if ( scale < 0 || scale > 1 ) warning("'scale' is outside [0, 1]")
    if ( scale != 0 ) lam <- lam^scale else lam <- 1
    if ( pc.biplot ) lam <- lam / sqrt(n) # back to sdev
    
    phases$rotation[,allv] <- t(t(phases$rotation[,allv]) * lam)
    phases$x[,alls] <- t(t(phases$x[,alls])/lam)

    ## transform data for better circular visibility?
    if ( arcsinh ) {
        ash <- function (x) 
            log(x + sqrt(x^2 + 1))
        phases$x[,alls] <-  apply(phases$x[,alls], 2, ash)
        phases$rotation[,allv] <- apply(phases$rotation[,allv], 2, ash)
    }  
    


    ## PC component coloring
    if ( !missing(z) ) {
        zs <- paste0('PC', z) # Rotated data
        if ( missing(col) )
            col <- num2col(phases$rotation[,zs], q=z.q)
    } else z.legend <- FALSE
    

    ## axis labels
    xlab <- xs
    ylab <- ys
    xvlab <- xv
    yvlab <- yv

    if ( arcsinh ) {
        xlab <- paste0("asinh(",xlab,")") # rotated data 
        ylab <- paste0("asinh(",ylab, ")")
        xvlab <- paste0("asinh(",xvlab,")") # eigenvector 
        yvlab <- paste0("asinh(",yvlab, ")")
    }

    if ( show.var ) {
        xlab <- paste0(xlab, " (", varp[xs], "%)") # rotated data 
        ylab <- paste0(ylab, " (", varp[ys], "%)")
        xvlab <- paste0(xvlab, " (", varp[xs], "%)") # eigenvector 
        yvlab <- paste0(yvlab, " (", varp[ys], "%)")
    }


    ## align xlim and ylim such that origins overlap
    ## copied from biplot.default

    ## TODO: fix for eigen|vector cases
    
    unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(phases$x[, xs])
    rangx2 <- unsigned.range(phases$x[, ys])
    rangy1 <- unsigned.range(phases$rotation[, xs])
    rangy2 <- unsigned.range(phases$rotation[, ys])

    ## NOTE: both PC dimensions are on the same scale!
    if(missing(xlim) && missing(ylim))
	xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand

    xlimv <- rangy1
    ylimv <- rangy2
    if ( scores ) {
        xlimv <- xlim*ratio
        ylimv <- ylim*ratio
    }

    ## TODO: allow to flip order of vectors and scores plots!
    
    if ( vectors ) { # plot eigenvectors
        ax <- c(3,4)
        if ( !scores ) ax <- c(1,2)
        monoplot(x=phases, type='rotation', #xy=phases$rotation,
                 xs=xs, ys=ys, lines=vlines, arrows=varrows,
                 labels=vlabels, labels.top=vlabels.top,
                 colf=colf, col=col, lwd=lwd, pch=pch, cex=cex,
                 txt.cex=txt.cex,
                 xlim=xlimv, ylim=ylimv, ax=ax, xlab=xvlab, ylab=yvlab,
                 axis=vaxis, ...)
    } else # empty plot
        plot(0, col=NA, axes=FALSE, xlab='', ylab='', xlim=xlimv, ylim=ylimv)
    
    if ( zero.axis ) {
        if ( zero.axis.label ) {
            axis(1, at=0, label=xlab)
            axis(2, at=0, label=ylab)
        }
        abline(h=0)
        abline(v=0)
    }
    ## draw arrows where length reflects % var
    ## TODO: draw higher vp to a set fraction of the axis
    if ( pc.arrows ) { 
        
        vp <- phases$summary['Proportion of Variance',]
        varx <- diff(par('usr')[1:2])*vp[x]
        vary <- diff(par('usr')[3:4])*vp[y]

        arrows(x0=par('usr')[1], y0=par('usr')[3],
               x1=par('usr')[1] + varx, lwd=pc.lwd, length=.05, xpd=TRUE)
        arrows(x0=par('usr')[1], y0=par('usr')[3],
               y1=par('usr')[3] + vary, lwd=pc.lwd, length=.05, xpd=TRUE)
        if ( FALSE ) { # TODO: plot at end of arrows
            shadowtext(0, par('usr')[3],
                       labels=paste(varp[1], '%'),
                       xpd=TRUE, pos=4, col=1, font=2)
            shadowtext(par('usr')[1], 0, 
                       labels=paste(varp[2], '%'),
                       xpd=TRUE, pos=4, srt=90, col=1, font=2)
        }
    }
    if ( scores ) { # plot scores/PCs
        if ( vectors )  par(new=TRUE)

        monoplot(x=phases, type='x', #xy=phases$x,
                 xs=xs, ys=ys, lines=slines, arrows=sarrows,
                 labels=slabels, labels.top=slabels.top,
                 colf=scolf, col=scol, lwd=slwd, pch=spch, cex=scex,
                 txt.cex=txt.cex,
                 xlim=xlim, ylim=ylim, ax=c(1,2), xlab=xlab, ylab=ylab,
                 axis=saxis, ...)
    }


    
    ## legend for PC-based coloring
    if ( z.legend ) 
        phcol.legend(leg.pos="bottomright",
                     legend=c("min","mid","max"), title=zs,
                     ##inset=c(-0.1,0),
                     xpd=TRUE, y.intersp=0.6, bg="#ffffffaa", box.col=NA)

    ## RETURN POTENTIALLY TRANSFORMED DATA!
    ## note: this should only be used for adding data to the canvas!
    ## indicate that this is scaled
    ## TODO: add scale and plot settings, e.g. axis limits!
    phases$SCALED <- TRUE 
    invisible(phases)
    
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
