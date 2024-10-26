
## CALCULATE PHASE SEGMENTS


#' Add a phase segmentation to a phase object.
#' @export
add_segments <- function(phases, ...) {
    segs <- get_segments(phases, ...)
    attr(phases, "segments") <- segs
    phases
}

#' Add inflection points to a phase object.
#' @export
add_inflections <- function(phases, ...) {
    segs <- get_inflections(phases, ...)
    attr(phases, "inflections") <- segs
    phases
}

#' Calculate inflection points.
#' @export
get_inflections <- function(phases, ma.win=ceiling(nrow(phases))*.02,
                            spar=.001, cut=TRUE, plot=FALSE, verb=0) {

    py <- phases$angle[phases$order]
    px <- phases$phase[phases$order]

    ## shift phases to remove circular jumps
    py <- remove_jumps(py, verb=verb)
              
    ## get maxima
    ## TODO: inflection is probably overkill and can likely be done
    ## without sm.spline, just moving average and manual diffs or linregs.
    maxl <- inflection(px, py, circular=TRUE,
                       ma.win=ma.win, spar=spar, plot=plot, n=1, cut=FALSE)

    if ( cut ) {
        rm <- maxl$x < -pi | maxl$x > pi
        maxl <- maxl[!rm, , drop=FALSE]
    }

    ## fast or slow change?
    type <- ifelse(maxl$ddd<0, "fast", "slow")
    
    ## report segments
    data.frame(x=maxl$x, y=maxl$y, dydx=maxl$dy, type=type)
}


## get segments from d(phi-phi')/dphi' maxima
#' Calculate phase segments.
#' @export
get_segments <- function(phases, ma.win=ceiling(nrow(phases))*.02,
                         spar=.001, cut=TRUE, plot=FALSE, verb=0) {

    py <- phases$angle[phases$order]
    px <- phases$phase[phases$order]

    ## shift phases to remove circular jumps
    py <- remove_jumps(py, verb=verb)
              
    ## get maxima
    ## TODO: inflection is probably overkill and can likely be done
    ## without sm.spline, just moving average and manual diffs or linregs.
    maxl <- inflection(px, py-px, circular=FALSE, circular.x=TRUE,
                       ma.win=ma.win, spar=spar, plot=plot, n=0, cut=FALSE)

    if ( cut ) {
        rm <- maxl$x < -pi | maxl$x > pi
        maxl <- maxl[!rm, , drop=FALSE]
    }

    ## next up or next down?
    type <- ifelse(maxl$ddd<0, "steady", "change")

    ## report segments
    data.frame(x=maxl$x, y=maxl$y, dydx=maxl$dy, type=type)
}

## primitive spline-based extrema detection for phase angle vs. rank order
## via chatGPT,
## TODO:
## * avoid pspline and just use base R or even just moving averages
inflection <- function(x, y, ma.win=length(x)/20, spar=.1, n=1, 
                       circular=TRUE, circular.x=FALSE, circular.y=FALSE,
                       cut=TRUE, plot=FALSE, verb=1, ...) {

    ## TODO: do all this on circular data, avoiding the phase angle
    ## pre-processing
    px <- x
    py <- y
    
    ## shift phases to remove circular jumps
    py <- remove_jumps(py, verb=verb)


    ## circularize and later cut central portion
    if ( circular ) {
        px <- c(px-2*pi, px, px+2*pi)
        py <- c(py-2*pi, py, py+2*pi)
        circular.x <- circular.y <-TRUE
    } else if ( circular.x ) {
        px <- c(px-2*pi, px, px+2*pi)
        py <- c(py, py, py)
    } else if ( circular.y ) {
        py <- c(py-2*pi, py, py+2*pi)
        px <- c(px, px, px)
    }

    if ( ma.win>0 ) {
        px <- as.numeric(ma(px, ma.win, circular=FALSE))
        py <- as.numeric(ma(py, ma.win, circular=FALSE))
        rm.ends <- !is.na(px)
        px <- px[rm.ends]
        py <- py[rm.ends]
    }
    
    ## smoothed function for derivatives!
    ## TODO: find something more stable! use gam!?
    pyf <- pspline::sm.spline(px, py, spar=spar, ...)

    ## get derivatives in high resolution
    pxh <- seq(min(px), max(px), length.out = length(px)*10)

    dpdt <-  predict(pyf, pxh, nderiv = n+0)
    dpddt <- predict(pyf, pxh, nderiv = n+1)

    ## second derivative sign changes (indicating potential roots)
    dsign <- which(diff(sign(dpddt)) != 0)
    
    ## interpolate all roots
    roots <- numeric(length(dsign))  
    for (i in seq_along(dsign)) {
        idx <- dsign[i]
        roots[i] <- approx(dpddt[idx:(idx+1)], pxh[idx:(idx+1)], xout = 0)$y
    }
    
    ## third derivatives at extrema
    ## NOTE/TODO: should this be dsign diff location +1,
    ## or the mean of dsign:(dsign+1) ?
    ddd <- predict(pyf, pxh[dsign + 1], nderiv = n+2)
    
    
    ## approximate original data at roots
    ys <- approx(py ~ px, xout=roots)$y
    ## approximate slope at roots
    dys <- approx(dpdt ~ pxh, xout=roots)$y

    if ( plot ) {

        ## get second deriv at roots
        ddys <- approx(dpddt ~ pxh, xout=roots)$y

        ## scale slopes of inflection points
        up <- ddd<0
        ddds <- ddd
        ddds[ up] <- ddd[ up]/max(abs(ddd[ up]))
        ddds[!up] <- ddd[!up]/max(abs(ddd[!up]))

        root.col <- ifelse(up, 2, 4)

        ## specific to implemented use with n=0 (maxima of first derivative,
        ## used for d(phi-phi')/dphi') or n=1 (used for inflection points
        ## of dphi/dphi', where the latter may soon be obsolete
        rlab <- ifelse(n==0,
                       expression(d*(phi[y]-phi[x])/d*phi[x]),
                       expression(d*phi[y]/d*phi[x]))
        rcol <- ifelse(n==0, 3, 2)
        xlab <- ifelse(n==0, expression(phi[x]), "")
        ylab <- ifelse(n==0, expression(phi[y]-phi[x]), expression(phi[y]))
        
        lm <- c(-pi, pi) #* 2
        ylm <- range(py[px>=lm[1] & px<=lm[2]])
        
        plot(px, py, xlim=lm, ylim=ylm,
             cex=.3,
             xlab=xlab,
             ylab=ylab, axes=FALSE)
        if ( n==1 ) 
            points(roots, ys, pch=3, col=root.col, cex=1.5*abs(ddds), lwd=3)

        ## axes for aligned plots between get_segments and get_inflections
        if ( circular.x ) circ.axis(ifelse(n==0, 1, 3))
        else axis(ifelse(n==0, 1, 3))
        
        if ( circular.y ) circ.axis(2) else axis(2)

        par(new=TRUE)
        plot(pxh, dpdt, xlim=lm, col=2, type="l", lwd=2,
             xlab=NA, ylab=NA, axes=FALSE)
        points(roots, dys, pch=1, col=root.col, cex=1.5*abs(ddds), lwd=2)

        axis(4, col=rcol, col.axis=rcol)
        mtext(rlab, 4, par("mgp")[1], col=rcol)

        ## axes for aligned plots between get_segments and get_inflections
        if ( n==1 )
            type <- ifelse(ddd<0, "fast", "slow")
        else if ( n==0 )
            type <- ifelse(ddd<0, "steady", "change")

        axis(ifelse(n==0, 3, 1),
             at=roots, labels=type, las=2, tcl=.25, mgp=c(0,.1,0))

        ## TODO: make optional?
        if ( n==0 ) {
            par(new=TRUE)
            plot(pxh, dpddt, xlim=lm, col=3, type="l",
                 xlab=NA, ylab=NA, axes=FALSE)
            abline(h=0, col=3)
            points(roots, ddys, col=root.col, pch=4)
        }
        if ( FALSE ) {
            par(new=TRUE)
            plot(pxh[dsign+1], ddd, xlim=lm, col=4, type="p", pch=19,
                 xlab=NA, ylab=NA, axes=FALSE, ylim=rev(range(ddd)))
            abline(h=0, col=4)
            axis(4, col=4, col.axis=4)
        }
    }
    
    ## collect data
    infl <- data.frame(x=roots,
                       y=ys,
                       dy=dys,
                       ddd=ddd,
                       type=ifelse(ddd<0, "max", "min"))

    ## cut padded data
    if ( cut ) 
        infl <- infl[infl$x >=-pi & infl$x<pi, ] 
    
    ## add first derivative as attribute
    if ( n==1 )
        attr(infl, "dp") <- approx(dpdt ~ pxh, xout=x)$y
    else if ( n==0 )
        attr(infl, "dp") <- approx(dpddt ~ pxh, xout=x)$y

    infl
}
