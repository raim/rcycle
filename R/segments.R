
## CALCULATE PHASE SEGMENTS

#' Classify phases into segments.
#' @export
phase_segments <- function(phi, breaks, lim=c(-pi, pi)) {

    breaks <- unique(c(lim[1], breaks, lim[2]))
    
    segs <- cut(phi, breaks=breaks, include.lowest = TRUE)
    
    ## numeric segments
    segs <- as.numeric(segs)
    segs <- segs-min(segs)+1

    ## fuse first and last
    segs[segs==max(segs)] <- min(segs)

    segs
}

## get derivatives from smoothed spline
find_maxima <- function(phi, theta, ND=1, spar=spar, lim=c(-pi,pi)) {
    
    dtpf <- pspline::sm.spline(phi, theta, spar=spar)  

    ## get roots of requested derivatives
    dtheta <- predict(dtpf, phi, nderiv = ND) 
    
    ## interpolate all roots
    roots <- find_roots(phi, dtheta)

    ## add all derivatives at roots
    dtheta   <- predict(dtpf, roots, nderiv = 1)
    ddtheta  <- predict(dtpf, roots, nderiv = 2)
    dddtheta <- predict(dtpf, roots, nderiv = 3)
    
    ## add min/max classification
    dmax <- ddtheta < 0
    ## add min/max classification
    ddmax <- dddtheta < 0
 
    roots <- data.frame(phi=roots,
                        dtheta  =  dtheta,
                        ddtheta = ddtheta,
                        dddtheta=dddtheta,
                        dmax = dmax,
                        ddmax=ddmax)

    ## remove roots from cyclization
    ## NOTE: include -pi and exclude pi
    rmc <- roots$phi <= lim[1] | roots$phi > lim[2]
    roots <- roots[!rmc,]

    ## segment ID
    roots$ID <- c(1:nrow(roots))
   
    list(breaks=roots, splinef=dtpf)
}

## primitive root finding
## TODO: better version? uniroot?
find_roots <- function(phi, dtheta) {

    dsign <- which(diff(sign(dtheta)) != 0)
    roots <- numeric(length(dsign))  
    for (i in seq_along(dsign)) {
        idx <- dsign[i]
        roots[i] <- approx(x = dtheta[idx:(idx+1)],
                           y =    phi[idx:(idx+1)],
                           xout = 0)$y
    }
    roots
}



#' @export
segments <- function(phases, 
                     method=c('shoulder', 'inflection', 'dpseg'),
                     spar=.001, #spari=100*spar,
                     P, Pscale=100, L=10, jumps=FALSE,
                     plot=TRUE, verb=0) {

    pca <- attr(phases, 'pca')

    ## TODO: make sure phi is within -pi:pi or otherwise adjust lim below
    
    ## original phi and theta
    ph <- pca$rotation$phi[pca$rotation$order]
    th <- pca$rotation$theta[pca$rotation$order]

    ## remove jumps from theta
    theta <- remove_jumps(th, shift=TRUE, verb=verb)
    
    ## circularize
    phi   <- c(   ph-2*pi,    ph,    ph+2*pi)
    theta <- c(theta-2*pi, theta, theta+2*pi)

    if ( verb>0 )
        cat(paste0('calculating segments for ', length(phi),  ' cells,\n',
                   'using method(s):', paste(method,collapse=';'), '\n'))


    ## SMOOTH SPLINE FITS
    ## theta = f(phi)
    
    ## NOTE: a low spar is required to correctly get segments,
    ## i.e. the shoulders of f(phi) or the extrema of
    ## d(theta-phi)/dphi. To get inflection points, points of max and
    ## min slope, we need the second and third derivatives, and
    ## require a stronger smoothing to avoid to get too many points.

    ## SHOULDERS
    if ( 'shoulder' %in% method ) {

        ## roots of first derivative
        ## d(theta-phi)/dphi = dtheta/dphi -1 = 0
        maxima <- find_maxima(phi, theta-phi, ND=1, spar=spar, lim=c(-pi, pi))

        shoulder <- maxima$breaks
        dtpf <- maxima$splinef
        
        ## classify phases by shoulder segments
        segs <- phase_segments(pca$rotation$phi,
                               breaks=shoulder$phi)
        
        x <- pca$rotation$phi
        pca$rotation$shoulder <- segs # shoulder segments
        pca$rotation$theta.s  <- predict(dtpf, x, nderiv = 0) + x
        pca$rotation$dtheta   <- predict(dtpf, x, nderiv = 1)
    
        pca$shoulder <- shoulder   
    }

    ## INFLECTIONS
    if ( 'inflection' %in% method ) {

        ## roots of second derivative
        ## d^2(theta-phi)/dphi^2 = 0
        infl <- find_maxima(phi, theta-phi, ND=2, spar=spar, lim=c(-pi, pi))

        inflection <- infl$breaks
        dtpf <- infl$splinef
        
        ## classify phases by max slope segments
        isegs <-  phase_segments(pca$rotation$phi,
                                 breaks=inflection$phi[inflection$ddmax])
    
        ## ADD TO PHASES OBJECT
        ## thetah  : smoothed theta, dtpf, deriv 0
        ## dtheta  : dtheta/dphi -1, dtpf, deriv 1
        ## ddtheta : d2theta/dphi2,  dtpi, deriv 2
        x <- pca$rotation$phi
        pca$rotation$inflection <- isegs # inflection segments
        pca$rotation$ddtheta    <- predict(dtpf, x, nderiv = 2)
    
        pca$inflection <- inflection
    }

    ## DPSEG: piecewise linear
    if ( 'dpseg' %in% method ) {

        ## estimate penalty
        if ( missing(P) )
            P <- dpseg::estimateP(x=ph, y=th-ph) * Pscale
        
        ## cut to make faster
        rmp <- phi < -1.5*pi | phi > 1.5*pi

        ## run dpseg
        dps <- dpseg::dpseg(x=phi[!rmp], y=theta[!rmp] - phi[!rmp],
                            jumps=jumps, P=P, minl=L, verb=verb, add.lm=TRUE)

        ## generate break table
        segs <- data.frame(phi  =dps$segments$x2,
                           slope=dps$segments$slope,
                           smax =dps$segments$slope >0)
        lim <- c(-pi, pi)
        rmc <- segs$phi <= lim[1] | segs$phi > lim[2]
        segs <- segs[!rmc,,drop=FALSE]

        ## only 1 segment
        if ( nrow(segs)==0 )
            segs <- data.frame(phi=pi, slope=NA, smax=TRUE)
        
        ## segment ID
        segs$ID <- c(1:nrow(segs))

        
        ## classify phases by max slope segments
        dsegs <-  phase_segments(pca$rotation$phi, breaks=segs$phi)
        
        
        ## ADD PHASES
        x <- pca$rotation$phi
        pca$rotation$dpseg   <- dsegs
        pca$rotation$theta.d <- predict(dps, xout=x)$y
    
        pca$dpseg <- segs
        pca$dpseg.object <- dps
    }
    
    pca$processing <- c(pca$processing,
                        paste0('segments:', paste0(method,collapse=';')))

    attr(phases, 'pca') <- pca

    if ( plot ) {
           if ( verb>0 )
               cat(paste('plotting segments\n'))
           plotSegments(phases, difference=TRUE, method=method[1])
    }
    
    phases
    ### extra functions, optionally already here:
    ## a) shift phases to max dtheta/dphi or max d(theta-phi)/dphi
    ## b) allow phase shifting, classification, etc. based on
    ##    segment means instead of just PC1/PC2

}


#' Plot the curve sketching of theta=f(phi).
#' @param difference plot the difference theta-phi instead of
#'     theta. This emphasizes the noise in the original data and the
#'     smoothing effects.
#' @export
plotSegments <- function(phases, difference=FALSE, shift=FALSE,
                         method='shoulder') {

    pca <- attr(phases, 'pca') 
    
    ord <- pca$rotation$order
    phi <- pca$rotation$phi[ord]
    theta <- pca$rotation$theta[ord]

    cl <- pca$rotation[,method][pca$rotation$order]
    theta <- remove_jumps(theta, shift=shift, verb=FALSE)

    deriv <- thetah <- NULL
    
    if ( 'shoulder' %in% method ) {
        thetah <- pca$rotation$theta.s[ord]

        deriv <- pca$rotation$dtheta[ord]
        breaks <- pca$shoulder
        ylab4 <- expression(d*theta/d*phi-1)
        magn <- "dtheta"
        maxn <- 'dmax'
        dcol <- 4
    }
    if ( 'inflection' %in% method ) {
        deriv <- pca$rotation$ddtheta[ord]
        breaks <- pca$inflection
        ylab4 <- expression(d^2*theta/d*phi^2)
        magn <- "ddtheta"
        maxn <- 'ddmax'
        dcol <- 3
    }
    if ( 'dpseg' %in% method ) {
        thetah <- pca$rotation$theta.d[ord] + phi
        breaks <- pca$dpseg
        ylab4 <- expression(piecewise~linear)
        magn <- "slope"
        maxn <- 'smax'
        dcol <- 5
    }


    if ( !is.null(deriv) ) {
        plot(phi, deriv, type='l', col=dcol,
             xlim=c(-pi, pi), xlab='', ylab='', axes=FALSE) 
        abline(h=0, col=dcol, lwd=.5)
        axis(4, col=dcol, col.axis=dcol)
        mtext(ylab4, 4, par('mgp')[1], col=dcol)
        par(new=TRUE)
    }

    ylab <- expression(angle~theta)
    yaxis <- circ.axis
    ylim <- c(-pi, pi)
    if ( !shift ) ylim <- range(theta)
    if ( difference ) {
        theta <- theta -phi
        if ( !is.null(thetah) ) thetah <- thetah -phi
        ylab <- expression(differene~theta-phi)
        yaxis <- axis
        ylim <- range(theta)
    }

    plot(1, col=NA, axes=FALSE,
         xlim=c(-pi, pi), ylim=ylim,
         xlab=expression(phase~phi), ylab=ylab)

    ## plot segments
    points(phi, theta, col=cl, cex=.5, pch=1)
    if ( !is.null(thetah) )
        lines(phi, thetah, col=1, lwd=1) # smoothed
       
    ## plot break points
    points(approx(phi, theta, breaks$phi),
           col=dcol, lwd=2,
           pch=ifelse(breaks[,maxn],24,25),
           cex=.5+1.5*abs(breaks[,magn])/max(abs(breaks[,magn])))
    circ.axis(1)
    yaxis(2)

    ## plot cohort phases
    ## NOTE: assuming prefix naming
    ## TODO: add official cohort name vector to object, or
    ## avoid long prefixes for separate PCA
    for ( i in 1:nrow(pca$x) )
        axis(3, at=pca$x$phi[i],
             labels=sub(".*_","",rownames(pca$x))[i], las=2)
} 


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

    theta <- phases$angle[phases$order]
    phi <- phases$phase[phases$order]

    ## shift phases to remove circular jumps
    theta <- remove_jumps(theta, verb=verb)
              
    ## get maxima
    ## TODO: inflection is probably overkill and can likely be done
    ## without sm.spline, just moving average and manual diffs or linregs.
    maxl <- inflection(phi, theta, circular=TRUE,
                       ma.win=ma.win, spar=spar, plot=plot, n=1, cut=FALSE)

    if ( cut ) {
        rm <- maxl$x < -pi | maxl$x > pi
        maxl <- maxl[!rm, , drop=FALSE]
    }

    ## fast or slow change?
    type <- ifelse(maxl$ddd<0, "fast", "slow")
    
    ## report segments
    segs <- data.frame(x=maxl$x, y=maxl$y, dydx=maxl$dy, type=type)

    ## report full smoothed data+derivatives
    smos <- attr(maxl, 'smos')
    
    segs
}


## get segments from d(phi-phi')/dphi' maxima
#' Calculate phase segments.
#' @export
get_segments <- function(phases, ma.win=ceiling(nrow(phases))*.02,
                         spar=.001, cut=TRUE, plot=FALSE, verb=0) {

    theta <- phases$angle[phases$order]
    phi <- phases$phase[phases$order]

    ## shift phases to remove circular jumps
    theta <- remove_jumps(theta, shift=TRUE, verb=verb)

    if ( plot ) {
        plot(phi, theta, col="gray", axes=FALSE, xlab='', ylab='', type='l')
        par(new=TRUE)
    }
              
    ## get maxima
    ## TODO: inflection is probably overkill and can likely be done
    ## without sm.spline, just moving average and manual diffs or linregs.
    maxl <- inflection(phi, theta-phi, circular=FALSE, circular.x=TRUE,
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
                       use.pspline=TRUE,
                       circular=TRUE, circular.x=FALSE, circular.y=FALSE,
                       cut=TRUE, plot=FALSE, verb=1, ...) {

    ## TODO: do all this on circular data, avoiding the phase angle
    ## pre-processing
    phi <- x
    theta <- y
    
    ## shift phases to remove circular jumps
    ## NOTE: only theta shifted, so phi coors stay the same
    theta <- remove_jumps(theta, shift=TRUE, verb=verb)


    ## circularize and later cut central portion
    if ( circular ) {
        phi <- c(phi-2*pi, phi, phi+2*pi)
        theta <- c(theta-2*pi, theta, theta+2*pi)
        circular.x <- circular.y <-TRUE
    } else if ( circular.x ) {
        phi <- c(phi-2*pi, phi, phi+2*pi)
        theta <- c(theta, theta, theta)
    } else if ( circular.y ) {
        theta <- c(theta-2*pi, theta, theta+2*pi)
        phi <- c(phi, phi, phi)
    }

    if ( ma.win>0 ) {
        phi <- as.numeric(ma(phi, ma.win, circular=FALSE))
        theta <- as.numeric(ma(theta, ma.win, circular=FALSE))
        rm.ends <- !is.na(phi)
        phi <- phi[rm.ends]
        theta <- theta[rm.ends]
    }
    
    ## smoothed function for derivatives!
    ## TODO: find something more stable! use gam!?
    if ( use.pspline )
        thetaf <- pspline::sm.spline(phi, theta, spar=spar, ...)
    else thetaf <- smooth.spline(phi, theta, spar=spar, ...)

    ## get derivatives in high resolution
    phih <- seq(min(phi), max(phi), length.out = length(phi)*10)

    if ( use.pspline ) {
        dpdt <-  predict(thetaf, phih, nderiv = n+0)
        dpddt <- predict(thetaf, phih, nderiv = n+1)
    } else {
        dpdt <-  predict(thetaf, phih, deriv = n+0)$y
        dpddt <- predict(thetaf, phih, deriv = n+1)$y
    }

    ## second derivative sign changes (indicating potential roots)
    dsign <- which(diff(sign(dpddt)) != 0)
    
    ## interpolate all roots
    roots <- numeric(length(dsign))  
    for (i in seq_along(dsign)) {
        idx <- dsign[i]
        roots[i] <- approx(dpddt[idx:(idx+1)], phih[idx:(idx+1)], xout = 0)$y
    }
    
    ## third derivatives at extrema
    ## NOTE/TODO: should this be dsign diff location +1,
    ## or the mean of dsign:(dsign+1) ?
    if ( use.pspline ) 
        ddd <- predict(thetaf, phih[dsign + 1], nderiv = n+2)
    else ddd <- predict(thetaf, phih[dsign + 1], deriv = n+2)$y

    
    ## approximate original data at roots
    ys <- approx(theta ~ phi, xout=roots)$y
    ## approximate slope at roots
    dys <- approx(dpdt ~ phih, xout=roots)$y

    if ( plot ) {

        ## get second deriv at roots
        ddys <- approx(dpddt ~ phih, xout=roots)$y

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
                       expression(d*(theta-phi)/d*phi),
                       expression(d*theta/d*phi))
        rcol <- ifelse(n==0, 3, 2)
        xlab <- ifelse(n==0, expression(phi), "")
        ylab <- ifelse(n==0, expression(theta-phi), expression(theta))
        
        lm <- c(-pi, pi) #* 2
        ylm <- range(theta[phi>=lm[1] & phi<=lm[2]])
        
        plot(phi, theta, xlim=lm, ylim=ylm,
             cex=.3, xlab=xlab, ylab=ylab, axes=FALSE)
        if ( n==1 ) 
            points(roots, ys, pch=3, col=root.col, cex=1.5*abs(ddds), lwd=3)

        ## axes for aligned plots between get_segments and get_inflections
        if ( circular.x ) circ.axis(ifelse(n==0, 1, 3))
        else axis(ifelse(n==0, 1, 3))
        
        if ( circular.y ) circ.axis(2) else axis(2)

        par(new=TRUE)
        plot(phih, dpdt, xlim=lm, col=2, type="l", lwd=2,
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
            plot(phih, dpddt, xlim=lm, col=3, type="l",
                 xlab=NA, ylab=NA, axes=FALSE)
            abline(h=0, col=3)
            points(roots, ddys, col=root.col, pch=4)
        }
        if ( FALSE ) {
            par(new=TRUE)
            plot(phih[dsign+1], ddd, xlim=lm, col=4, type="p", pch=19,
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
    
    if ( use.pspline ) {
        smos <- cbind(predict(thetaf, x, nderiv = 0),
                      predict(thetaf, x, nderiv = 1),
                      predict(thetaf, x, nderiv = 2))
    } else {
        smos <- cbind(predict(thetaf, x, deriv = 0)$y,
                      predict(thetaf, x, deriv = 1)$y,
                      predict(thetaf, x, deriv = 2)$y)
    }

    ## all smoothed data as table
    attr(infl, 'smos') <- smos

    ## add first derivative as attribute
    if ( n==1 )
        attr(infl, "dp") <- approx(dpdt ~ phih, xout=x)$y
    else if ( n==0 )
        attr(infl, "dp") <- approx(dpddt ~ phih, xout=x)$y

    infl
}


## primitive spline-based extrema detection for phase angle vs. rank order
## via chatGPT,
## TODO:
## * avoid pspline and just use base R or even just moving averages
curves <- function(x, y, ma.win=length(x)/20, spar=.1, n=1,
                        use.pspline=TRUE,
                        circular=TRUE, circular.x=FALSE, circular.y=FALSE,
                        cut=TRUE, plot=FALSE, verb=1, ...) {

    ## TODO: do all this on circular data, avoiding the phase angle
    ## pre-processing
    phi <- x
    theta <- y
    
    ## shift phases to remove circular jumps
    ## NOTE: only theta shifted, so phi coors stay the same
    theta <- remove_jumps(theta, shift=TRUE, verb=verb)


    ## circularize and later cut central portion
    if ( circular ) {
        phi <- c(phi-2*pi, phi, phi+2*pi)
        theta <- c(theta-2*pi, theta, theta+2*pi)
        circular.x <- circular.y <-TRUE
    } else if ( circular.x ) {
        phi <- c(phi-2*pi, phi, phi+2*pi)
        theta <- c(theta, theta, theta)
    } else if ( circular.y ) {
        theta <- c(theta-2*pi, theta, theta+2*pi)
        phi <- c(phi, phi, phi)
    }

    if ( ma.win>0 ) {
        phi <- as.numeric(ma(phi, ma.win, circular=FALSE))
        theta <- as.numeric(ma(theta, ma.win, circular=FALSE))
        rm.ends <- !is.na(phi)
        phi <- phi[rm.ends]
        theta <- theta[rm.ends]
    }
    
    ## smoothed function for derivatives!
    ## TODO: find something more stable! use gam!?
    if ( use.pspline )
        thetaf <- pspline::sm.spline(phi, theta, spar=spar, ...)
    else thetaf <- smooth.spline(phi, theta, spar=spar, ...)

    ## get derivatives in high resolution
    phih <- seq(min(phi), max(phi), length.out = length(phi)*10)

    if ( use.pspline ) {
        dpdt <-  predict(thetaf, phih, nderiv = n+0)
        dpddt <- predict(thetaf, phih, nderiv = n+1)
    } else {
        dpdt <-  predict(thetaf, phih, deriv = n+0)$y
        dpddt <- predict(thetaf, phih, deriv = n+1)$y
    }

    ## second derivative sign changes (indicating potential roots)
    dsign <- which(diff(sign(dpddt)) != 0)
    
    ## interpolate all roots
    roots <- numeric(length(dsign))  
    for (i in seq_along(dsign)) {
        idx <- dsign[i]
        roots[i] <- approx(dpddt[idx:(idx+1)], phih[idx:(idx+1)], xout = 0)$y
    }
    
    ## third derivatives at extrema
    ## NOTE/TODO: should this be dsign diff location +1,
    ## or the mean of dsign:(dsign+1) ?
    if ( use.pspline ) 
        ddd <- predict(thetaf, phih[dsign + 1], nderiv = n+2)
    else ddd <- predict(thetaf, phih[dsign + 1], deriv = n+2)$y

    
    ## approximate original data at roots
    ys <- approx(theta ~ phi, xout=roots)$y
    ## approximate slope at roots
    dys <- approx(dpdt ~ phih, xout=roots)$y

    if ( plot ) {

        ## get second deriv at roots
        ddys <- approx(dpddt ~ phih, xout=roots)$y

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
                       expression(d*(theta-phi)/d*phi),
                       expression(d*theta/d*phi))
        rcol <- ifelse(n==0, 3, 2)
        xlab <- ifelse(n==0, expression(phi), "")
        ylab <- ifelse(n==0, expression(theta-phi), expression(theta))
        
        lm <- c(-pi, pi) #* 2
        ylm <- range(theta[phi>=lm[1] & phi<=lm[2]])
        
        plot(phi, theta, xlim=lm, ylim=ylm,
             cex=.3, xlab=xlab, ylab=ylab, axes=FALSE)
        if ( n==1 ) 
            points(roots, ys, pch=3, col=root.col, cex=1.5*abs(ddds), lwd=3)

        ## axes for aligned plots between get_segments and get_inflections
        if ( circular.x ) circ.axis(ifelse(n==0, 1, 3))
        else axis(ifelse(n==0, 1, 3))
        
        if ( circular.y ) circ.axis(2) else axis(2)

        par(new=TRUE)
        plot(phih, dpdt, xlim=lm, col=2, type="l", lwd=2,
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
            plot(phih, dpddt, xlim=lm, col=3, type="l",
                 xlab=NA, ylab=NA, axes=FALSE)
            abline(h=0, col=3)
            points(roots, ddys, col=root.col, pch=4)
        }
        if ( FALSE ) {
            par(new=TRUE)
            plot(phih[dsign+1], ddd, xlim=lm, col=4, type="p", pch=19,
                 xlab=NA, ylab=NA, axes=FALSE, ylim=rev(range(ddd)))
            abline(h=0, col=4)
            axis(4, col=4, col.axis=4)
        }
    }
    
    ## 
    roots <- data.frame(x=roots,
                        type=ifelse(ddd<0, "max", "min"))

    ## cut padded data
    if ( cut ) 
        infl <- infl[infl$x >=-pi & infl$x<pi, ] 
    
    if ( use.pspline ) {
        curves <- cbind(y=predict(thetaf, x, nderiv = 0),
                        dy=predict(thetaf, x, nderiv = 1),
                        ddy=predict(thetaf, x, nderiv = 2),
                        dddy=predict(thetaf, x, nderiv = 3))
    } else {
        curves <- cbind(y=predict(thetaf, x, deriv = 0)$y,
                        dy=predict(thetaf, x, deriv = 1)$y,
                        ddy=predict(thetaf, x, deriv = 2)$y,
                        dddy=predict(thetaf, x, deriv = 3)$y)
    }


 
    list(roots=roots, curves=curves)
}
