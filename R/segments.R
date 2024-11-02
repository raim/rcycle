
## CALCULATE PHASE SEGMENTS


#' @export
segments <- function(phases, spar=.001, plot=TRUE, verb=1) {

    pca <- attr(phases, 'pca')
    
    phi   <- pca$rotation$phase[pca$rotation$order]
    theta <- pca$rotation$angle[pca$rotation$order]

    ## remove jumps from theta
    theta <- remove_jumps(theta, shift=TRUE, verb=verb)
    
    ## circularize
    phi   <- c(  phi-2*pi,   phi,   phi+2*pi)
    theta <- c(theta-2*pi, theta, theta+2*pi)


    ## SMOOTH SPLINE FITS
    ## theta = f(phi)
    dtpf <- pspline::sm.spline(phi, theta, spar=spar)
         
    ## get smoothed theta
    thetah <- predict(dtpf, phi)
    ## theta-phi
    dtph <- thetah - phi
        
    ## FIND ROOTS:
    ## d(theta-phi)/dphi = dtheta/dphi -1
    ## first derivative -> roots, 2nd derivative min/max
    ddtph <- predict(dtpf, phi, nderiv=1) -1
    
    dsign <- which(diff(sign(ddtph)) != 0)
    ## interpolate all roots
    roots <- numeric(length(dsign))  
    for (i in seq_along(dsign)) {
        idx <- dsign[i]
        roots[i] <- approx(x = ddtph[idx:(idx+1)],
                           y = phi[idx:(idx+1)],
                           xout = 0)$y
    }
    ## is it a maximum or minimum?
    ## TODO: filter based on ddy? e.g. if smaller than 5% of max?
    ddy <- predict(dtpf, roots, nderiv = 2) 
    ismax <- ddy < 0


    ## classify phases by segments
    ## TODO: cyclize: fuse first and last
    segs <- cut(pca$rotation$phase, breaks=roots, include.lowest = TRUE)

    ## TODO: also get dtheta/dphi max/min as alignment anchors

    if ( plot ) {
        plot(phi, theta-phi, pch=20, cex=.3, col="gray",
             xlim=c(-pi, pi), xlab='', ylab='', axes=FALSE)
        lines(phi, dtph, col=2, lwd=2)
        abline(v=roots, lwd=.5, col=2)
        axis(4, col=2, col.axis=2)
        mtext(expression(theta-phi), 4, par('mgp')[1], col=2)
        par(new=TRUE)
        plot(phi, theta,
             col=cut(phi, breaks=roots, include.lowest = TRUE),
             cex=.1, pch=20,
             xlab=expression(phase~phi), xlim=c(-pi, pi),ylim=c(-pi, pi),
             ylab=expression(angle~theta), axes=FALSE)
        circ.axis(1:2)
    } 

    x <- pca$rotation$phase
    pca$rotation <- cbind(pca$rotation,
                          segment=segs,
                          y=predict(dtpf, x, nderiv = 0) -x,
                          dy=predict(dtpf, x, nderiv = 1) -x,
                          ddy=predict(dtpf, x, nderiv = 2) -x,
                          dddy=predict(dtpf, x, nderiv = 3) -x)
    pca$roots <- roots

    ## TODO: add roots or segment table
    pca
    
    ## 1. smoothed theta(phi)
    ## 2. add dtheta/dphi and d(theta-phi)/dphi to pca$rotation
    ## 3. define segments and add as extra table,
    ##    TODO: correct segment coors when phase shifting
    ## * (i) inflection segments d(theta-phi)/dphi,
    ##   (ii) maximal slopes,
    ##   (iii) fuse (i) and (ii), take max (ii) between segments from (i)
    ## 4. classify segments by comparison with cohort phase angles

    ### extra functions, optionally already here:
    ## a) shift phases to max dtheta/dphi or max d(theta-phi)/dphi
    ## b) allow phase shifting, classification, etc. based on
    ##    segment means instead of just PC1/PC2

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
