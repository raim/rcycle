
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
find_maxima <- function(phi, theta, ND=1, spar=spar, lim=c(-pi,pi),
                        col, colf=segmenTools::arno) {

    ## calculate penalized spline 
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
 
    ## color by min/max derivative
    colv <- predict(dtpf, roots, nderiv = ND+1)
    dcol <- segmenTools::num2col(colv, colf=colf)

    roots <- data.frame(phi=roots,
                        dtheta  =  dtheta,
                        ddtheta = ddtheta,
                        dddtheta=dddtheta,
                        dmax = dmax,
                        ddmax=ddmax,
                        dcol=dcol)

    ## remove roots from cyclization
    ## NOTE: include -pi and exclude pi
    rmc <- roots$phi <= lim[1] | roots$phi > lim[2]
    roots <- roots[!rmc,]

    ## segment ID
    roots$ID <- c(1:nrow(roots))

    ## colors
    roots$col <- 1:nrow(roots)

   
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
                     names=method, 
                     spar=.001, #spari=100*spar,
                     P, Pscale=.05, L=10, jumps=FALSE,
                     plot=TRUE, verb=0, ...) {

    ## TODO: make sure phi is within -pi:pi or otherwise adjust lim below
    
    ## original phi and theta
    ph <- phases$rotation$phi[phases$rotation$order]
    th <- phases$rotation$theta[phases$rotation$order]

    ## remove jumps from theta
    theta <- remove_jumps(th, center=FALSE, verb=verb)
    
    ## circularize
    phi   <- c(   ph-2*pi,    ph,    ph+2*pi)
    theta <- c(theta-2*pi, theta, theta+2*pi)

    if ( verb>0 )
        cat(paste0('calculating segments for ', length(phi),  ' cells,\n',
                   'using method(s): ', paste(method,collapse=';'), '\n'))


    ## SMOOTH SPLINE FITS
    ## theta = f(phi)
    
    ## NOTE: a low spar is required to correctly get segments,
    ## i.e. the shoulders of f(phi) or the extrema of
    ## d(theta-phi)/dphi. To get inflection points, points of max and
    ## min slope, we need the second and third derivatives, and
    ## require a stronger smoothing to avoid to get too many points.

    ## TODO: every segmentation yields different results to be
    ## added to phases$rotation; make this function return those,
    ## and only add in separate function; 

    ## SHOULDERS
    if ( 'shoulder' %in% method ) {

        ## roots of first derivative
        ## d(theta-phi)/dphi = dtheta/dphi -1 = 0
        maxima <- find_maxima(phi, theta-phi, ND=1, spar=spar, lim=c(-pi, pi))

        ## TODO: correct derivative for theta-phi!
        ##maxima$breaks$dtheta <- maxima$breaks$dtheta 

        shoulder <- maxima$breaks
        dtpf <- maxima$splinef
        
        ## classify phases by shoulder segments
        segs <- phase_segments(phases$rotation$phi,
                               breaks=shoulder$phi)

        x <- phases$rotation$phi
        df <- data.frame(segment=segs,
                         theta  = predict(dtpf, x, nderiv = 0) + x,
                         dtheta = predict(dtpf, x, nderiv = 1),
                         ddtheta= predict(dtpf, x, nderiv = 2))
        nm <- names[method=='shoulder']
        colnames(df) <- paste0(nm, "_", colnames(df))

        ## replace previus
        old <- which(colnames(phases$rotation)%in%colnames(df))
        if ( length(old)>0 )
            phases$rotation <-phases$rotation[,-old] 
        phases$rotation <- cbind(phases$rotation,
                                 df)
        
        ## add breaks
        phases[[names[method=='shoulder']]] <- shoulder   
    }

    ## INFLECTIONS
    if ( 'inflection' %in% method ) {

        ## roots of second derivative
        ## d^2(theta-phi)/dphi^2 = 0
        infl <- find_maxima(phi, theta-phi, ND=2, spar=spar, lim=c(-pi, pi))

        
        inflection <- infl$breaks
        dtpf <- infl$splinef
        
        ## only take maxima! segments: max-to-max
        inflection <- inflection[inflection$ddmax,]

        ## ADD ID AND COLORS
        inflection$ID <- inflection$col <- 1:nrow(inflection)
     
        ## classify phases by max slope segments
        segs <-  phase_segments(phases$rotation$phi,
                                breaks=inflection$phi)
    
        ## ADD TO PHASES OBJECT
        ## thetah  : smoothed theta, dtpf, deriv 0
        ## dtheta  : dtheta/dphi -1, dtpf, deriv 1
        ## ddtheta : d2theta/dphi2,  dtpi, deriv 2
        x <- phases$rotation$phi
        df <- data.frame(segment= segs,
                         theta  = predict(dtpf, x, nderiv = 0) + x,
                         dtheta = predict(dtpf, x, nderiv = 1),
                         ddtheta= predict(dtpf, x, nderiv = 2))
        nm <- names[method=='inflection']
        colnames(df) <- paste0(nm, "_", colnames(df))

        ## replace previus
        old <- which(colnames(phases$rotation)%in%colnames(df))
        if ( length(old)>0 )
            phases$rotation <-phases$rotation[,-old] 
        phases$rotation <- cbind(phases$rotation,
                                 df)
        
        ## add breaks
        phases[[names[method=='inflection']]] <- inflection 

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
        brks <- data.frame(phi  =dps$segments$x2,
                           slope=dps$segments$slope,
                           smax =dps$segments$slope >0)
        lim <- c(-pi, pi)
        rmc <- brks$phi <= lim[1] | brks$phi > lim[2]
        brks <- brks[!rmc,,drop=FALSE]

        ## only 1 segment
        if ( nrow(brks)==0 )
            brks <- data.frame(phi=pi, slope=NA, smax=TRUE)
        
        ## segment ID
        brks$ID <- c(1:nrow(brks))

        
        ## classify phases by max slope segments
        dsegs <-  phase_segments(phases$rotation$phi, breaks=brks$phi)
        
        ## ADD TO PHASES OBJECT
        ## thetah  : smoothed theta, dtpf, deriv 0
        ## dtheta  : dtheta/dphi -1, dtpf, deriv 1
        ## ddtheta : d2theta/dphi2,  dtpi, deriv 2
        x <- phases$rotation$phi
        df <- data.frame(segment= dsegs,
                         theta  = predict(dps, xout=x)$y + x)
        nm <- names[method=='dpseg']
        colnames(df) <- paste0(nm, "_", colnames(df))

        ## replace previus
        old <- which(colnames(phases$rotation)%in%colnames(df))
        if ( length(old)>0 )
            phases$rotation <-phases$rotation[,-old] 
        phases$rotation <- cbind(phases$rotation,
                              df)
        
        ## add breaks
        phases[[names[method=='dpseg']]] <- brks 

        ## add full dpseg object; todo: required?
        phases$dpseg.object <- dps
    }
    
    phases$processing <- c(phases$processing,
                        paste0('segments:', paste0(method,collapse=';')))

    if ( plot ) {
           if ( verb>0 )
               cat(paste('plotting segments\n'))
           plotSegments(phases, 
                        method=names[1], ...)
    }
    
    invisible(phases)
    
}


#' Plot the curve sketching of theta=f(phi).
#' @param difference plot the difference theta-phi instead of
#'     theta. This emphasizes the noise in the original data and the
#'     smoothing effects.
#' @export
plotSegments <- function(phases, difference=FALSE, center=TRUE,
                         method='shoulder') {

    
    ord <- phases$rotation$order
    phi <- phases$rotation$phi[ord]
    theta.orig <- phases$rotation$theta[ord]

    ## remove jumps from all theta
    ## TODO: remove_jumps doesn't work for thetah in dpseg and shoulder
    ## and the current approach via approx_phase has errors
    theta <- remove_jumps(theta.orig, center=center, verb=FALSE)

    ## shift cohort theta accordingly
    ## TODO: is this correct?
    thetac <- approx_phase(x=theta.orig, y=theta, xout=phases$x$theta)$y
    phic <- approx_phase(x=theta, y=phi, xout=thetac)$y

    deriv <- thetah <- breaks <- brks <- cl <- NULL
    
    ## get common structure 
    if ( method %in% names(phases) ) {
        
        cln <- paste0(method,'_segment')
        cl <- phases$rotation[,cln][phases$rotation$order]
        base.cex <- 1
        
        thn <- paste0(method,'_theta')
        thetah <- phases$rotation[[thn]][ord]
        
        breaks <- brks <- phases[[method]]
    }

    if (  method %in% c('shoulder','slope') ) {
        
        dthn <- paste0(method,'_dtheta')
        deriv <- phases$rotation[[dthn]][ord]
        ylab4 <- expression(d*theta/d*phi-1)

        if ( method=='slope' ) {
            brks <- NULL
            cl <- num2col(deriv)
            base.cex <- 0
        }
        
        magn <- "dtheta"
        maxn <- 'dmax'
        dcol <- bcol <- 4
    }
    if ( method %in% c('inflection') ) {


        dthn <- paste0(method,'_ddtheta')
        deriv <- phases$rotation[[dthn]][ord]
        ylab4 <- expression(d^2*theta/d*phi^2)
           
        brks <- breaks[breaks$ddmax,] ## TODO: handle upstream
        if ( method=='slope' )
            brks <- NULL
        
        magn <- "ddtheta"
        maxn <- 'ddmax'
        dcol <- bcol <- 3
    }
    if ( 'dpseg' %in% method ) {

        deriv <- rep(0, length(thetah))
        
        ylab4 <- expression(piecewise~linear)
        magn <- "slope"
        maxn <- 'smax'
        dcol <- NA
        bcol <- 5
    }

    ## align smoothed theta
    ## TODO: is there a better solution?
    ##if ( !is.null(thetah) ) 
    ##    thetah <- align_phase(thetah, target=theta, center=center)
        ##thetah <- approx_phase(x=theta.orig, y=theta, xout=thetah)$y
    
    if ( !is.null(deriv) ) {
        
        plot(phi, deriv, type='l', col=dcol,
             xlim=c(-pi, pi), xlab='', ylab='', axes=FALSE) 

        abline(h=0, col=dcol, lwd=.25)
        
        if ( !is.null(brks) )
            abline(v=brks$phi, lwd=.25, col=bcol) 

        axis(4, col=dcol, col.axis=dcol)
        mtext(ylab4, 4, par('mgp')[1], col=dcol)
        par(new=TRUE)
    }

    ylab <- expression(angle~theta)
    yaxis <- circ.axis
    ylim <- c(-pi, pi)
    if ( !center ) ylim <- range(theta)
    
    if ( difference ) {
        theta <- theta - phi
        if ( !is.null(thetah) ) thetah <- thetah -phi
        ylab <- expression(difference~theta-phi)
        yaxis <- axis
        ylim <- range(theta)
    }

    plot(1, col=NA, axes=FALSE,
         xlim=c(-pi, pi), ylim=ylim,
         xlab=expression(phase~phi), ylab=ylab)

    ## plot segments
    if ( !is.null(cl) ) {
        
        points(phi, theta, col=cl, cex=.5, pch=1)
        if ( !is.null(thetah) ) 
            lines(phi, thetah, col=1, lwd=1) # smoothed
       
        if ( method=='slope' & !center ) {
            abline(v=0, col=1, lwd=.25)
            abline(h=0, col=1, lwd=.25)
        }
        ## plot break points
        points(approx(phi, theta, breaks$phi),
               col=dcol, lwd=2,
               pch=ifelse(breaks[,maxn],24,25),
               cex=base.cex+1.5*abs(breaks[,magn])/max(abs(breaks[,magn])))
    } else points(phi, theta, col=1, cex=.5, pch=1)
    circ.axis(1)
    yaxis(2)

    ## plot cohort phases
    ## NOTE: assuming prefix naming
    ## TODO: add official cohort name vector to object, or
    ## avoid long prefixes for separate PCA
    for ( i in 1:nrow(phases$x) )
        axis(3, at=phases$x$phi[i],
             labels=sub(".*_","",rownames(phases$x))[i], las=2)
    ## if right y-axis is available, plot cohort theta
    if ( is.null(deriv) ) 
        for ( i in 1:nrow(phases$x) )
            axis(4, at=thetac[i] - ifelse(difference, phic[i], 0),
                 labels=sub(".*_","",rownames(phases$x))[i], las=2)
    
} 

