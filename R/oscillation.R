
## DO WAVELET ANALYSIS OF CELLS SORTED BY PHASE

#' Wavelet analysis of counts along a circular pseudophase.
#' @param phi a circular pseudophase assignement for the counts vector
#'     or matrix.
#' @param counts a vector or matrix of counts of the same dimension
#'     and order as the pseudophase phi.
#' @param ID if counts is a matrix, ID specifies the row number or row
#'     name for which the analysis is to be performed.
#' @param verb verbosity level, set to >2 to also get verbose output
#'     from the wavelet function.
#' @param ... arguments to \code{\link[WavletComp]{analyze.wavelet}}.
#'@export
get_wavelet <- function(phi, counts, ID, verb=1, ...) {

    ## default: expect a single counts vector,
    ## if ID is provided get data from full count matrix
    
    if ( !missing(ID) ) counts <- counts[ID,]

    ## ORDER
    if ( any(diff(phi)<0) ) {
        counts <- counts[order(phi)]
        phi <- sort(phi)
        if ( verb>0 )
            cat(paste("NOTE: re-ordered data along phi\n"))
    }

    ## time in fraction of cycle in radian
    dt <- mean(diff(phi)) # 
    
    ## pad circular data to remove edge effects
    
    n <- length(counts)
    n2 <- floor(n/2)
    counts <- c(counts[(n2+1):n],
                counts,
                counts[1:n2])
    ctme <- c(phi[(n2+1):n]-2*pi,
              phi,
              phi[1:n2]+2*pi)
    
    ## remove rare duplicated time points
    dups <- which(duplicated(ctme))
    if ( length(dups)>0 ) {
        warning(paste('removing', length(dups),
                      'duplicated time points\n'))
        counts <- counts[-dups]
        ctme <- ctme[-dups]
    }

    ## minimal period
    ##if ( missing(lowerPeriod) )
    lowerPeriod <- 2*dt
    ## maximal period: take all
    ##if ( missing(upperPeriod) )
    upperPeriod <- length(counts)*dt
    
    clet <- WaveletComp::analyze.wavelet(data.frame(x=counts), "x",
                                         lowerPeriod = lowerPeriod, 
                                         upperPeriod = upperPeriod, 
                                         dt = dt, 
                                         ##loess.span = LOESS.SPAN,
                                         ##make.pval = NPERM>0,
                                         ##method = "white.noise",
                                         ##n.sim = NPERM, 
                                         verbose = verb>1,
                                         ...)
    ## RESULT: expanded/filtered  data 
    list(wlet=clet, data=data.frame(time=ctme, counts=counts))
}

## TODO: finish custom plotter
plot_wavelet <- function(wlet, type = "Power", col, p.min) {

    clet <- wlet$wlet
    ctme <- wlet$data[,1] # time in radian
    
    ## period in fraction of a full cycle
    period <- clet$Period/(2*pi)
    ylab <- paste0('cycles')
    ylim <- log2(c(5e-4,3))

    ## values to plot
    values <- clet[[type]]

    if ( type%in%c("Power") ) {
    } else if ( type%in%c("Power.pval","Power.xy.pval") ) {

        ## estimate minimal p-value cutoff from data
        if ( missing(p.min) )
            p.min <- min(values[values>0])

    } else if ( type%in%c("Coherence") ) {
        ## values from 0 to 1
    }

    ## colorscale
    
    
}

## Wavelet coherence analysis between two counts vectors along a
## circular pseudophase.
## @param phi a circular pseudophase assignment for the counts matrix.
## @param counts a vector or matrix of counts where
##     and order as the pseudophase phi.
## @param ID1 if counts is a matrix, ID specifies the row number or
##     row name for which the analysis is to be performed.
## @param ID2 if counts is a matrix, ID specifies the row number or
##     row name for which the analysis is to be performed.
## @param verb verbosity level, set to >2 to also get verbose output
##     from the wavelet function.
## @param ... arguments to \code{\link[WavletComp]{analyze.wavelet}}.
##@export - TODO: implement for ID vector
get_coherence <- function(phi, counts, ID1, ID2, verb=0, ...) {


    clet <- analyze.coherency(df,
                              my.pair = c(iID, jID),
                              lowerPeriod = 2*dt, 
                              upperPeriod = upperPeriod, 
                              dt = dt, dj = 1/20,
                              loess.span = LOESS.SPAN,
                              make.pval = CNPERM>0,
                              method = "white.noise",
                              n.sim = CNPERM,
                              verb = verb>1)

}


#' Morlet wavelet function (via chatGPT)
#' @param t time vector in sec
#' @param central frequency in Hz
morlet_wavelets <- function(t, f0 = 6, s0 = 1, dj = 1/12, j = 0) {
    
    ## Convert scale index j to scale value
    scale <- s0 * 2^(j * dj)
  
    ## Morlet parameters
    pi <- base::pi
    sigma <- scale / (2 * pi * f0)  # Gaussian width
    norm_factor <- (pi^(-1/4)) / sqrt(sigma)
  
    ## Wavelet
    wavelet <- norm_factor *
        exp(1i * 2 * pi * f0 * t / scale) *
        exp(-t^2 / (2 * sigma^2))
    
    return(list(
        wavelet = wavelet,
        scale = scale,
        period = (4 * pi * scale) / (f0 + sqrt(2 + f0^2))
    ))
}

show_morlet <- function(t = seq(-5, 5, length.out = n), n = 1000,
                        f0 = 6, s0 = 1, dj = 1/12, j = 0,
                        xlab = "Time", ylab = "Amplitude",
                        main = "Morlet Wavelet", legend = TRUE, ...) {

    ## Generate Morlet wavelet
    w <- morlet_wavelets(t=t, f0 = f0, s0 = s0, dj = dj, j = j)

    ## Plot real and imaginary parts
    plot(t, Re(w$wavelet),
         type = "l", col = "blue",
         ylim = range(c(Re(w$wavelet), Im(w$wavelet))),
         ylab = ylab, xlab = xlab, main = main, ...)
    lines(t, Im(w$wavelet), col = "red")
    lines(t, Mod(w$wavelet), col = "black", lty = 2)
    if ( legend )
        legend("topright", legend = c("Re(ψ)", "Im(ψ)", "|ψ|"),
               col = c("blue", "red", "black"), lty = c(1, 1, 2))
    return(w)
}
