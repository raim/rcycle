
library(rcycle) # this library
library(segmenTools) # for num2col only
library(RSpectra) # for eigs_sym
load.rcycle <- function() {
    srcs <- list.files(path='~/programs/rcycle/R', pattern='*.R$',
                       full.names = TRUE)
    for ( src in srcs )
        source(src)
}
load.rcycle()

## Compares SVD/PCA on simulated data (sine waves, with noise, just noise)

## 1. generate data genes X samples,
## 2. norm by total counts (single cells),
## 3. row-centering or standardization:
##     - analyze effect of row normalization (genes across samples),
##       followed by column scaling (within sample across genes).
##     - compare to normalization problem as e.g. in @Lehmann2013,
## 4. "Eigengenes": do not transpose
##     - raw data: svd, eigen(cov()), prcomp(scale.=FALSE, center=FALSE)
##     - row-centered: same same,
##     - column centering: svd, eigen(cor), prcom(scale.=FALSE)
##     - column scaling: svd, eigen(cor), prcom(scale.=FALSE)
## 5. single cell: transpose and do SVD=PCA,

## TODO: get sine wave simulation from cohorts_simulate.R and
## PWM simulation from models.R and pwmode.R in ChemostatData/src/


ftyp <- "png"
H <- W <- 3
res <- 200


tau <- 2
ncells <- 100
ncycles <- 1

## cohorts
NY <- 5
cls <- paste0("c",1:NY)
# cohort colors
cls.col <- setNames(1:NY, cls)


## time vector
xtime <- seq(0,ncycles*tau, length.out=ncells)
tcol <- num2col(xtime) 


## STATE MATRIX
ys <- matrix(NA, nrow=NY, ncol=length(xtime))
rownames(ys) <- cls

## mean count, amplitude and phase shift
mcnt <- 300 # mean level of transcripts
ramp <- 0.1 # relative amplitude
shft <- 1/NY # phase shift

## uniform mean counts, amplitudes and phase shifts
mcnts <- rep(mcnt, NY)
ramps <- rep(ramp, NY)
shfts <- (1:NY -1)*shft

## loop over standard deviations for abundance and amplitudes,
## with and without variable phase shifts

## generate all combinations for options

## simulated data options: from regular sine waves to realistic oscillators
sims <- expand.grid(amp=0:1, phs=0:1, mean=0:1)
## data processing options
proc <- expand.grid(cols=0:1, rows=0:1) #, tt=0:1)

for ( v in 1:nrow(sims)  ) {

    ## simulated time series
    mean <- sims[v,"mean"] # randomize mean abundances
    phs <- sims[v,"phs"] # randomize shift phases
    amp <- sims[v,"amp"] # randomize amplitudes

 
    rtyp <- paste(paste0(colnames(sims), sims[v,]),collapse="_")
    
    set.seed(2)

    ## variable counts, amplitudes and phase shifts
    ##shfts <- sort(minmax(rnorm(NY, mean=shft, sd=sd)))
    ##
    shfts <- (1:NY -1)*1/NY
    if ( phs==1 ) ## PHASE SHIFTS
      shfts <- sort(sample(1:100, NY)/100)  
    mcnts <- rep(mcnt, NY)
    if ( mean==1) ## MEANS - ONLY THING THAT REQUIRES ROW CENTERING
        mcnts <- sample(10:5000, NY)
    ramps <- rep(ramp, NY)
    if ( amp==1 ) ## AMPLITUDES
        ramps <- sort(sample(1:50, NY)/100)  


    ## sine wave of transcript levels
    for ( i in 1:NY ) {
    
        ## regular phase shifts between states
        phase <- tau * shfts[i]
        
        ## calculate osci expression
        ## TODO: use PWM model instead of sinus
        ys[i,] <- mcnts[i]*( 1 + ramps[i] * sin(2*pi/(tau) * xtime
                                                - 2*pi*phase/tau))
    }

    
    ## PCA doesn't care about order, can be confirmed by sampling columns
    ## but NOTE, currently sampling breaks the plots below
    ys <- ys#[,sample(1:ncol(ys))]

    ## TEST NORM. STEPS

    ## raw simulated data
    ##  total normalized
    yn <- normalize_counts(ys)
    ##  row-norm
    yr <- yn - apply(yn, 1, mean)
    ## column scaling as for prcomp(scale.=TRUE, center=TRUE)
    yc <- apply(yr, 2, scale)

    par(mfcol=c(1,4), mai=c(.5,.5,.05,.05), mgp=c(1.3,.3,0), tcl=-.25)
    matplot(xtime, t(ys), type="l", lty=1)
    legend("topright", rtyp)
    matplot(xtime, t(yn), type="l", lty=1)
    legend("topright", "fraction of total")
    matplot(xtime, t(yr), type="l", lty=1)
    legend("topright", "transcript centering")
    matplot(xtime, t(yc), type="l", lty=1)
    legend("topright", "sample scaling")
        
    if ( FALSE ) {
    
        plot(ys[1,], yc[1,], col=NA, xlim=range(c(ys)), ylim=range(c(yc)),
             xlab="simulated data", "normalized data")
        for ( j in 1:nrow(ys) )
            points(ys[j,], yc[j,], col=j, type='l')
        legend("topright", rtyp)
    }

    for ( p in 1:nrow(proc)  ) {

        total <- 1 #proc[p,"total"] # randomize amplitudes
        rows <- proc[p,"rows"] # row-center before scaling
        cols <- proc[p,"cols"] # column scaling

        styp <- paste(paste0(colnames(proc), proc[p,]),collapse="_")

        ## APPLY SELECTED PROCESSING STEPS
        ## normalize each time point to total counts
        ## NOTE: LESS WARPING w/o TOTAL COUNT NORM
        yn <- ys

        ## fraction of total count/intensity
        if ( total==1 )
            yn <- normalize_counts(ys)
        ## row-centering! see @Lehmann2013
        if ( rows==1 )
            yn <- yn - apply(yn, 1, mean)
        ## column scaling as for prcomp(scale.=TRUE, center=TRUE)
        if ( cols==1 )
            yn <- apply(yn, 2, scale)

        par(mfcol=c(1,5), mai=c(.5,.5,.05,.05), mgp=c(1.3,.3,0), tcl=-.25)
        matplot(xtime, t(ys), type="l", lty=1, xlab='time/h')
        legend("topright", c(rtyp))
        matplot(xtime, t(yn), type="l", lty=1, xlab='time/h')
        legend("topright", c(styp))

        plot(ys[1,], yc[1,], col=NA, xlim=range(c(ys)), ylim=range(c(yc)),
             xlab="simulated data", ylab="normalized data")
        for ( j in 1:nrow(ys) )
            points(ys[j,], yc[j,], col=j, type='l')
          
        pca <- prcomp(yn, scale.=FALSE, center=FALSE)
        plotPC(pca, sarrows=TRUE, vlines=TRUE, scol=1:NY, slwd=2, scale=1,
               zero.axis=TRUE, saxis=FALSE, vaxis=FALSE)
        legend("bottomright", c(styp))

        plot(xtime, atan2(pca$rotation[,2], pca$rotation[,1]),
             xlab='time/h', ylab=expression(atan2(PC2,PC1)), axes=FALSE)
        axis(1);circ.axis(2)
    }
}
