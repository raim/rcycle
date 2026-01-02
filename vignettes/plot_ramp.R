#!/usr/bin/Rscript --vanilla

library(rcycle)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'


## TODO:
## * move plots to vignette,
## * generate more tests of model results,
## * add transcriptome osci conditions!

W <- H <- 2.5

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 2*1.7
mu <- 0
gamma <- dr+mu
k0 <- 10

## osci and growth params
phi <- .5
tau <- 2

## vary over tau and phi
taus <- seq(0,7.5, .1)
phis <- 0:100/100

## amplitudes dependence on period tau
models <- c('k', 'k_dr', 'dr', 'k_dr_k0')
tmns <- matrix(NA, nrow=length(taus), ncol=length(models))
colnames(tmns) <- models
for ( mod in models ) {
    tmns[,mod] <- get_ramp(gamma=gamma, phi=phi, tau=taus,
                           k=k, k0=k0, # only required with basal expression!
                           model = mod, relative = FALSE)
}

plotdev(file.path(out.path, 'pwm_ramp_tau'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, tmns, type='l', lty=1:ncol(tmns), col=1:ncol(tmns),
        ##ylim=c(0,3e3), #range(tmns, na.rm=TRUE),        
        xlab=expression(period~tau/h), ylab=axis_labels['ramp'])
legend('topleft', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(phi==.(phi)), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_ramp_tau_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, log10(tmns), type='l', lty=1, col=1:ncol(tmns),
        axes=FALSE,
        xlab=expression(period~tau/h), ylab=axis_labels['ramp'])
logaxis(2)
axis(1)
logaxis(4, labels=FALSE)
box()
legend('topright', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(phi==.(phi)), 3, 0)
dev.off()

## amplitudes dependence on duty cycle phi
pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
for ( mod in models ) {
    pmns[,mod] <- get_ramp(gamma=gamma, phi=phis, tau=tau,
                           k=k, k0=k0, # only required with basal expression!
                           model = mod, relative=FALSE)
}

plotdev(file.path(out.path, 'pwm_ramp_phi'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
ppmns <- pmns
ppmns[ppmns<0] <- NA
matplot(phis, ppmns, type='l', lty=1:ncol(pmns), col=1:ncol(pmns),
        ##ylim=c(0,1e3), # max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['ramp'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_ramp_phi_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(log10(0.03), log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['ramp'])
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
box()
legend('bottomleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()


### 3D - CALCULATE ALL PHI/TAU COMBINATIONS

taus <- seq(0,7.5, .5)
phis <- 1:10/10
phta <- matrix(NA, nrow=length(phis), ncol=length(taus))

PHT <- list()

mod <- 'k'
for ( mod in models ) {

    for ( i in 1:length(phis) ) 
        for ( j in 1:length(taus) )
            phta[i,j] <- get_ramp(gamma=gamma, phi=phis[i], tau=taus[j],
                                  k=k, k0=k0, 
                                  model = mod, relative=FALSE)

    PHT[[mod]] <- phta
    phta[is.infinite(phta)] <- NA
    plot(phis, phta[,ncol(phta)])

    plotdev(file.path(out.path, paste0('pwm_ramp_3d_', mod)),
            type='pdf', width=W, height=H)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
    image_matrix(x = taus, y = phis,
                 z = phta[nrow(phta):1,],
                 axis = 1:2, col = viridis::viridis(100),
                 xlab = axis_labels['tau'], ylab = axis_labels['phi'],
                 main = axis_labels['ramp'])
    figlabel(paste(mod), pos = 'topright', cex=1.2)
    dev.off()
    plotdev(file.path(out.path, paste0('pwm_ramp_3d_', mod, '_log')),
            type='pdf', width=W, height=H)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
    image(x = phis, y = taus, z = log10(phta),
          col = viridis::viridis(100),
          ylab = axis_labels['tau'], xlab = axis_labels['phi'],
          main = axis_labels['ramp'])
    figlabel(paste(mod), pos = 'topright', cex=1.2)
    dev.off()

    z <- phta
    nrz <- nrow(z)
    ncz <- ncol(z)
    color <- viridis::viridis(100)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, 100)
    plotdev(file.path(out.path, paste0('pwm_ramp_3d_', mod, '_persp')),
            type='pdf', width=W, height=H)
    par(mai=c(.0,.0,.0,.0), mgp=c(2,0.3,0), tcl=-.25)
    persp(x=phis, y=taus, z= phta, theta = -60, phi = 0,
          col=color[facetcol],
          zlab = 'RNA amplitude',
          xlab = 'duty cycle',
          ylab = 'period',
          ticktype = 'simple')
    ## TODO: draw ticks for z-axis!
    figlabel(mod, pos = 'topleft')
    dev.off()
    
    ## TODO: useful 3D plots?
    if ( interactive() ) {

   
        library('rgl')

        grid <- expand.grid(x = phis, y = taus)
        plot3d(x = grid$x, y = grid$y, z = phta,
               ylab = axis_labels['tau'], xlab = axis_labels['phi'],
               zlab = axis_labels['ramp'])
        surface3d(x = phis, y = taus, z = phta,
                  color = "lightblue",
                  back = "lines")
        
        persp(x=phis, y=taus, z= phta)
        contour(x=phis, y=taus, z= phta)

    }
}
