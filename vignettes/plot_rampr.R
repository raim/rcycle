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
k0 <- 100

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
                           model = mod, relative = TRUE, force.relative = FALSE)
}

plotdev(file.path(out.path, 'pwm_rampr_tau'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0,6), #range(tmns, na.rm=TRUE),        
        xlab=expression(period~tau/h), ylab=axis_labels['rampr'])
legend('topleft', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(phi==.(phi)), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_rampr_tau_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, log10(tmns), type='l', lty=1, col=1:ncol(tmns),
        axes=FALSE,
        xlab=expression(period~tau/h), ylab=axis_labels['rampr'])
logaxis(2)
axis(1)
logaxis(4, labels=FALSE)
box()
legend('topright', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(phi==.(phi)), 3, 0)
dev.off()

## amplitudes dependence on duty cycle phi
## NOTE: this also tests the version
pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
mns <- amps <- pmns
for ( mod in models ) 
    for ( i in 1:length(phis) ) {
        pmns[i,mod] <- get_ramp(gamma=gamma, phi=phis[i], tau=tau,
                               k=k, k0=k0, 
                               model = mod, relative = TRUE,
                               force.relative = TRUE)
        amps[i,mod] <- get_ramp(gamma=gamma, phi=phis[i], tau=tau,
                                k=k, k0=k0, 
                                model = mod, relative = FALSE,
                                force.relative = FALSE)
        mns[i,mod] <- get_rmean(k=k, gamma=gamma, k0=k0, phi=phis[i], tau=tau,
                                model = mod)
    }

## TODO: analyze mean vs. amplitudea along phis
## NOTE: u-shaped for model k, and
##       shifted hyperbolic for models dr
if ( FALSE ) 
    for ( mod in models ) {
        lplot(amps[,mod], mns[,mod])#, log='y')
        lplot(amps[,mod], mns[,mod], log='xy')
        Sys.sleep(1)
    }

plotdev(file.path(out.path, 'pwm_rampr_phi'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
ppmns <- pmns
##ppmns[ppmns<0] <- NA
matplot(phis, ppmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0,6), # max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=axis_labels['phi'], ylab=axis_labels['rampr'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_rampr_phi_all'),
        type='pdf', width=W, height=2.5*H, bg=NA)
par(mfrow=c(3,1), mai=c(.35,.5,.05,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, amps, type='l', lty=1, col=1:ncol(pmns),
        xlab=axis_labels['phi'], ylab=axis_labels['ramp'])
matplot(phis, mns, type='l', lty=1, col=1:ncol(pmns),
        ylim = c(0,500),
        xlab=axis_labels['phi'], ylab=axis_labels['rmean'])
matplot(phis, amps/mns, type='l', lty=1, col=1:ncol(pmns),
        xlab=axis_labels['phi'], ylab=axis_labels['rampr'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
#mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_rampr_phi_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(log10(0.03), log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=axis_labels['phi'], ylab=axis_labels['rampr'])
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
phis <- seq(0,1,.1)
phta <- matrix(NA, nrow=length(phis), ncol=length(taus))

PHT <- list()

mod <- 'k'
for ( mod in models ) {

    for ( i in 1:length(phis) ) 
        for ( j in 1:length(taus) )
            phta[i,j] <- get_ramp(gamma=gamma, phi=phis[i], tau=taus[j],
                                  k=k, k0=k0, 
                                  model = mod, relative = TRUE)

    PHT[[mod]] <- phta
    phta[is.infinite(phta)] <- NA
    plot(phis, phta[,ncol(phta)])

    plotdev(file.path(out.path, paste0('pwm_rampr_3d_', mod)),
            type='pdf', width=W, height=H)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
    image_matrix(x = taus, y = phis,
                 z = phta[nrow(phta):1,],
                 axis = 1:2, col = viridis::viridis(100),
                 xlab = axis_labels['tau'], ylab = axis_labels['phi'],
                 main = axis_labels['rampr'])
    figlabel(paste(mod), pos = 'topright', cex=1.2)
    dev.off()
    plotdev(file.path(out.path, paste0('pwm_rampr_3d_', mod, '_log')),
            type='pdf', width=W, height=H)
    par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
    image(x = phis, y = taus, z = log10(phta),
          col = viridis::viridis(100),
          ylab = axis_labels['tau'], xlab = axis_labels['phi'],
          main = axis_labels['rampr'])
    figlabel(paste(mod), pos = 'topright', cex=1.2)
    dev.off()

    z <- phta
    nrz <- nrow(z)
    ncz <- ncol(z)
    color <- viridis::viridis(100)
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    facetcol <- cut(zfacet, 100)
    plotdev(file.path(out.path, paste0('pwm_rampr_3d_', mod, '_persp')),
            type='pdf', width=W, height=H)
    par(mai=c(.1,.0,.0,.0), mgp=c(2,0.3,0), tcl=-.25)
    persp(x=phis, y=taus, z= phta, theta = -60, phi = 0,
          col=color[facetcol],
          zlab = 'rel. RNA ampl.',
          xlab = 'duty cycle',
          ylab = 'period', zlim = c(0,10))
    figlabel(mod, pos = 'topleft')
    dev.off()
    
    ## TODO: useful 3D plots?
    if ( interactive() ) {

        image(x = phis, y = taus, z = phta,
              col = viridis::viridis(100),
              ylab = axis_labels['tau'], xlab = axis_labels['phi'],
              main = axis_labels['rampr'])

        library('rgl')

        grid <- expand.grid(x = phis, y = taus)
        plot3d(x = grid$x, y = grid$y, z = phta,
               ylab = axis_labels['tau'], xlab = axis_labels['phi'],
               zlab = axis_labels['rampr'])
        surface3d(x = phis, y = taus, z = phta,
                  color = "lightblue",
                  back = "lines")
        
        persp(x=phis, y=taus, z= phta)
        contour(x=phis, y=taus, z= phta)

    }
}


## vary gamma

gammas <- seq(0, 10, length.out=100)

## absolute amplitude
rampg <- list()
for ( mod in models ) 
    rampg[[mod]] <- get_ramp(gamma=gammas, phi=phi, tau=tau,
                             k=k, k0=k0, 
                             model = mod, relative = FALSE)
rampg <- do.call(cbind, rampg)
plotdev(file.path(out.path, 'pwm_gamma_ramp'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(gammas, rampg, type ='l', lty=1, col=1:ncol(rampg),
        xlab = axis_labels['gamma'], ylab = axis_labels['ramp'])
legend('right',
       colnames(rampg), col=1:ncol(rampg), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
dev.off()

## relative amplitude
ramprg <- list()
for ( mod in models ) 
    ramprg[[mod]] <- get_ramp(gamma=gammas, phi=phi, tau=tau,
                             k=k, k0=k0, 
                             model = mod, relative = TRUE)
ramprg <- do.call(cbind, ramprg)
plotdev(file.path(out.path, 'pwm_gamma_rampr'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(gammas, ramprg, type ='l', lty=1, col=1:ncol(ramprg),
        xlab = axis_labels['gamma'], ylab = axis_labels['rampr'])
legend('bottomright',
       colnames(ramprg), col=1:ncol(ramprg), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
dev.off()

## mean abundance
rmeang <- list()
for ( mod in models ) 
    rmeang[[mod]] <- get_rmean(gamma=gammas, phi=phi, tau=tau,
                               k=k, k0=k0, 
                               model = mod)
rmeang <- do.call(cbind, rmeang)
plotdev(file.path(out.path, 'pwm_gamma_rmean'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(gammas, log10(rmeang), type ='l', lty=1, col=1:ncol(rmeang),
        xlab = axis_labels['gamma'], ylab = axis_labels['rmean'],
        axes = FALSE)
logaxis(2)
axis(1)
box()
legend('topright',
       colnames(rmeang), col=1:ncol(rmeang), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
dev.off()

## vary k and calculate ramp and rmean
ks <- seq(0, 1000, length.out=100)

## absolute amplitude
rampg <- list()
for ( mod in models ) 
    rampg[[mod]] <- get_ramp(gamma=gamma, phi=phi, tau=tau,
                             k=ks, k0=k0, 
                             model = mod, relative = FALSE)
rmeang <- list()
for ( mod in models ) 
    rmeang[[mod]] <- get_rmean(gamma=gamma, phi=phi, tau=tau,
                               k=ks, k0=k0, 
                               model = mod)
rmeang <- do.call(cbind, rmeang)
rampg <- do.call(cbind, rampg)

plotdev(file.path(out.path, 'pwm_k_ramp'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.2,0.3,0), tcl=-.25)
matplot(ks, (rampg), type ='l', lty=1, col=1:ncol(rampg),
        xlab = axis_labels['k'], ylab = axis_labels['ramp'],
        axes = FALSE)
axis(2)
axis(1)
box()
legend('topleft',
       colnames(rampg), col=1:ncol(rampg), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
dev.off()


plotdev(file.path(out.path, 'pwm_k_rmean'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.2,0.3,0), tcl=-.25)
matplot(ks, (rmeang), type ='l', lty=1, col=1:ncol(rmeang),
        xlab = axis_labels['k'], ylab = axis_labels['rmean'],
        axes = FALSE)
axis(2)
axis(1)
box()
legend('topleft',
       colnames(rmeang), col=1:ncol(rmeang), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_k_rmean_ramp'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.2,0.3,0), tcl=-.25)
matplot(log10(rmeang), log10(rampg), log='', type ='l', lty=1, col=1:ncol(rmeang),
        xlab = axis_labels['rmean'], ylab = axis_labels['ramp'],
        axes = FALSE)
logaxis(1)
logaxis(2)
box()
arrows(y0=log10(20), y1=log10(100),
       x0=log10(50), x1=log10(300), length=.05, col=8, lwd=2, lty=2)
text(log10(200), log10(30), 'k', col=8, font=2)
legend('topleft',
       colnames(rmeang), col=1:ncol(rmeang), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h*';'~phi==.(phi)), 3, 0)
##abline(a=0, b=1, lty=2, col=8)
dev.off()
