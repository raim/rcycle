#!/usr/bin/Rscript --vanilla

library(rcycle)


#. temp requirements
source('/home/raim/programs/rcycle/R/models.R')
minmax <- function(x, ...) (x-min(x, ...))/(max(x, ...)-min(x, ...))
library(segmenTools)
library(deSolve)
out.path <- '/home/raim/programs/rcycle/vignettes'

## TODO:
## * move to vignettes and instead test pw_fourier, pw_simple, pw_sawtooth,


## osci and growth params
stime <- seq(0, 3, .01)
phocs <- c(.01, seq(.2, .8, .2), .99)
tosc <- 1
mu <- 0

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 1.7
gamma <- dr+mu
k0 <- 10


plotdev(file.path(out.path, 'pwm_cartoon'), type='pdf', width=5, height=2)
par(mfcol=c(1,2), mai=c(.35,.5,.35,.1), mgp=c(.5,0,0))
par(mai=c(.35,.5,.35,.4))
plot(1, col=NA, xlim=range(stime), ylim=c(0, length(phocs)), axes=FALSE,
     xlab='time', ylab=NA)
text(mean(stime), length(phocs)+1,
     labels='Alternating pulse waves\nof gene expression:', font=2, xpd=TRUE)
mtext('growth rate', 2, .6)
arrows(x0=-.2, y0=.5, y1=length(phocs)-.5, xpd=TRUE, length=.1)
arrows(x0=.5, x1=max(stime)-.5, y0=-.7, xpd=TRUE, length=.1)
for ( i in seq_along(phocs) ) {
    phoc <- phocs[i]
    ploc <- 1 - phoc
    
    hocon <- as.numeric(pw_fourier(t=stime, k=1,
                                   phi=phoc, tau=tosc, theta=0*pi,
                                   N=1e3)>.5)
    locon <- as.numeric(pw_fourier(t=stime, k=1,
                                   phi=ploc, tau=tosc, theta=pi,
                                   N=1e3)>.5)
    
    lines(stime, .8*hocon+i-1, type='l', col=2)
    lines(stime, .8*locon+i-1, col=4, lty=3)
}
text(x=rep(max(stime), length(phocs)), y=seq_along(phocs)-.5,
     labels=phocs, col=2, xpd=TRUE, pos=4)
##text(x=max(stime), y=-.75, labels=expression(phi), col=2, xpd=TRUE, cex=1.2,
##     pos=4)
text(x=max(stime)+1, y=-1.5, labels=expression(atop(duty), cycle~varphi),
     col=2, xpd=TRUE, cex=1, pos=2)

par(mai=c(.35,.35,.35,.1))
plot(phocs, col=NA, xlim=c(0, length(phocs)), ylim=c(0, 1), axes=FALSE,
     ylab=NA, xlab='growth rate')
mtext('proteome fraction', 2, 0)
arrows(x0=0, y0=.05, y1=.95, xpd=TRUE, length=.1)
arrows(x0=.5, x1=length(phocs)-.5, y0=-.1, xpd=TRUE, length=.1)
lines(phocs, col=2, type='b')
lines(1-phocs, col=4, type='b', lty=3)
text(length(phocs)/2, 1.15,
     labels='Continuous output:', font=2, xpd=TRUE)
dev.off()


###  CARTOON FOR PIECEWISE ODE SOLUTIONS

state <- c(R=0)
tosc <- 3
thoc <- 1.5
mu <- .1

## NOTE: parameters for average RP from pwmavg.R/pwmode.R
params <- c(k=263.9,
            k0=100,
            dr=1.7,
            nb=3,
            dp=0.069,
            ell=250,
            mu=mu)
times <- seq(0,50*tosc, tosc/200)
## start pulse at time 0
pulses <- pwm_simple(time=times, tau=tosc, phi=thoc/tosc, theta = 0)
hocf <- approxfun(x = times, y = pulses,
                  method = "constant", rule = 2) # HOC pulse

models <- c('k', 'k_dr', 'dr', 'k_dr_k0')
outl <- list()
for ( mod in models ) {

    state <- c(R=get_rmean(dr=params['dr'],
                           k=params['k'],
                           k0=params['k0'],
                           phi=thoc/tosc,
                           tau=tosc,
                           mu=params['mu'], model = mod))
    outl[[mod]] <- ode(y = state,
                       times = times,
                       func = get(paste0('pwmode_', mod), mode = 'function'),
                       parms = params, hocf = hocf)
}
## CUT ENDS for proper relative ylims
outl <- lapply(outl, function(x)  {
    x[x[,1]<10*tosc, 2] <- mean(x[,2]); x})

## show pulse
pulse <- cbind(time=times, kappa=pulses)


plotdev(file.path(out.path, 'pwm_cartoon_piecewise'),
        type='pdf', width=3, height=2)
par(mai=c(.35,.35,.1,.35), mgp=c(1.4,0.3,0), tcl=-.25)
plot(0, ylim=c(0,1), col=NA,
     xlim=c(45*tosc-.05*tosc,
            46*tosc+.05*tosc),
     type='l', axes=FALSE,
     xlab='', ylab='')
for ( i in seq_along(models) ) {
    dat <- outl[[models[i]]]
    lines(dat[,1], minmax(dat[,2]), col=i, lty=i, lwd=.5)
    mn <- mean(minmax(dat[,2]))
    abline(h=mn, col=i, lty=i, lwd=.1)
    if (i==1) 
        axis(4, at=mn, label=axis_labels['rmeanau'], las=2, col=i)
}
axis(2, at=0:1, labels=expression(R[0], R[1]), las=2, col=NA)
arrows(x0=par('usr')[1]-.2, y0=0.1, y1=.9, code=3, length=.1, xpd=TRUE)
axis(4, at=0:1, labels=expression(R[tau], R[1]), las=2, col=NA)
mtext(axis_labels['ramp'], 2, .6)
axis(1, at=45*tosc, labels=0, col=NA) 
axis(1, at=45*tosc+thoc, labels=expression(tau[on]), col=NA) 
axis(1, at=46*tosc, labels=expression(tau), col=NA)
abline(h=0:1, lty=3, col=8)
abline(v=seq(0:100)*tosc, lty=3, col=8)
abline(v=seq(0:100)*tosc+thoc, lty=3, col=8)

arrows(y0=0, x0=45*tosc, x1=45*tosc+thoc,
       code=3, length=.1, col=8, lwd=1.5)
arrows(y0=0, x0=45*tosc+thoc, x1=46*tosc,
       code=3, length=.1, col=8, lwd=1.5)
if ( FALSE )
    legend(x=45*tosc+1.1*thoc, y=1, names(outl),
           col=1:length(outl), lty=1:length(outl), bg='#ffffff77', box.col=NA,
           seg.len=.75, y.intersp=.75, inset=c(-.25,0), xpd=TRUE)         
dev.off()


plotdev(file.path(out.path, 'pwm_cartoon_piecewise_abs'),
        type='pdf', width=3, height=2)
par(mai=c(.5,.75,.1,.35), mgp=c(1.4,0.3,0), tcl=-.25)
ylm <- range(unlist(lapply(outl, function(x) range(x[,2], na.rm=TRUE))))
plot(0, col = NA, ylim=ylm,
     xlim=c(45*tosc-.05*tosc,
            46*tosc+.05*tosc),
     type='l', axes=FALSE,
     xlab='time/h', ylab='R(t)/(n/cell)')
for ( i in seq_along(models) ) {
    dat <- outl[[models[i]]]
    lines(dat[,1], dat[,2], col=i, lty=i)
}
axis(1)
axis(2)
abline(v=seq(0:100)*tosc, lty=3, col=8)
abline(v=seq(0:100)*tosc+thoc, lty=3, col=8)

legend('topright', names(outl),
       col=1:length(outl), lty=1:length(outl), bg='#ffffff77', box.col=NA,
       seg.len=.75, y.intersp=.75, inset=c(-.2,0), xpd=TRUE)         
dev.off()
