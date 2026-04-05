
library(rcycle)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 2*1.7
mu <- 0
gamma <- dr+mu
k0 <- 0

## osci and growth params
phi <- .25
tau <- 5

## Simulate PWM with the analytical model
## TODO:
## * simulate all with ODE models, and compare,
## * use analytic vs. ODE for tests!

mod='k'

ravg <- get_rmean(k=k, gamma=gamma, k0=k0, phi=phi, tau=tau,
                  model = mod, use.coth = TRUE)

time <- seq(0,5*tau,.01)
y1 <- pwm_k(t=time, R0=ravg, k=k, dr=dr, k0=k0, mu=mu, phi=1-phi, tau=tau,
            theta=pi)$R
y2 <- pwm_k(t=time, R0=ravg, k=k, dr=dr, k0=k0, mu=mu, phi=phi, tau=tau)$R

## add noise:
n1 <- jitter(y1, 100)
n2 <- jitter(y2, 100)

plot(time, n1, type='l', ylim=c(0,80), ylab=axis_labels["r"])
lines(time, n2, col=2)

## NOTE: anti-correlation in bins: always -1?
bins <- seq(min(time), max(time), length.out=100) # resolution
bintme <- bincor <- c()
for ( k in 2:length(bins) ) {
    idx <- time>=bins[k-1] & time<bins[k]
    cr <- cor.test(n1[idx], n2[idx])
    bincor[k-1] <- cr$estimate
    bintme[k-1] <- mean(time[idx])
    if ( cr$estimate> -.9) {
        lines(time[idx], y1[idx], col=1, lwd=5)
        lines(time[idx], y2[idx], col=2, lwd=5)
    }
        
}

## correlation becomes about 0 in steady state
## regions, ONLY when noise is added.
par(mfcol=c(2,1))
plot(time, n1, type='l', ylim=c(0,80), ylab=axis_labels["r"], xlab='time/h')
lines(time, n2, col=2)
plot(bintme, bincor, type='l', xlab='time/h', ylab="binned correlation")
abline(h=0)
