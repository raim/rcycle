library(rcycle)


#. temp requirements
source('/home/raim/programs/rcycle/R/models.R')

## Testin root finding functions

## TODO:
## * write unit tests,
## * understand all root finding functions, and find
##   best solution, ideally just stats::uniroot,
##   or informed selection from uniroot.all.
## * trace NA cases?
## * better understand difficulties for get_tau for model k.

a <- 1.32991318770238
gamma <- 0.411938922118397
phi <- 0.5
A <- 2.44109012400914e-05
k <- 3.5480961104784e-05

R <- A/a
Rmin <- A - A/2

## get_tau roots

x <- seq(-10,10,.01)
tau_k <- root_tau_k(x, a = a, gamma = gamma, phi = phi, A = NA, k = NA)
roots_k <- rootSolve::uniroot.all(f=root_tau_k,
                                  a = a, gamma = gamma, phi = phi,
                                  lower = min(x), upper = max(x))

plot(x, tau_k, type='l', xlab='period/h', ylab='root_tau_k_dr')
abline(h = 0, col=2, lty=2)
abline(v = roots_k, col=2, lty=2)

tau_k_dr <- root_tau_k_dr(x, a = a, gamma = gamma, phi = NA, A = A, k = k)
roots_k_dr <- rootSolve::uniroot.all(f=root_tau_k_dr,
                                     a = a, gamma = gamma, A = A, k = k,
                                     lower = min(x), upper = max(x))

plot(x, tau_k_dr, type='l', xlab='period/h', ylab='root_tau_k_dr',
     ylim = c(-10,10))
abline(h = 0, col=2, lty=2)
abline(v = roots_k_dr, col=2, lty=2)



## get_degradation roots

tau <- seq(-100,100,.01)
gt <- gamma*tau

for ( model in c('k', 'dr', 'k_dr', 'k_dr_k0') ) {

    rootf <- get(paste0('root_', model), mode = 'function')
    y <- rootf(gt, a = a, phi = .1, R = R, Rmin = Rmin)

    plot(gt, y, type='l', xlab=expression(gamma*tau), ylab='root_k')
    abline(h = 0, col=2, lty=2)
    ##Sys.sleep(2)
}


gt_k <- root_k(gt, a = a, phi = phi)
gt_k_dr <- root_k_dr(gt, a = a, phi = phi)
gt_dr <- root_dr(gt, a = a, phi = .1)
gt_k_dr_k0 <- root_k_dr_k0(gt, a = a, phi = .1, R = R, Rmin = Rmin)

plot(gt, gt_k, type='l', xlab=expression(gamma*tau), ylab='root_k')
abline(h = 0, col=2, lty=2)

plot(gt, gt_k_dr, type='l', xlab=expression(gamma*tau), ylab='root_k_dr',
     ylim=c(-10,10))
abline(h = 0, col=2, lty=2)

plot(gt, gt_dr, type='l', xlab=expression(gamma*tau), ylab='root_dr',
     ylim=c(-10,10))
abline(h = 0, col=2, lty=2)

plot(gt, gt_k_dr_k0, type='l', xlab=expression(gamma*tau), ylab='root_k_dr_k0')
abline(h = 0, col=2, lty=2)

source('/home/raim/programs/rcycle/R/models.R')

gt_k_dr_k0 <- root_k_dr_k0(gt, a = a, phi = .1, R = R, Rmin = Rmin)
gt_k_dr_k0_old <- root_k_dr_k0_old(gt, a = a, phi = .1, R = R, Rmin = Rmin)
gt_k_dr_k0_coth <- root_k_dr_k0_coth(gt, a = a, phi = .1, R = R, Rmin = Rmin)

## relation of old the new root finding function: a*R
plot(gt_k_dr_k0, gt_k_dr_k0_old/(a*R))
abline(a=0, b=1)

plot(gt_k_dr_k0, gt_k_dr_k0_coth)
abline(a=0, b=1)

roots_k_dr_k0 <- rootSolve::uniroot.all(f=root_k_dr_k0,
                                        a = a, phi = .1, R = R, Rmin = Rmin,
                                        lower = min(x), upper = max(x))
roots_k_dr_k0_old <- rootSolve::uniroot.all(f=root_k_dr_k0_old,
                                            a = a, phi = .1, R = R, Rmin = Rmin,
                                            lower = min(x), upper = max(x))
roots_k_dr_k0_coth <- rootSolve::uniroot.all(f=root_k_dr_k0_coth,
                                             a = a, phi = .1, R = R, Rmin = Rmin,
                                             lower = min(x), upper = max(x))

plot(gt, gt_k_dr_k0, type='l', xlab=expression(gamma*tau), ylab='root_k',
     ylim = c(-1,10))
lines(gt, gt_k_dr_k0_old/(a*R), col =2, lty=3 , lwd=3)
abline(h=0)
abline(v = roots_k_dr_k0_old, col=4, lty=1, lwd=4)
abline(v = roots_k_dr_k0_coth, col=3, lty=3, lwd=3)
abline(v = roots_k_dr_k0, col=2, lty=3, lwd=2)
