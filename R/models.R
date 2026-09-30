
## MODELS for PERIODIC GENE EXPRESSION

## PULSE WAVE GENERATING FUNCTIONS

pwm_simple <-  function(time, tau, phi, theta = 0) {
  # Returns integer 0/1 vector: 1 when pulse is on.
  # theta in radians; tau period in same time units as time.
  omega <- 2 * pi / tau
  # convert phase shift theta (radians) to time shift
  t_shifted <- (time - theta / omega) %% tau
  as.integer(t_shifted < (phi * tau))
}


## pulse wave function - direct and very SLOW calculation
## wikipedia: "Note that, for symmetry, the starting time (t = 0)
## in this expansion is halfway through the first pulse."
## TODO: why do we need to scale?
#' Pulse wave, Fourier series.
#'
#' A pulse wave between 0 and \code{k} with duty cycle \code{phi} and period
#' \code{tau}, summed directly from its Fourier series (slow). The series
#' starts halfway through the first pulse (\code{t = 0} is the pulse centre).
#' @param t time points.
#' @param k height of the pulse.
#' @param phi duty cycle, the fraction of the period in the ON phase.
#' @param tau period.
#' @param N number of harmonics summed.
#' @param shift phase shift, in units of time.
#' @param theta phase shift, in radians.
#' @param start.on shift time by half a pulse, to start with the ON phase.
#' @param alpha smoothing, by exponentially damping higher harmonics:
#' 0, no smoothing; up to 0.5, mild; above 0.5, noticeable; above 2, towards
#' a sine.
#' @return numeric vector of the pulse wave at \code{t}.
#' @examples
#' \dontrun{
#' t <- seq(0,10, length=100)
#' plot(t, pw_fourier(t=t, tau=1.5, k=1, phi=.4), type='l', col=2)
#' lines(t, pw_fourier(t=t, tau=1.5, k=1, phi=.6, theta=pi), type='l', col=4)
#' }
#' @keywords internal
pw_fourier <- function(t=0, k, phi, tau, N=1e4, shift=0, theta=0,
                       start.on=FALSE, alpha=0) {

    ## shift to start with on phase
    ## NOTE: to start at on we'd shift in the cos function:
    ## cos(n*omega*t - n*pi*phi)
    if ( start.on )
        t <- t - phi*tau/2
 
    
    ## shift time to get phase shift (passed in unit of time!)
    t <- t - shift

    ## sum expression
    ## TODO: faster vectorization instead of loop here but we have
    ## both n and t vectors?
    omega <- 2*pi/tau
    sincos <- 0
    for ( n in 1:N )  
        sincos <- sincos +
            1/n * sin(pi*n*phi) * cos(n*omega*t - n*theta)* exp(-alpha*n)
    x <- 2/pi * sincos
    
    ## shift to 0 and k
    x <- k*(phi + x) 
    x
}

## pulse wave via sawtooth: subtract phase-shifted
## NOTE: sawtooth starts with off phase, but we can shift time
## here to ensure it starts with an on phase
## TODO: different formulation to avoid min/max scaling?
pw_sawtooth <- function(t=0, k=1, tau, thoc=.5, alpha=0, shift=0,
                        start.on=FALSE, start.like.pw=TRUE, scale=TRUE) {

    if ( start.like.pw ) # shift to start like pwf with half pulse!
        t <- t - thoc/2
    else if ( start.on ) # shift to start with on phase
        t <- t + (tau-thoc)
        
    ## shift time to get phase shift (passed in unit of time!)
    t <- t - shift
    
    st <-  atan(pracma::cot(pi* t      /tau)) 
    st2 <- atan(pracma::cot(pi*(t+thoc)/tau))

    ## TODO: implement smoothing of saw tooth function!
    if ( alpha!=0 ) {
        stop('smoothing not implemented yet for saw tooth function')
        sig2 <- (t+thoc)/alpha
        sig2 <- 1/(1-exp(-sig2))
        sig1 <- t/alpha
        sig1 <- 1/(1-exp(-sig1))
        
        st2 <- st2*sig2
        st <- st*sig1
        ##return(sig)
    } 
    
    x <- 1/pi * (st2 - st)
    
    ## scale between 0 and k
    ## TODO: use different formula to avoid min/max scaling?
    if ( scale )
        x <- k*(x-min(x))/(max(x)-min(x)) # x <- k*thoc/tau + x
    x
}

### ODE Models


#' ODE of pulse width-modulated transcription.
#'
#' Right-hand side of \code{dR/dt = kappa*k - (dr+mu)*R}, for
#' \code{\link[deSolve]{ode}}, with \code{kappa} the ON (1) or OFF (0)
#' state at \code{time}.
#' @param time time point.
#' @param state named state vector, \code{c(R = ...)}.
#' @param parameters named parameter vector: \code{k}, \code{dr}, \code{mu};
#' and \code{k0} for \code{pwmode_k_dr_k0}.
#' @param hocf function of time returning the pulse state \code{kappa}, 0 or
#' 1 (\emph{not} scaled by \code{k}).
#' @return list with the derivative \code{dR}, as required by
#' \code{\link[deSolve]{ode}}.
#' @seealso \code{\link{get_rcycle}} for the periodic steady state.
#' @export
pwmode_k <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        ## NOTE: compared to older functions, here
        ## kappa is 0 or 1!
        dR = kappa*k - (dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated degradation.
#'
#' \code{dR/dt = k - ((1-kappa)*dr + mu)*R}: constant transcription,
#' degradation in the OFF phase only.
#' @inheritParams pwmode_k
#' @export
pwmode_dr <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = k - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated transcription and anti-phasic
#' degradation.
#'
#' \code{dR/dt = kappa*k - ((1-kappa)*dr + mu)*R}.
#' @inheritParams pwmode_k
#' @export
pwmode_k_dr <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = kappa*k - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated transcription and anti-phasic degradation
#' and basal transcription.
#'
#' \code{dR/dt = kappa*k + (1-kappa)*k0 - ((1-kappa)*dr + mu)*R}.
#' @inheritParams pwmode_k
#' @export
pwmode_k_dr_k0 <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = kappa*k + (1-kappa)*k0 - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}


### ANALYTIC


#' Mean abundance, closed form without dilution in the ON phase.
#'
#' The cycle mean of the pulse width-modulated models without any loss in
#' the ON phase (\code{gamma = dr + mu} acts in the OFF phase only); exact for
#' model \code{"k"} and, for all models, for \code{mu = 0}. The closed forms
#' of the slides (\code{pwm_equ.md}); use \code{\link{get_rmean}} for the
#' exact mean with dilution in both phases.
#' @inheritParams get_rcycle
#' @param gamma total loss rate, \code{dr + mu}; if missing, calculated from
#' \code{dr} and \code{mu}.
#' @param mu growth rate, used only to calculate a missing \code{gamma}.
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"} (or \code{"k_dr_k0_coth"}); the first is used.
#' @param use.coth evaluate the OFF-phase term via \code{coth}, otherwise via
#' \code{expm1}.
#' @return numeric vector of cycle means.
#'@export
get_rmean_nogrowth <- function(k, gamma, k0, dr, mu, phi, tau,
                      model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      use.coth = TRUE) {

    if ( length(model)>1 ) model <- model[1]

    if ( missing(gamma) )
        gamma <- dr+mu

    ## base model: phi*k/gamma
    rmn <- phi*k/gamma


    ## expand, if the model is requested for a vector of periods or duty cycles
    if ( !missing(tau) )
        if ( length(tau) > length(rmn) ) 
            rmn <- rep(rmn, length(tau))

    ## add terms for all other models
    ## phi^2*k*tau/2 * coth(gamma*tau*(1-phi)/2)
    if ( model %in% c('dr', 'k_dr', 'k_dr_k0', 'k_dr_k0_coth') ) {

        ## exponent
        y <- gamma*tau*(1-phi)
        if ( use.coth ) 
            term2 <- pracma::coth(y/2)
        else {
            ye <- expm1(y) # exp(y)-1
            term2 <- (ye+2)/ye
        }

        rmn <- rmn + (phi^2*k*tau/2)*term2
    }
    if ( model %in% c('dr') ) {
        ##cat(paste('adding k/gamma\n'))
        rmn <- rmn + k/gamma
    }
    if ( model %in% c('k_dr_k0', 'k_dr_k0_coth') ) {
        ##cat(paste('adding k0/gamma\n'))
        rmn <- rmn + k0/gamma
    }
    unname(rmn)
    
}

#' Exact periodic steady state of the pulse width-modulated ODE models.
#'
#' Start and end of the ON phase, minimum, maximum and cycle mean of the
#' periodic steady state of the ODEs \code{pwmode_k}, \code{pwmode_dr},
#' \code{pwmode_k_dr} and \code{pwmode_k_dr_k0}, with dilution \code{mu}
#' acting in both phases. Within each phase \code{R} relaxes monotonically
#' (towards \code{k/mu} in the ON phase, towards \code{k0/gamma} in the OFF
#' phase), so the extremes lie at the phase boundaries, \code{R0} and
#' \code{R1}. For the models with phase-switched degradation, \code{R} rises
#' during the ON phase if \code{k/mu > k0/gamma} (the normal regime) and falls
#' otherwise (the reversed regime, only possible for \code{"k_dr_k0"} with a
#' high basal rate); \code{R0} is then the maximum.
#'
#' Derivation: ON phase (duration \code{a = phi*tau}),
#' \code{dR/dt = k - lon*R}, with \code{lon = mu} (\code{lon = gamma} for model
#' \code{"k"}); OFF phase (\code{b = (1-phi)*tau}), \code{dR/dt = k0 - gamma*R},
#' with \code{k0 = k} for \code{"dr"}, \code{k0 = 0} for \code{"k"} and
#' \code{"k_dr"}. With \code{eA = exp(-lon*a)}, \code{eB = exp(-gamma*b)},
#' \code{iA = (1-eA)/lon}, \code{iB = (1-eB)/gamma}:
#' \code{R0 = (k0*iB + k*iA*eB)/(1 - eA*eB)}, \code{R1 = R0*eA + k*iA}, and
#' the mean \code{(R0*iA + k*(a-iA)/lon + R1*iB + k0*(b-iB)/gamma)/tau},
#' evaluated stably for small \code{lon} (limit: \code{iA = a},
#' \code{(a-iA)/lon = a^2/2}, the ON phase of \code{\link{get_rmean_nogrowth}}).
#' @param k transcription rate in the ON phase.
#' @param gamma total loss rate in the OFF phase, \code{dr + mu}; if missing,
#' calculated from \code{dr} and \code{mu}.
#' @param k0 basal transcription rate in the OFF phase (model
#' \code{"k_dr_k0"} only).
#' @param dr degradation rate; if missing, \code{gamma - mu}.
#' @param mu growth rate (dilution, both phases); if missing or \code{NA}, 0
#' (no growth: \code{gamma} or \code{dr} is then the total loss rate in the OFF
#' phase, and the result equals the closed forms, \code{\link{get_rmean_nogrowth}}).
#' @param phi duty cycle, the fraction of the period in the ON phase.
#' @param tau period.
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"}.
#' @return data frame with columns \code{R0} (start of the ON phase),
#' \code{R1} (end of the ON phase), \code{Rmin}, \code{Rmax}, \code{mean}.
#' @seealso \code{\link{get_rmean}}, \code{\link{get_ramp}},
#' \code{\link{get_rates}}
#'@export
get_rcycle <- function(k, gamma, k0 = 0, dr, mu, phi, tau,
                             model = c('k', 'dr', 'k_dr', 'k_dr_k0')) {

    model <- match.arg(model)
    if ( missing(mu) ) mu <- 0
    mu[is.na(mu)] <- 0                          # no growth
    if ( missing(gamma) ) gamma <- dr + mu
    cy <- .pwm_cycle(k = k, k0 = k0, gamma = gamma, mu = mu, phi = phi, tau = tau,
                     model = model)
    res <- as.data.frame(cy)
    rownames(res) <- NULL
    res
}

## internal: get_rcycle without argument handling and data frame, for the
## root finding in get_rates and get_times (called thousands of times per gene)
.pwm_cycle <- function(k, k0, gamma, mu, phi, tau, model) {
    lon <- if ( model == 'k' ) gamma else mu   # loss rate in the ON phase
    if ( model == 'dr' ) k0 <- k
    if ( model %in% c('k', 'k_dr') ) k0 <- 0

    a <- phi*tau
    b <- (1-phi)*tau
    x <- lon*a
    small <- abs(x) < 1e-8
    lon1 <- ifelse(small, 1, lon)
    ## (1-exp(-lon*a))/lon and (a - iA)/lon, stable for lon -> 0
    iA <- ifelse(small, a - lon*a^2/2, -expm1(-x)/lon1)
    jA <- ifelse(small, a^2/2 - lon*a^3/6, (a - iA)/lon1)
    eA <- exp(-x)
    eB <- exp(-gamma*b)
    iB <- -expm1(-gamma*b)/gamma

    R0 <- (k0*iB + k*iA*eB)/(1 - eA*eB)   # start of ON phase
    R1 <- R0*eA + k*iA                    # end of ON phase
    ion <- R0*iA + k*jA                   # integral over ON phase
    ioff <- R1*iB + k0*(b - iB)/gamma     # integral over OFF phase
    list(R0 = R0, R1 = R1, Rmin = pmin(R0, R1), Rmax = pmax(R0, R1),
         mean = (ion + ioff)/tau)
}

#' Exact cycle mean of the pulse width-modulated ODE models.
#'
#' Periodic steady-state mean abundance of the ODEs \code{pwmode_k},
#' \code{pwmode_dr}, \code{pwmode_k_dr} and \code{pwmode_k_dr_k0}, with
#' dilution \code{mu} acting in both phases (see
#' \code{\link{get_rcycle}} for the derivation). For \code{model="k"}
#' this is \code{phi*k/gamma}. For the models with phase-switched
#' degradation, \code{\link{get_rmean_nogrowth}} (the closed forms of the
#' slides) is the exact mean of a model without dilution in the ON phase
#' (\code{gamma = dr + mu} only in the OFF phase) and is therefore higher; the
#' two agree for \code{mu = 0} or missing.
#' @inheritParams get_rcycle
#' @return numeric vector of cycle means.
#'@export
get_rmean <- function(k, gamma, k0 = 0, dr, mu, phi, tau,
                            model = c('k', 'dr', 'k_dr', 'k_dr_k0')) {

    model <- match.arg(model)
    if ( missing(mu) ) mu <- 0
    mu[is.na(mu)] <- 0                          # no growth
    if ( missing(gamma) )
        gamma <- dr + mu
    if ( model == 'k' ) # recycled over tau
        return(unname(rep_len(phi*k/gamma,
                              max(length(k), length(gamma), length(phi), length(tau)))))
    get_rcycle(k = k, gamma = gamma, k0 = k0, mu = mu, phi = phi, tau = tau,
                     model = model)$mean
}

#' Exact amplitude of the pulse width-modulated ODE models.
#'
#' Absolute (\code{Rmax - Rmin}) or relative (\code{(Rmax - Rmin)/mean})
#' amplitude of the periodic steady state, with dilution \code{mu} in both
#' phases (see \code{\link{get_rcycle}}). For model \code{"k"} this is
#' the same as \code{\link{get_ramp_nogrowth}}. For the models with phase-switched
#' degradation, \code{\link{get_ramp_nogrowth}} uses \code{k*phi*tau}, the rise in an ON
#' phase without dilution; with dilution the rise is
#' \code{(k/mu - R0)*(1 - exp(-mu*phi*tau))}. The relative amplitude does not
#' depend on the scale of \code{k} (with \code{k0} given relative to it), so
#' \code{k} defaults to 1 for it.
#' @inheritParams get_rcycle
#' @param relative relative amplitude, \code{(Rmax - Rmin)/mean}; otherwise
#' absolute.
#' @return numeric vector of amplitudes.
#'@export
get_ramp <- function(k = 1, gamma, k0 = 0, dr, mu, phi, tau,
                           relative = TRUE,
                           model = c('k', 'dr', 'k_dr', 'k_dr_k0')) {

    model <- match.arg(model)
    if ( missing(mu) ) mu <- 0
    mu[is.na(mu)] <- 0                          # no growth
    if ( missing(gamma) )
        gamma <- dr + mu
    cy <- get_rcycle(k = k, gamma = gamma, k0 = k0, mu = mu, phi = phi,
                           tau = tau, model = model)
    ramp <- cy$Rmax - cy$Rmin
    if ( relative ) ramp <- ramp/cy$mean
    unname(ramp)
}

#' Exact rates of the pulse width-modulated ODE models from abundance data.
#'
#' Recovers the transcription rate \code{k}, the degradation rate \code{dr}
#' and, for model \code{"k_dr_k0"}, the basal rate \code{k0}, from the mean
#' abundance \code{R}, the relative amplitude \code{a} (or the absolute
#' amplitude \code{A}) and, for \code{"k_dr_k0"}, the minimum \code{Rmin} (or
#' maximum \code{Rmax}), given the duty cycle, the period and the growth rate
#' \code{mu}, with dilution in both phases (see \code{\link{get_rcycle}}).
#' Without growth (\code{mu} missing, \code{NA} or 0 for all inputs), the
#' closed forms, \code{\link{get_rates_nogrowth}}, are exact and are used;
#' \code{dr} is then the total loss rate in the OFF phase (degradation and
#' any dilution), and \code{lower}, \code{upper} default to its range for
#' \code{gamma*tau}. Model \code{"k"} is exact in
#' \code{\link{get_rates_nogrowth}} for any \code{mu} and is passed on to it.
#'
#' The relative quantities \code{a} and \code{Rmin/R} do not depend on the
#' scale of \code{k}: \code{dr} (and \code{q = k0/k}) are found from them by
#' root finding, then \code{k = R/mean(k = 1)}. For \code{"dr"} and
#' \code{"k_dr"}, \code{a} is monotone in \code{dr}. For \code{"k_dr_k0"},
#' \code{Rmin/R} is monotone in \code{q} within each regime (normal: \code{R}
#' rises during the ON phase, \code{q < gamma/mu}; reversed: \code{R} falls,
#' \code{q > gamma/mu}), but not across them, so the regime must be chosen;
#' for each \code{dr}, \code{q} is matched to \code{Rmin/R}, and \code{dr} to
#' \code{a} (checked numerically to have a single root over typical YRO
#' conditions).
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"}.
#' @param a relative amplitude, \code{(Rmax - Rmin)/R}.
#' @param A absolute amplitude, used if \code{a} is missing.
#' @param R mean abundance.
#' @param Rmin minimal abundance (model \code{"k_dr_k0"}).
#' @param Rmax maximal abundance, used if \code{Rmin} is missing.
#' @param phi duty cycle.
#' @param tau period.
#' @param mu growth rate; if missing, \code{NA} or 0, no growth (see above);
#' an element that is \code{NA} counts as 0.
#' @param regime for \code{"k_dr_k0"}: \code{"normal"} (R rises during the ON
#' phase) or \code{"reversed"}.
#' @param lower,upper range of \code{dr} searched (without growth: of
#' \code{gamma*tau}, as in \code{\link{get_rates_nogrowth}}, default 1e-6 and 20).
#' @param n number of grid points on a log scale used to bracket the roots.
#' @param verb verbosity.
#' @param ... passed to \code{\link{get_rates_nogrowth}} (model \code{"k"}
#' or no growth).
#' @return data frame with columns \code{k}, \code{dr} and, for
#' \code{"k_dr_k0"}, \code{k0}; \code{NA} where no solution exists.
#'@export
get_rates <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                            a = NA, A = NA, R = NA, Rmin = NA, Rmax = NA,
                            phi, tau, mu,
                            regime = c('normal', 'reversed'),
                            lower = 1e-4, upper = 1e3, n = 200,
                            verb = 0, ...) {

    ## NOTE: model is not matched against the choices here, so that the
    ## variants of the closed forms (e.g. 'k_dr_coth', 'k_dr_k0_coth') can be
    ## passed on to get_rates_nogrowth without growth
    model <- model[1]
    regime <- match.arg(regime)
    if ( all(is.na(a)) ) a <- A/R
    if ( missing(mu) ) mu <- NA
    nogrowth <- all(is.na(mu) | mu == 0)
    if ( model == 'k' | nogrowth ) {
        args <- list(model = model, a = a, R = R, Rmin = Rmin, Rmax = Rmax,
                     phi = phi, tau = tau, mu = mu, verb = verb, ...)
        if ( !missing(lower) ) args$lower <- lower
        if ( !missing(upper) ) args$upper <- upper
        return(do.call(get_rates_nogrowth, args))
    }
    ## model k (exact in the closed form for any mu) and all models without
    ## growth were passed on above; what is left with growth are the models
    ## with phase-switched degradation. Other names (the _coth variants of the
    ## closed forms) exist only without growth.
    if ( !model %in% c('dr', 'k_dr', 'k_dr_k0') )
        stop("model '", model, "' is only available without growth, see get_rates_nogrowth")
    mu[is.na(mu)] <- 0

    m <- Rmin/R
    if ( model == 'k_dr_k0' ) m <- ifelse(is.na(m), Rmax/R - a, m)

    one <- function(a, R, m, phi, tau, mu) {
        na <- if ( model == 'k_dr_k0' ) c(k = NA, dr = NA, k0 = NA) else c(k = NA, dr = NA)
        if ( any(!is.finite(c(a, R, phi, tau, mu))) ) return(na)
        cyc <- function(dr, q) .pwm_cycle(k = 1, k0 = if ( model == 'dr' ) 1 else q,
                                          gamma = dr + mu, mu = mu, phi = phi, tau = tau,
                                          model = model)
        ## k0/k matching Rmin/R for a given dr, in the chosen regime
        qfit <- function(dr) {
            if ( !is.finite(m) ) return(NA)
            qb <- if ( mu > 0 ) (dr + mu)/mu else Inf        # regime boundary
            rng <- if ( regime == 'normal' ) c(0, min(qb*(1 - 1e-9), 1e8))
                   else c(qb*(1 + 1e-9), qb*1e6)
            if ( !all(is.finite(rng)) ) return(NA)
            h <- function(q) { cy <- cyc(dr, q); cy$Rmin/cy$mean - m }
            hr <- c(h(rng[1]), h(rng[2]))
            if ( any(!is.finite(hr)) || sign(hr[1]) == sign(hr[2]) ) return(NA)
            stats::uniroot(h, rng, tol = 1e-12)$root
        }
        f <- function(dr) {
            q <- if ( model == 'k_dr_k0' ) qfit(dr) else 0
            if ( is.na(q) ) return(NA)
            cy <- cyc(dr, q)
            (cy$Rmax - cy$Rmin)/cy$mean - a
        }
        grid <- exp(seq(log(lower), log(upper), length.out = n))
        fg <- sapply(grid, f)
        ok <- which(is.finite(fg[-n]) & is.finite(fg[-1]) & sign(fg[-n]) != sign(fg[-1]))
        if ( length(ok) == 0 ) {
            if ( verb > 0 ) cat('get_rates: no solution for dr\n')
            return(na)
        }
        if ( length(ok) > 1 & verb > 0 )
            cat(paste('get_rates:', length(ok), 'roots for dr, taking the largest\n'))
        j <- max(ok)
        dr <- stats::uniroot(f, grid[c(j, j + 1)], tol = 1e-12)$root
        q <- if ( model == 'k_dr_k0' ) qfit(dr) else 0
        k <- R/cyc(dr, q)$mean
        if ( model == 'k_dr_k0' ) c(k = k, dr = dr, k0 = q*k) else c(k = k, dr = dr)
    }
    res <- do.call(rbind, Map(one, a = a, R = R, m = m, phi = phi, tau = tau, mu = mu))
    res <- as.data.frame(res)
    rownames(res) <- NULL
    res
}

#' Abundance amplitudes, closed form without dilution in the ON phase.
#'
#' As \code{\link{get_rmean_nogrowth}}: exact for model \code{"k"} and for
#' \code{mu = 0}; the amplitude of the models with phase-switched degradation
#' is \code{k*phi*tau}. Use \code{\link{get_ramp}} for the exact amplitude
#' with dilution in both phases.
#'
#' For the models with phase-switched degradation, the absolute amplitude is
#' calculated first if \code{k} is given (and \code{force.relative} is
#' \code{FALSE}), and divided by the mean for the relative amplitude;
#' otherwise the relative amplitude is calculated directly, which requires
#' \code{k} only for model \code{"k_dr_k0"}.
#' @inheritParams get_rmean_nogrowth
#' @param relative relative amplitude, \code{(Rmax - Rmin)/mean}; otherwise
#' absolute.
#' @param force.relative calculate the relative amplitude directly, also if
#' \code{k} is given.
#' @param ... unused.
#' @return numeric vector of amplitudes.
#'@export
get_ramp_nogrowth <- function(gamma, dr, mu, phi, tau, relative = TRUE,
                     k, k0, force.relative = FALSE, use.coth = FALSE,
                     model = c('k', 'dr', 'k_dr', 'k_dr_k0'), ...) {

    if ( length(model)>1 ) model <- model[1]

    if ( model %in% c('dr', 'k_dr', 'k_dr_k0', 'k_dr_k0_coth') ) {

        if ( !missing(k) & !force.relative ) {


            ## ABS. AMPLITUDE FIRST
            ## DIRECTLY, VIA K, PHI AND TAU
            ramp <- phi*tau*k

            ## for relative amplitude we need to divide by the mean
            if ( relative ) {
                ## TODO: numerical instabilities via mean?

                if ( missing(gamma) )
                    gamma <- dr+mu

                rmean <- get_rmean_nogrowth(
                    k = k,
                    gamma = gamma,
                    phi = phi,
                    tau = tau,
                    k0 = k0,
                    model = model, use.coth = use.coth) 
                ramp <- ramp/rmean
            }
        } else {
            
            ## REL. AMPLITUDE FIRST
            ## NOTE: numerical instabilities via both
            ##       relative amplitude and mean?
            ## NOTE: k for relative amplitude is only
            ##       required for model with basal expression
            
            if ( missing(gamma) )
                gamma <- dr+mu

            gt <- gamma*tau

            ## TODO: test numerics coth vs. expm1
            ## term2 <- phi/2 * pracma::coth(gt*(1-phi)/2)
            y <- gt*(1-phi)
            if ( use.coth ) 
                term2 <-  phi/2 * pracma::coth(y/2)
            else {
                ye <- expm1(y) # exp(y)-1
                term2 <- phi/2 * (ye+2)/ye
            }
            
            ##if ( model == 'k_dr' )
            ##    term1 <- 1/gt
            ##else if ( model == 'dr' ) 
            ##    term1 <- (1+1/phi)/gt
            ##else if ( model == 'k_dr_k0' )
            ##    term1 <- (1+k0/(k*phi))/gt
            term1 <- switch (model,
                             k_dr = 1/gt,
                             dr = (1+1/phi)/gt,
                             k_dr_k0 = (1+k0/(k*phi))/gt)
            
            ramp <- 1/(term1 + term2)
            
            if ( missing(k0) ) k0 <- NA

            ## for absolute amplitude we need to multiple by the mean
            if ( !relative ) {
                rmean <- get_rmean_nogrowth(
                    k = k,
                    gamma = gamma,
                    phi = phi,
                    tau = tau,
                    k0 = k0,
                    model = model, use.coth = use.coth)
                ramp <- ramp*rmean 
            }
        }
    } else if ( model %in% c('k') ) {

        if ( missing(gamma) )
            gamma <- dr+mu

        gt <- gamma*tau
        term1 <- 1-exp(-gt*phi)
        term2 <- 1-exp(gt*(phi-1))
        term3 <- 1-exp(-gt) 

        term <- term1*term2/term3

        if ( relative ) ramp <- term/phi
        else ramp <- term*k/gamma
    }

    unname(ramp)
    
}

get_rna <- function() {

    ## TODO: wrapper for get_rmean and get_ramp
}


#' Mean protein abundance.
#'
#' Steady-state protein abundance \code{P = R*rho*l/(mu + dp)}, from the mean
#' transcript abundance \code{R}, which is calculated with
#' \code{\link{get_rmean}} if missing.
#' @param R mean transcript abundance; if missing, from
#' \code{\link{get_rmean}} with \code{mu}, \code{phi}, \code{r.model} and
#' \code{...}.
#' @param rho translating ribosomes per transcript.
#' @param l translation elongation rate (per ribosome, protein per time).
#' @param dp protein degradation rate.
#' @param mu growth rate.
#' @param phip duty cycle of translation (\code{p.model = "phi"}).
#' @param phi duty cycle of transcription, passed to \code{\link{get_rmean}}.
#' @param r.model transcript model, passed to \code{\link{get_rmean}} as
#' \code{model}.
#' @param p.model \code{"const"}: constant translation; \code{"phi"}:
#' translation in the ON phase only, \code{P} multiplied by \code{phip},
#' which assumes that translation phase and \code{R} are uncorrelated.
#' @param ... further arguments to \code{\link{get_rmean}}: \code{k},
#' \code{dr} or \code{gamma}, \code{tau}, \code{k0}.
#' @return numeric vector of mean protein abundances.
#' @export
get_pmean <- function(R, rho, l, dp, mu, phip, phi,
                      r.model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      p.model = c('const', 'phi'), ...)  { # P(mu)

    if ( length(r.model)>1 ) r.model <- r.model[1]
    if ( length(p.model)>1 ) p.model <- p.model[1]

    ## NOTE: the RNA duty cycle phi and the rates (k, dr or gamma, tau, k0)
    ## are passed on to get_rmean; they used to be referenced as free
    ## variables here, which silently took them from the calling environment.
    ## phi is a formal argument, since phi= in ... would partially match phip
    if ( missing(R) )
        R <- get_rmean(mu = mu, phi = phi, model = r.model, ...)
    p <- R*rho*l/(mu+dp)
    
    ## NOTE: translation in HOC only; multiplying by phip assumes that the
    ## translation phase and R are uncorrelated, an approximation when R
    ## oscillates with the same phases
    if ( p.model %in% c('phi') ) # translation in HOC only
        p <- p*phip

    unname(p)
}

#' Exact duty cycle and period of the pulse width-modulated ODE models.
#'
#' Recovers the duty cycle \code{phi} and the period \code{tau} from the mean
#' abundance \code{R} and the relative amplitude \code{a} (or the absolute
#' amplitude \code{A}), given the rates \code{k}, \code{dr} (and \code{k0}) and
#' the growth rate \code{mu}, with dilution in both phases (see
#' \code{\link{get_rcycle}}). The counterpart of \code{\link{get_rates}}, which
#' recovers the rates given the times. Without growth (\code{mu} missing,
#' \code{NA} or 0 for all inputs) and for model \code{"k"} (exact for any
#' \code{mu}), the closed forms of \code{\link{get_times_nogrowth}} are used,
#' with \code{gamma = dr + mu}; \code{lower}, \code{upper} are then passed on.
#'
#' For a given \code{tau}, the mean is monotone in \code{phi} (increasing, or
#' decreasing where the basal rate \code{k0} dominates), so \code{phi(tau)} is
#' found from \code{R}; then \code{tau} from \code{a(phi(tau), tau)}. Roots
#' with \code{phi} at its bounds are discarded; where several remain, the
#' largest \code{tau} is taken (as \code{get_times_nogrowth}). This happens
#' when \code{R} hardly oscillates (e.g. \code{k0/gamma} close to
#' \code{k/mu} for \code{"k_dr_k0"}), and the amplitude then carries no
#' information about \code{tau}.
#'
#' For \code{"k_dr_k0"}, \code{k0} is used if given; otherwise \code{Rmin} is
#' used instead: for given \code{phi} and \code{tau} the mean is linear in
#' \code{k0}, \code{k0 = (R - R(k0=0))/R(k=0, k0=1)}, and \code{phi} is then
#' found from \code{Rmin}.
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"}.
#' @param a relative amplitude, \code{(Rmax - Rmin)/R}.
#' @param A absolute amplitude, used if \code{a} is missing.
#' @param R mean abundance, in the units of \code{k}.
#' @param Rmin minimal abundance (model \code{"k_dr_k0"} without \code{k0}).
#' @param k transcription rate.
#' @param k0 basal transcription rate (model \code{"k_dr_k0"}).
#' @param dr degradation rate.
#' @param mu growth rate; if missing, \code{NA} or 0, no growth (see above);
#' an element that is \code{NA} counts as 0.
#' @param gamma total loss rate in the OFF phase, \code{dr + mu}; used for
#' \code{dr} if that is missing.
#' @param lower,upper range of \code{tau} searched (without growth: as in
#' \code{\link{get_times_nogrowth}}, default 1e-6 and 100).
#' @param n number of grid points on a log scale used to bracket the roots.
#' @param verb verbosity.
#' @param ... passed to \code{\link{get_times_nogrowth}} (model \code{"k"} or
#' no growth).
#' @return data frame with columns \code{phi} and \code{tau}; \code{NA} where
#' no solution exists.
#'@export
get_times <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      a = NA, A = NA, R = NA, Rmin = NA,
                      k = NA, k0 = NA, dr = NA, mu = NA, gamma = NA,
                      lower = 1e-3, upper = 100, n = 200, verb = 0, ...) {

    model <- model[1]
    if ( all(is.na(a)) ) a <- A/R
    if ( missing(mu) ) mu <- NA
    nogrowth <- all(is.na(mu) | mu == 0)
    mu0 <- ifelse(is.na(mu), 0, mu)
    if ( all(is.na(dr)) ) dr <- gamma - mu0
    if ( model == 'k' | nogrowth ) {
        args <- list(model = model, a = a, R = R, Rmin = Rmin, k = k,
                     gamma = dr + mu0, k0 = k0, verb = verb, ...)
        if ( !missing(lower) ) args$lower <- lower
        if ( !missing(upper) ) args$upper <- upper
        return(do.call(get_times_nogrowth, args))
    }
    ## model k (exact in the closed form for any mu) and all models without
    ## growth were passed on above; other names exist only without growth
    if ( !model %in% c('dr', 'k_dr', 'k_dr_k0') )
        stop("model '", model, "' is only available without growth, see get_times_nogrowth")

    one <- function(a, R, Rmin, k, k0, dr, mu) {
        na <- c(phi = NA, tau = NA)
        if ( any(!is.finite(c(a, R, k, dr, mu))) ) return(na)
        usek0 <- model != 'k_dr_k0' || is.finite(k0)
        if ( !usek0 && !is.finite(Rmin) ) return(na)
        cyc <- function(phi, tau, k, k0) .pwm_cycle(k = k, k0 = if ( model == 'dr' ) k else
                                                       if ( model == 'k_dr' ) 0 else k0,
                                                    gamma = dr + mu, mu = mu, phi = phi,
                                                    tau = tau, model = model)
        eps <- 1e-6
        pgrid <- c(eps, seq(0.01, 0.99, length.out = 50), 1 - eps)
        k0of <- function(phi, tau) {
            if ( usek0 ) return(ifelse(is.finite(k0), k0, 0))
            ## the mean is linear in k0 (NA where k0 would be negative)
            k0 <- (R - cyc(phi, tau, k, 0)$mean)/cyc(phi, tau, 0, 1)$mean
            ifelse(k0 < 0, NA, k0)
        }
        ## phi(tau) from R (k0 given) or from R and Rmin (k0 unknown); h is
        ## evaluated on a phi grid in one vectorised call, then refined in the
        ## bracket
        phifit <- function(tau) {
            h <- function(phi) { k0x <- k0of(phi, tau)
                cy <- cyc(phi, tau, k, k0x)
                if ( usek0 ) cy$mean - R else cy$Rmin - Rmin }
            hg <- h(pgrid)
            j <- which(is.finite(hg[-length(hg)]) & is.finite(hg[-1]) &
                       sign(hg[-length(hg)]) != sign(hg[-1]))
            if ( length(j) == 0 ) return(NA)
            j <- j[1]
            stats::uniroot(h, pgrid[c(j, j + 1)], tol = 1e-12)$root
        }
        f <- function(tau) {
            phi <- phifit(tau)
            if ( is.na(phi) ) return(NA)
            cy <- cyc(phi, tau, k, k0of(phi, tau))
            (cy$Rmax - cy$Rmin)/cy$mean - a
        }
        grid <- exp(seq(log(lower), log(upper), length.out = n))
        fg <- sapply(grid, f)
        ok <- which(is.finite(fg[-n]) & is.finite(fg[-1]) & sign(fg[-n]) != sign(fg[-1]))
        if ( length(ok) == 0 ) {
            if ( verb > 0 ) cat('get_times: no solution for tau\n')
            return(na)
        }
        taus <- sapply(ok, function(j) stats::uniroot(f, grid[c(j, j + 1)], tol = 1e-12)$root)
        phis <- sapply(taus, phifit)
        keep <- is.finite(phis) & phis > 1e-4 & phis < 1 - 1e-4
        if ( !any(keep) ) return(na)
        taus <- taus[keep]; phis <- phis[keep]
        if ( length(taus) > 1 & verb > 0 )
            cat(paste('get_times:', length(taus), 'roots for tau, taking the largest\n'))
        j <- which.max(taus)
        c(phi = phis[j], tau = taus[j])
    }
    res <- do.call(rbind, Map(one, a = a, R = R, Rmin = Rmin, k = k, k0 = k0,
                              dr = dr, mu = mu0))
    res <- as.data.frame(res)
    rownames(res) <- NULL
    res
}

#' Duty cycle and period, closed form without dilution in the ON phase.
#'
#' The inverse of \code{\link{get_rmean_nogrowth}} and
#' \code{\link{get_ramp_nogrowth}} for the duty cycle and period, given the
#' rates: exact for model \code{"k"} and without growth; for the models with
#' phase-switched degradation, \code{phi = A/(k*tau)} and \code{tau} from the
#' relative amplitude (\code{root_tau_*}), where \code{gamma} is the total loss
#' rate in the OFF phase. Use \code{\link{get_times}} with the growth rate for
#' the exact duty cycle and period with dilution in both phases.
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"}; a root function \code{root_tau_<model>} must exist.
#' @param a relative amplitude, \code{(Rmax - Rmin)/R}.
#' @param R mean abundance.
#' @param Rmin minimal abundance (model \code{"k_dr_k0"}).
#' @param k transcription rate.
#' @param gamma total loss rate in the OFF phase; if \code{NA}, \code{dr + mu}.
#' @param dr degradation rate, used if \code{gamma} is \code{NA}.
#' @param mu growth rate, used if \code{gamma} is \code{NA}; \code{NA}
#' counts as 0.
#' @param k0 unused; the closed form of \code{"k_dr_k0"} uses \code{Rmin}.
#' @param lower,upper range of \code{tau} searched.
#' @param tol tolerance of the root finding.
#' @param verb verbosity.
#' @param ... unused.
#' @return data frame with columns \code{phi} and \code{tau}.
#' @seealso \code{\link{get_tau}}
#' @export
get_times_nogrowth <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      a = NA, R = NA, Rmin = NA, 
                      k = NA, gamma = NA, dr = NA, mu = NA, k0 = NA, 
                      lower = 1e-6, upper = 100, tol = 1e-9,
                      verb = 0, ...) {

    if ( length(model)>1 ) model <- model[1]

    if ( all(is.na(gamma)) )
        gamma <- dr + ifelse(is.na(mu), 0, mu) # mu = NA: no growth


    ## model k:
    ## 1.) get phi from R/(k/gamma)
    ## 2.) get tau from relative amplitude via uniroot
    phi <- NA
    if ( model %in% c('k') )
        phi <- R*gamma/k

    ## TODO: get ton=phi*tau for models dr ..
    ## and use in root functions?


    ## NOTE: using Map allows vectorization of input
    tau <- unlist(Map(get_tau,
                      model = model,
                      a = a,  R = R, Rmin = Rmin, 
                      k = k, gamma = gamma,
                      phi = phi,
                      lower = lower, upper = upper, tol = tol,
                      verb = verb))
    
    ## TODO: omit A from arguments and solve this
    ## in root functions?
    A <- a*R
    if ( model %in% c('k_dr', 'dr', 'k_dr_k0') )
        phi <- A/(k*tau)

    if ( length(phi)==1 )
        phi <- rep(phi, length(tau))

    if ( all(is.na(tau)) )
        warning('no periods calculated')
    if ( all(is.na(phi)) )
        warning('no duty cycles calculated')
    
    res <- cbind.data.frame(phi=phi, tau=tau)
    rownames(res) <- NULL
    res

    ## models dr:
    ## fit ton=phi*tau and untangel via absolute amplitude
}

#' Period from abundance data, closed form without dilution in the ON phase.
#'
#' Finds the period \code{tau} as the root of the model's relative amplitude
#' equation, \code{root_tau_<model>}, with
#' \code{\link[rootSolve]{uniroot.all}}; where several roots exist, the
#' largest is taken. For the models with phase-switched degradation,
#' \code{phi = A/(k*tau)}, so only \code{tau > A/k} is searched. Called by
#' \code{\link{get_times_nogrowth}} for one set of values.
#' @param a relative amplitude, \code{(Rmax - Rmin)/R}.
#' @param R mean abundance.
#' @param Rmin minimal abundance (model \code{"k_dr_k0"}).
#' @param k transcription rate.
#' @param gamma total loss rate in the OFF phase.
#' @param phi duty cycle (model \code{"k"}); for the other models calculated
#' from \code{A/(k*tau)}.
#' @param model model name; a function \code{root_tau_<model>} must exist.
#' @param lower,upper range of \code{tau} searched.
#' @param tol tolerance of the root finding.
#' @param verb verbosity; with 0, errors of the root finding are silent.
#' @param ... unused.
#' @return the period, or \code{NA} if no root was found.
#' @export
get_tau <- function(a, R = NA, Rmin = NA,  k, gamma, phi,
                    model, ## must exist as root finding function
                    lower = 1e-6, upper = 1e4, tol = 1e-9,
                    verb = 0, ...) {

    ## Solve f(x) = 0 for x in a reasonable range
    ## where x = tau
    
    ## get model-specific function for root-finding,
    ## each returns x = tau
    rootf <- get(paste0('root_tau_', model), mode = 'function')

    ## TODO: omit A from arguments and solve this
    ## in root functions?
    A <- a*R

    ## NOTE: for the models with phase-switched degradation, phi = A/(k*tau)
    ## so tau > A/k; the root functions have a pole at tau = A/k (phi = 1)
    ## and are not defined below it, and uniroot.all's grid could step over
    ## the root; search above the pole
    if ( model %in% c('dr', 'k_dr', 'k_dr_k0') && is.finite(A/k) && A/k > 0 )
        lower <- max(lower, A/k*(1 + 1e-6))
    
    ## TODO: find better solution or test whether taking the
    ## the highest root is always appropriate. 
    solution <- try(rootSolve::uniroot.all(rootf,
                                           a = a, gamma = gamma, phi = phi,
                                           A = A, k = k, R = R, Rmin = Rmin,
                                           lower = lower, upper = upper,
                                           tol = tol),
                    silent = verb==0)
    if ( FALSE ) {
        solution <- try(stats::uniroot(rootf, a = a, gamma = gamma, phi = phi,
                                       A = A, k = k, R = R, Rmin = Rmin,
                                       lower = lower, upper = upper, tol = tol),
                        silent = verb==0)
    }
    
    if ( class(solution)=="try-error" | length(solution)==0 ) {
        if ( verb>0 )
            cat(paste0('Model <', model, '> failed with:',
                       '\n\ta= ', a,
                       ';\n\tgamma= ', gamma,
                       ';\n\tphi= ', phi,
                       ';\n\tA= ', A,
                       ';\n\tR= ', R,
                       ';\n\tRmin= ', Rmin,
                       ';\n\tk= ', k,
                       ';\n\t#solutions= ', length(solution),
                       ';\n\tclass= ', class(solution),
                       ';\n'))
        return(NA)
    }
    

    ## return tau
    if ( length(solution)>1 ) {
        if ( verb>0 ) 
            cat(paste0('taking max of ', length(solution), ' roots\n'))
        solution <- max(solution)
    }
    solution
}

## where x=tau
root_tau_k <- function(x, a, gamma, phi, A = NA, k = NA, R = NA, Rmin = NA) {

    GT <- gamma*x

    lhs <- a * phi
    term1 <- (- expm1(-phi * GT)) / (- expm1(- GT))
    term2 <- ( - expm1(GT * (phi - 1)))
    rhs <- term1 * term2
    return(lhs - rhs)
}


## where x=tau
root_tau_k_dr <- function(x, a, gamma, phi = NA, A, k, R = NA, Rmin = NA) {

    GT <- gamma*x
    phi <- A/(k*x)

    lhs <- 1/a
    term1 <- phi/expm1(GT * (1-phi)) 
    term2 <- phi/2
    term3 <- 1/GT
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

## where x=tau
root_tau_dr <- function(x, a, gamma, phi = NA, A, k, R = NA, Rmin = NA) {
    
    GT <- gamma*x
    phi <- A/(k*x)

    lhs <- 1/a
    term1 <- phi/expm1(GT * (1-phi)) 
    term2 <- phi/2
    term3 <- (1+1/phi)/GT
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

root_tau_k_dr_k0 <- function(x, a, gamma, phi = NA, A, k, R, Rmin) {

    GT <- gamma*x
    phi <- A/(k*x)

    lhs <- (1 - Rmin/R)/a
    ##term1 <- (1-phi)/expm1(x*(1-phi)) # NOTE: wrong but worked e.g. for ATO3 in chin12 data
    term1 <- (phi-1)/expm1(GT*(1-phi))
    term2 <- phi/2
    term3 <- 1/GT
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

#' Rates from abundance data, closed form without dilution in the ON phase.
#'
#' The inverse of \code{\link{get_rmean_nogrowth}} and
#' \code{\link{get_ramp_nogrowth}}: exact for model \code{"k"} and for
#' \code{mu = 0} or \code{mu = NA}, where \code{dr} is the total loss rate in
#' the OFF phase (degradation and dilution). \code{lower}, \code{upper} bound
#' \code{gamma*tau}. Use \code{\link{get_rates}} with the growth rate for the
#' exact rates with dilution in both phases.
#'
#' For the models with phase-switched degradation, \code{k = A/(phi*tau)}
#' (\code{\link{get_transcription}}) and \code{gamma*tau} from the relative
#' amplitude (\code{\link{get_degradation}}); for model \code{"k"},
#' \code{gamma*tau} first and then \code{k = R*gamma/phi}. For
#' \code{"k_dr_k0"}, \code{k0} from \code{Rmin} or \code{Rmax}.
#' @param model one of \code{"k"}, \code{"dr"}, \code{"k_dr"},
#' \code{"k_dr_k0"}, or the variants \code{"k_dr_coth"}, \code{"dr_coth"},
#' \code{"k_dr_k0_coth"}; a function \code{root_<model>} must exist.
#' @param a relative amplitude, \code{(Rmax - Rmin)/R}.
#' @param A absolute amplitude, used if \code{a} is missing.
#' @param R mean abundance.
#' @param Rmin minimal abundance (model \code{"k_dr_k0"}).
#' @param Rmax maximal abundance, used for \code{Rmin} if that is missing.
#' @param phi duty cycle.
#' @param tau period.
#' @param mu growth rate, subtracted from \code{gamma} to give \code{dr};
#' \code{NA} counts as 0. A \code{dr} below 0 is returned as \code{NA}.
#' @param k transcription rate; if \code{NA}, calculated.
#' @param k0 unused.
#' @param gamma unused; overwritten by the fitted loss rate.
#' @param lower,upper range of \code{gamma*tau} searched.
#' @param tol tolerance of the root finding.
#' @param verb verbosity.
#' @param ... unused.
#' @return data frame with columns \code{k}, \code{dr} and, for
#' \code{"k_dr_k0"}, \code{k0}.
#'@export
get_rates_nogrowth <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      a = NA, A = NA, R = NA, Rmin = NA, Rmax =NA,
                      phi = NA, tau = NA, mu = NA,
                      k = NA, k0 = NA, gamma = NA, 
                      lower = 1e-6, upper = 20, tol = 1e-9,
                      verb = 0, ...) {

    if ( length(model)>1 ) model <- model[1]

    ## REQUIRED:
    ## * period tau and duty cycle phi,
    ## * relative amplitude a and mean abundance R,
    ## * additionally, the model with basal transcription requires
    ##   Rmin or Rmax
    
    ## generate required values
    if ( all(is.na(a)) ) a <- A/R

    if ( !model %in% c('k')  ) {
        
        ## get transcription rate from abs. amplitude
        if ( all(is.na(k)) ) {
            if ( all(is.na(a*R)) ) 
                warning('transcription rate k requires absolute amplitude A')
            
            k <- get_transcription(A = a*R, phi = phi, tau = tau, model = model)
        }
        ## only required for k_dr_k0
        ## get Rmin from Rmax and k
        if ( all(is.na(Rmin)) & !any(is.na(Rmax)) ) Rmin <- Rmax - k*phi*tau
    }

    ## NOTE: using Map allows vectorization of input
    dr <- unname(unlist(Map(get_degradation,
                     model = model,
                     a=a, R=R, Rmin=Rmin, 
                     phi = phi, tau = tau, mu = mu,
                     lower = lower, upper = upper, tol = tol,
                     verb = verb)))
    ## a total loss rate below the dilution rate has no solution: dilution
    ## alone gives gamma >= mu (model k with growth, where dr = gamma - mu)
    dr[!is.na(dr) & dr < 0] <- NA

    gamma <- dr + ifelse(is.na(mu), 0, mu) # element-wise; mu may be a vector or a scalar
    if ( model %in% c('k') ) {
        k <- get_transcription(R = R, gamma = gamma, phi = phi, model = model)
        if ( length(k)==1 )
            k <- rep(k, length(dr))
    }

    res <- cbind.data.frame(k=k, dr=dr)
    rownames(res) <- NULL

    ## add k0
    if ( model %in% c('k_dr_k0', 'k_dr_k0_coth') ) {
        res$k0 <- get_basal(k = k, gamma = gamma, tau = tau, phi = phi,
                            Rmin = Rmin, Rmax = Rmax, verb = 0)
        if ( all(is.na(res$k0)) )
            warning('no basal transcription rates calculated')
    }
    
    if ( all(is.na(k)) )
        warning('no transcription rates calculated')
    if ( all(is.na(dr)) )
        warning('no degradation rates calculated')
    
    res
}

## TODO: map fluorescence to a rough transcript/cell number
#' Transcription rate from abundance data, closed form.
#'
#' For model \code{"k"}, from the mean abundance, \code{k = R*gamma/phi};
#' for the models with phase-switched degradation, from the absolute
#' amplitude, \code{k = A/(phi*tau)}.
#' @param R mean abundance (model \code{"k"}).
#' @param A absolute amplitude (models with phase-switched degradation).
#' @param phi duty cycle.
#' @param dr degradation rate (model \code{"k"}, with \code{mu} replaces
#' \code{gamma}).
#' @param mu growth rate (model \code{"k"}, see \code{dr}).
#' @param gamma total loss rate (model \code{"k"}).
#' @param tau period (models with phase-switched degradation).
#' @param model model name.
#' @return the transcription rate.
#' @keywords internal
get_transcription <- function(R = NA, A = NA, phi = NA,
                              dr = NA, mu = NA, gamma = NA, tau = NA,
                              model = c('k', 'dr','k_dr', 'k_dr_k0')) {

    if ( model %in% c('k') ) {
        if ( !is.na(dr) & !is.na(mu) )
            gamma <- mu+dr
        k <- R*gamma/phi
    } else 
        k <- A/(phi*tau)
    
    k
}

## TODO: use this!
get_basal <- function(k, gamma, tau, phi, Rmin, Rmax, verb = 1) {

    ##beta <- exp(gamma*tau*(1-phi))
    betam1 <- expm1(gamma*tau*(1-phi)) 
    
    ## NOTE: an Rmin or Rmax that is NA counts as not given (get_rates
    ## passes both, with NA defaults); where both are available, their
    ## basal rates are averaged, element-wise
    k0min <- k0max <- NA
    if ( !missing(Rmin) ) 
        k0min <- (Rmin - k*phi*tau/betam1)*gamma
    if ( !missing(Rmax) ) 
        k0max <- (Rmax - k*phi*tau*(1 + 1/betam1))*gamma
    both <- !is.na(k0min) & !is.na(k0max)
    if ( verb>0 & any(both) )
        cat(paste('mean of rmin- and rmax-based basal rates\n'))
    k0 <- ifelse(both, (k0min + k0max)/2, ifelse(is.na(k0min), k0max, k0min))
    k0
}


#' Fit a transcript degradation rate from relative amplitude and
#' oscillation parameters, using the \code{\link[stats]{uniroot}}  function
#'
#' @param a relative abundance amplitude, (max(x)-min(x))/mean(x).
#' @param R mean abundance, required for model with basal transcription.
#' @param Rmin minimal abundance, required for model with basal transcription.
#' @param phi duty cycle.
#' @param tau period.
#' @param mu growth rate, used only for correction (gamma=growth+degradation).
#' @param lower the lower end point of the interval to be  searched by \code{\link[stats]{uniroot}}.
#' @param upper the upper end point of the interval to be  searched by \code{\link[stats]{uniroot}}.
#' @param model model name; a function \code{root_<model>} must exist.
#' @param tol error tolerance of the \code{\link[stats]{uniroot}}  call.
#' @param verb verbosity; with 0, errors of the root finding are silent.
#' @param ... unused.
#' @return the degradation rate \code{gamma - mu} (or \code{gamma} if
#' \code{mu} is missing or \code{NA}), \code{NA} if no root was found.
#' @keywords internal
get_degradation <- function(a,  R, Rmin, 
                            phi, tau, mu,
                            model, ## must exist as root finding function
                            lower = 1e-6, upper = 1e4, tol = 1e-9,
                            verb = 0, ...) {
    

    
    ## get model-specific function for root-finding,
    ## each returns x = gamma * tau
    rootf <- get(paste0('root_', model), mode = 'function')
    

    ## Solve f(x) = 0 for x in a reasonable range
    ## where x = gamma*tau
    solution <- try(stats::uniroot(rootf, a=a, phi=phi, R=R, Rmin=Rmin,
                                   lower = lower, upper = upper, tol = tol),
                    silent = verb==0)
    if ( class(solution)=="try-error" ) {
        if ( verb>0 )
            cat(paste0('Model <', model, '> failed with:',
                       '\n\ta=', a,
                       '\n\tR=', R,
                       '\n\tRmin=', Rmin,
                       '\n\tphi=', phi,
                       '\n\ttau=', tau,
                       '\n\tmu=',  mu, '\n'))
        return(NA)
    }
        
    ## Recover gamma from x = gamma * tau
    GT <- solution$root
    gamma <- GT / tau

    ## correct for growth
    ## gamma = degradation + growth
    if ( !missing(mu) && !is.na(mu) )
        gamma <- gamma - mu

    gamma
}

## ROOT FINDING FUNCTIONS USED IN get_degradation

## transcription with const. degradation,
## where x = gamma * tau
root_k <- function(x, a, phi, R = NA, Rmin = NA) {
    lhs <- a * phi
    ##term1 <- (1 - exp(-phi * x)) / (1 - exp(-x))
    ##term2 <- (1 - exp(x * (phi - 1)))
    term1 <- (- expm1(-phi * x)) / (- expm1(-x))
    term2 <- ( - expm1(x * (phi - 1)))
    rhs <- term1 * term2
    return(lhs - rhs)
}
    

## anti-phasic transcription/degradation,
## where x = gamma * tau
root_k_dr <- function(x, a, phi, R = NA, Rmin = NA) {
        
    lhs <- 1/a
    term1 <- phi/expm1(x * (1-phi)) 
    term2 <- phi/2
    term3 <- 1/x
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

root_dr <- function(x, a, phi, R = NA, Rmin = NA) {
    lhs <- 1/a
    term1 <- phi/expm1(x * (1-phi)) 
    term2 <- phi/2
    term3 <- (1+1/phi)/x
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

root_k_dr_k0 <- function(x, a, phi, R, Rmin) {
    lhs <- (1 - Rmin/R)/a
    ##term1 <- (1-phi)/expm1(x*(1-phi)) # NOTE: wrong but worked e.g. for ATO3 in chin12 data
    term1 <- (phi-1)/expm1(x*(1-phi))
    term2 <- phi/2
    term3 <- 1/x 
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

## NOTE: keep for testing?
root_k_dr_k0_old <- function(x, a, phi, R, Rmin) {
    
    lhs <- R - Rmin
    betam1 <- expm1(x*(1-phi)) # exp() +1
    ##term1 <- (1-phi) /betam1 # TODO: should this be (phi-1)  !?
    term1 <- (phi-1) /betam1 # TODO: correct now?
    term2 <- phi/2
    term3 <- 1/x 
    rhs <- a*R*(term1 + term2 + term3)
    return(lhs - rhs)
}

root_k_dr_coth <- function(x, a, phi, R = NA, Rmin = NA) {
    
    lhs <- 1/a
    term1 <- 1/x
    term2 <- phi/2 * pracma::coth((x-x*phi)/2)
    rhs <- term1 + term2
    return(lhs - rhs)
}
root_dr_coth <- function(x, a, phi, R = NA, Rmin = NA) {
    lhs <- 1/a
    term1 <- (1+1/phi)/x  
    term2 <- phi/2 * pracma::coth((x-x*phi)/2)
    rhs <- term1 + term2
    return(lhs - rhs)
}

root_k_dr_k0_coth <- function(x, a, phi, R, Rmin) {
    lhs <- (1 - Rmin/R)/a
    term1 <- 1/x
    term2 <- phi/2 * pracma::coth((x-x*phi)/2)
    term3 <- - 1/expm1(x*(1-phi)) 
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}



### ANALYTIC MODEL

#' PWM of transcription, analytic solution.
#'
#' Analytic solution of the pulse-wave ODE for transcript abundance,
#' optionally extended to proteins via numeric integration.
#'
#' @param t Numeric vector of time points.
#' @param R0 Initial RNA abundance at t=0.
#' @param k Transcription rate during pulse (ON state).
#' @param gamma total turnover rate, gamma=dr+mu.
#' @param dr RNA degradation rate, required if gamma is missing.
#' @param mu Growth/dilution rate, required if gamma is missing.
#' @param phi Duty cycle (fraction of period ON, between 0 and 1).
#' @param tau Oscillation period.
#' @param k0 Basal transcription rate (default 0).
#' @param alpha Fourier damping factor (default 0).
#' @param theta Phase shift in radians (default 0).
#' @param N Number of Fourier terms to use (default 500).
#' @param P0 protein abundance at time t[1] (!).
#' @param dp protein degradation rate.
#' @param ell transcript elongation rate.
#' @param rho translating ribosomes per mRNA.
#' @param shift phase shift, in units of time; this currently also shifts
#' \code{R0} (with a warning).
#'
#' @return Data frame with columns:
#'   - `time`: original time points
#'   - `R`: RNA abundance
#'   - `pulse`: binary 0/1 pulse state
#'   - `P`: protein abundance (if protein params given)
#'@export
pwm_k <- function(t, R0, k, gamma, dr,  mu, k0=0, 
                  phi, tau,
                  alpha=0, shift=0, theta=0, N=5e2,
                  P0, dp, ell, rho) {


    otime <- t # store original time before shifting
    
    ## phase shift by shifting time
    ## TODO: fix this, this simple solution shifts R0 to t=shift
    if ( shift != 0 ) {
        t <- t - shift
        warning("phase shift also shifts R0")
    }
    
    ## total turnover
    if ( missing(gamma) )
        gamma <- mu+dr

    ## angular frequency
    omega <- 2*pi/tau

    ## exponential decay term
    EGT <- exp(-gamma*t)
    
    ## contribution from initial concentration
    r0 <- EGT*R0

    ## contribution from constant term
    rc <- (k*phi + k0)/gamma * (1-EGT)

    ## contribution from periodic term

    ## phase shift theta present?
    shifted <- !missing(theta)
    if ( shifted ) shifted <- shifted & theta!=0
    
    sm <- 0
    if ( !shifted ) 
        for ( n in 1:N ) { # NO PHASE SHIFT
            
            ra <- (1/n) * sin(n*pi*phi) * exp(-alpha*n)
            nom <- t*n*omega 
            ra <- ra * (cos(nom) + n*omega*sin(nom)/gamma - EGT)
            ra <- ra / (gamma + n^2*omega^2/gamma)
            
            sm <- sm + ra
        }
    else
        for ( n in 1:N ) { # WITH PHASE SHIFT
            
            ra <- (1/n) * sin(n*pi*phi) * exp(-alpha*n)
            tnom <- n*(t*omega - theta) 
            ra <- ra * (gamma*cos(tnom) + n*omega*sin(tnom) +
                        EGT *(-gamma*cos(n*theta) + n*omega*sin(n*theta)))
            ra <- ra / (gamma^2 + n^2*omega^2)
            
            sm <- sm + ra
        }
    rp <- 2*k/pi *sm
        
    rt <- data.frame(time=otime, R=r0+rc+rp)

    ## protein if parameters are present!
    ## TODO:
    ## * test and fix this, esp. for times >> 0!!??
    ## * for t>0 do we need to still calculate all R from t==0 and integrate?
    ## * get proper analytic solution instead of just integrating over RNA?
    if ( !missing(P0) ) {

        if ( min(t)>0 )
            warning("analytic solution for proteins is currently untested",
                    "and appears wrong for min(t)>0!",
                    "Use ODE model for correct solutions.")

        ## fix for t>>0 probem: set initial time to 0
        ## TODO: is this correct?
        tp <- t-min(t)
        
        ## total turnover
        gammap <- mu+dp

        ## exponential decay term
        EGP <- exp(-gammap*tp)
    
        ## RNA contributions

        ## beta*int(R(t)):
        ## integrate calculated RNA values
        intr <- rt[,"R"] * exp(gammap * tp)
        ## Trapezoidal rule
        intr <- cumsum((intr[-1] + intr[-nrow(rt)]) / 2) * diff(tp) 

        ## translation rate
        beta <- ell*rho
        intr <- beta * intr
        
        ## protein abundance
        ## TODO: beta *rt[1,"R"] instead of 0
        pt <- EGP*(P0 + c(0, intr))

        ## bind RNA and protein abundances
        rt <- cbind.data.frame(rt, P=pt)
   
    }

    rt
}

### SOME PLOT UTILS

## plot labels
## TODO: get axis/unit functions used for stan model script
axis_labels <- c(rmean=expression('\u27E8'*R*'\u27E9'/(n/cell)),
                 rmeanau=expression('\u27E8'*R*'\u27E9'),
                 ramp=expression(tilde(R)/(n/cell)),
                 rampr=expression(tilde(r)),
                 r=expression(R(t)/(n/cell)),
                 phi=expression(duty~cycle~varphi),
                 tau=expression(period~tau/h),
                 mu=expression(growth~rate~mu/h^-1),
                 k=expression(k/(n/h)),
                 dr=expression(delta[R]/h^-1),
                 hl=expression(tau[1/2]/h),
                 bi=expression(budding~index~varphi[bud]/'%'),
                 gamma=expression(gamma/h^-1))
