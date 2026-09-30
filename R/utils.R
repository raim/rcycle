
## SIMPLE UTILITY FUNCTIONS


## norm. between 0 and 1
minmax <- function(x, ...) (x-min(x, ...))/(max(x, ...)-min(x, ...))

## moving average, from segmenTools
## TODO: also calculate quantiles, sd, etc, for each point,
## see segmenTools::clusterAverages
ma <- function (x, n = 5, circular = FALSE) {
    stats::filter(x, rep(1/n, n), sides = 2, circular = circular)
}

### CELL CYCLE TIMING
## TODO: fuse with ChemostatData/models.R

## budding time<->fraction after @Slater1977/@Tyson1979
## \tau_{bud} = \frac{\log(i_\{bud} + 1)}{\mu}
##
#' Budding index, budded time and growth rate.
#'
#' Conversions between the budding index \code{BI} (fraction of budded
#' cells), the budded time \code{BT} and the growth rate \code{mu}, after
#' Slater et al. (1977) and Tyson et al. (1979):
#' \code{BT = log(BI + 1)/mu}.
#' @param BI budding index, a fraction (not in percent).
#' @param BT duration of the budded phase.
#' @param mu growth rate.
#' @return \code{fraction2time}: the budded time; \code{time2fraction}: the
#' budding index; \code{fraction2growth}: the growth rate.
#' @name budding
NULL
#' @rdname budding
#' @export
fraction2time <- function(BI, mu) log(BI+1)/mu
#' @rdname budding
#' @export
time2fraction <- function(BT, mu) exp(BT*mu)-1
#' @rdname budding
#' @export
fraction2growth <- function(BI, BT) log(BI+1)/BT

