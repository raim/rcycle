
## SIMPLE UTILITY FUNCTIONS


## norm. between 0 and 1
minmax <- function(x, ...) (x-min(x, ...))/(max(x, ...)-min(x, ...))

## moving average, from segmenTools
## TODO: also calculate quantiles, sd, etc, for each point,
## see segmenTools::clusterAverages
ma <- function (x, n = 5, circular = FALSE) 
{
    stats::filter(x, rep(1/n, n), sides = 2, circular = circular)
}

### CELL CYCLE TIMING
## TODO: fuse with ChemostatData/models.R

## budding time<->fraction after @Slater1977/@Tyson1979
## \tau_{bud} = \frac{\log(i_\{bud} + 1)}{\mu}
##
#' @export
fraction2time <- function(BI, mu) log(BI+1)/mu
#' @export
time2fraction <- function(BT, mu) exp(BT*mu)-1
#' @export
fraction2growth <- function(BI, BT) log(BI+1)/BT

