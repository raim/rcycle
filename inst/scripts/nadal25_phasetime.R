### DIVISION-CYCLE PSEUDOPHASE IN TIME, NADAL-RIBELLES ET AL. 2025

## Cell age by ergodic rate analysis (ERA), landmark durations and
## per-cell quality for the two deep Nadal-Ribelles et al. 2025 sets:
## ctrl (exponential growth in YPD) and nacl (15 min after 0.4 M NaCl,
## not a steady state, so ERA does not strictly apply). The growth rate
## is not known: ages are in units of the cycle. Cohort states are
## computed from the counts with the gene lists of the Jackson et al.
## 2020 extract (same 6717 genes); the published cell cycle classes
## (G1, S, G2M) serve as an external check.
## Usage: Rscript nadal25_phasetime.R [nadal25.rda] [jackson20.rda]
##        [classes directory]
## Output: the list R and table tab, saved to nadal25_phasetime.rda.

library(rcycle)
library(Matrix)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
infiles <- c("~/work/teaching/pwm/scmw/nadal25.rda",
             "~/work/teaching/pwm/scmw/jackson20.rda",
             "~/data/yeast/processedData/scRNAseq")
if ( length(args) ) infiles[seq_along(args)] <- args
load(infiles[1])
load(infiles[2])

## shared functions
utils.file <- system.file("scripts", "phasetime_utils.R", package="rcycle")
if ( utils.file=="" )
    utils.file <- "~/programs/rcycle/inst/scripts/phasetime_utils.R"
source(utils.file)

## cohort gene lists, and the states from counts (mean fraction per gene)
need <- c(phase.cohorts, mark.cohorts, "daughter")
cohorts <- J[[1]]$cohorts[need]
get_st <- function(counts, total) {
    w <- matrix(0, nrow=length(cohorts), ncol=nrow(counts),
                dimnames=list(names(cohorts), rownames(counts)))
    for ( k in names(cohorts) )
        w[k, intersect(cohorts[[k]], rownames(counts))] <-
            1/length(cohorts[[k]])
    st <- as.matrix(w %*% counts)
    t(t(st)/total)
}

R <- list()
cls <- list()
for ( d in names(N) ) {
    cat(paste("data set", d, "\n"))
    x <- N[[d]]
    x$states <- get_st(x$counts, x$total)
    x$cohorts <- cohorts
    x$mu <- NA
    id <- sub("nadal-ribelles25_", "", d)
    R[[id]] <- one_medium(x, boot=20)

    ## published classes
    cf <- file.path(infiles[3], paste0(d, "_classes.tsv"))
    if ( file.exists(cf) ) {
        cl <- read.delim(cf)
        cls[[id]] <- cl$CDC[match(colnames(x$counts), cl$ID)]
    }
}
tab <- phase_table(R)
if ( !exists("outfile") ) outfile <- "nadal25_phasetime.rda"
save(R, tab, cls, file=outfile)
print(tab, digits=2)
