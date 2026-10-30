### DIVISION-CYCLE PSEUDOPHASE IN TIME, JACKSON ET AL. 2020

## Cell age by ergodic rate analysis (ERA), landmark and group
## durations, per-cell quality, for the 11 media of Jackson et al. 2020.
## Input: a light extract of the owner's pseudophase pipeline, list J
## with per medium counts, total, states, cohorts, cdc.srt and mu
## (growth rate, natural log per hour).
## The phase is from all periodic Cyclebase genes (cdc_bins_*, five peak
## time bins); the landmarks are the peaks of the full Spellman 1998
## cohorts (sp_*) along it. The pipeline's cdc_* cohorts, Spellman
## without growth-cohort and HOC/LOC genes, lose most cycling genes and
## give the wrong cohort order in several media (ROADMAP.md).
## Usage: Rscript jackson20_phasetime.R [jackson20.rda]
## Output: the list R, saved to jackson20_phasetime.rda.

library(rcycle)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
infile <- "~/work/teaching/pwm/scmw/jackson20.rda"
if ( length(args) ) infile <- args[1]
load(infile)

## shared functions
utils.file <- system.file("scripts", "phasetime_utils.R", package="rcycle")
if ( utils.file=="" )
    utils.file <- "~/programs/rcycle/inst/scripts/phasetime_utils.R"
source(utils.file)

if ( !exists("boot.medium") ) boot.medium <- "jackson20_YPD"
R <- list()
for ( m in names(J) ) {
    cat(paste("medium", m, "\n"))
    R[[sub("jackson20_", "", m)]] <-
        one_medium(J[[m]], boot=ifelse(m %in% boot.medium, 20, 0))
}

tab <- phase_table(R)
## durations in units of T where T is unknown (no growth)
tab[is.na(tab$T), grep("^G|^S|^w", colnames(tab))] <- NA
if ( !exists("outfile") ) outfile <- "jackson20_phasetime.rda"
save(R, tab, file=outfile)
print(tab, digits=2)
