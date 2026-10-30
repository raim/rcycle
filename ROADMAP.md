# rcycle roadmap: division-cycle pseudophase in time

Agreed with the owner on 2026-10-30. Goal: estimate how long cells
spend in each phase of the division cycle, and how variable that is,
from one snapshot of single-cell transcriptomes. This builds on what
`get_pseudophase` already returns (the PC1/PC2 angle `theta`, the rank
phase `phi`, cohort phases) and on the derivative and root machinery
in `R/segments.R` (points of maximal change, used to align phases).

## Why

* With `scale=TRUE` every cell is scaled to unit variance across
  cohorts, so all cells sit on a circle: amplitude is lost and noise
  cells still get an angle.
* The rank phase `phi` spreads cells evenly over the circle, so it
  hides where cells spend longer. `calibrate()` turns `phi` into time
  linearly, which ignores that an exponentially growing population has
  more young than old cells.
* For the Jackson et al. 2020 media, the stored per-cell classes
  (`cdc_class`, cohort of maximal fold change) give a G1 share that
  does not grow with the doubling time (Spearman -0.45 over 9 growing
  media).

## Ergodic rate analysis (ERA)

Kafri et al. 2013 (Box 2): in an asynchronous population growing
exponentially at rate `a = log(2)/T` (doubling time `T`), the speed
along the cycle path `l` is `v(l) = a (2 - F(l))/f(l)`, with `f` the
density of cells along `l` and `F` its cumulative fraction from
division. Integrating `dt = dl/v` gives the time since division

    age = T * (1 - log2(2 - F))

so `F = 0` is age 0 and `F = 1` is age `T`. Only the order of the cells
enters, so `theta` and `phi` give the same age; the density is in `F`.
The same balance for a sub-population that receives a share `s` of the
newborn cells (both cells of each division count, `s = 1/2` for
daughters or mothers in budding yeast), and that holds a fraction `f`
of all cells before it leaves a phase, gives its time in that phase

    tau = -T * log2(1 - f/(2*s))

which is the age formula for `s = 1`.

Assumptions: steady state (balanced, asynchronous exponential growth),
a single path that every cell follows in the same direction, and
division as the origin. Budding yeast divides asymmetrically: mothers
and daughters both restart in G1, daughters stay longer; the age above
is then the population mean, and the group formula splits it.

Pilot (2026-10-30, stored `cdc_phi`, origin at the M/G1 cohort peak,
G1 ends at the S cohort peak): the ERA G1 share rises with doubling
time, Spearman +0.75 over the 9 growing media (YPD 0.22 of T, Proline
0.43), against -0.45 for the `cdc_class` shares.

## Plan

All four parts were approved. Parts 1-3 and the principal curve of
part 4 are implemented; count-aware fitting is open.

1. **Cell age** (`R/phasetime.R`).
   `era_age()` and `era_fraction()` (the formula and its inverse);
   `phase_time()` adds the age of each cell and cohort to a phases
   object, origin at a phase or a cohort peak (`state_phase`), cells
   counted in the direction of increasing phase; `phase_landmarks()`
   gives cohort peaks and points of maximal change (maximum of the
   smoothed `d(theta-phi)/dphi` between two peaks, from
   `get_segments`); `phase_durations()` turns landmarks into ages and
   durations. G1/S is reported both at the S peak and at the point of
   maximal change between the G1 and S peaks. Validation: G1 duration
   against doubling time over the 11 Jackson et al. 2020 media.
2. **Per-cell quality** (`R/quality.R`).
   `project_cells()` puts new cohort state columns into an existing
   PCA (exact for the cells it was built from); `cell_quality()` adds
   the unscaled amplitude (`amplitude * phases$scale`), the share of a
   cell's variance across cohorts in the PC1/PC2 plane (`r2`), and
   flags cells that are not above a Poisson noise reference simulated
   from their own totals (`noise_states()`: cohort counts drawn from
   the population mean fractions, without a cycle); `phase_bootstrap()`
   gives a per-cell phase uncertainty by resampling genes within
   cohorts.
3. **Variability of durations** (`R/phasetime.R`).
   `group_durations()` splits a phase by a cell grouping (G1 by the
   `daughter` and `aged` cohort states) with the group formula above;
   `age_profile()` bins any per-cell value along age (mean, sd, cv),
   for the spread of log total counts as a size proxy, which should
   narrow towards G1 exit if size is controlled there; cells off the
   trajectory (`r2`, radius, `pphi`) against the noise reference.
4. **Geometry**.
   `principal_phase()` (`R/principal.R`) fits a closed principal curve
   (Fourier series in each coordinate, Hastie-Stuetzle iterations)
   through the cells in the first `k` PCs, so more than two PCs, and
   so more cohorts than five cohort means, can define the path; the
   arc-length phase `lambda` and the distance to the curve are
   returned, and `phase_time()` takes `phase="lambda"`.
   Open: count-aware (Poisson) phase fitting with uncertainty, i.e.
   maximum likelihood of each cell's counts along the fitted cohort
   profiles, which would replace both the angle and the bootstrap.

Example: `vignettes/phase_time.Rmd` (Jackson et al. 2020, YPD, and the
11-media check; Nadal-Ribelles et al. 2025), compiled to
`vignettes/phase_time.html`. Its Methods section gives the math of all
steps (states, PCA, ERA derivation, group durations, landmarks,
projection, noise reference, bootstrap, principal curve, daughter
calls). Step-by-step figures: `plotPhasePlane`, `plotPhaseOrder`,
`plotEra`, `plotAgeStates` (`R/plotphasetime.R`). Scripts:
`inst/scripts/phasetime_utils.R` (shared), `jackson20_phasetime.R`,
`nadal25_phasetime.R`.

## First results (2026-10-30, `inst/scripts/jackson20_phasetime.R`)

Phase from all periodic Cyclebase genes (`cdc_bins_*`), landmarks from
the full Spellman cohorts (`sp_*`); see the next section for why.

* Cohort peaks in cycle order in 10 of 11 media (not YPDDiauxic).
* G1 (M/G1 peak to G1S) does not take a larger share of T in slower
  media: G1/T 0.11-0.25 without noise cells, Spearman +0.19 with T over
  the 9 growing media; YPD 0.13 of T = 1.05 h. G1 lengthens about in
  proportion to T. (With the pipeline's YRO-filtered `cdc_*` cohorts
  the share rose with T, Spearman +0.78, carried by the media in wrong
  cohort order, whose G1 included M; +0.43 for the 6 in order.)
* S takes about the same time in all media (G1S to SG2, 0.2-0.4 h;
  Spearman with T +0.27 in hours, -0.49 as share of T, -0.83 for the 8
  media in order), so it is a smaller share of the cycle in slower
  media. The histone (S) peak is 0.17-0.29 h wide at half maximum, a
  share of T that falls with T (-0.87); the G2 peak likewise (-0.58).
  Widths include the per-cell phase uncertainty: upper bounds, not the
  cell-to-cell spread of S. G2 between points of maximal change jumps
  between media (0.00-0.37 of T); the S-to-G2 peak interval is steadier,
  0.11-0.21 of T, and scales with T (share vs T +0.24, -0.08 in order).
  9 media, one snapshot each; Proline's G1 and S peaks coincide
  without noise cells.
* Noise cells (amplitude not above Poisson cells of the same totals,
  B = 20): 3-21% in growing media, 42% in YPDDiauxic, 85% in CStarve.
* Per-cell uncertainty from gene resampling (YPD): median 0.19 of the
  cycle (0.20 h). Single cells are placed only roughly; durations come
  from the population order.
* Size proxy (YPD): mean log total is lowest from G1 to S; its spread
  is narrowest around the G2 peak, not at G1 exit.
* The `daughter` cohort is high in ~90% of cells right after division
  and again in late M: it marks M/G1, not daughters, so mother and
  daughter G1 (`group_durations`) need another marker, or the daughter
  counts relative to age-matched cells.
* Nadal-Ribelles et al. 2025 (`nadal25_phasetime.R`, growth rate
  unknown, ages in units of the cycle): control (YPD) cohorts in order,
  G1 0.12, S 0.21 of the cycle, as in Jackson YPD; 36% noise cells.
  After 15 min of 0.4 M NaCl (not a steady state) G1 takes 0.33. Their
  published S class agrees with our S peak (76-93% of cells at ages
  0.15-0.35); their G1 class also covers much of the second half of
  the cycle, perhaps cells without a clear S or G2M signal.
* The last interval, from the G2/M point to the next M/G1 peak, is long
  (0.3-0.65 of T in Jackson, 0.54 in Nadal-Ribelles control): M and
  early G1 are poorly resolved by five cohorts, or the M/G1 peak lies
  after division. Open: finer bins there, and an origin check.
* The principal curve (k = 3 or 4) did not repair wrong cohort orders
  with the filtered cohorts (it fixed YPDDiauxic, broke YPDRapa and
  AmmoniumSulfate); the cause was the gene filter.

## Which cohorts for the division-cycle phase (2026-10-30)

The B8 filter removes from the Spellman cohorts every gene of the
growth cohorts (`cgc_` Ribi, RP, AA, mRP, S/C) and of HOC/LOC: 3383
genes, half the genome and 97% of the top 1000 Cyclebase genes. It
leaves 376 of 790 Spellman genes (S: 25 of 71). Compared over the 11
media (scratch analysis, cohort states recomputed from counts):

| cohorts | genes | order wrong | split-half R | noise cells (YPD) |
|---|---|---|---|---|
| `cdc_` Spellman, filtered | 376 | 6 media | 0.17-0.41 | 30% |
| `sp_` Spellman, full | 790 | 4 | 0.41-0.77 | 10% |
| `cdc_bins_` Cyclebase p<.05, 5 bins | 3080 | 2 (YPDDiauxic, CStarve) | 0.33-0.77 | 21% |
| Cyclebase top 300, 5 bins, filtered | 300 | 6 | 0.13-0.44 | 68% |
| Cyclebase top 1500, 12 bins | 1500 | (neighbour swaps) | 0.31-0.84 | 12% |

Split-half R: phases from two random halves of each cohort's genes,
circular agreement (1 = identical). The filter halves the
reproducibility and triples the noise cells, and it explains the wrong
orders: the strongest cycling genes are also in the growth cohorts.
Adopted (owner, 2026-10-30): the phase from all periodic Cyclebase
genes, landmarks from the full Spellman cohorts' peaks along it. Open:
more bins; where a growth cohort is to be read against the phase,
leave only that cohort's genes out of the phase, not all of HOC/LOC.

## Angle or principal curve? (tested 2026-10-30)

`curve_vs_angle.R` (session notes): the 9 growing Jackson media and
Nadal-Ribelles ctrl, phase from all periodic Cyclebase genes in 5, 12
or 20 peak-time bins, by the PC1/PC2 angle (`get_pseudophase`) or the
closed principal curve in min(bins-1, 4) PCs (`principal_phase`),
noise cells removed. Means over the 10 data sets:

| bins | method | order ok | curve collapsed | split-half R | Cyclebase R | G1 CI width / T | S CI width / T |
|---|---|---|---|---|---|---|---|
| 5 | angle | 0.9 | - | 0.71 | 0.76 | 0.13 | 0.05 |
| 12 | angle | 1.0 | - | 0.68 | 0.72 | 0.13 | 0.04 |
| 20 | angle | 0.9 | - | 0.67 | 0.67 | 0.12 | 0.06 |
| 5 | curve | 0.8 | 0.4 | 0.49 | 0.72 | 0.11 | 0.10 |
| 12 | curve | 0.7 | 0.3 | 0.51 | 0.70 | 0.09 | 0.08 |
| 20 | curve | 0.8 | 0.6 | 0.51 | 0.65 | 0.14 | 0.14 |

Split-half R: agreement of the cell order from two halves of the genes.
Cyclebase R: circular agreement of the peaks of 20 bins of held-out
genes along the age (phase from the other half) with their Cyclebase
peak times, an independent time axis. The angle is better on every
measure; the curve collapses onto a few points in 30-60% of the runs
(now caught: `principal_phase` stops and keeps the previous phase),
and doubles the S intervals. More bins do not help the angle either:
fewer genes per bin and the plane cannot use the extra structure.
This is our simple curve (Fourier series on the scaled loadings); a
curve in cohort-state space, or a smoothing-spline curve, might do
better, but the evidence does not support a switch.

**Decision: keep `get_pseudophase` (the angle, 5 bins) as the core.**
Redesign around it rather than replacing it:

* `division_phase` as the front end, with `method=` choosing the
  ordering ("angle" default; "curve" experimental; later "poisson");
  ERA, landmarks, CIs, quality and plots stay the same for all;
* the phase tables apart from the `prcomp` object (README TODO);
* finer resolution (M to early G1) through count-aware fitting (step 4
  below), which can use many cohorts or genes without the PCA plane,
  rather than through more bins in the angle.

## Benchmark against other tools (planned, 2026-10-30)

Few cell cycle inference tools have been tested on yeast; most were
built and validated on mammalian cells. Candidates, on the same yeast
data (Jackson 2020 media, Nadal-Ribelles 2025 ctrl, Gasch 2017):

* discrete scores: Seurat `CellCycleScoring`, scran `cyclone` (classes
  only; Seurat calls G1 when neither score is positive);
* tricycle (Zheng et al. 2022): projection onto a reference PCA,
  continuous angle; needs a yeast reference (orthologs of the mouse
  reference, or one built from Cyclebase genes);
* Revelio, reCAT, Cyclum: continuous orderings from markers, graphs or
  autoencoders; peco needs FUCCI training data (none for yeast);
* RNA velocity (DeepCycle, VeloCycle, which estimates cycle speed
  directly): expected to fail in yeast, ~5% of genes have introns;
  check the unspliced fraction first.

Measures: cohort peak order, split-half reproducibility (phase from two
halves of the genes), agreement with published classes (Nadal-Ribelles)
and orderings (the Jackson ordering files read by the owner's
pipeline), noise-cell fraction, and the ERA durations on top of each
ordering (ERA is independent of the ordering method). External
reference for durations: G1 or unbudded fractions by imaging or flow
cytometry in the same media, from the literature.

## Next (before returning to the growth law / YRO paper)

Division point (the long M to M/G1 interval):

1. Finer cohorts in M and early G1: Cyclebase bins at 60-100% of the
   cycle, and known clusters in order, CLB2 (G2/M), Mcm1/ECB (late M,
   early G1: CLN3, SWI4, CDC6, MCM genes), Ace2/Swi5 (cytokinesis and
   daughter G1: SIC1, ASH1, CTS1, DSE1-4, EGT2), SBF/MBF (late G1:
   CLN1/2, histones in S).
2. Size as a division clock: total mRNA should halve at division (mother
   and daughter share it). The steepest drop of mean log total along
   the phase (we see its maximum before and its minimum after the M/G1
   peak) marks division independently of the cohorts; capture
   efficiency adds noise, so check against the daughter cohort onset.
3. Signal amplitude: the unscaled amplitude and the noise fraction along
   phase; a dip where cohorts change over may mark the transition.
4. Count-aware phase: fit each cell's counts as Poisson along smooth
   cohort (or gene) profiles of phase, maximum likelihood with an
   uncertainty per cell; more cohorts resolve short phases without the
   PCA plane.

First tests of 1 and 2 (2026-10-30, Jackson YPD and Nadal-Ribelles
ctrl; script `division.R` in the session notes, not in the package):

* Cyclebase genes in 20 peak-time bins: their peak ages along our ERA
  age follow the Cyclebase peak times (from synchronized cultures)
  with slope 1.03 (YPD) and 1.00 (Nadal-Ribelles), up to the wrap at
  the last bin; Spearman 0.99 in Nadal-Ribelles. So the ERA age agrees
  with an independent time axis, and the long interval after the G2/M
  point is real: Cyclebase also places 35-40% of the cycle after the
  CLB2 wave (peak times 60-100%).
* Marker waves along age, YPD: ECB (Mcm1) 0.04 and Ace2/Swi5 0.06
  of T, right after the M/G1 origin; SBF/MBF 0.18, histones 0.24, CLB2
  0.55. Nadal-Ribelles alike, except a flat ECB cohort. The M/G1 peak
  is a good proxy for division.
* Size clock: mean log total does not halve at any point (range about
  0.2 in ln, not ln 2); it is highest right after the origin, falls
  through G1 to a minimum in S, then rises. Its steepest drop lies
  0.06 (YPD) to 0.17 (Nadal-Ribelles) of T after the M/G1 peak; if
  mother-bud pairs are captured as one cell until separation, that is
  where cells separate. Too weak and noisy to set the origin; total
  counts mix size with capture.

Tools for the paper:

Steps 5-8 implemented (2026-10-30, `R/division.R`): `division_phase`
(5, with `exclude` for 7 and `noise.rm`), `cohort_states`,
`duration_ci` (6), `time_weights`/`time_mean` (8); the scripts use
them. First results (vignette):

* 95% intervals from resampling cells: S is tight in the large sets
  (YPD about +-0.01 h) and still a falling share of T; G1 intervals
  are wide in the small sets (MinimalGlucose about 0.03-0.20 of T,
  Proline down to 0.05 of T), so a G1 trend over the 9 media cannot be
  read from these data.
* Growth cohorts in YPD, each against a phase without its genes (G1
  0.11-0.15 h, as with all genes): S/C peaks right after division
  (+14%) and is lowest at G1/S and in S (-12%); Ribi is low from late
  G1 through S (-9%) and rises through G2/M; RP, AA and mRP vary by
  less than 5%.
* Snapshot vs time averages differ little for these cohorts (share of
  cells above the median 0.50, of the cycle 0.50-0.52): the cohorts
  vary little along age, so the ERA weights matter for durations more
  than for cohort averages.

5. A package function for the whole chain (states, phase from Cyclebase
   genes, landmarks from Spellman cohorts, orientation, ERA age, noise
   flags), replacing `inst/scripts/phasetime_utils.R`.
6. Confidence intervals for durations by resampling cells (and media
   replicates where present); needed before reading trends over 9 media.
7. Leave-cohort-out phases: a growth cohort (Ribi, RP, S/C) read
   against a phase that does not contain its genes, rather than
   removing all HOC/LOC genes.
8. Growth cohorts along ERA age, not rank phase: time shares of YRO
   states within the division cycle (the time-share balance of the
   paper) need time weights, which the ERA age gives; noise cells
   excluded or down-weighted.
9. Regenerate the pipeline's division pseudophase
   (genomeBrowser `clusterTranscriptomes_pseudophase.R`, B8 fix) with
   the new cohorts, and run Gasch 2017.

## Not here

The Heineike data sets in the pseudophase cache are unpublished and
stay out of public outputs, including this package's vignettes.
